"""Independent per-input partitions, named rules and authorized alert review."""
import asyncio
import csv
import io
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import pandas as pd
from wrappers.admet import ADMETWrapper
from biomolexplorer.artifact_naming import assign_block_ids,publish
from biomolexplorer.artifact_choices import produced_files
from biomolexplorer.artifacts import records
from biomolexplorer.catalog import new_stage
from biomolexplorer.workspace import AccessDenied
import test_results_display as fixtures
from test_frontend_navigation import Page

class StandardADMETOutputsTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.source=self.root/'source';self.source.mkdir()
        self.data=pd.DataFrame([dict(molecule_chembl_id='NITRO',canonical_smiles='O=[N+]([O-])c1ccccc1'),
            dict(molecule_chembl_id='ETHANOL',canonical_smiles='CCO'),dict(molecule_chembl_id='HEAVY',canonical_smiles='C'*45)])
        self.data.to_csv(self.source/'B01__chembl__mol.csv',index=False)

    def run_admet(self,**kwargs):
        wrapper=ADMETWrapper(str(self.root/'output'),str(self.source),admet_profile='custom',**kwargs)
        with patch.object(wrapper,'generate_plot'):wrapper.run_pipeline()
        return wrapper

    def test_each_input_has_only_four_barrier_partitions_one_named_rule_and_one_alert_table(self):
        for name in ('B01__chembl__sim','pubchem_BACE1','zinc_BCHE'):
            self.data.to_csv(self.source/(name+'.csv'),index=False)
        wrapper=self.run_admet(bbb_target='BBB+',hia_target='HIA+',admet_mode='filter',rule_set='lipinski',alert_catalogs=['PAINS','BRENK'])
        files=list((self.root/'output').glob('*.csv'));self.assertEqual(len(files),24)
        for source in self.source.glob('*.csv'):
            base=source.stem
            self.assertEqual({p.name for p in files if p.name.startswith(base+'_')},
                {base+'_'+role+'.csv' for role in ('BBB+','BBB-','HIA+','HIA-','rules_lipinski','alerts')})
            for endpoint in ('BBB','HIA'):
                positive=pd.read_csv(self.root/'output'/(base+'_'+endpoint+'+.csv'))
                negative=pd.read_csv(self.root/'output'/(base+'_'+endpoint+'-.csv'))
                plus=set(positive.molecule_chembl_id);minus=set(negative.molecule_chembl_id)
                self.assertFalse(plus&minus);self.assertEqual(plus|minus,set(self.data.molecule_chembl_id))
            rules=pd.read_csv(self.root/'output'/(base+'_rules_lipinski.csv'))
            self.assertTrue((rules.Lipinski_violations<=0).all())
            self.assertEqual(set(rules.Rule_set),{'lipinski'})
            with (self.root/'output'/(base+'_alerts.csv')).open() as stream:
                alerts=list(csv.DictReader(stream))
            self.assertEqual({r['molecule_chembl_id'] for r in alerts},set(self.data.molecule_chembl_id))
            self.assertTrue(all(r['Alert_PAINS'] in ('true','false') and r['Alert_BRENK'] in ('true','false') for r in alerts))
            self.assertEqual(next(r for r in alerts if r['molecule_chembl_id']=='NITRO')['Alert_BRENK'],'true')
        self.assertTrue(any(p.name.endswith('BBB-.csv') and len(pd.read_csv(p)) for p in files))

    def test_rule_names_and_publication_survive_short_names_and_reopening(self):
        self.run_admet(rule_set='oral')
        stage=new_stage('admet');assign_block_ids([stage])
        output=self.root/'output';files=[str(p) for p in output.iterdir()]
        published=publish(stage,files,output)
        names={Path(p).name for p in published}
        self.assertIn('B01__admet__rules-lipinski-veber.csv',names)
        self.assertIn('B01__admet__alerts.csv',names)
        item=dict(stage,configuration=stage,artifacts=published)
        self.assertEqual(len([p for p in produced_files(item) if p.suffix=='.csv']),6)
        labels=[r['output'] for r in records(item)]
        self.assertIn('rules-lipinski-veber',labels)
        report=json.loads(next(Path(p) for p in published if p.endswith('.json') and 'report' in p).read_text())
        self.assertEqual(report['schema_version'],3)
        self.assertTrue(all(Path(p).exists() for p in report['outputs'].values()))

    def test_disabled_rules_and_alert_catalogs_do_not_generate_extra_csvs(self):
        self.run_admet(rule_set='none',alert_catalogs=[])
        self.assertEqual(len(list((self.root/'output').glob('*.csv'))),4)

    def test_each_supported_rule_uses_one_named_file_and_its_own_violation_limit(self):
        for rule,label,columns in [('lipinski','lipinski',['Lipinski']),('veber','veber',['Veber']),
                                  ('oral','lipinski-veber',['Lipinski','Veber']),('fragment','fragment',['Fragment'])]:
            with self.subTest(rule=rule):
                output=self.root/rule
                wrapper=ADMETWrapper(str(output),str(self.source),admet_profile='custom',rule_set=rule,
                    max_rule_violations=1,alert_catalogs=[])
                with patch.object(wrapper,'generate_plot'):wrapper.run_pipeline()
                files=list(output.glob('*.csv'));self.assertEqual(len(files),5)
                frame=pd.read_csv(output/('B01__chembl__mol_rules_'+label+'.csv'))
                self.assertFalse(frame.empty)
                for column in columns:
                    allowed=0 if rule=='oral' and column=='Veber' else 1
                    self.assertTrue((frame[column+'_violations']<=allowed).all())

class AlertReviewTests(unittest.TestCase):
    setUp=fixtures.ResultsDisplayTests.setUp
    tearDown=fixtures.ResultsDisplayTests.tearDown

    def add_alerts(self):
        path=self.admet_csv.with_name('admet_alerts.csv')
        path.write_text('molecule_chembl_id,canonical_smiles,TPSA,WLOGP,BBB,HIA,Alert_PAINS,Alert_BRENK,StructuralAlerts\n'+
            ''.join(f'MOL{i},CCO,20,1,BBB+,HIA+,false,'+('true,BRENK: fixture' if i%2 else 'false,')+'\n' for i in range(31)))
        with self.store.connect() as db:
            run=self.store.get_run(self.token,self.rid)
            run['stages'][1]['artifacts'].append(str(path))
            db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(run['stages']),self.rid))
        return path

    def test_page_is_per_input_authorized_boolean_and_paginated(self):
        path=self.add_alerts();args=(self.token,self.pid,self.rid,self.admet['id'],str(self.admet_csv))
        first=self.service.admet_alert_page(*args)
        self.assertEqual(first['path'],str(path));self.assertEqual(first['total'],31);self.assertEqual(len(first['rows']),25)
        self.assertFalse(first['rows'][0]['flags']['BRENK']);self.assertTrue(first['rows'][1]['flags']['BRENK'])
        second=self.service.admet_alert_page(*args,25);self.assertEqual(len(second['rows']),6)
        search=self.service.admet_alert_page(*args,query='MOL30');self.assertEqual(search['total'],1)
        with self.assertRaises(AccessDenied):self.service.admet_alert_page(self.guest,self.pid,self.rid,self.admet['id'],str(path))
        with self.assertRaises(AccessDenied):self.service.admet_alert_page(*args[:-1],str(self.compounds))

    def test_popup_pages_and_opens_selected_molecule_without_losing_the_table(self):
        import flet as ft
        from biomolexplorer.ui.admet_alerts import ADMETAlertTable
        async def scenario():
            self.add_alerts();page=Page()
            async def call(function,*args,**kwargs):return function(*args,**kwargs)
            async def guard(action):await action()
            ui=SimpleNamespace(page=page,call=call,guard=guard)
            result=SimpleNamespace(ui=ui,service=self.service,token=self.token,project_id=self.pid,run_id=self.rid,stage=self.admet,valid=lambda:True)
            popup=ADMETAlertTable(result,str(self.admet_csv));await popup.open()
            self.assertTrue(popup.dialog.open);self.assertEqual(len(popup.table.rows),25)
            await popup.next.on_click(None);self.assertEqual(len(popup.table.rows),6)
            await popup.previous.on_click(None)
            await popup.table.rows[0].cells[-1].content.content.on_click(None)
            self.assertTrue(popup.dialog.open)
            preview=page.dialogs[-1];self.assertTrue(isinstance(preview.content,ft.Image))
            self.assertEqual(preview.title.value,'MOL0')
        asyncio.run(scenario())
