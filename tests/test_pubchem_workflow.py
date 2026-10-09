"""Independent provider outputs, curated references and shared downstream inputs."""
import asyncio
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock,patch

import pandas as pd
from crawlers.pubchem import PubChemSimilarMols
from biomolexplorer.catalog import new_stage
from biomolexplorer.flow import connect,input_ports
from biomolexplorer.operations import validate_operation
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.pubchem_retrieval import retrieve_pubchem,reference_choices
from biomolexplorer.ui.guided import GuidedForm
from biomolexplorer.ui.file_selection import FileSelection
import test_pipeline_execution as execution


class PubChemRetrievalTests(unittest.TestCase):
    def extra(self):
        return pd.DataFrame([{'molecule_chembl_id':'PUBCHEM3','canonical_smiles':'CCC','source':'PubChem'}])

    def test_uploaded_aliases_all_or_single_produce_only_pubchem_compounds(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);(root/'mine.csv').write_text('name,smiles\nA,CCO\nB,CCN\n')
            self.assertEqual(reference_choices([root/'mine.csv']),{'A':'CCO','B':'CCN'})
            for mode,expected in (('all',['A','B']),('single',['B'])):
                with patch.object(PubChemSimilarMols,'search',return_value=self.extra()) as search:
                    result=retrieve_pubchem(str(root/mode),base_input_path=str(root),input_file='mine.csv',selection_mode=mode,compound_id='B')
                self.assertEqual(search.call_args.args[0]['molecule_chembl_id'].tolist(),expected)
                self.assertEqual(result['molecule_chembl_id'].tolist(),['PUBCHEM3'])
                output=pd.read_csv(root/mode/'compounds/PubChem/compounds.csv')
                self.assertEqual(output['source'].tolist(),['PubChem'])
                report=json.loads((root/mode/'retrieval_report.json').read_text())
                self.assertEqual([r['molecule_chembl_id'] for r in report['references']],expected)
            with self.assertRaisesRegex(ValueError,'não está'):
                retrieve_pubchem(str(root/'missing'),base_input_path=str(root),input_file='mine.csv',selection_mode='single',compound_id='MISSING')
            with self.assertRaises(ValueError):retrieve_pubchem(str(root/'bad'),base_input_path=str(root),input_file='../mine.csv')

    def test_manual_smiles_needs_no_upstream_block(self):
        with tempfile.TemporaryDirectory() as folder,patch.object(PubChemSimilarMols,'search',return_value=self.extra()) as search:
            result=retrieve_pubchem(folder,reference_source='manual',reference='OCC')
            self.assertEqual(search.call_args.args[0]['canonical_smiles'].tolist(),['CCO'])
            self.assertEqual(result['source'].tolist(),['PubChem'])
        stage=new_stage('retrieve_pubchem');stage['parameters'].update(reference_source='manual',reference='CCO')
        self.assertEqual(input_ports(stage),[])
        validate_operation(stage['operation'],stage['parameters'])
        stage['parameters']['reference_type']='cid';stage['parameters']['reference']='invalid'
        with self.assertRaisesRegex(ValueError,'CID PubChem'):validate_operation(stage['operation'],stage['parameters'])

    def test_cid_and_unambiguous_name_resolve_without_assuming_first_hit(self):
        with tempfile.TemporaryDirectory() as folder:
            crawler=PubChemSimilarMols(folder)
            try:
                def request(endpoint,data,params=None):
                    if endpoint=='name/cids/JSON':return {'IdentifierList':{'CID':[7]}}
                    self.assertEqual(data,{'cid':'7'})
                    return {'PropertyTable':{'Properties':[{'CID':7,'SMILES':'CCO'}]}}
                crawler._request=Mock(side_effect=request)
                for kind,value in (('cid','7'),('name','ethanol')):
                    result=crawler.reference(kind,value)
                    self.assertEqual(result['molecule_chembl_id'].tolist(),['PUBCHEM7'])
                crawler._request=Mock(return_value={'IdentifierList':{'CID':[7,8]}})
                with self.assertRaisesRegex(ValueError,'único composto'):crawler.reference('name','ambiguous')
            finally:crawler.session.close()

    def test_reference_controls_block_inactive_modes_and_drop_unused_connections(self):
        source=new_stage('retrieve_compounds');stage=new_stage('retrieve_pubchem')
        connect([source,stage],source['id'],stage['id'],'base_input_path')
        ui=SimpleNamespace(current={'pipeline':[source,stage]},page=SimpleNamespace(update=lambda:None),
            pubchem_references={'stage:'+source['id']:{'compounds.csv':{'A':'CCO','B':'CCN'}}})
        form=GuidedForm(ui,stage,[],True);form.layout()
        mode=form.field_controls['selection_mode'];mode.value='single';mode.on_select(SimpleNamespace(control=mode))
        self.assertTrue(form.field_controls['reference'].disabled)
        self.assertFalse(form.field_controls['compound_id'].disabled)
        self.assertEqual({o.key for o in form.field_controls['compound_id'].options},{'A','B'})
        form.field_controls['compound_id'].value='B'
        self.assertEqual(form.read()['input_processing'],'merge')
        origin=form.field_controls['reference_source'];origin.value='manual';origin.on_select(SimpleNamespace(control=origin))
        self.assertTrue(form.field_controls['base_input_path'].disabled)
        self.assertFalse(form.field_controls['reference'].disabled)
        self.assertTrue(form.processing.disabled)
        form.field_controls['reference'].value='C/C=C\\C'
        saved=form.read()
        self.assertFalse(saved['bindings']);self.assertNotIn('base_input_path',saved['parameters'])
        self.assertEqual(saved['parameters']['reference'],'C/C=C\\C')

    def test_popup_only_lists_compounds_in_checked_files_and_preserves_selection(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);a=root/'first.csv';b=root/'second.csv'
            a.write_text('molecule_chembl_id,canonical_smiles\nA,CCO\n');b.write_text('name,smiles\nB,CCN\n')
            source=new_stage('retrieve_compounds');stage=new_stage('retrieve_pubchem')
            connect([source,stage],source['id'],stage['id'],'base_input_path')
            pending={'name':stage['name'],'configuration':stage}
            run={'stages':[dict(source,status='succeeded',artifacts=[str(a),str(b)])]}
            form=FileSelection(run,pending)
            check=form.rows['base_input_path'][1][0];check.value=True;check.on_change(SimpleNamespace(control=check))
            form.selection_mode.value='single';form.selection_mode.on_select(SimpleNamespace(control=form.selection_mode))
            self.assertEqual({o.key for o in form.compound_id.options},{'B'})
            form.compound_id.value='B'
            selected=form.read();self.assertEqual(selected['parameters']['compound_id'],'B')
            self.assertEqual(selected['input_processing'],'merge')
            reopened=FileSelection(run,pending,state=form.state())
            self.assertEqual(reopened.compound_id.value,'B')
            check.value=False;check.on_change(SimpleNamespace(control=check))
            self.assertEqual(form.compound_id.value,'')
            with self.assertRaisesRegex(ValueError,'composto de referência'):form.read()


class PubChemPipelineTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    upload=execution.PipelineExecutionTests.upload
    imported=execution.PipelineExecutionTests.imported
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait

    def launch_fixture(self,args,**kwargs):
        fixture=Path(__file__).parent/'fixtures/chembl_worker.py'
        return self.original_launch([args[0],str(fixture),*args[3:]],**kwargs)

    def confirm(self,service,run,compound_id=None,processing=None):
        pending=next(s for s in run['stages'] if s['status']=='awaiting_input')
        form=FileSelection(run,pending)
        for rows in form.rows.values():
            for check,_ in rows:check.value=True
        if compound_id:
            form.refresh_references();form.selection_mode.value='single'
            form.selection_mode.on_select(SimpleNamespace(control=form.selection_mode));form.compound_id.value=compound_id
        if processing:form.processing.value=processing
        return self.wait(service.resume(self.token,run['id'],form.read()))

    def test_real_worker_manual_reference_produces_shared_compound_contract(self):
        stage=new_stage('retrieve_pubchem');stage['parameters'].update(reference_source='manual',reference='CCO')
        self.save([stage]);self.original_launch=subprocess.Popen
        with patch('biomolexplorer.jobs.subprocess.Popen',side_effect=self.launch_fixture):
            service=PipelineService(self.store,worker_python=sys.executable,cpu_workers=1)
            try:
                run=self.wait(service.submit(self.token,self.project_id))
                self.assertEqual(run['status'],'succeeded',run['error'])
                dataset=next(Path(p) for p in run['stages'][0]['artifacts'] if '/compounds/' in p and p.endswith('compounds.csv'))
                with dataset.open() as stream:self.assertEqual([r['molecule_chembl_id'] for r in csv.DictReader(stream)],['PUBCHEM3'])
            finally:service.close()

    def test_real_worker_accepts_an_uploaded_csv_without_a_retrieval_block(self):
        asset=self.upload('my_references.csv','name,smiles\nA,CCO\nB,CCN\n')
        stage=new_stage('retrieve_pubchem');stage['bindings']['base_input_path']={'asset':asset}
        stage['parameters'].update(selection_mode='single',compound_id='B');stage['input_processing']='merge'
        self.save([stage]);self.original_launch=subprocess.Popen
        with patch('biomolexplorer.jobs.subprocess.Popen',side_effect=self.launch_fixture):
            service=PipelineService(self.store,worker_python=sys.executable,cpu_workers=1)
            try:
                run=self.wait(service.submit(self.token,self.project_id))
                self.assertEqual(run['status'],'succeeded',run['error'])
                report=next(Path(p) for p in run['stages'][0]['artifacts'] if p.endswith('retrieval_report.json'))
                self.assertEqual([r['molecule_chembl_id'] for r in json.loads(report.read_text())['references']],['B'])
            finally:service.close()

    def test_real_workers_curate_one_reference_and_analyze_sources_separately_or_merged(self):
        asset=self.upload('mine.csv','name,smiles\nA,CCO\nB,CCN\n');source=self.imported(asset)
        pubchem=new_stage('retrieve_pubchem');connect([source,pubchem],source['id'],pubchem['id'],'base_input_path')
        downstream=new_stage('fingerprints')
        stages=[source,pubchem,downstream]
        for origin in (source,pubchem):connect(stages,origin['id'],downstream['id'],'base_input_path')
        self.original_launch=subprocess.Popen
        with patch('biomolexplorer.jobs.subprocess.Popen',side_effect=self.launch_fixture):
            service=PipelineService(self.store,worker_python=sys.executable,cpu_workers=1)
            try:
                for processing,count in (('individual',2),('merge',1)):
                    self.save(stages)
                    run=self.wait(service.submit(self.token,self.project_id,reuse_results=False))
                    self.assertEqual(run['status'],'awaiting_input',run['error'])
                    run=self.confirm(service,run,compound_id='B')
                    self.assertEqual(run['status'],'awaiting_input',run['error'])
                    report=next(Path(p) for p in run['stages'][1]['artifacts'] if p.endswith('retrieval_report.json'))
                    self.assertEqual([r['molecule_chembl_id'] for r in json.loads(report.read_text())['references']],['B'])
                    run=self.confirm(service,run,processing=processing)
                    self.assertEqual(run['status'],'succeeded',run['error'])
                    self.assertEqual(len(run['stages'][2].get('batches',[run['stages'][2]])),count)
                    data=[Path(p) for p in run['stages'][2]['artifacts'] if p.endswith('.csv')]
                    records=[]
                    for p in data:
                        with p.open() as stream:records.extend(r['molecule_chembl_id'] for r in csv.DictReader(stream))
                    self.assertEqual(set(records),{'A','B','PUBCHEM3'})
            finally:service.close()
