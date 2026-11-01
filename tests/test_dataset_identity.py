"""Target provenance and dataset names survive batching and downstream reads."""
import asyncio
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd
import test_redocking_results as fixtures
from biomolexplorer.artifacts import index_run
from biomolexplorer.catalog import new_stage
from biomolexplorer.dataset_identity import admet_origin,output_stem,redocking_targets
from biomolexplorer.operations import execute_operation
from biomolexplorer.redocking_results import RedockingResults
from biomolexplorer.result_files import ResultFiles
from biomolexplorer.result_origins import ResultOrigins
from biomolexplorer.ui.redocking_results import RedockingResultsTable


def source(operation,name,target,paths):
    config=new_stage(operation);config['name']=name
    if target:config['parameters']['target']=target
    return dict(id=config['id'],name=name,operation=operation,configuration=config,
                status='succeeded',artifacts=list(map(str,paths)))


class IdentityTests(unittest.TestCase):
    def test_known_file_tags_and_conservative_unknown_inputs(self):
        for file,tag in [('molecules.csv','MOL'),('similars.csv','SIM'),('ACHE_MOL.csv','MOL'),
            ('BCHE_SIM_BBB+.csv','SIM'),('/data/PubChem/similars/ref.csv','PUBCHEM'),
            ('/data/ZINC/compounds.csv','ZINC'),('compounds.csv','INPUT')]:
            self.assertTrue(output_stem(file).endswith('_'+tag))

    def test_transitive_admet_origin_distinguishes_chembl_pubchem_and_zinc(self):
        for operation,file,tag in [('retrieve_compounds','molecules.csv','MOL'),
            ('retrieve_compounds','similars.csv','SIM'),('retrieve_pubchem','compounds.csv','PUBCHEM'),
            ('retrieve_zinc','compounds.csv','ZINC')]:
            producer=source(operation,'ACHE',None,['/data/'+file])
            origins=ResultOrigins([producer],'run')
            identity=admet_origin(['/data/'+file],origins)
            self.assertEqual(identity,dict(label='ACHE',tag=tag))
            self.assertTrue(output_stem(file,identity).startswith('ACHE_'))

    def test_pubchem_branch_does_not_inherit_chembl_molecule_tag(self):
        chembl=source('retrieve_compounds','ACHE',None,['/molecules.csv'])
        pubchem=source('retrieve_pubchem','PubChem expansion',None,['/pubchem/compounds.csv'])
        pubchem['configuration']['bindings']={'base_input_path':{'stage':chembl['id'],'selector':'molecules.csv'}}
        pubchem['input_files']={'base_input_path':['/molecules.csv']}
        index_run([chembl,pubchem],'run')
        self.assertEqual(admet_origin(['/pubchem/compounds.csv'],ResultOrigins([chembl,pubchem],'run')),
            dict(tag='PUBCHEM',label='ACHE'))

    def test_merge_keeps_all_labels_and_declares_mixed_collections(self):
        a=source('retrieve_compounds','ACHE',None,['/a/molecules.csv'])
        b=source('retrieve_compounds','BCHE',None,['/b/similars.csv'])
        identity=admet_origin(['/a/molecules.csv','/b/similars.csv'],ResultOrigins([a,b],'r'))
        self.assertEqual(identity,dict(tag='MIX',label='ACHE + BCHE'))
        self.assertIn('ACHE_BCHE',output_stem('selected_compounds.csv',identity).replace('_+_','_'))

    def test_ambiguous_metadata_does_not_assign_all_targets_to_each_pdb(self):
        routes=[[dict(operation='retrieve_structures',target=t,path='/'+t+'/pdb_codes.csv'),dict(operation='redocking')]
            for t in ('ACHE','BCHE')]
        self.assertEqual(redocking_targets(routes,'1ABC'),[])

    def test_names_are_safe_bounded_and_do_not_duplicate_trailing_tags(self):
        self.assertEqual(output_stem('molecules_MOL.csv'),'molecules_MOL')
        self.assertNotEqual(output_stem('A B.csv'),output_stem('A_B.csv'))
        name=output_stem('compounds.csv',dict(label='../ACHE / '+('long '*100),tag='PUBCHEM'))
        self.assertNotIn('/',name);self.assertFalse(name.startswith('.'));self.assertLess(len(name),150)


class RedockingTargetTests(unittest.TestCase):
    def setUp(self):
        self.f=fixtures.RedockingResultsTests();self.f.setUp();self.addCleanup(self.f.tearDown)

    def test_merged_metadata_resolves_each_pdb_using_its_retrieval_block(self):
        stage=self.f.stage();root=self.f.store.project_dir(self.f.pid)/'source';root.mkdir()
        paths=[]
        for code in ('1ABC','2DEF'):
            path=root/(code+'.pdb');path.write_text('ATOM');paths.append(path)
        a=source('retrieve_structures','PDB ACHE','ACHE',[paths[0]])
        b=source('retrieve_structures','PDB BCHE','BCHE',[paths[1]])
        metadata=Path(next(p for p in stage['artifacts'] if p.endswith('pdb_codes.csv')))
        metadata.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RMSD\n1ABC,LIG,1,A,0.2\n2DEF,LIG,1,A,0.3\n')
        stage['artifacts']=[str(metadata)]
        stage['configuration']['bindings']={'base_input_path':{'sources':[{'stage':a['id']},{'stage':b['id']}]}}
        stage['input_files']={'base_input_path':list(map(str,paths))}
        stages=[a,b,stage];index_run(stages,self.f.rid)
        with self.f.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),self.f.rid))
        rows=RedockingResults(self.f.store).simulations(self.f.token,self.f.pid,self.f.rid,stage['id'])
        self.assertEqual([(r['pdb'],r['target']) for r in rows],[('1ABC','ACHE'),('2DEF','BCHE')])
        import csv,io
        service=RedockingResults(self.f.store)
        exported=list(csv.DictReader(io.StringIO(service.rmsd_csv(self.f.token,self.f.pid,self.f.rid,stage['id']).decode())))
        self.assertEqual([r['Target'] for r in exported],['ACHE','BCHE'])
        self.assertEqual([r['best_rmsd_for'] for r in rows],[['ACHE'],['BCHE']])

    def test_imported_explicit_target_is_supported_without_invented_provenance(self):
        stage=self.f.stage();metadata=Path(next(p for p in stage['artifacts'] if p.endswith('pdb_codes.csv')))
        metadata.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RMSD,Target\n1ABC,LIG,1,A,0.2,BACE1\n')
        with self.f.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps([stage]),self.f.rid))
        rows=RedockingResults(self.f.store).simulations(self.f.token,self.f.pid,self.f.rid,stage['id'])
        self.assertEqual(rows[0]['target'],'BACE1')

    def test_target_filter_search_best_ties_and_global_pagination(self):
        ui=self.f.ui();result=SimpleNamespace(ui=ui,valid=lambda:True,writable=False)
        rows=[dict(id=str(i),pdb=f'{i:04}',ligand='LIG',residue='1',chain='A',rmsd=score,
            target=target,targets=[target],status='succeeded')
            for i,(target,score) in enumerate([('ACHE',0.2),('ACHE',0.1),('ACHE',0.1),('BCHE',0.3),('BCHE',0.5)])]
        rows.append(dict(rows[0],id='failed',status='preparation_failed',rmsd=None))
        table=RedockingResultsTable(result,rows)
        self.assertEqual(table.table.columns[0].label,'Target')
        self.assertEqual(sum(r.color=='#DCFCE7' for r in table.table.rows),3)
        self.assertEqual(sum(r.color=='#FEF2F2' for r in table.table.rows),1)
        table.target_filter.value='target:BCHE';table.filter_changed()
        self.assertEqual(len(table.table.rows),2)
        self.assertEqual(table.table.rows[0].cells[6].content.value,'0.300')
        table.highlight.value=False;table.filter_changed()
        self.assertTrue(all(r.color is None for r in table.table.rows))
        table.search.value='0004';table.filter_changed();self.assertEqual(len(table.table.rows),1)
        table.target_filter.value='all';table.search.value='';table.limit=1;table.filter_changed()
        table.move(1);self.assertEqual(len(table.table.rows),1)
        self.assertEqual(table.table.rows[0].color,'#DCFCE7' if table.highlight.value else None)


class AdmetNameIntegrationTests(unittest.TestCase):
    def test_every_published_admet_artifact_has_a_dataset_tag_and_preserves_the_reader_contract(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);input=root/'input';input.mkdir()
            pd.DataFrame([dict(molecule_chembl_id='A',canonical_smiles='CCO'),
                dict(molecule_chembl_id='BAD',canonical_smiles='invalid')]).to_csv(input/'molecules.csv',index=False)
            result=execute_operation('admet',dict(base_input_path=str(input),input_file='molecules.csv',
                admet_profile='exploratory',dataset_origin=dict(tag='MOL',label='ACHE')),root/'out')
            names=[Path(p).name for p in result.artifacts]
            self.assertTrue(names);self.assertTrue(all('_MOL' in name for name in names),names)
            self.assertTrue({'ACHE_molecules_MOL_evaluated.csv','ACHE_molecules_MOL_egg.png',
                'ACHE_molecules_MOL_BBB+.csv','ACHE_molecules_MOL_admet_report.json',
                'ACHE_molecules_MOL_molecule_exclusions.json'}<=set(names))
            self.assertEqual(result.details['excluded_records'],1)
            self.assertEqual(json.loads((root/'out/ACHE_molecules_MOL_admet_report.json').read_text())['dataset_origin']['tag'],'MOL')


class AdmetPipelineIdentityTests(unittest.TestCase):
    def test_individual_admet_datasets_reach_fingerprints_with_distinct_tagged_names(self):
        import sys
        import test_pipeline_execution as execution
        from biomolexplorer.pipeline import PipelineService
        from biomolexplorer.ui.file_selection import FileSelection
        fixture=execution.PipelineExecutionTests();fixture.setUp();self.addCleanup(fixture.tearDown)
        mol=fixture.upload('molecules.csv','name,smiles\nMOL_A,CCO\n')
        sim=fixture.upload('similars.csv','name,smiles\nSIM_A,CCC\n')
        imported=fixture.imported(mol);imported['name']='ACHE';imported['parameters']['target']='ACHE'
        imported['parameters']['asset_ids']=[mol,sim]
        admet=new_stage('admet');admet['bindings']['base_input_path']={'stage':imported['id']}
        fingerprint=new_stage('fingerprints');fingerprint['bindings']['base_input_path']={'stage':admet['id']}
        fixture.save([imported,admet,fingerprint])
        service=PipelineService(fixture.store,worker_python=sys.executable,cpu_workers=1)
        self.addCleanup(service.close)
        run=fixture.wait(service.submit(fixture.token,fixture.project_id))
        for operation in ('admet','fingerprints'):
            self.assertEqual(run['status'],'awaiting_input',run.get('error'))
            pending=next(s for s in run['stages'] if s['status']=='awaiting_input')
            self.assertEqual(pending['operation'],operation)
            form=FileSelection(run,pending);form.processing.value='individual'
            for rows in form.rows.values():
                for check,reference in rows:
                    name=Path(reference.get('selector','')).name
                    check.value=operation=='admet' or name.endswith('_evaluated.csv')
            run=fixture.wait(service.resume(fixture.token,run['id'],form.read()))
        self.assertEqual(run['status'],'succeeded',run.get('error'))
        admet_stage=next(s for s in run['stages'] if s['id']==admet['id'])
        self.assertEqual(len(admet_stage['batches']),2)
        names=[Path(p).name for p in admet_stage['artifacts']]
        self.assertEqual(len(names),len(set(names)),names)
        self.assertTrue(any('_MOL_evaluated.csv' in name for name in names))
        self.assertTrue(any('_SIM_evaluated.csv' in name for name in names))
        fingerprint_stage=next(s for s in run['stages'] if s['id']==fingerprint['id'])
        self.assertEqual(len(fingerprint_stage['batches']),2)
        tables=[pd.read_csv(p) for p in fingerprint_stage['artifacts'] if p.endswith('.csv')]
        self.assertEqual(set().union(*(set(t['molecule_chembl_id']) for t in tables)),{'MOL_A','SIM_A'})
