"""Clear All preserves source data/edges and atomically resets computed state."""
import copy
import json
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch
from uuid import uuid4

from biomolexplorer.catalog import new_stage
from biomolexplorer.result_files import ResultFiles
from biomolexplorer.workspace import WorkspaceStore,AccessDenied
from biomolexplorer.clear_results import reset_selections


class ClearAllResultsTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.store=WorkspaceStore(self.temp.name);self.token=self.store.register('Owner','clear@example.org','test-password')
        self.project=self.store.create_project(self.token,'Clear test');self.pid=self.project['id'];self.root=self.store.project_dir(self.pid)
        ticket=self.store.prepare_upload(self.token,self.pid,'original.csv','compounds');(self.store.staging/ticket).write_text('molecule_chembl_id,canonical_smiles\nA,CCO\n')
        self.asset=self.store.finish_upload(self.token,ticket);self.original=self.store.asset_path(self.token,self.pid,self.asset)
        imported=new_stage('import_results');imported['parameters'].update(kind='compounds',asset_ids=[self.asset])
        structures=new_stage('retrieve_structures');structures['parameters']['pdb_query']='ACHE'
        redocking=new_stage('redocking');redocking['parameters'].update(pdb_codes=[['1ABC','LIG',1,'A']],preparation_options={'ligand':{'charge_type':'am1'}})
        redocking['bindings']['base_input_path']={'stage':structures['id'],'selector':'1ABC.pdb','artifact_id':'a'*64}
        fingerprints=new_stage('fingerprints');fingerprints['bindings']['base_input_path']={'sources':[
            {'stage':imported['id'],'selector':'one.csv','artifact_id':'b'*64},
            {'stage':imported['id'],'selector':'two.csv','artifact_id':'c'*64}]}
        fingerprints['parameters']['files']=['one.csv','two.csv'];self.pipeline=[imported,structures,redocking,fingerprints]
        self.store.save_pipeline(self.token,self.pid,self.pipeline,0)
        self.rid=uuid4().hex;self.generated=self.root/'runs'/self.rid/redocking['id']/'artifacts'/'pose.pdbqt'
        self.generated.parent.mkdir(parents=True);self.generated.write_text('computed pose')
        self.partial=self.generated.parent/'unpublished.mol2';self.partial.write_text('partial result')
        self.log=self.generated.parent.parent/'execution.log';self.log.write_text('native log')
        self.cached=self.root/'.input-cache/old/data.csv';self.cached.parent.mkdir(parents=True);self.cached.write_text('old inputs')
        self.ccd=self.root/'.biomolexplorer/cache/calculations/ccd/ABC.cif';self.ccd.parent.mkdir(parents=True);self.ccd.write_text('source chemistry')
        self.stages=[dict(s,status='succeeded',configuration=copy.deepcopy(s),artifacts=[str(self.original)] if s['operation']=='import_results' else [str(self.generated)] if s['operation']=='redocking' else [],
            reused=True,redocking_failures=2,cache_key='cached',artifact_records=[],artifact_manifest={}) for s in self.pipeline]
        with self.store.connect() as db:
            uid=self.store.user(self.token)['id']
            db.execute('INSERT INTO runs VALUES(?,?,?,?,?,?,?,?)',(self.rid,self.pid,uid,'succeeded',json.dumps(self.stages),time.time(),time.time(),None))

    def test_results_selections_and_cache_reset_without_deleting_sources_or_edges(self):
        before=self.store.project(self.token,self.pid);run_before=self.store.get_run(self.token,self.rid)
        result=ResultFiles(self.store).clear_project_results(self.token,self.pid)
        after=self.store.project(self.token,self.pid)
        self.assertGreater(result['files_removed'],0)
        self.assertFalse(self.generated.exists());self.assertFalse(self.partial.exists());self.assertFalse(self.cached.exists())
        self.assertTrue(self.original.is_file());self.assertTrue(self.log.is_file());self.assertTrue(self.ccd.is_file())
        self.assertEqual(after['pipeline'][2]['bindings']['base_input_path'],{'stage':self.pipeline[1]['id'],'selector':'auto'})
        self.assertNotIn('pdb_codes',after['pipeline'][2]['parameters'])
        self.assertEqual(after['pipeline'][2]['parameters'],new_stage('redocking')['parameters'])
        self.assertEqual(after['pipeline'][3]['bindings']['base_input_path'],{'stage':self.pipeline[0]['id'],'selector':'auto'})
        self.assertNotIn('files',after['pipeline'][3]['parameters']);self.assertGreater(after['revision'],before['revision'])
        run=self.store.get_run(self.token,self.rid)
        self.assertEqual(run['updated'],run_before['updated'])
        self.assertEqual(run['stages'][0]['artifacts'],[])
        self.assertTrue(run['stages'][0]['outputs_removed']);self.assertTrue(self.original.exists())
        for stage in run['stages'][1:]:
            self.assertEqual(stage['artifacts'],[]);self.assertTrue(stage['outputs_removed']);self.assertFalse(stage['reused']);self.assertNotIn('cache_key',stage)
        self.assertEqual(self.store.history(self.token,self.pid)[0]['action'],'clear_all_results')

    def test_failure_restores_files_configuration_and_execution_state(self):
        before=self.store.project(self.token,self.pid)
        with patch('biomolexplorer.project_state.record',side_effect=RuntimeError('audit failure')):
            with self.assertRaisesRegex(RuntimeError,'audit failure'):ResultFiles(self.store).clear_project_results(self.token,self.pid)
        self.assertTrue(self.generated.is_file());self.assertTrue(self.partial.is_file());self.assertTrue(self.cached.is_file())
        self.assertEqual(self.store.project(self.token,self.pid)['pipeline'],before['pipeline'])
        self.assertEqual(self.store.get_run(self.token,self.rid)['stages'],self.stages)

    def test_active_run_blocks_cleanup(self):
        with self.store.connect() as db:db.execute("UPDATE runs SET status='awaiting_input' WHERE id=?",(self.rid,))
        with self.assertRaisesRegex(ValueError,'cancele'):ResultFiles(self.store).clear_project_results(self.token,self.pid)
        self.assertTrue(self.generated.exists());self.assertTrue(self.cached.exists())

    def test_viewer_cannot_clear_and_provided_results_are_preserved(self):
        viewer=self.store.register('Viewer','viewer-clear@example.org','test-password')
        self.store.invite(self.token,self.pid,'viewer-clear@example.org','viewer');self.store.accept_invitation(viewer,self.pid)
        with self.assertRaises(AccessDenied):ResultFiles(self.store).clear_project_results(viewer,self.pid)
        self.stages[2]['provided']=True
        with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(self.stages),self.rid))
        ResultFiles(self.store).clear_project_results(self.token,self.pid)
        self.assertTrue(self.generated.exists());self.assertTrue(self.original.exists())

    def test_reset_is_idempotent_and_preserves_parameters_and_direct_sources(self):
        once,_=reset_selections(self.pipeline);twice,count=reset_selections(once)
        self.assertEqual(once,twice);self.assertEqual(count,0)
        self.assertEqual(once[0]['parameters']['asset_ids'],[self.asset])
        self.assertEqual(once[1]['parameters']['pdb_query'],'ACHE')

    def test_standalone_pair_selections_are_reset_without_losing_shared_preparation(self):
        stage=new_stage('redocking');stage['bindings']={}
        stage['parameters'].update(pdb_codes=[['1ABC','LIG',1,'A']],preparation_pairs={'old':{}},
            preparation_options={'ligand':{'charge_type':'am1'}},base_input_path='originals')
        pipeline,_=reset_selections([stage]);params=pipeline[0]['parameters']
        self.assertNotIn('pdb_codes',params);self.assertNotIn('preparation_pairs',params)
        self.assertEqual(params['preparation_options'],stage['parameters']['preparation_options'])
        self.assertEqual(params['base_input_path'],'originals')


if __name__=='__main__':unittest.main()
