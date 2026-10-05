"""Project folders, collaboration, rollback and migration preserve experiments."""
import copy
import json
import sys
import tempfile
import time
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

from biomolexplorer.catalog import new_stage
from biomolexplorer.flow import connect
from biomolexplorer.pipeline import PipelineService,validate_pipeline
from biomolexplorer.workspace import WorkspaceStore,AccessDenied
from biomolexplorer.input_validation import merge_csv,validate_bundle


class ProjectPortabilityTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='biomol portable ')
        self.root=Path(self.temp.name)
        self.store=WorkspaceStore(self.root/'workspace')
        self.owner=self.store.register('Owner','owner@example.org','owner-password-2026')
        self.editor=self.store.register('Editor','editor@example.org','editor-password-2026')
        self.viewer=self.store.register('Viewer','viewer@example.org','viewer-password-2026')
        self.folder=self.root/'experiments'/'Study with spaces'
        self.project=self.store.create_project(self.owner,'Portable study',directory=self.folder)
        self.pid=self.project['id']

    def tearDown(self):self.temp.cleanup()

    def upload(self,name='compounds.csv',content='name,smiles\nMOL1,CCO\nMOL2,CCN\n',kind='compounds'):
        ticket=self.store.prepare_upload(self.owner,self.pid,name,kind)
        (self.store.staging/ticket).write_text(content)
        return self.store.finish_upload(self.owner,ticket)

    def save(self,stages):
        p=self.store.project(self.owner,self.pid)
        return self.store.save_pipeline(self.owner,self.pid,stages,p['revision'])

    def wait(self,service,run):
        for _ in range(1600):
            current=self.store.get_run(self.owner,run['id'])
            if current['status'] not in ('queued','running'):return current
            time.sleep(.025)
        self.fail('execution timeout')

    def share(self,token,role):
        self.store.invite(self.owner,self.pid,self.store.user(token)['email'],role)
        self.store.accept_invitation(token,self.pid,expected_role=role)

    def test_selected_folder_contains_pipeline_assets_results_and_versions(self):
        asset=self.upload()
        stage=new_stage('retrieve_compounds')
        stage['provided_results']={'kind':'compounds','asset_ids':[asset]}
        self.save([stage])
        service=PipelineService(self.store)
        try:run=self.wait(service,service.submit(self.owner,self.pid))
        finally:service.close()
        self.assertEqual(run['status'],'succeeded',run['error'])
        self.assertTrue(all(Path(p).is_relative_to(self.folder) for p in run['stages'][0]['artifacts']))
        self.assertTrue(self.store.asset_path(self.owner,self.pid,asset).is_relative_to(self.folder))
        mirror=json.loads((self.folder/'project.json').read_text())
        self.assertEqual(mirror['pipeline'][0]['id'],stage['id'])
        self.assertTrue((self.folder/'.history'/'blobs').is_dir())
        self.assertFalse((self.store.projects_root/self.pid).exists())
        reopened=WorkspaceStore(self.store.root)
        self.assertEqual(reopened.project(self.owner,self.pid)['directory'],str(self.folder))
        self.assertEqual(reopened.asset_path(self.owner,self.pid,asset),self.store.asset_path(self.owner,self.pid,asset))
        self.assertEqual(reopened.get_run(self.owner,run['id'])['stages'][0]['artifacts'],run['stages'][0]['artifacts'])

    def test_explicit_empty_folder_does_not_fall_back_to_workspace_storage(self):
        for directory in ('','   '):
            with self.subTest(directory=directory),self.assertRaisesRegex(ValueError,'Informe a pasta'):
                self.store.create_project(self.owner,'Missing folder',directory=directory)
        self.assertEqual(len(self.store.list_projects(self.owner)),1)

    def test_project_folder_cannot_overlap_another_project_or_overwrite_files(self):
        for directory in (self.folder,self.folder/'nested',self.root/'experiments',self.store.root,self.store.staging/'oops'):
            with self.subTest(directory=directory),self.assertRaises(ValueError):
                self.store.create_project(self.owner,'Other',directory=directory)

    def test_rollback_restores_before_change_including_files_and_configuration(self):
        original=self.upload()
        stage=new_stage('admet')
        stage['bindings']={'base_input_path':{'asset':original}}
        self.save([stage])
        second=self.upload('extra.csv','name,smiles\nMOL3,CCC\n')
        event=next(e for e in self.store.history(self.owner,self.pid) if e['summary']=='Arquivo fornecido: extra.csv')
        stage['parameters']['target']='Changed'
        self.save([stage])
        self.store.update_project(self.owner,self.pid,'Renamed','changed','#6366F1',['tag'])
        self.store.rollback(self.owner,self.pid,event['id'])
        project=self.store.project(self.owner,self.pid)
        self.assertEqual(project['name'],'Portable study')
        self.assertNotEqual(project['pipeline'][0]['parameters'].get('target'),'Changed')
        self.assertEqual([a['id'] for a in self.store.assets(self.owner,self.pid)],[original])
        self.assertFalse((self.folder/'assets'/second).exists())
        self.assertEqual(self.store.history(self.owner,self.pid)[0]['action'],'rollback')

    def test_history_authorization_and_membership_rollback(self):
        self.share(self.editor,'editor')
        self.share(self.viewer,'viewer')
        self.store.revoke(self.owner,self.pid,self.store.user(self.editor)['id'])
        event=self.store.history(self.owner,self.pid)[0]
        with self.assertRaises(AccessDenied):self.store.history(self.editor,self.pid)
        with self.assertRaises(AccessDenied):self.store.rollback(self.viewer,self.pid,event['id'])
        self.store.rollback(self.owner,self.pid,event['id'])
        self.assertEqual(self.store.project(self.editor,self.pid)['role'],'editor')
        self.assertEqual(event['actor_name'],'Owner')
        events=self.store.history(self.viewer,self.pid)
        self.assertEqual([e['created'] for e in events],sorted((e['created'] for e in events),reverse=True))

    def test_collaborative_merge_preserves_independent_changes_and_rejects_same_field(self):
        self.share(self.editor,'editor')
        first,second=new_stage('retrieve_compounds'),new_stage('admet')
        initial=self.save([first,second])
        local,remote=copy.deepcopy(initial['pipeline']),copy.deepcopy(initial['pipeline'])
        local[0]['name']='Owner edit';remote[1]['name']='Editor edit'
        self.store.save_pipeline(self.owner,self.pid,local,initial['revision'])
        saved=self.store.save_pipeline(self.editor,self.pid,remote,initial['revision'],base_pipeline=initial['pipeline'])
        self.assertEqual([s['name'] for s in saved['pipeline']],['Owner edit','Editor edit'])
        remote[0]['name']='Conflicting edit'
        with self.assertRaisesRegex(ValueError,'mesmo campo'):
            self.store.save_pipeline(self.editor,self.pid,remote,initial['revision'],base_pipeline=initial['pipeline'])

    def test_multiple_inputs_normalize_aliases_deduplicate_and_are_recorded(self):
        first=self.upload()
        second=self.upload('additional.csv','molecule_chembl_id,canonical_smiles\nMOL1,OCC\nMOL3,CCC\n')
        stage=new_stage('admet')
        stage['bindings']={'base_input_path':{'sources':[{'asset':first},{'asset':second}]}}
        validate_pipeline([stage]);self.save([stage])
        service=PipelineService(self.store,worker_python=sys.executable)
        try:run=self.wait(service,service.submit(self.owner,self.pid))
        finally:service.close()
        self.assertEqual(run['status'],'succeeded',run['error'])
        item=run['stages'][0]
        self.assertEqual(len(item['input_files']['base_input_path']),2)
        import pandas as pd
        file=next(p for p in item['artifacts'] if p.endswith('selected_compounds.csv'))
        self.assertEqual(set(pd.read_csv(file).molecule_chembl_id),{'MOL1','MOL2','MOL3'})

    def test_mouse_connections_can_have_multiple_origins(self):
        first,second,target=new_stage('retrieve_compounds'),new_stage('retrieve_compounds'),new_stage('admet')
        stages=[first,second,target]
        connect(stages,first['id'],target['id'],'base_input_path')
        connect(stages,second['id'],target['id'],'base_input_path')
        self.assertEqual(len(target['bindings']['base_input_path']['sources']),2)
        self.assertEqual(validate_pipeline(stages)[-1],target['id'])

    def test_supplied_admet_results_skip_worker_and_generate_interactive_egg(self):
        asset=self.upload('admet.csv','molecule_chembl_id,canonical_smiles,TPSA,WLOGP\nMOL1,CCO,20.2,-0.01\n')
        stage=new_stage('admet')
        stage['provided_results']={'kind':'compounds','asset_ids':[asset]}
        self.save([stage])
        with patch('biomolexplorer.pipeline.JobManager',side_effect=AssertionError('scientific work must be skipped')):
            service=PipelineService(self.store)
            try:run=self.wait(service,service.submit(self.owner,self.pid))
            finally:service.close()
        self.assertEqual(run['status'],'succeeded',run['error'])
        self.assertTrue(run['stages'][0]['provided'])
        view=json.loads(Path(next(p for p in run['stages'][0]['artifacts'] if p.endswith('.biomol-view.json'))).read_text())
        self.assertEqual(view['nodes'][0]['id'],'MOL1')
        self.assertEqual(view['nodes'][0]['properties']['canonical_smiles'],'CCO')

    def test_supplied_result_contract_explains_missing_admet_columns(self):
        path=self.store.asset_path(self.owner,self.pid,self.upload())
        with self.assertRaisesRegex(ValueError,'TPSA.*WLOGP'):
            validate_bundle([path],'compounds','admet')

    def test_export_import_rebases_and_reuses_completed_experiments_on_new_computer(self):
        asset=self.upload()
        source=new_stage('retrieve_compounds');source['provided_results']={'kind':'compounds','asset_ids':[asset]}
        analysis=new_stage('admet');analysis['bindings']={'base_input_path':{'stage':source['id'],'selector':'compounds.csv'}}
        self.save([source,analysis])
        service=PipelineService(self.store,worker_python=sys.executable)
        try:
            original=self.wait(service,service.submit(self.owner,self.pid))
            self.assertEqual(original['status'],'awaiting_input')
            original=self.wait(service,service.resume(self.owner,original['id'],analysis))
        finally:service.close()
        self.assertEqual(original['status'],'succeeded',original['error'])
        archive=self.store.export_project(self.owner,self.pid)
        other=WorkspaceStore(self.root/'computer-two')
        token=other.register('New owner','new@example.org','new-owner-password')
        imported=other.import_project(token,archive,self.root/'migrated-study')
        self.assertEqual(len(other.list_runs(token,imported['id'])),1)
        self.assertEqual(other.project(token,imported['id'])['owner'],other.user(token)['id'])
        with patch('biomolexplorer.pipeline.JobManager',side_effect=AssertionError('experiment unexpectedly repeated')):
            service=PipelineService(other,worker_python=None)
            try:
                run=service.submit(token,imported['id'])
                for _ in range(400):
                    run=other.get_run(token,run['id'])
                    if run['status'] not in ('queued','running'):break
                    time.sleep(.025)
                self.assertEqual(run['status'],'succeeded',run['error'])
            finally:service.close()
        self.assertEqual(run['status'],'succeeded',run['error'])
        self.assertTrue(all(s.get('reused') for s in run['stages']))
        self.assertTrue(all(Path(p).is_relative_to(self.root/'migrated-study') for s in run['stages'] for p in s['artifacts']))
        with zipfile.ZipFile(archive) as package:
            self.assertNotIn('workspace.sqlite3',package.namelist())

    def test_import_on_same_workspace_handles_identifier_collisions(self):
        self.upload()
        archive=self.store.export_project(self.owner,self.pid)
        imported=self.store.import_project(self.owner,archive,self.root/'copy-study')
        self.assertNotEqual(self.store.assets(self.owner,self.pid)[0]['id'],self.store.assets(self.owner,imported['id'])[0]['id'])
        self.assertEqual(self.store.read_file(self.owner,self.pid,self.store.asset_path(self.owner,self.pid,self.store.assets(self.owner,self.pid)[0]['id'])),
                         self.store.read_file(self.owner,imported['id'],self.store.asset_path(self.owner,imported['id'],self.store.assets(self.owner,imported['id'])[0]['id'])))

    def test_archive_rejects_path_traversal_without_creating_project(self):
        archive=self.root/'malicious.zip'
        with zipfile.ZipFile(archive,'w') as z:z.writestr('../escape.txt','bad')
        with self.assertRaisesRegex(ValueError,'caminho inválido'):
            self.store.import_project(self.owner,archive,self.root/'invalid-import')
        self.assertFalse((self.root/'invalid-import').exists())

    def test_failed_csv_union_is_atomic_and_cannot_be_reused_as_partial_input(self):
        first=self.store.asset_path(self.owner,self.pid,self.upload())
        second=self.store.asset_path(self.owner,self.pid,self.upload('conflict.csv','name,smiles\nMOL1,CCC\n'))
        destination=self.folder/'inputs'/'merged.csv'
        for _ in range(2):
            with self.assertRaisesRegex(ValueError,'conflitantes'):merge_csv([first,second],destination)
            self.assertFalse(destination.exists())

    def test_different_fingerprint_lengths_are_rejected_across_files(self):
        first=self.store.asset_path(self.owner,self.pid,self.upload('first.csv','molecule_chembl_id,fingerprint\nMOL1,"[1, 0]"\n','fingerprints'))
        second=self.store.asset_path(self.owner,self.pid,self.upload('second.csv','molecule_chembl_id,fingerprint\nMOL2,"[1, 0, 1]"\n','fingerprints'))
        with self.assertRaisesRegex(ValueError,'tamanhos diferentes'):
            merge_csv([first,second],self.folder/'merged.csv','fingerprints')

    def test_blank_compound_identifiers_are_generated_consistently(self):
        asset=self.upload('blank.csv','name,smiles\n,CCO\n')
        source=self.store.asset_path(self.owner,self.pid,asset)
        destination=self.folder/'merged.csv';merge_csv([source],destination)
        import csv
        with destination.open() as stream:identifier=next(csv.DictReader(stream))['molecule_chembl_id']
        copied=self.folder/'copy.csv';copied.write_bytes(source.read_bytes())
        PipelineService._normalize_compounds(copied)
        with copied.open() as stream:self.assertEqual(next(csv.DictReader(stream))['molecule_chembl_id'],identifier)

    def test_corrupted_archive_is_rejected_without_partial_project(self):
        self.upload()
        exported=self.store.export_project(self.owner,self.pid)
        damaged=self.root/'damaged.zip'
        with zipfile.ZipFile(exported) as source,zipfile.ZipFile(damaged,'w') as target:
            filename=next(n for n in source.namelist() if n.startswith('data/assets/'))
            for name in source.namelist():target.writestr(name,b'corrupted' if name==filename else source.read(name))
        before=len(self.store.list_projects(self.owner))
        with self.assertRaisesRegex(ValueError,'alterado|incompleto'):
            self.store.import_project(self.owner,damaged,self.root/'damaged-project')
        self.assertEqual(len(self.store.list_projects(self.owner)),before)
        self.assertFalse((self.root/'damaged-project').exists())

    def test_migrated_completed_dock6_does_not_require_the_old_machine_installation(self):
        import test_pipeline_execution as execution
        manager,calls=execution.PipelineExecutionTests.immediate_manager(self)
        old_install=self.root/'old-dock6-installation';old_install.mkdir()
        stage=new_stage('docking_dock6')
        for key in ('base_input_path','base_selected_mols','base_vina_path'):
            directory=self.folder/key;directory.mkdir();stage['parameters'][key]=str(directory)
        stage['parameters'].update(dock6_app_path=str(old_install),pdb_code=['1ABC','LIG',1,'A'],charge_type='gas')
        self.save([stage])
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,dock6_path=old_install)
            try:run=self.wait(service,service.submit(self.owner,self.pid))
            finally:service.close()
        self.assertEqual(run['status'],'succeeded',run['error'])
        archive=self.store.export_project(self.owner,self.pid)
        other=WorkspaceStore(self.root/'other-computer')
        token=other.register('New owner','new@example.org','new-owner-password-2026')
        imported=other.import_project(token,archive,self.root/'imported-dock6')
        with patch('biomolexplorer.pipeline.JobManager',side_effect=AssertionError('completed docking should be reused')):
            service=PipelineService(other,dock6_path=None)
            try:
                copied=service.submit(token,imported['id'])
                for _ in range(400):
                    copied=other.get_run(token,copied['id'])
                    if copied['status'] not in ('queued','running'):break
                    time.sleep(.025)
            finally:service.close()
        self.assertEqual(copied['status'],'succeeded',copied['error'])
        self.assertTrue(copied['stages'][0]['reused'])


if __name__=='__main__':unittest.main()
