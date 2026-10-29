"""Confirmed state transitions and incremental scientific work units."""
import asyncio
import copy
import json
from pathlib import Path
import sys
import time
import unittest
from unittest.mock import patch
from uuid import uuid4

import test_pipeline_execution as execution
import test_pipeline_editor_regressions as editor
from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.result_state import signature


class IncrementalActionsTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    finish_selected=execution.PipelineExecutionTests.finish_selected
    upload=execution.PipelineExecutionTests.upload
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def test_add_algorithm_then_change_only_its_parameters_preserves_other_algorithm(self):
        asset=self.upload();stage=new_stage('fingerprints');stage['input_processing']='individual'
        stage['bindings']['base_input_path']={'asset':asset};self.save([stage])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:
                def run():
                    value=self.finish_selected(service,service.submit(self.token,self.project_id))
                    self.assertEqual(value['status'],'succeeded',value['error']);return value
                first=run();self.assertEqual(len(calls),1)
                old=first['stages'][0]['artifacts']
                stage['parameters']['maccs']=True
                self.assertEqual(self.store.result_removal_plan(self.token,self.project_id,[stage])['affected'],0)
                self.save([stage]);second=run();self.assertEqual(len(calls),2)
                self.assertTrue(second['stages'][0]['batches'][0]['reused'])
                self.assertTrue(all(Path(p).exists() for p in old))
                maccs=second['stages'][0]['batches'][1]['artifacts']
                stage['parameters']['radius']=3
                proposal=self.store.result_removal_plan(self.token,self.project_id,[stage])
                self.assertEqual([b['id'] for b in proposal['blocks']],[stage['id']])
                with self.assertRaisesRegex(ValueError,'Confirme'):self.save([stage])
                project=self.store.project(self.token,self.project_id)
                self.store.save_pipeline(self.token,self.project_id,[stage],project['revision'],prune_results=True,confirmation=proposal['confirmation'])
                self.assertTrue(all(Path(p).exists() for p in maccs))
                third=run();self.assertEqual(len(calls),3)
                self.assertTrue(third['stages'][0]['batches'][1]['reused'])
            finally:service.close()

    def test_presentation_and_inactive_engine_changes_have_equivalent_signatures(self):
        stage=new_stage('redocking');other=copy.deepcopy(stage)
        other.update(name='Renamed',position={'x':500,'y':300});other['parameters']['density']=.9
        self.assertEqual(signature(stage),signature(other))
        other['parameters']['exhaustiveness']=99
        self.assertNotEqual(signature(stage),signature(other))

    def test_joint_preparation_keeps_one_job_for_both_native_formats(self):
        from test_redocking_incremental import ATOM,LIG
        source=self.store.project_dir(self.project_id)/'input/Target';source.mkdir(parents=True)
        (source/'1ABC.pdb').write_text(ATOM+LIG)
        stage=new_stage('prepare_structures')
        stage['parameters'].update(base_input_path=str(source.parent),target='Target',pdb_codes=[['1ABC','LIG',1,'A']],docking_engines='both')
        self.save([stage]);manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:run=self.finish_selected(service,service.submit(self.token,self.project_id))
            finally:service.close()
        self.assertEqual(run['status'],'succeeded',run['error'])
        self.assertEqual(len(calls),1);self.assertEqual(calls[0][1]['docking_engines'],'both')

    def test_explicit_preparation_defaults_and_admet_profile_labels_are_equivalent(self):
        from biomolexplorer.redocking_config import DEFAULTS
        from biomolexplorer.admet_config import resolve_config
        stage=new_stage('redocking');other=copy.deepcopy(stage)
        other['parameters']['preparation_options']={'receptor':dict(DEFAULTS),'ligand':dict(DEFAULTS),'cofactors':[]}
        self.assertEqual(signature(stage),signature(other))
        stage=new_stage('admet');other=copy.deepcopy(stage)
        other['parameters'].update(resolve_config(stage['parameters']));other['parameters']['admet_profile']='custom'
        other['parameters']['alert_catalogs'].reverse()
        self.assertEqual(signature(stage),signature(other))

    def test_omitted_options_use_worker_defaults_instead_of_ui_presets(self):
        from biomolexplorer.result_state import request
        stage=new_stage('fingerprints')
        for action in ('maccs','pharmacophore'):stage['parameters'].pop(action)
        self.assertEqual(request(stage)['axes']['actions'],['morgan','maccs','pharmacophore'])
        old=new_stage('admet');old['parameters']={}
        self.assertNotEqual(signature(old),signature(new_stage('admet')))

    def test_adding_action_reuses_algorithms_from_an_older_combined_job(self):
        self.exercise_combined_job(remove_existing=False)

    def test_removing_action_from_an_older_combined_job_preserves_other_actions(self):
        self.exercise_combined_job(remove_existing=True)

    def exercise_combined_job(self,remove_existing):
        asset=self.upload();stage=new_stage('fingerprints');stage['input_processing']='merge'
        stage['bindings']['base_input_path']={'asset':asset};stage['parameters']['maccs']=True;self.save([stage])
        base,calls=self.immediate_manager()
        class NamedOutputs(base):
            def submit(manager,operation,parameters):
                job=super().submit(operation,parameters);outputs=[]
                for action in ('morgan','maccs','pharmacophore'):
                    if not parameters.get(action):continue
                    path=manager.config.workspace/(action+'_compounds.csv');path.write_text(action+' scientific result');outputs.append(str(path))
                job['result']['artifacts']=outputs
                return job
        with patch('biomolexplorer.pipeline.JobManager',NamedOutputs):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:
                def combined(service,pid,config,results):yield config,config['name']
                with patch.object(PipelineService,'_work_units',combined):
                    first=self.finish_selected(service,service.submit(self.token,self.project_id))
                self.assertEqual(first['status'],'succeeded',first['error']);self.assertEqual(len(calls),1)
                original=first['stages'][0]['artifacts']
                if remove_existing:
                    stage['parameters']['morgan']=False
                    project=self.store.project(self.token,self.project_id)
                    self.store.save_pipeline(self.token,self.project_id,[stage],project['revision'],prune_results=True)
                    removed=next(path for path in original if Path(path).name.startswith('morgan_'))
                    self.assertFalse(Path(removed).exists())
                    original=[path for path in original if path!=removed]
                    self.assertTrue(all(Path(path).exists() for path in original))
                stage['parameters']['pharmacophore']=True;self.save([stage])
                second=self.finish_selected(service,service.submit(self.token,self.project_id))
                self.assertEqual(second['status'],'succeeded',second['error']);self.assertEqual(len(calls),2)
                self.assertEqual([b['reused'] for b in second['stages'][0]['batches']],[True,False] if remove_existing else [True,True,False])
                self.assertTrue(set(original)<=set(second['stages'][0]['artifacts']))
            finally:service.close()

    def test_removing_one_pair_preserves_shared_receptor_and_other_pose(self):
        from biomolexplorer.stage_cache import artifact_manifest
        stage=new_stage('redocking');stage['parameters']['pdb_codes']=[['1ABC','LIG',1,'A'],['1ABC','OTH',2,'A']]
        self.save([stage]);directory=self.store.project_dir(self.project_id)/'artifacts';directory.mkdir()
        receptor=directory/'receptor.pdbqt';one=directory/'one.pdbqt';two=directory/'two.pdbqt';table=directory/'pdb_codes.csv';index=directory/'redocking_cache.json'
        for path in (receptor,one,two,table):path.write_text('valid scientific result')
        entries=[dict(record=record,files={path.name:'digest' for path in (receptor,pose)}) for record,pose in zip(stage['parameters']['pdb_codes'],(one,two))]
        index.write_text(json.dumps(dict(version=1,entries=entries)))
        artifacts=list(map(str,(receptor,one,two,table,index)))
        item=dict(stage,configuration=copy.deepcopy(stage),artifacts=artifacts,artifact_manifest=artifact_manifest(artifacts),status='succeeded')
        rid=uuid4().hex
        with self.store.connect() as db:db.execute('INSERT INTO runs VALUES(?,?,?,?,?,?,?,?)',(rid,self.project_id,self.store.user(self.token)['id'],'succeeded',json.dumps([item]),time.time(),time.time(),None))
        stage['parameters']['pdb_codes']=stage['parameters']['pdb_codes'][1:]
        proposal=self.store.result_removal_plan(self.token,self.project_id,[stage]);self.assertGreater(proposal['affected'],0)
        project=self.store.project(self.token,self.project_id)
        self.store.save_pipeline(self.token,self.project_id,[stage],project['revision'],prune_results=True,confirmation=proposal['confirmation'])
        self.assertFalse(one.exists());self.assertFalse(table.exists());self.assertTrue(two.exists());self.assertTrue(receptor.exists())
        self.assertTrue(index.exists());self.assertTrue(self.store.get_run(self.token,rid)['stages'][0]['partial_results'])
        self.assertEqual(self.store.result_removal_plan(self.token,self.project_id,[stage])['affected'],0)


class AutomaticConfirmationTests(unittest.TestCase):
    setUp=editor.PipelineEditorRegressionTests.setUp

    def seed(self):
        self.stage['parameters']['pdb_query']='ACHE'
        dependent=new_stage('redocking');dependent['bindings']['base_input_path']={'stage':self.stage['id']}
        independent=new_stage('retrieve_structures');independent['parameters']['pdb_query']='BCHE'
        pipeline=[self.stage,dependent,independent]
        self.store.save_pipeline(self.token,self.pid,pipeline,0)
        self.paths=[];items=[]
        for stage in pipeline:
            path=self.store.project_dir(self.pid)/('result-'+stage['id']+'.csv');path.write_text('valid')
            self.paths.append(path)
            items.append(dict(stage,configuration=copy.deepcopy(stage),status='succeeded',artifacts=[str(path)]))
        # Legacy results without file lineage still honor explicit DAG edges.
        self.rid=uuid4().hex
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES(?,?,?,?,?,?,?,?)',(self.rid,self.pid,self.store.user(self.token)['id'],'succeeded',json.dumps(items),time.time(),time.time(),None))
        self.ui.current=self.store.project(self.token,self.pid)
        self.ui.base_pipeline=copy.deepcopy(self.ui.current['pipeline'])
        self.ui.current['pipeline'][0]['parameters']['pdb_query']='BACE1';self.ui.dirty=True
        return pipeline

    def test_cancel_restores_configuration_in_memory_and_database(self):
        pipeline=self.seed()
        async def scenario():
            pending=asyncio.create_task(self.ui.save_draft(silent=True))
            await asyncio.sleep(0)
            dialog=self.page.dialogs[-1]
            dialog.actions[0].on_click(None)
            self.assertIs(await pending,False)
        asyncio.run(scenario())
        self.assertEqual(self.ui.current['pipeline'],pipeline)
        self.assertEqual(self.store.project(self.token,self.pid)['pipeline'],pipeline)
        self.assertTrue(all(p.exists() for p in self.paths))

    def test_confirm_lists_and_deletes_cascade_but_preserves_independent_branch(self):
        from test_frontend_navigation import walk
        import flet as ft
        pipeline=self.seed()
        async def scenario():
            pending=asyncio.create_task(self.ui.save_draft(silent=True));await asyncio.sleep(0)
            dialog=self.page.dialogs[-1]
            labels=[c.value for c in walk(dialog.content) if isinstance(c,ft.Text)]
            self.assertIn('• '+pipeline[0]['name'],labels)
            self.assertIn('• '+pipeline[1]['name'],labels)
            dialog.actions[1].on_click(None);await pending
        asyncio.run(scenario())
        self.assertFalse(self.paths[0].exists());self.assertFalse(self.paths[1].exists());self.assertTrue(self.paths[2].exists())
        self.assertEqual(self.store.project(self.token,self.pid)['pipeline'][0]['parameters']['pdb_query'],'BACE1')

    def test_confirmation_is_bound_to_exact_result_state(self):
        self.seed();project=self.store.project(self.token,self.pid);candidate=self.ui.current['pipeline']
        proposal=self.store.result_removal_plan(self.token,self.pid,candidate)
        candidate=copy.deepcopy(candidate);candidate[0]['parameters']['pdb_query']='different'
        with self.assertRaisesRegex(ValueError,'Revise'):
            self.store.save_pipeline(self.token,self.pid,candidate,project['revision'],prune_results=True,confirmation=proposal['confirmation'])
        self.assertTrue(all(p.exists() for p in self.paths))
