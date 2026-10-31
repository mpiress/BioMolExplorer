"""Completed experiments survive downstream additions and code-only cache drift."""
import copy
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.workspace import WorkspaceStore
from biomolexplorer.published_results import catalog,verify
from biomolexplorer.stage_cache import file_digest
import test_pipeline_execution as execution


class PublishedResultTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    upload=execution.PipelineExecutionTests.upload
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    finish_selected=execution.PipelineExecutionTests.finish_selected
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def execute(self,service):
        run=self.finish_selected(service,service.submit(self.token,self.project_id))
        self.assertEqual(run['status'],'succeeded',run['error']);return run

    def test_completed_experiments_keep_original_version_after_restart_and_new_md(self):
        source=new_stage('retrieve_compounds');admet=new_stage('admet')
        admet['bindings']['base_input_path']={'stage':source['id']};self.save([source,admet])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:first=self.execute(service)
            finally:service.close()
            records=copy.deepcopy(first['stages'])
            for item in records:item['processing_version']='original-producer-v1';item.pop('adopted_processing_version',None)
            with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(records),first['id']))
            files={p:(file_digest(p),Path(p).stat().st_mtime_ns) for item in records for p in item['artifacts']}
            root=self.store.project_dir(self.project_id);rec=root/'4M0E_A.noH.pdb';lig=root/'CHEMBL1179697.mol2'
            rec.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n')
            lig.write_text('@<TRIPOS>MOLECULE\nCHEMBL1179697\n1 0 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 LIG 0\n@<TRIPOS>BOND\n')
            md=new_stage('dynamics_dock6');md['depends_on']=[admet['id']]
            md['parameters'].update(base_input_path=str(rec),base_selected_mols=str(lig))
            self.save([source,admet,md]);service=PipelineService(WorkspaceStore(self.store.root),worker_python=sys.executable)
            try:second=self.wait(service.submit(self.token,self.project_id))
            finally:service.close()
        self.assertEqual(second['status'],'awaiting_input',second['error'])
        self.assertEqual(len(calls),2) # Retrieval and ADMET from the original run only.
        for item in second['stages'][:2]:
            self.assertTrue(item['reused']);self.assertEqual(item['reuse_kind'],'published_result')
            self.assertEqual(item['processing_version'],'original-producer-v1')
            self.assertEqual(item['reuse_provenance']['producer_version'],'original-producer-v1')
        self.assertEqual(files,{p:(file_digest(p),Path(p).stat().st_mtime_ns) for p in files})
        available=self.store.available_results(self.token,self.project_id)
        self.assertEqual({item['id'] for item in available},{source['id'],admet['id']})

    def test_changed_configuration_or_tampered_output_is_not_a_published_result(self):
        stage=new_stage('retrieve_compounds');self.save([stage]);manager,_=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:run=self.execute(service)
            finally:service.close()
        self.assertEqual(len(catalog(self.store,self.project_id,[stage])),1)
        changed=copy.deepcopy(stage);changed['parameters']['search_term']='DIFFERENT'
        self.assertEqual(catalog(self.store,self.project_id,[changed]),{})
        output=Path(run['stages'][0]['artifacts'][0]);output.write_text('tampered')
        self.assertEqual(catalog(self.store,self.project_id,[stage]),{})

    def test_current_input_changes_invalidate_only_affected_units(self):
        first=self.upload('a.csv');second=self.upload('b.csv','name,smiles\nB,CCC\n')
        stage=new_stage('fingerprints');stage['input_processing']='individual'
        stage['bindings']['base_input_path']={'sources':[{'asset':first},{'asset':second}]};self.save([stage])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:run=self.execute(service)
            finally:service.close()
            records=copy.deepcopy(run['stages'])
            records[0]['processing_version']='old-science'
            for batch in records[0]['batches']:batch['processing_version']='old-science';batch['cache_key']='old-'+batch['cache_key']
            with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(records),run['id']))
            third=self.upload('c.csv','name,smiles\nC,CCCC\n')
            stage['bindings']['base_input_path']['sources'].append({'asset':third});self.save([stage])
            before=len(calls);service=PipelineService(self.store,worker_python=sys.executable)
            try:result=self.execute(service)
            finally:service.close()
            self.assertEqual(len(calls),before+1)
            self.assertEqual(sum(bool(batch['reused']) for batch in result['stages'][0]['batches']),2)

    def test_invalid_new_downstream_block_does_not_hide_completed_upstream(self):
        source=new_stage('retrieve_compounds');self.save([source]);manager,_=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:self.execute(service)
            finally:service.close()
        md=new_stage('dynamics_dock6');md['bindings']['base_input_path']={'stage':source['id']}
        self.assertEqual(set(catalog(self.store,self.project_id,[source,md])),{source['id']})

    def test_explicit_no_reuse_still_recalculates_without_deleting_previous_files(self):
        source=new_stage('retrieve_compounds');self.save([source]);manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:
                first=self.execute(service)
                second=self.finish_selected(service,service.submit(self.token,self.project_id,reuse_results=False))
            finally:service.close()
        self.assertEqual(second['status'],'succeeded',second['error']);self.assertEqual(len(calls),2)
        self.assertTrue(all(Path(p).is_file() for p in first['stages'][0]['artifacts']))

    def test_cancelled_latest_run_keeps_older_results_accessible_in_runs(self):
        import asyncio,time
        from uuid import uuid4
        from types import SimpleNamespace
        from unittest.mock import AsyncMock
        import flet as ft
        from biomolexplorer.ui.app import WorkspaceUI
        from test_frontend_navigation import Page,walk
        source=new_stage('retrieve_compounds');self.save([source]);manager,_=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:first=self.execute(service)
            finally:service.close()
        cancelled=uuid4().hex
        item=dict(id=source['id'],operation=source['operation'],name=source['name'],configuration=source,status='cancelled',artifacts=[])
        with self.store.connect() as db:
            user=self.store.user(self.token)['id']
            for index in range(101):
                rid=cancelled if index==100 else uuid4().hex
                db.execute('INSERT INTO runs VALUES(?,?,?,?,?,?,?,?)',(rid,self.project_id,user,
                    'cancelled',json.dumps([item]),time.time()+index+1,time.time()+index+1,'Execução cancelada.'))
        self.assertNotIn(first['id'],{run['id'] for run in self.store.list_runs(self.token,self.project_id)})
        async def scenario():
            ui=WorkspaceUI(Page(),self.store,SimpleNamespace(dock6_path=None));ui.token=self.token
            self.addCleanup(ui.structure_viewer.close)
            ui.current=self.store.project(self.token,self.project_id);ui.tab='Execuções'
            async def call(fn,*args,**kwargs):return fn(*args,**kwargs)
            ui.call=call;ui.draw_project=AsyncMock()
            controls=await ui.runs_view(True)
            self.assertFalse(any(isinstance(c,ft.Text) and c.value=='Resultados disponíveis do pipeline' for c in walk(controls)))
            self.assertEqual(sum(isinstance(c,ft.Text) and c.value==ui.message('run.label.latest') for c in walk(controls)),1)
            self.assertEqual(set(ui.stage_results),{(first['id'],source['id'])})
            result=ui.stage_results[(first['id'],source['id'])]
            await result.load()
            self.assertTrue(result.loaded)
            self.assertEqual(result.run_id,first['id'])
            self.assertEqual(result.stage['artifacts'],first['stages'][0]['artifacts'])
            self.assertEqual(self.store.get_run(self.token,cancelled)['status'],'cancelled')
            self.assertTrue(all(Path(p).is_file() for p in first['stages'][0]['artifacts']))
        asyncio.run(scenario())

    def test_automatic_source_growth_cannot_restore_only_the_old_selection(self):
        from biomolexplorer.published_results import retained_configuration
        source=new_stage('retrieve_compounds');admet=new_stage('admet');admet['bindings']['base_input_path']={'stage':source['id']}
        self.save([source,admet]);manager,_=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:run=self.execute(service)
            finally:service.close()
        old=next(item for item in run['stages'] if item['id']==admet['id'])
        files=next(item['artifacts'] for item in run['stages'] if item['id']==source['id'])
        extra=self.store.project_dir(self.project_id)/'new_MOLS.csv';extra.write_text('molecule_chembl_id,canonical_smiles\nNEW,CCC\n')
        configuration=retained_configuration(self.store,self.project_id,admet,old,run['stages'],{source['id']:[*files,str(extra)]})
        self.assertIsNone(configuration)

    def test_invalid_legacy_request_does_not_hide_earlier_intact_results(self):
        from uuid import uuid4
        import time
        source=new_stage('retrieve_compounds');self.save([source]);manager,_=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:run=self.execute(service)
            finally:service.close()
        broken=copy.deepcopy(run['stages']);broken[0]['configuration']['parameters']=None
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES(?,?,?,?,?,?,?,?)',(uuid4().hex,self.project_id,self.store.user(self.token)['id'],
                'succeeded',json.dumps(broken),time.time()+1,time.time()+1,None))
        available=catalog(self.store,self.project_id,[source])
        self.assertEqual(available[source['id']]['publication_run_id'],run['id'])
