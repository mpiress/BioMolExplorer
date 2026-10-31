"""DAG ordering, durable priorities and non-destructive user pause/resume."""
import asyncio
import copy
import threading
import time
import unittest
from unittest.mock import patch,AsyncMock,Mock
import test_pipeline_execution as fixture
import test_runs_pipeline_controls as ui_fixture
from biomolexplorer.catalog import new_stage
from biomolexplorer.execution_order import ordered,move
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.workspace import WorkspaceStore


class OrderTests(unittest.TestCase):
    setUp=fixture.PipelineExecutionTests.setUp
    tearDown=fixture.PipelineExecutionTests.tearDown
    save=fixture.PipelineExecutionTests.save
    wait=fixture.PipelineExecutionTests.wait
    immediate_manager=fixture.PipelineExecutionTests.immediate_manager

    def test_move_carries_transitive_shared_dependents_and_rejects_crossing_parents(self):
        stages=[dict(id=k,depends_on=d) for k,d in [('a',[]),('b',[]),('c',['a']),('d',['c','b']),('e',[])]]
        result=move(stages,['a','b','c','d','e'],'a',1)
        self.assertEqual(result,['b','a','c','d','e'])
        self.assertEqual(move(stages,result,'c',-1),result)
        result=move(stages,result,'a',1)
        self.assertEqual(result,['b','e','a','c','d'])
        self.assertEqual(ordered(stages,result),result)

    def test_order_persists_reconciles_new_nodes_and_does_not_edit_pipeline(self):
        a,b=new_stage('retrieve_compounds'),new_stage('retrieve_structures');b['parameters']['pdb_ids']='4M0E';self.save([a,b])
        self.assertFalse(self.store.execution_order(self.token,self.project_id)['confirmed'])
        before=self.store.project(self.token,self.project_id)
        self.store.save_execution_order(self.token,self.project_id,[b['id'],a['id']])
        reopened=WorkspaceStore(self.temporary.name)
        self.assertEqual(reopened.execution_order(self.token,self.project_id),dict(order=[b['id'],a['id']],confirmed=True))
        after=reopened.project(self.token,self.project_id)
        self.assertEqual(before['pipeline'],after['pipeline'])
        self.assertEqual(before['revision'],after['revision'])
        self.assertEqual(self.store.history(self.token,self.project_id)[0]['action'],'execution_order')
        c=new_stage('retrieve_structures');self.save([a,b,c])
        self.assertEqual(self.store.execution_order(self.token,self.project_id)['order'],[b['id'],a['id'],c['id']])
        self.save([a,c]);self.assertEqual(self.store.execution_order(self.token,self.project_id)['order'],[a['id'],c['id']])

    def test_scheduler_obeys_priority(self):
        a,b=new_stage('retrieve_compounds'),new_stage('retrieve_structures');b['parameters']['pdb_ids']='4M0E';self.save([a,b])
        self.store.save_execution_order(self.token,self.project_id,[b['id'],a['id']])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store)
            try:
                run=self.wait(service.submit(self.token,self.project_id,reuse_results=False))
                self.assertEqual(run['status'],'succeeded',run.get('error'))
                self.assertEqual([x[0] for x in calls],[b['operation'],a['operation']])
                previous={item['id']:item['artifacts'][:] for item in run['stages']}
                self.store.save_execution_order(self.token,self.project_id,[a['id'],b['id']])
                again=self.wait(service.submit(self.token,self.project_id,reuse_results=True))
                self.assertEqual(again['status'],'succeeded',again.get('error'))
                self.assertEqual(len(calls),2)
                self.assertTrue(all(item.get('reused') for item in again['stages']))
                self.assertEqual(previous,{item['id']:item['artifacts'] for item in again['stages']})
            finally:service.close()

    def test_pause_finishes_current_task_survives_restart_and_resumes_without_repeating(self):
        a,b=new_stage('retrieve_compounds'),new_stage('retrieve_structures');b['parameters']['pdb_ids']='4M0E';self.save([a,b])
        original,calls=self.immediate_manager();entered,release=threading.Event(),threading.Event()
        class Blocking(original):
            def submit(manager,operation,parameters):
                manager.operation=operation
                return super().submit(operation,parameters)
            def get(manager,job_id):
                if manager.operation==a['operation']:entered.set();release.wait(8)
                return super().get(job_id)
        with patch('biomolexplorer.pipeline.JobManager',Blocking):
            service=PipelineService(self.store)
            try:
                run=service.submit(self.token,self.project_id,reuse_results=False)
                self.assertTrue(entered.wait(8));service.pause(self.token,run['id']);release.set()
                paused=self.wait(run)
                self.assertEqual(paused['status'],'awaiting_input',paused.get('error'))
                self.assertEqual(paused['stages'][0]['status'],'succeeded')
                self.assertEqual(paused['stages'][1]['pause_reason'],'user')
                files=copy.deepcopy(paused['stages'][0]['artifact_manifest'])
            finally:release.set();service.close()
            service=PipelineService(self.store)
            try:
                run=self.wait(service.continue_run(self.token,run['id']))
                self.assertEqual(run['status'],'succeeded',run.get('error'))
                self.assertEqual(len(calls),2)
                self.assertEqual(files,run['stages'][0]['artifact_manifest'])
                self.assertEqual(run['stages'][0]['artifacts'],paused['stages'][0]['artifacts'])
            finally:service.close()


class OrderDialogTests(unittest.IsolatedAsyncioTestCase):
    setUp=ui_fixture.SaveAndTemplatesTests.setUp

    async def test_first_run_asks_once_and_manual_dialog_restores_order(self):
        a,b=new_stage('retrieve_compounds'),new_stage('retrieve_structures')
        self.ui.current['pipeline']=[a,b];self.ui.dirty=True
        self.ui.service.submit=AsyncMock();self.ui.service.existing_results=lambda *args:False
        await self.ui._submit_pipeline()
        dialog=self.page.dialogs[-1]
        self.assertEqual(dialog.title.value,'Ordem de execução')
        self.ui.service.submit.assert_not_called()
        # Cancel preserves the unconfirmed state and does not start execution.
        dialog.actions[0].on_click(None)
        self.assertFalse(self.store.execution_order(self.token,self.pid)['confirmed'])
        await self.ui.execution_order_dialog()
        dialog=self.page.dialogs[-1]
        await dialog.actions[-1].on_click(None)
        self.assertTrue(self.store.execution_order(self.token,self.pid)['confirmed'])
        self.assertFalse(dialog.open)
        await self.ui.execution_order_dialog()
        self.assertEqual(self.page.dialogs[-1].title.value,'Ordem de execução')


    async def test_confirming_first_order_starts_once_and_skips_future_automatic_prompt(self):
        stage=new_stage('retrieve_compounds');self.ui.current['pipeline']=[stage];self.ui.dirty=True
        self.ui.service.existing_results=Mock(return_value=False)
        self.ui.service.submit=Mock(return_value=dict(id='run',project_id=self.pid,status='queued',stages=[]))
        self.ui.draw_project=AsyncMock()
        with patch('biomolexplorer.ui.run_progress.RunProgress'):
            await self.ui._submit_pipeline()
            dialog=self.page.dialogs[-1]
            await dialog.actions[-1].on_click(None)
            await dialog.actions[-1].on_click(None)
            self.ui.service.submit.assert_called_once()
            self.assertFalse(dialog.open)
            self.assertFalse(self.ui._configuration_loading.dialog.open)
            await self.ui._submit_pipeline()
            self.assertEqual(self.ui.service.submit.call_count,2)


class BatchCheckpointTests(unittest.TestCase):
    setUp=fixture.PipelineExecutionTests.setUp
    tearDown=fixture.PipelineExecutionTests.tearDown
    save=fixture.PipelineExecutionTests.save

    def test_repeated_pause_keeps_completed_batches_even_with_reuse_disabled(self):
        from biomolexplorer.pipeline import RunPaused
        from pathlib import Path
        stage=new_stage('fingerprints');stage['input_processing']='individual'
        first,second=copy.deepcopy(stage),copy.deepcopy(stage)
        first['parameters']['unit']='one';second['parameters']['unit']='two'
        item=dict(configuration=stage,input_mode='curated',batches=[])
        state=dict(event=threading.Event(),pause_event=threading.Event(),reuse_results=False)
        calls=[]
        def execute(*args):
            variant=args[-2];key=variant['parameters']['unit'];calls.append(key)
            path=Path(self.temporary.name)/(key+'.csv');path.write_text('id\n'+key+'\n')
            item['input_manifest']={key:key}
            if key=='one':state['pause_event'].set()
            return [str(path)]
        service=PipelineService(self.store)
        try:
            with patch.object(service,'_work_units',return_value=[(first,'one'),(second,'two')]),patch.object(service,'_cache_key',side_effect=lambda *a,**k:a[2]['parameters']['unit']),patch.object(service,'_execute_single',side_effect=execute),patch.object(service,'_update'):
                args=('run',self.project_id,'user',item,[item],Path(self.temporary.name),{},state)
                with self.assertRaises(RunPaused):service._execute_stage(*args)
                self.assertEqual(len(item['batches']),1)
                with self.assertRaises(RunPaused):service._execute_stage(*args)
                self.assertEqual(len(item['batches']),1)
                state['pause_event'].clear()
                result=service._execute_stage(*args)
                self.assertEqual(calls,['one','two'])
                self.assertEqual(len(result),2)
        finally:service.close()

    def test_late_pause_never_resurrects_completed_run(self):
        stage=new_stage('retrieve_compounds');self.save([stage])
        service=PipelineService(self.store)
        try:
            with patch.object(service._executor,'submit'):
                run=service.submit(self.token,self.project_id)
                state=service._active[run['id']]
                state['stages'][0]['status']='succeeded'
                service._update(run['id'],'succeeded',state['stages'])
                stale=dict(run,status='running')
                original=self.store.get_run;reads=iter([stale])
                with patch.object(self.store,'get_run',side_effect=lambda *args:next(reads,None) or original(*args)):
                    result=service.pause(self.token,run['id'])
                self.assertEqual(result['status'],'succeeded')
                self.assertEqual(original(self.token,run['id'])['status'],'succeeded')
        finally:service.close()

    def test_cancel_paused_execution_does_not_start_jobs(self):
        stage=new_stage('retrieve_compounds');self.save([stage])
        service=PipelineService(self.store)
        try:
            with patch.object(service._executor,'submit'),patch('biomolexplorer.pipeline.JobManager') as manager:
                run=service.submit(self.token,self.project_id)
                service.pause(self.token,run['id'])
                state=service._active[run['id']]
                service._run(run['id'],self.project_id,self.store.user(self.token)['id'],state['stages'])
                paused=self.store.get_run(self.token,run['id'])
                self.assertEqual(paused['status'],'awaiting_input')
                service.cancel(self.token,run['id'])
                self.assertEqual(self.store.get_run(self.token,run['id'])['status'],'cancelled')
                manager.assert_not_called()
        finally:service.close()


class PauseFeedbackTests(unittest.TestCase):
    def test_progress_exposes_pause_then_resume_and_blocks_viewers(self):
        import test_run_progress as feedback
        fixture=feedback.RunFeedbackTests();fixture.setUp()
        ui,progress,run=fixture.ui,fixture.progress,fixture.run
        ui.pause_run=AsyncMock();ui.continue_paused_run=AsyncMock()
        self.assertTrue(progress.pause.visible);self.assertFalse(progress.resume.visible)
        asyncio.run(progress.request_pause(None));ui.pause_run.assert_awaited_once_with('run')
        run['status']='awaiting_input';run['stages'][0].update(status='awaiting_input',pause_reason='user')
        progress.update(run)
        self.assertFalse(progress.pause.visible);self.assertTrue(progress.resume.visible)
        self.assertTrue(progress.results.visible)
        self.assertIn('pausada',progress.heading.value.lower())
        asyncio.run(progress.request_resume(None));ui.continue_paused_run.assert_awaited_once_with('run')
        ui.current['role']='viewer';progress.update(run)
        self.assertFalse(progress.resume.visible)
        asyncio.run(progress.request_resume(None));self.assertEqual(ui.continue_paused_run.await_count,1)
