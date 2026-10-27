"""Cancellation buttons stop the run and clear pending configuration views."""
import asyncio
import copy
import threading
import unittest
from types import SimpleNamespace
from unittest.mock import AsyncMock,Mock,patch

from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.ui.app import WorkspaceUI
from biomolexplorer.ui.run_progress import RunProgress
from test_frontend_navigation import Page
import test_pipeline_execution as execution


class CancellationUITests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def pending_run(self,service):
        source=new_stage('retrieve_compounds');pending=new_stage('admet')
        pending['bindings']={'base_input_path':{'stage':source['id']}}
        self.save([source,pending])
        run=self.wait(service.submit(self.token,self.project_id))
        self.assertEqual(run['status'],'awaiting_input',run.get('error'))
        return run,pending

    def ui(self,service,run):
        ui=WorkspaceUI(Page(),self.store,service)
        self.addCleanup(ui.structure_viewer.close)
        ui.token=self.token;ui.current=self.store.project(self.token,self.project_id)
        ui.base_pipeline=copy.deepcopy(ui.current['pipeline'])
        ui.current_run=run['id'];ui.polling=None
        ui.draw_project=AsyncMock();ui.notify=Mock()
        async def call(function,*args,**kwargs):return function(*args,**kwargs)
        ui.call=call
        ui.run_progress=RunProgress(ui,run)
        return ui

    def assert_cancelled(self,ui,run,dialog,calls):
        updated=self.store.get_run(self.token,run['id'])
        self.assertEqual(updated['status'],'cancelled')
        self.assertEqual([s['status'] for s in updated['stages']],['succeeded','cancelled'])
        self.assertEqual(ui.run_progress.run['status'],'cancelled')
        self.assertFalse(ui.run_progress.cancel.visible)
        self.assertFalse(dialog.open)
        self.assertFalse(ui.editing_stage)
        self.assertIsNone(ui.current_run)
        ui.draw_project.assert_awaited()
        ui.notify.assert_called_with('Execução cancelada.')
        self.assertEqual([op for op,_ in calls],['retrieve_compounds'])

    def test_settings_cancel_button_cancels_pending_execution(self):
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store);self.addCleanup(service.close)
            run,pending=self.pending_run(service);ui=self.ui(service,run)
            async def scenario():
                await ui.open_stage_dialog(pending['id'],waiting_run_id=run['id'])
                dialog=ui.page.dialogs[-1]
                self.assertEqual(dialog.actions[0].content,'Cancelar execução')
                await dialog.actions[0].on_click(SimpleNamespace(control=dialog.actions[0]))
                self.assert_cancelled(ui,run,dialog,calls)
                self.assertIsNone(ui.settings_body)
            asyncio.run(scenario())

    def test_file_selection_cancel_updates_without_polling_and_discards_draft(self):
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store);self.addCleanup(service.close)
            run,pending=self.pending_run(service);ui=self.ui(service,run)
            async def scenario():
                await ui.configure_pending(run['id'])
                dialog=ui.selection_dialog;key=ui.selection_dialog_key
                cancel=next(action for action in dialog.actions if action.content=='Cancelar execução')
                await cancel.on_click(SimpleNamespace(control=cancel))
                self.assert_cancelled(ui,run,dialog,calls)
                self.assertNotIn(key,ui.pending_file_selections)
                # A delayed Flutter dismissal must not resurrect cancelled drafts.
                dialog.on_dismiss(None)
                self.assertNotIn(key,ui.pending_file_selections)
            asyncio.run(scenario())

    def test_running_cancel_restarts_monitoring_until_terminal_state(self):
        manager,calls=self.immediate_manager();entered=threading.Event()
        class Running(manager):
            def submit(self,operation,parameters):
                job=super().submit(operation,parameters);job['status']='running';return job
            def get(self,job_id):
                entered.set();return super().get(job_id)
        self.save([new_stage('retrieve_compounds'),new_stage('retrieve_compounds')])
        with patch('biomolexplorer.pipeline.JobManager',Running):
            service=PipelineService(self.store);self.addCleanup(service.close)
            run=service.submit(self.token,self.project_id)
            self.assertTrue(entered.wait(5))
            ui=self.ui(service,run)
            async def scenario():
                ui.page.run_task=lambda handler,*args:asyncio.create_task(handler(*args))
                await ui.cancel_run(run['id'])
                if isinstance(ui.polling,asyncio.Task):await asyncio.wait_for(ui.polling,4)
                updated=self.store.get_run(self.token,run['id'])
                self.assertEqual(updated['status'],'cancelled')
                self.assertEqual([s['status'] for s in updated['stages']],['cancelled','skipped'])
                self.assertEqual(ui.run_progress.run['status'],'cancelled')
                self.assertFalse(ui.run_progress.cancel.visible)
                self.assertIsNone(ui.current_run)
                self.assertEqual(len(calls),1)
            asyncio.run(scenario())


if __name__=='__main__':unittest.main()
