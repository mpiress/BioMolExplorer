"""Actual UI callbacks, persisted state, single-flight feedback and presets."""
import asyncio
import copy
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import AsyncMock,Mock,patch
import flet as ft
from biomolexplorer.catalog import new_stage,PRESETS
from biomolexplorer.pipeline import validate_pipeline
from biomolexplorer.ui.app import WorkspaceUI
from biomolexplorer.workspace import WorkspaceStore
from test_frontend_navigation import Page,walk
import test_clear_all_results as clear_fixture


class SaveAndTemplatesTests(unittest.IsolatedAsyncioTestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup)
        self.store=WorkspaceStore(temporary.name);self.token=self.store.register('Editor','controls@example.org','controls-password')
        self.project=self.store.create_project(self.token,'Controls');self.pid=self.project['id']
        self.page=Page();self.ui=WorkspaceUI(self.page,self.store,SimpleNamespace(dock6_path=None))
        self.addCleanup(self.ui.structure_viewer.close)
        self.ui.token=self.token;self.ui.current=self.project
        self.ui.notify=Mock()
        async def call(function,*args,**kwargs):return function(*args,**kwargs)
        self.ui.call=call

    async def test_sidebar_save_persists_even_if_dirty_flag_was_missing_and_blocks_repeats(self):
        self.ui.current['pipeline']=[new_stage('retrieve_structures')];self.ui.dirty=True
        await self.ui.draw_project()
        button=next(c for c in walk(self.page.controls) if isinstance(c,ft.Button) and getattr(c.content,'value',None)=='Salvar')
        self.ui.current['pipeline'][0]['name']='Saved manually'
        self.ui.dirty=False
        entered,release=asyncio.Event(),asyncio.Event();writes=[]
        async def call(function,*args,**kwargs):
            if function==self.store.save_pipeline:
                writes.append(copy.deepcopy(args));entered.set();await release.wait()
            return function(*args,**kwargs)
        self.ui.call=call
        pending=asyncio.create_task(button.on_click(None))
        await asyncio.wait_for(entered.wait(),3)
        self.assertTrue(self.ui._configuration_loading.dialog.open)
        await button.on_click(None);self.assertEqual(len(writes),1)
        release.set();await pending
        reopened=WorkspaceStore(self.store.root)
        self.assertEqual(reopened.project(self.token,self.pid)['pipeline'][0]['name'],'Saved manually')
        self.assertFalse(self.ui._configuration_loading.dialog.open)
        self.ui.notify.assert_called_once_with('Pipeline salvo.')

    async def test_save_acknowledges_autosaved_state_and_recovers_after_failure(self):
        await self.ui.event(self.ui.save_pipeline_action)(None)
        self.ui.notify.assert_called_with('Pipeline salvo.')
        stage=new_stage('retrieve_compounds');self.ui.current['pipeline']=[stage]
        original=self.ui.call
        async def failed(function,*args,**kwargs):
            if function==self.store.save_pipeline:raise OSError('save failed')
            return await original(function,*args,**kwargs)
        self.ui.call=failed
        await self.ui.event(self.ui.save_pipeline_action)(None)
        self.assertFalse(self.ui._manual_save_running)
        self.assertFalse(self.ui._configuration_loading.dialog.open)
        self.assertEqual(self.store.project(self.token,self.pid)['pipeline'],[])
        self.ui.call=original
        await self.ui.event(self.ui.save_pipeline_action)(None)
        self.assertEqual(self.store.project(self.token,self.pid)['pipeline'][0]['id'],stage['id'])

    async def test_every_template_applies_persists_and_cannot_apply_twice(self):
        self.ui.draw_project=AsyncMock()
        for name,operations in PRESETS.items():
            with self.subTest(template=name):
                count=len(self.ui.current['pipeline'])
                await self.ui.preset_dialog();dialog=self.page.dialogs[-1]
                choice=next(c for c in walk(dialog) if isinstance(c,ft.Dropdown))
                choice.value=name
                await dialog.actions[-1].on_click(None)
                await dialog.actions[-1].on_click(None)
                pipeline=self.store.project(self.token,self.pid)['pipeline']
                self.assertEqual(len(pipeline),count+len(operations))
                self.assertEqual([s['operation'] for s in pipeline[count:]],operations)
                validate_pipeline(pipeline)
                self.assertFalse(dialog.open)
                self.assertFalse(self.ui._configuration_loading.dialog.open)
                self.assertEqual(self.ui.current['pipeline'],pipeline)

    async def test_selected_run_requires_an_explicit_selected_stage(self):
        self.ui.current['pipeline']=[new_stage('retrieve_compounds')]
        self.ui.service.submit=Mock()
        with self.assertRaisesRegex(ValueError,'Selecione um bloco'):
            await self.ui._submit_pipeline(selected=True)
        self.ui.service.submit.assert_not_called()

    async def test_selected_stage_toolbar_passes_only_target_to_scheduler(self):
        source=new_stage('retrieve_compounds');selected=new_stage('admet')
        selected['bindings']['base_input_path']={'stage':source['id']}
        self.ui.current['pipeline']=[source,selected];self.ui.dirty=True
        await self.ui.draw_project();self.ui.selected=selected['id']
        button=next(c for c in walk(self.page.controls) if isinstance(c,ft.Button) and getattr(c.content,'value',None)=='Executar até seleção')
        self.ui.service.existing_results=Mock(return_value=0)
        self.ui.service.submit=Mock(return_value=dict(id='run',project_id=self.pid,status='queued',stages=[]))
        self.ui.draw_project=AsyncMock()
        await self.ui.save_draft()
        self.store.save_execution_order(self.token,self.pid,[source['id'],selected['id']])
        with patch('biomolexplorer.ui.run_progress.RunProgress'):
            await button.on_click(None)
        self.ui.service.submit.assert_called_once_with(self.token,self.pid,[selected['id']],reuse_results=True)
        self.assertFalse(self.ui.submitting)



class CleanupCallbackTests(unittest.IsolatedAsyncioTestCase):
    def setUp(self):
        clear_fixture.ClearAllResultsTests.setUp(self)
        self.page=Page();self.ui=WorkspaceUI(self.page,self.store,SimpleNamespace(dock6_path=None))
        self.addCleanup(self.ui.structure_viewer.close)
        self.ui.token=self.token;self.ui.current=self.store.project(self.token,self.pid)
        self.ui.base_pipeline=copy.deepcopy(self.ui.current['pipeline'])
        self.ui.draw_project=AsyncMock();self.ui.notify=Mock()
        async def call(function,*args,**kwargs):return function(*args,**kwargs)
        self.ui.call=call

    async def exercise(self,results_only):
        before=copy.deepcopy(self.ui.current['pipeline']);original_bytes=self.original.read_bytes()
        await self.ui.clear_all_results(results_only)
        dialog=self.page.dialogs[-1]
        self.assertTrue(self.generated.exists())
        await self.ui.clear_all_results(results_only)
        self.assertIs(self.page.dialogs[-1],dialog)
        dialog.actions[0].on_click(None)
        self.assertTrue(self.generated.exists())
        await self.ui.clear_all_results(results_only);dialog=self.page.dialogs[-1]
        entered,release=asyncio.Event(),asyncio.Event();deletions=[]
        async def call(function,*args,**kwargs):
            if getattr(function,'__name__','') in ('remove_project_results','clear_project_results'):
                deletions.append(function.__name__);entered.set();await release.wait()
            return function(*args,**kwargs)
        self.ui.call=call
        pending=asyncio.create_task(dialog.actions[-1].on_click(None))
        await asyncio.wait_for(entered.wait(),3)
        self.assertTrue(self.ui._configuration_loading.dialog.open)
        self.assertTrue(all(c.disabled for c in dialog.actions))
        await dialog.actions[-1].on_click(None)
        release.set();await pending
        await dialog.actions[-1].on_click(None)
        self.assertEqual(len(deletions),1)
        self.assertFalse(self.generated.exists());self.assertEqual(self.original.read_bytes(),original_bytes)
        pipeline=self.store.project(self.token,self.pid)['pipeline']
        if results_only:self.assertEqual(pipeline,before)
        else:
            self.assertEqual([s['id'] for s in pipeline],[s['id'] for s in before])
            self.assertNotIn('pdb_codes',pipeline[2]['parameters'])
            self.assertEqual(pipeline[2]['bindings']['base_input_path']['stage'],before[2]['bindings']['base_input_path']['stage'])
        self.assertFalse(self.ui._configuration_loading.dialog.open)

    async def test_remove_all_confirmation_and_real_persistence(self):await self.exercise(True)
    async def test_clear_all_confirmation_and_real_persistence(self):await self.exercise(False)

    async def test_failed_removal_leaves_data_and_restores_retry(self):
        await self.ui.clear_all_results(True);dialog=self.page.dialogs[-1]
        original=self.ui.call
        async def failed(function,*args,**kwargs):
            if getattr(function,'__name__','')=='remove_project_results':raise OSError('removal failed')
            return await original(function,*args,**kwargs)
        self.ui.call=failed
        await dialog.actions[-1].on_click(None)
        self.assertTrue(self.generated.exists())
        self.assertTrue(all(not c.disabled for c in dialog.actions))
        self.assertFalse(self.ui._configuration_loading.dialog.open)
        self.ui.call=original
        await dialog.actions[-1].on_click(None)
        self.assertFalse(self.generated.exists())


class DirectedExecutionTests(unittest.TestCase):
    from test_pipeline_execution import PipelineExecutionTests as _Fixture
    setUp=_Fixture.setUp
    tearDown=_Fixture.tearDown
    upload=_Fixture.upload
    imported=_Fixture.imported
    save=_Fixture.save
    wait=_Fixture.wait
    finish_selected=_Fixture.finish_selected
    immediate_manager=_Fixture.immediate_manager

    def test_selected_stage_excludes_both_descendants_and_independent_branches(self):
        import sys
        from biomolexplorer.pipeline import PipelineService
        source=self.imported(self.upload());selected=new_stage('admet')
        selected['bindings']['base_input_path']={'stage':source['id']}
        downstream=new_stage('fingerprints');downstream['bindings']['base_input_path']={'stage':selected['id']}
        independent=new_stage('retrieve_compounds')
        self.save([downstream,independent,selected,source]);manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:run=self.finish_selected(service,service.submit(self.token,self.project_id,[selected['id']]))
            finally:service.close()
        self.assertEqual(run['status'],'succeeded',run['error'])
        self.assertEqual([item['id'] for item in run['stages']],[source['id'],selected['id']])
        self.assertEqual([operation for operation,_ in calls],['admet'])
