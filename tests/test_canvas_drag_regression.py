"""Drag identity survives redraws; saved results remain usable by new blocks."""
import asyncio
import copy
import unittest
from types import SimpleNamespace
from unittest.mock import Mock,patch

import flet as ft
from biomolexplorer.catalog import new_stage,operation_fields
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.ui.flow_canvas import FlowCanvas
from biomolexplorer.ui.input_editor import InputEditor
from biomolexplorer.ui.app import WorkspaceUI
from test_frontend_navigation import Page,walk
import test_pipeline_execution as execution


def output(editor,stage_id):
    return next(c for c in walk(editor.stack) if isinstance(c,ft.Draggable) and c.data==stage_id)


def port(editor,stage):
    return next(c for c in walk(editor.node(stage,editor.stages.index(stage))) if isinstance(c,ft.DragTarget))


class CanvasDragRegressionTests(unittest.TestCase):
    def setUp(self):
        self.source=new_stage('retrieve_compounds');self.target=new_stage('admet')
        self.ui=SimpleNamespace(current={'pipeline':[self.source,self.target]},selected=None,
            dirty=False,page=SimpleNamespace(update=Mock()),notify=Mock())
        self.editor=FlowCanvas(self.ui,True)

    def test_source_control_identity_survives_redraw_and_autosave_state_replacement(self):
        original=output(self.editor,self.source['id'])
        self.ui.current=copy.deepcopy(self.ui.current)
        self.editor.refresh()
        self.assertIs(output(self.editor,self.source['id']),original)
        target=self.editor.stages[1]
        asyncio.run(port(self.editor,target).on_accept(SimpleNamespace(src=original)))
        self.assertEqual(target['bindings']['base_input_path']['stage'],self.source['id'])

    def test_delayed_drop_resolves_only_the_exact_captured_drag_source(self):
        original=output(self.editor,self.source['id'])
        original.on_drag_start(SimpleNamespace(control=original))
        self.editor.refresh()
        asyncio.run(port(self.editor,self.target).on_accept(SimpleNamespace(src=None,src_id=original._i)))
        self.assertEqual(self.target['bindings']['base_input_path']['stage'],self.source['id'])
        self.assertIsNone(self.editor.active_drag)
        self.ui.notify.assert_not_called()

    def test_unresolvable_or_removed_source_reports_error_without_mutation(self):
        original=output(self.editor,self.source['id'])
        original.on_drag_start(SimpleNamespace(control=original))
        before=copy.deepcopy(self.editor.stages)
        asyncio.run(port(self.editor,self.target).on_accept(SimpleNamespace(src=None,src_id=-1)))
        self.assertEqual(self.editor.stages,before)
        self.assertIn('origem do arraste',self.ui.notify.call_args.args[0])
        original.on_drag_start(SimpleNamespace(control=original))
        self.editor.stages.remove(self.source);self.editor.refresh()
        asyncio.run(port(self.editor,self.target).on_accept(SimpleNamespace(src=None,src_id=original._i)))
        self.assertEqual(self.target['bindings'],{})
        self.assertEqual(self.ui.notify.call_args.args[0],'O bloco da conexão não existe.')

    def test_library_drop_also_survives_missing_live_source(self):
        draggable=next(c for c in walk(self.editor.build()) if isinstance(c,ft.Draggable) and c.group=='stage-library')
        draggable.on_drag_start(SimpleNamespace(control=draggable))
        asyncio.run(self.editor.drop(SimpleNamespace(src=None,src_id=draggable._i,
            local_position=SimpleNamespace(x=900,y=200))))
        self.assertEqual(self.editor.stages[-1]['operation'],draggable.data)
        self.assertEqual(len(self.editor.stages),3)

    def test_unexecuted_source_connects_without_requiring_result_metadata(self):
        self.ui.artifact_items={};self.ui.artifact_choices={}
        asyncio.run(port(self.editor,self.target).on_accept(SimpleNamespace(src=output(self.editor,self.source['id']))))
        field=next(f for f in operation_fields('admet') if f['name']=='base_input_path')
        editor=InputEditor(self.ui,self.target,field,[],True)
        self.assertIsNone(editor.selection)
        self.assertEqual(editor.read(),{'stage':self.source['id'],'selector':'auto'})


class ExecutedPipelineEditingTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def test_new_consumer_recovers_persisted_outputs_and_reuses_producer(self):
        producer=new_stage('retrieve_compounds');self.save([producer])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store)
            try:
                original=self.wait(service.submit(self.token,self.project_id,reuse_results=False))
                self.assertEqual(original['status'],'succeeded',original.get('error'))
                metadata=copy.deepcopy(original['stages'][0]['artifact_records'])
                ui=WorkspaceUI(Page(),self.store,service)
                self.addCleanup(ui.structure_viewer.close)
                ui.token=self.token;ui.current=self.store.project(self.token,self.project_id)
                ui.schedule_save=lambda:None
                async def call(function,*args,**kwargs):return function(*args,**kwargs)
                ui.call=call
                canvas=FlowCanvas(ui,True)
                asyncio.run(canvas.add('admet'))
                consumer=canvas.stages[-1]
                draggable=output(canvas,producer['id'])
                draggable.on_drag_start(SimpleNamespace(control=draggable))
                canvas.refresh()
                asyncio.run(port(canvas,consumer).on_accept(SimpleNamespace(src=None,src_id=draggable._i)))
                # Opening the real dialog discovers outputs from SQLite, not an old UI session.
                asyncio.run(ui.open_stage_dialog(consumer['id']))
                field=next(f for f in operation_fields('admet') if f['name']=='base_input_path')
                editor=InputEditor(ui,consumer,field,[],True)
                picker=editor.selection.pickers['base_input_path'];picker.add_all_sources()
                consumer['bindings']['base_input_path']=editor.read()
                self.assertEqual(consumer['bindings']['base_input_path']['artifact_id'],metadata[0]['artifact_id'])
                self.save(canvas.stages)
                final=self.wait(service.submit(self.token,self.project_id,reuse_results=True))
                self.assertEqual(final['status'],'succeeded',final.get('error'))
                self.assertTrue(final['stages'][0]['reused'])
                self.assertEqual([op for op,_ in calls],['retrieve_compounds','admet'])
                self.assertEqual(self.store.get_run(self.token,original['id'])['stages'][0]['artifact_records'],metadata)
                asyncio.run(canvas.link(producer['id'],consumer['id'],'base_input_path'))
                from biomolexplorer.flow import disconnect
                disconnect(consumer,'base_input_path')
                refreshed=InputEditor(ui,consumer,field,[],True)
                self.assertIsNone(refreshed.selection)
            finally:service.close()

