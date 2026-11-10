"""Canvas work scales with affected visual connections, not selected files."""
import asyncio
import copy
from types import SimpleNamespace
import unittest
from unittest.mock import Mock,patch

import flet as ft
from biomolexplorer.catalog import new_stage
from biomolexplorer.ui.flow_canvas import FlowCanvas
from biomolexplorer.ui.localization import LocalizedPage
import test_pipeline_editor_regressions as fixtures


def pipeline(count=100,files_per_connection=1):
    stages=[new_stage('retrieve_compounds')]+[new_stage('admet') for _ in range(count-1)]
    for index,stage in enumerate(stages):
        stage['name']=f'Block {index+1:03d}'
        if index:
            refs=[dict(stage=stages[index-1]['id'],selector=f'file-{i}.csv') for i in range(files_per_connection)]
            stage['bindings']['base_input_path']=refs[0] if len(refs)==1 else {'sources':refs}
    return stages


class CanvasScalingTests(unittest.TestCase):
    def editor(self,stages):
        ui=SimpleNamespace(current={'pipeline':stages},selected=None,dirty=False,
            page=SimpleNamespace(update=Mock()),notify=Mock(),schedule_save=Mock())
        return FlowCanvas(ui,True)

    def test_file_selections_publish_one_visual_edge_without_altering_bindings(self):
        stages=pipeline(100,50);original=copy.deepcopy(stages)
        editor=self.editor(stages)
        self.assertEqual(len(editor.edges.shapes),99)
        self.assertEqual(stages,original)
        self.assertEqual(len(editor._node_states[stages[-1]['id']][2]),1)
        self.assertIn('(+49)',editor._node_states[stages[-1]['id']][2][0])

    def test_drag_touches_only_incident_edges_and_preserves_unrelated_geometry(self):
        stages=pipeline();editor=self.editor(stages)
        node=editor.nodes[stages[40]['id']];header=node.content.controls[0]
        unrelated=(stages[70]['id'],stages[71]['id'],'base_input_path')
        shape=editor._edge_shapes[unrelated];elements=shape.elements
        with patch.object(editor,'_edge_elements',wraps=editor._edge_elements) as geometry, \
             patch('biomolexplorer.flow.position',side_effect=AssertionError('Full graph scan during drag')):
            for _ in range(10):header.on_pan_update(SimpleNamespace(local_delta=SimpleNamespace(x=1,y=1)))
        self.assertEqual(geometry.call_count,20)
        self.assertIs(editor._edge_shapes[unrelated],shape);self.assertIs(shape.elements,elements)
        with patch.object(editor,'node',wraps=editor.node) as rebuild:
            editor.refresh()
        self.assertEqual(rebuild.call_count,0)
        self.assertEqual(editor.ui.schedule_save.call_count,0)
        header.on_pan_end(None)
        editor.ui.schedule_save.assert_called_once()

    def test_geometry_layers_do_not_walk_thousands_of_localization_attributes(self):
        editor=self.editor(pipeline(100,50))
        page=LocalizedPage(SimpleNamespace(update=Mock(),controls=[],overlay=[]))
        with patch.object(page.translator,'__call__',side_effect=AssertionError('Geometry has no labels')):
            page.localize(editor.grid);page.localize(editor.edges)
        self.assertFalse(hasattr(editor.edges,'_biomol_messages'))

    def test_grid_growth_keeps_edges_attached_and_layout_is_built_once(self):
        editor=self.editor(pipeline(2))
        view=editor.build(False)
        self.assertIs(view,editor.build(False))
        node=editor.nodes[editor.stages[-1]['id']]
        node.content.controls[0].on_pan_update(SimpleNamespace(local_delta=SimpleNamespace(x=4500,y=3500)))
        self.assertIn(editor.edges,editor.stack.controls)
        self.assertGreater(editor.grid.width,node.left+290)
        self.assertGreater(editor.grid.height,node.top+202)
        self.assertEqual(len(editor.edges.shapes),1)

    def test_palette_search_updates_only_the_library_and_preserves_category_state(self):
        editor=self.editor(pipeline(2));editor.build()
        section=next(control for control in editor.library.controls if isinstance(control,ft.ExpansionTile))
        section.on_change(SimpleNamespace(data='true'))
        search=next(control for control in editor.library.controls if isinstance(control,ft.TextField))
        search.value='ADMET';search.on_change(None)
        editor.ui.page.update.assert_called_with(editor.library)
        search.value='';search.on_change(None)
        self.assertTrue(section.expanded)


class CanvasOpeningTests(unittest.TestCase):
    setUp=fixtures.PipelineEditorRegressionTests.setUp

    def test_open_and_reenter_canvas_do_not_read_scientific_history_or_reset_viewport(self):
        async def scenario():
            with patch.object(self.store,'list_runs',side_effect=AssertionError('Canvas must not load execution history')):
                first=await self.ui.pipeline_view(True)
                editor=self.ui.flow_editor
                editor.viewport_scale=.4;editor.viewport_offset=ft.Offset(-1200,-800)
                second=await self.ui.pipeline_view(True)
                self.assertIs(self.ui.flow_editor,editor);self.assertIs(first,second)
                self.assertEqual(editor.viewport_scale,.4)
                self.assertEqual(editor.viewport_offset,ft.Offset(-1200,-800))
                await self.ui.pipeline_view(False)
                self.assertIsNot(self.ui.flow_editor,editor)
                self.assertFalse(self.ui.flow_editor.writable)
        asyncio.run(scenario())

    def test_remote_edit_keeps_the_canvas_but_discards_stale_undo_snapshots(self):
        from unittest.mock import AsyncMock
        async def scenario():
            self.ui.dirty=False
            await self.ui.pipeline_view(True)
            editor=self.ui.flow_editor;editor.remember()
            fresh=copy.deepcopy(self.ui.current)
            fresh['revision']+=1;fresh['updated']+=1
            fresh['pipeline'][0]['name']='Collaborator edit'
            self.ui.call=AsyncMock(return_value=fresh)
            self.ui.draw_project=AsyncMock()
            with patch('biomolexplorer.ui.project_tools.asyncio.sleep',side_effect=[None,asyncio.CancelledError()]):
                with self.assertRaises(asyncio.CancelledError):await self.ui.poll_project(self.pid)
            self.assertIs(self.ui.flow_editor,editor)
            self.assertEqual(editor.history,[])
            self.assertEqual(self.ui.current['pipeline'][0]['name'],'Collaborator edit')
        asyncio.run(scenario())
