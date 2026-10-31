"""Whole-flow PNG export is independent of viewport and editing permissions."""
import asyncio
import copy
import io
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import AsyncMock, Mock

from PIL import Image
import flet as ft
from biomolexplorer import flow
from biomolexplorer.catalog import new_stage
from biomolexplorer.ui.flow_canvas import FlowCanvas
from biomolexplorer.ui.pipeline_png import MAX_PIXELS, render_pipeline
from test_frontend_navigation import walk


class PipelinePngTests(unittest.TestCase):
    def stages(self):
        source=new_stage('retrieve_structures');source['name']='Estruturas PDB · α'
        target=new_stage('redocking');target['name']='Redocking incremental'
        source['position']={'x':50,'y':60};target['position']={'x':1500,'y':900}
        flow.connect([source,target],source['id'],target['id'],'base_input_path')
        return [source,target]

    def test_entire_flow_including_far_block_is_rendered_without_mutation(self):
        stages=self.stages();before=copy.deepcopy(stages)
        data=render_pipeline(stages);image=Image.open(io.BytesIO(data))
        self.assertEqual(image.format,'PNG');self.assertEqual(stages,before)
        self.assertEqual(image.size,((1500-50+flow.NODE_WIDTH+260)*2,(900-60+flow.node_height(stages[-1])+260)*2))
        self.assertGreater(image.info['dpi'][0],299)
        # Both headers appear, including the block outside the default viewport.
        for stage in stages:
            x=round((stage['position']['x']-50+130+flow.NODE_WIDTH-18)*2)
            y=round((stage['position']['y']-60+130+20)*2)
            self.assertEqual(image.getpixel((x,y)),(13,148,136) if stage['operation']=='retrieve_structures' else (217,119,6))
        # The connection is drawn in the large empty region between nodes.
        self.assertTrue(any(image.getpixel((image.width//2,y))==(13,148,136) for y in range(image.height)))

    def test_molecular_dynamics_category_exports_with_shared_palette(self):
        from biomolexplorer.catalog import TITLES
        from biomolexplorer.ui.pipeline_palette import COLORS
        self.assertTrue({v[2] for v in TITLES.values()}<=set(COLORS))
        md=next(k for k,v in TITLES.items() if v[2]=='Dinâmica Molecular')
        image=Image.open(io.BytesIO(render_pipeline([new_stage(md)])))
        self.assertEqual(image.format,'PNG')

    def test_large_flow_has_bounded_memory_and_empty_flow_exports(self):
        stages=self.stages();stages[-1]['position']={'x':11000,'y':39000}
        image=Image.open(io.BytesIO(render_pipeline(stages)))
        self.assertLessEqual(image.width*image.height,MAX_PIXELS)
        self.assertGreater(image.height,5000)
        self.assertEqual(Image.open(io.BytesIO(render_pipeline([]))).format,'PNG')

    def ui(self):
        async def call(fn,*args):return fn(*args)
        async def guard(action):await action()
        return SimpleNamespace(current={'pipeline':self.stages()},selected=None,token='session',tr=lambda s:s,
            page=SimpleNamespace(update=Mock()),call=call,guard=guard,picker=SimpleNamespace(save_file=AsyncMock()))

    def test_camera_is_available_to_viewers_and_saves_png(self):
        ui=self.ui();editor=FlowCanvas(ui,False)
        camera=next(c for c in walk(editor.build(False)) if isinstance(c,ft.IconButton) and c.icon==ft.Icons.PHOTO_CAMERA_OUTLINED)
        self.assertFalse(camera.disabled)
        asyncio.run(camera.on_click(None))
        kwargs=ui.picker.save_file.await_args.kwargs
        self.assertEqual(kwargs['file_name'],'pipeline.png')
        self.assertEqual(Image.open(io.BytesIO(kwargs['src_bytes'])).format,'PNG')

    def test_switching_project_during_export_does_not_save_old_project(self):
        ui=self.ui();editor=FlowCanvas(ui,False)
        async def call(fn,*args):
            result=fn(*args);ui.current={'pipeline':[]};return result
        ui.call=call
        asyncio.run(editor.export_png(None))
        ui.picker.save_file.assert_not_awaited()


if __name__=='__main__':unittest.main()
