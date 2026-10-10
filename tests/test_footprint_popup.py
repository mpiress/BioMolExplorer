"""The integrated footprint preview keeps downloads authorized and unchanged."""
import asyncio
import io
import shutil
import unittest
from types import SimpleNamespace
from unittest.mock import AsyncMock,patch

import flet as ft
from test_docking_scene import SceneCapabilityTests
from biomolexplorer.docking_results import DockingResults
from biomolexplorer.pdf_preview import pdf_preview
from biomolexplorer.ui.footprint import footprint_dialog
from biomolexplorer.workspace import AccessDenied

from PIL import Image
_preview=io.BytesIO();Image.new('RGB',(1100,1000),'white').save(_preview,format='PNG')
PREVIEW=_preview.getvalue()


class FootprintPopupTests(unittest.TestCase):
    setUp=SceneCapabilityTests.setUp
    save=SceneCapabilityTests.save
    docking=SceneCapabilityTests.docking

    def ui(self):
        async def call(fn,*args,**kwargs):return fn(*args,**kwargs)
        async def guard(fn):return await fn()
        def show(dialog):dialog.open=True;page.dialogs.append(dialog)
        page=SimpleNamespace(width=1200,height=900,on_resize=None,dialogs=[],update=lambda:None,show_dialog=show)
        return SimpleNamespace(store=self.store,token=self.token,page=page,call=call,guard=guard,
            picker=SimpleNamespace(save_file=AsyncMock()))

    def selection(self,table):
        data=DockingResults(self.store).page(self.token,self.pid,self.rid,self.sid,str(table))
        return dict(table=str(table),index=0,version=data['version'])

    def test_popup_renders_graph_downloads_original_and_restores_resize(self):
        async def scenario():
            table,receptor,pose,pdf=self.docking('dock6');ui=self.ui()
            previous=AsyncMock();ui.page.on_resize=previous
            with patch('biomolexplorer.ui.footprint.pdf_preview',return_value=PREVIEW):
                await footprint_dialog(ui,self.pid,self.rid,self.sid,self.selection(table),'M1',lambda:True,origin='docked_pose')
            dialog=ui.page.dialogs[-1]
            controls=dialog.content.content.controls
            self.assertTrue(controls[2].visible)
            self.assertIsInstance(controls[3].content,ft.InteractiveViewer)
            self.assertEqual(controls[3].content.content.src,PREVIEW)
            ui.page.width=400;ui.page.height=650
            await ui.page.on_resize(None)
            previous.assert_awaited_once()
            self.assertLessEqual(dialog.content.width,ui.page.width)
            self.assertAlmostEqual(controls[3].width/controls[3].height,1.1)
            self.assertLessEqual(dialog.content.height,ui.page.height-170)
            await dialog.actions[0].on_click(None)
            ui.picker.save_file.assert_awaited_once_with(file_name='M1_footprint.pdf',src_bytes=pdf.read_bytes())
            table.write_text(table.read_text().replace('M1','M2'))
            with self.assertRaises(ValueError):await dialog.actions[0].on_click(None)
            dialog.actions[1].on_click(None)
            self.assertFalse(dialog.open)
            self.assertIs(ui.page.on_resize,previous)
        asyncio.run(scenario())

    def test_preview_failure_still_allows_pdf_download(self):
        async def scenario():
            table,receptor,pose,pdf=self.docking('dock6');ui=self.ui()
            with patch('biomolexplorer.ui.footprint.pdf_preview',side_effect=ValueError('unavailable')):
                await footprint_dialog(ui,self.pid,self.rid,self.sid,self.selection(table),'M1',lambda:True)
            dialog=ui.page.dialogs[-1]
            self.assertFalse(dialog.content.content.controls[2].visible)
            await dialog.actions[0].on_click(None)
            self.assertEqual(ui.picker.save_file.await_args.kwargs['src_bytes'],pdf.read_bytes())
        asyncio.run(scenario())

    def test_logout_revokes_pdf_download(self):
        async def scenario():
            table,receptor,pose,pdf=self.docking('dock6');ui=self.ui()
            with patch('biomolexplorer.ui.footprint.pdf_preview',return_value=PREVIEW):
                await footprint_dialog(ui,self.pid,self.rid,self.sid,self.selection(table),'M1',lambda:True)
            self.store.logout(self.token)
            with self.assertRaises(AccessDenied):await ui.page.dialogs[-1].actions[0].on_click(None)
            ui.picker.save_file.assert_not_awaited()
        asyncio.run(scenario())


class PDFPreviewTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which('pdftoppm'),'Poppler is required for PDF rendering')
    def test_existing_pdf_is_rendered_as_a_bounded_png(self):
        from matplotlib.figure import Figure
        from PIL import Image
        figure=Figure(figsize=(6,4));figure.subplots().plot([0,1,2],[-1,-4,-2])
        pdf=io.BytesIO();figure.savefig(pdf,format='pdf')
        image=pdf_preview(pdf.getvalue())
        with Image.open(io.BytesIO(image)) as preview:
            self.assertEqual(preview.format,'PNG')
            self.assertLessEqual(max(preview.size),2400)
            self.assertGreater(max(preview.size),1000)

    def test_invalid_document_is_rejected(self):
        with self.assertRaises(ValueError):pdf_preview(b'not a pdf')


if __name__=='__main__':unittest.main()
