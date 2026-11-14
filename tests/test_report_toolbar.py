"""Report access and existing export formats from the canvas toolbar."""
import asyncio
import io
import json
import unittest
import zipfile
from types import SimpleNamespace
from unittest.mock import AsyncMock, patch

import flet as ft
from biomolexplorer.ui.flow_canvas import FlowCanvas
from biomolexplorer.experimental_report import build_report
from test_frontend_navigation import walk
import test_pipeline_editor_regressions as fixtures


class ReportToolbarTests(unittest.TestCase):
    setUp=fixtures.PipelineEditorRegressionTests.setUp

    def test_report_and_screenshot_share_wrapping_toolbar_in_both_editor_modes(self):
        for writable in (True,False):
            for library in (True,False):
                editor=FlowCanvas(self.ui,writable)
                toolbar=editor.build(include_library=library).controls[0]
                self.assertTrue(toolbar.wrap)
                camera=next(i for i,c in enumerate(toolbar.controls) if isinstance(c,ft.PopupMenuButton))
                report=toolbar.controls[camera+1].content
                self.assertEqual(report.icon,ft.Icons.SCIENCE_OUTLINED)
                self.assertFalse(report.disabled)
                self.assertEqual(report.on_click,editor.export_experimental_report)

    def test_toolbar_opens_persisted_report_and_exports_existing_formats(self):
        async def scenario():
            editor=FlowCanvas(self.ui,True)
            self.store.save_pipeline(self.token,self.pid,self.ui.current['pipeline'],self.ui.current['revision'])
            pipeline_before=self.store.project(self.token,self.pid)['pipeline']
            button=next(c for c in walk(editor.build(False)) if isinstance(c,ft.IconButton) and c.icon==ft.Icons.SCIENCE_OUTLINED)
            await button.on_click(None)
            dialog=self.page.dialogs[-1]
            items=[c for c in walk(dialog) if isinstance(c,ft.PopupMenuItem)]
            self.assertEqual({c.data for c in items},{'pdf','docx','json'})
            calls=[];original_call=self.ui.call
            async def counted(function,*args,**kwargs):
                calls.append(function)
                return await original_call(function,*args,**kwargs)
            self.ui.call=counted
            with patch('biomolexplorer.ui.flow_canvas.save_file',new_callable=AsyncMock) as save:
                for item in items:
                    await item.on_click(SimpleNamespace(control=item))
                    data=save.call_args.kwargs['src_bytes']
                    self.assertEqual(save.call_args.kwargs['file_name'],'experimental-report.'+item.data)
                    if item.data=='pdf':self.assertTrue(data.startswith(b'%PDF-'))
                    elif item.data=='docx':self.assertIn('word/document.xml',zipfile.ZipFile(io.BytesIO(data)).namelist())
                    else:self.assertEqual(json.loads(data)['schema'],'biomolexplorer.experimental-report')
                self.assertEqual(calls.count(build_report),1)
            self.assertEqual(self.store.project(self.token,self.pid)['pipeline'],pipeline_before)
            self.assertFalse(self.ui._experimental_export_busy)
            printer=next(c for c in walk(dialog) if isinstance(c,ft.Button) and c.icon==ft.Icons.PRINT)
            for web in (False,True):
                self.page.web=web
                with patch('flet.UrlLauncher.launch_url',new_callable=AsyncMock) as launch:
                    await printer.on_click(None)
                    url=launch.call_args.args[0]
                    self.assertTrue(url.startswith('/experimental-reports/' if web else 'file:'))
                    preview=next((self.store.staging/'report-previews').glob('*.html'))
                    self.assertIn('window.print()',preview.read_text())
        asyncio.run(scenario())

    def test_preview_does_not_read_results_and_old_project_cannot_export(self):
        async def scenario():
            editor=FlowCanvas(self.ui,True)
            with patch('biomolexplorer.experimental_report.build_report',side_effect=AssertionError('Preview must not build documents')) as build:
                await editor.export_experimental_report(None)
                build.assert_not_called()
                dialog=self.page.dialogs[-1]
                self.ui.current=self.store.create_project(self.token,'Other project')
                item=next(c for c in walk(dialog) if isinstance(c,ft.PopupMenuItem))
                with patch('biomolexplorer.ui.flow_canvas.save_file',new_callable=AsyncMock) as save:
                    await item.on_click(SimpleNamespace(control=item))
                    save.assert_not_called();build.assert_not_called()
            self.assertFalse(self.ui._experimental_export_busy)
        asyncio.run(scenario())


if __name__=='__main__':unittest.main()
