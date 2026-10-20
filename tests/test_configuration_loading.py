"""Visible feedback survives concurrent reads and always releases its overlay."""
import asyncio
import unittest

import flet as ft

from biomolexplorer.catalog import new_stage
from biomolexplorer.ui.branding import image_bytes
import test_pipeline_editor_regressions as fixtures
from test_frontend_navigation import walk


class ConfigurationLoadingTests(unittest.TestCase):
    setUp = fixtures.PipelineEditorRegressionTests.setUp

    def test_branded_spinner_is_published_before_reads_and_removed_on_success(self):
        async def scenario():
            entered=asyncio.Event();release=asyncio.Event()
            async def call(function,*args,**kwargs):
                if function==self.store.project:
                    entered.set();await release.wait()
                return function(*args,**kwargs)
            self.ui.call=call
            pending=asyncio.create_task(self.ui.open_stage_dialog(self.stage['id']))
            await asyncio.wait_for(entered.wait(),1)
            self.assertEqual(len(self.page.overlay),1)
            overlay=self.page.overlay[0]
            self.assertEqual((overlay.left,overlay.top,overlay.right,overlay.bottom),(0,0,0,0))
            self.assertEqual(overlay.alignment,ft.Alignment.CENTER)
            controls=list(walk(overlay))
            spinner=next(c for c in controls if isinstance(c,ft.ProgressRing))
            self.assertIsNone(spinner.value)
            emblem=next(c for c in controls if isinstance(c,ft.Image))
            self.assertEqual(emblem.src,image_bytes('logo.png'))
            release.set();await pending
            self.assertEqual(self.page.overlay,[])
            self.assertTrue(self.page.dialogs[-1].open)
        asyncio.run(scenario())

    def test_exception_and_cancellation_remove_feedback_and_allow_retry(self):
        async def scenario():
            async def fail(*args,**kwargs):raise RuntimeError('Failed to read')
            self.ui.call=fail
            with self.assertRaisesRegex(RuntimeError,'Failed to read'):
                await self.ui.open_stage_dialog(self.stage['id'])
            self.assertEqual(self.page.overlay,[])
            entered=asyncio.Event();release=asyncio.Event()
            async def wait(*args,**kwargs):
                entered.set();await release.wait()
            self.ui.call=wait
            pending=asyncio.create_task(self.ui.open_stage_dialog(self.stage['id']))
            await asyncio.wait_for(entered.wait(),1)
            pending.cancel()
            with self.assertRaises(asyncio.CancelledError):await pending
            self.assertEqual(self.page.overlay,[])
            self.assertEqual(self.ui._configuration_loading.requests,set())
            async def call(function,*args,**kwargs):return function(*args,**kwargs)
            self.ui.call=call
            await self.ui.open_stage_dialog(self.stage['id'])
            self.assertTrue(self.page.dialogs[-1].open)
            self.assertEqual(self.page.overlay,[])
        asyncio.run(scenario())

    def test_old_request_cannot_hide_indicator_while_new_request_is_loading(self):
        async def scenario():
            second=new_stage('docking_vina');self.ui.current['pipeline'].append(second)
            first_entered=asyncio.Event();first_release=asyncio.Event()
            second_entered=asyncio.Event();second_release=asyncio.Event();calls=0
            async def call(function,*args,**kwargs):
                nonlocal calls
                if function==self.store.project:
                    calls+=1
                    if calls==1:first_entered.set();await first_release.wait()
                    else:second_entered.set();await second_release.wait()
                return function(*args,**kwargs)
            self.ui.call=call
            old=asyncio.create_task(self.ui.open_stage_dialog(self.stage['id']))
            await asyncio.wait_for(first_entered.wait(),1)
            indicator=self.page.overlay[0]
            latest=asyncio.create_task(self.ui.open_stage_dialog(second['id']))
            await asyncio.sleep(0)
            first_release.set();await old
            await asyncio.wait_for(second_entered.wait(),1)
            self.assertEqual(len(self.page.overlay),1)
            self.assertIs(self.page.overlay[0],indicator)
            self.assertEqual(self.page.dialogs,[])
            second_release.set();await latest
            self.assertEqual(self.page.overlay,[])
            self.assertEqual(len(self.page.dialogs),1)
        asyncio.run(scenario())

    def test_leaving_project_clears_only_our_overlay_and_invalidates_pending_form(self):
        async def scenario():
            other=ft.Text('Existing overlay');self.page.overlay=[other]
            entered=asyncio.Event();release=asyncio.Event()
            async def call(function,*args,**kwargs):
                entered.set();await release.wait();return function(*args,**kwargs)
            self.ui.call=call
            pending=asyncio.create_task(self.ui.open_stage_dialog(self.stage['id']))
            await asyncio.wait_for(entered.wait(),1)
            self.assertEqual(len(self.page.overlay),2)
            self.ui.close_project_views()
            self.assertEqual(len(self.page.overlay),1)
            self.assertIs(self.page.overlay[0],other)
            release.set();await pending
            self.assertEqual(self.page.dialogs,[])
            self.assertEqual(len(self.page.overlay),1)
        asyncio.run(scenario())

    def test_loading_feedback_is_translated(self):
        async def scenario():
            self.ui.page.set_language('en')
            entered=asyncio.Event();release=asyncio.Event()
            async def call(function,*args,**kwargs):
                entered.set();await release.wait();return function(*args,**kwargs)
            self.ui.call=call
            pending=asyncio.create_task(self.ui.open_stage_dialog(self.stage['id']))
            await asyncio.wait_for(entered.wait(),1)
            values=[c.value for c in walk(self.page.overlay[0]) if isinstance(c,ft.Text)]
            self.assertIn('Loading settings…',values)
            self.assertIn('Reading block information. Please wait a moment.',values)
            release.set();await pending
        asyncio.run(scenario())


if __name__=='__main__':unittest.main()
