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


class ConfigurationSavingTests(unittest.TestCase):
    def setUp(self):
        fixtures.PipelineEditorRegressionTests.setUp(self)
        self.stage['parameters']['pdb_query']='ACHE'

    def test_apply_publishes_modal_splash_before_reads_and_saves_only_once(self):
        from types import SimpleNamespace
        async def scenario():
            await self.ui.open_stage_dialog(self.stage['id'])
            dialog=self.page.dialogs[-1];button=dialog.actions[-1]
            entered=asyncio.Event();release=asyncio.Event();writes=[]
            async def call(function,*args,**kwargs):
                if function==self.store.project and len(args)>2:
                    # This is the first permission read in Apply, before form validation.
                    self.assertTrue(button.disabled)
                    self.assertTrue(all(c.disabled for c in dialog.actions))
                    splash=self.page.dialogs[-1]
                    self.assertIsNot(splash,dialog);self.assertTrue(splash.open);self.assertTrue(splash.modal)
                    entered.set();await release.wait()
                if function==self.store.save_pipeline:writes.append(args)
                return function(*args,**kwargs)
            self.ui.call=call
            click=SimpleNamespace(control=button)
            pending=asyncio.create_task(button.on_click(click))
            await asyncio.wait_for(entered.wait(),1)
            splash=self.page.dialogs[-1]
            await button.on_click(click)
            self.assertEqual(writes,[])
            release.set();await pending
            await button.on_click(click)  # An already queued click after completion.
            self.assertEqual(len(writes),1)
            self.assertFalse(splash.open);self.assertFalse(dialog.open)
        asyncio.run(scenario())

    def test_failed_apply_hides_only_splash_and_restores_buttons_for_retry(self):
        from types import SimpleNamespace
        async def scenario():
            await self.ui.open_stage_dialog(self.stage['id'])
            dialog=self.page.dialogs[-1];button=dialog.actions[-1];attempts=0
            async def call(function,*args,**kwargs):
                nonlocal attempts
                if function==self.store.project and len(args)>2:
                    attempts+=1
                    if attempts==1:raise ValueError('test saving failure')
                return function(*args,**kwargs)
            self.ui.call=call
            before=len(self.page.dialogs)
            await button.on_click(SimpleNamespace(control=button))
            self.assertTrue(dialog.open)
            self.assertTrue(all(not c.disabled for c in dialog.actions))
            self.assertFalse(self.page.dialogs[before].open)
            self.assertTrue(self.page.dialogs[-1].open)  # Keep the validation/error notification.
            await button.on_click(SimpleNamespace(control=button))
            self.assertEqual(attempts,2);self.assertFalse(dialog.open)
        asyncio.run(scenario())

    def test_cancelled_handler_releases_splash_and_preserves_disabled_state(self):
        from types import SimpleNamespace
        async def scenario():
            dialog=ft.AlertDialog(actions=[ft.TextButton('Cancel',disabled=True),ft.TextButton('Apply')])
            self.ui.page.show_dialog(dialog)
            entered=asyncio.Event();release=asyncio.Event()
            async def save(e):entered.set();await release.wait()
            action=self.ui.event_action(save,configuration_dialog=lambda:dialog)
            pending=asyncio.create_task(action(SimpleNamespace(control=dialog.actions[-1])))
            await asyncio.wait_for(entered.wait(),1);splash=self.page.dialogs[-1]
            pending.cancel()
            with self.assertRaises(asyncio.CancelledError):await pending
            self.assertFalse(splash.open);self.assertTrue(dialog.open)
            self.assertTrue(dialog.actions[0].disabled);self.assertFalse(dialog.actions[1].disabled)
        asyncio.run(scenario())

    def test_save_failure_keeps_form_open_and_allows_a_new_attempt(self):
        from types import SimpleNamespace
        async def scenario():
            await self.ui.open_stage_dialog(self.stage['id'])
            dialog=self.page.dialogs[-1];button=dialog.actions[-1]
            entered=asyncio.Event();release=asyncio.Event();attempts=0
            async def call(function,*args,**kwargs):
                nonlocal attempts
                if function==self.store.save_pipeline:
                    attempts+=1
                    if attempts==1:
                        entered.set();await release.wait()
                        raise OSError('test write failure')
                return function(*args,**kwargs)
            self.ui.call=call
            pending=asyncio.create_task(button.on_click(SimpleNamespace(control=button)))
            await asyncio.wait_for(entered.wait(),1)
            splash=self.page.dialogs[-1]
            self.assertTrue(splash.open);self.assertTrue(button.disabled)
            release.set();await pending
            self.assertFalse(splash.open);self.assertTrue(dialog.open)
            self.assertFalse(button.disabled)
            await button.on_click(SimpleNamespace(control=button))
            self.assertEqual(attempts,2);self.assertFalse(dialog.open)
        asyncio.run(scenario())

    def test_saving_splash_follows_session_language(self):
        from types import SimpleNamespace
        async def scenario(language):
            self.ui.page.set_language(language)
            dialog=ft.AlertDialog(actions=[ft.TextButton('Apply')]);self.ui.page.show_dialog(dialog)
            async def save(e):
                values=[c.value for c in walk(self.page.dialogs[-1]) if isinstance(c,ft.Text)]
                self.assertIn('Saving settings…' if language=='en' else 'Salvando configurações…',values)
                self.assertTrue(self.page.dialogs[-1].open)
            await self.ui.event_action(save,configuration_dialog=lambda:dialog)(SimpleNamespace(control=dialog.actions[0]))
            self.assertFalse(self.page.dialogs[-1].open)
        for language in ('en','pt'):asyncio.run(scenario(language))


if __name__=='__main__':unittest.main()
