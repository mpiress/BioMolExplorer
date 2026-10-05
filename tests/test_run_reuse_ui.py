"""Existing results require an explicit reuse decision before submitting a run."""
import asyncio
import unittest
from types import SimpleNamespace
from unittest.mock import AsyncMock

try:
    from biomolexplorer.ui.app import WorkspaceUI
except ImportError:
    WorkspaceUI=None


@unittest.skipIf(WorkspaceUI is None,'Install the UI extra')
class ReusePromptTests(unittest.TestCase):
    def setUp(self):
        self.dialogs=[]
        self.ui=object.__new__(WorkspaceUI)
        ui=self.ui
        ui.token='token';ui.current={'id':'project'};ui.selected='stage';ui.submitting=False
        ui.save_draft=AsyncMock();ui.force_execution=False
        def show(dialog):dialog.open=True;self.dialogs.append(dialog)
        ui.page=SimpleNamespace(show_dialog=show,update=lambda:None,pop_dialog=lambda:None)
        ui.service=SimpleNamespace(existing_results=lambda *args:2,submit=AsyncMock())
        async def call(function,*args,**kwargs):return function(*args,**kwargs)
        ui.call=call

    def test_submission_waits_for_choice_and_both_answers_reach_the_backend(self):
        for reuse,action in ((True,2),(False,1)):
            with self.subTest(reuse=reuse):
                asyncio.run(self.ui.run_pipeline(selected=True))
                self.ui.service.submit.assert_not_called()
                dialog=self.dialogs[-1]
                submit=AsyncMock()
                self.ui._submit_pipeline=submit
                asyncio.run(dialog.actions[action].on_click(None))
                submit.assert_awaited_once_with(True,reuse,'stage')
                self.assertFalse(dialog.open)
                del self.ui._submit_pipeline

    def test_cancel_or_session_switch_never_starts_execution(self):
        asyncio.run(self.ui.run_pipeline())
        dialog=self.dialogs[-1]
        dialog.actions[0].on_click(None)
        self.ui.service.submit.assert_not_called()
        self.ui.token='new-session'
        self.ui._submit_pipeline=AsyncMock()
        asyncio.run(dialog.actions[-1].on_click(None))
        self.ui._submit_pipeline.assert_not_called()

    def test_repeated_execute_clicks_keep_only_one_prompt(self):
        asyncio.run(self.ui.run_pipeline())
        asyncio.run(self.ui.run_pipeline())
        self.assertEqual(len(self.dialogs),1)
        self.ui.save_draft.assert_awaited_once()
        self.ui.service.submit.assert_not_called()

    def test_cancel_and_dismiss_release_prompt_and_ignore_old_responses(self):
        for close in ('cancel','dismiss'):
            with self.subTest(close=close):
                asyncio.run(self.ui.run_pipeline())
                previous=self.dialogs[-1]
                if close=='cancel':previous.actions[0].on_click(None)
                else:
                    previous.open=False
                    previous.on_dismiss(None)
                self.assertIsNone(self.ui.reuse_dialog)
                self.assertIsNone(self.ui.reuse_dialog_key)
                asyncio.run(self.ui.run_pipeline())
                current=self.dialogs[-1]
                self.assertIsNot(current,previous)
                asyncio.run(previous.actions[-1].on_click(None))
                self.assertIs(self.ui.reuse_dialog,current)
                self.assertTrue(current.open)
                self.ui.service.submit.assert_not_called()
                current.actions[0].on_click(None)

    def test_leaving_project_invalidates_prompt(self):
        asyncio.run(self.ui.run_pipeline())
        dialog=self.dialogs[-1]
        self.ui.close_project_views()
        self.assertFalse(dialog.open)
        self.assertIsNone(self.ui.reuse_dialog)
        self.assertIsNone(self.ui.reuse_dialog_key)
        asyncio.run(dialog.actions[-1].on_click(None))
        self.ui.service.submit.assert_not_called()

    def test_opening_other_project_closes_prompt(self):
        asyncio.run(self.ui.run_pipeline())
        dialog=self.dialogs[-1]
        self.ui.store=SimpleNamespace(project=lambda token,project_id:{'id':project_id,'pipeline':[]})
        self.ui.draw_project=AsyncMock()
        self.ui.page.run_task=lambda *args:None
        asyncio.run(self.ui.open_project('other-project'))
        self.assertFalse(dialog.open)
        self.assertIsNone(self.ui.reuse_dialog)
        asyncio.run(dialog.actions[-1].on_click(None))
        self.ui.service.submit.assert_not_called()

    def test_response_after_session_or_project_switch_is_ignored(self):
        for switch in ('session','project'):
            with self.subTest(switch=switch):
                self.ui.token='token';self.ui.current={'id':'project'}
                asyncio.run(self.ui.run_pipeline())
                dialog=self.dialogs[-1]
                if switch=='session':self.ui.token='new-session'
                else:self.ui.current={'id':'other-project'}
                asyncio.run(dialog.actions[-1].on_click(None))
                self.assertFalse(dialog.open)
                self.assertIsNone(self.ui.reuse_dialog)
                self.ui.service.submit.assert_not_called()

    def test_late_result_lookup_cannot_open_prompt_after_leaving(self):
        for switch in ('session','project'):
            with self.subTest(switch=switch):
                self.ui.token='token';self.ui.current={'id':'project'}
                async def scenario():
                    started=asyncio.Event();finish=asyncio.Event()
                    async def lookup(function,*args,**kwargs):
                        started.set()
                        await finish.wait()
                        return 2
                    self.ui.call=lookup
                    pending=asyncio.create_task(self.ui.run_pipeline())
                    await started.wait()
                    if switch=='session':self.ui.token='new-session'
                    else:self.ui.current={'id':'other-project'}
                    finish.set()
                    await pending
                asyncio.run(scenario())
                self.assertEqual(self.dialogs,[])
                self.assertFalse(self.ui.submitting)
                self.ui.service.submit.assert_not_called()
