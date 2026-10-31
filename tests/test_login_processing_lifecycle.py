"""Real authentication followed by workspace reads before dialog dismissal."""
import tempfile
import unittest
from types import SimpleNamespace
import flet as ft

from biomolexplorer.workspace import WorkspaceStore,AccessDenied
from biomolexplorer.ui.app import WorkspaceUI
from biomolexplorer.ui.configuration_loading import ProcessingSplash


from test_frontend_navigation import Page,walk


class RetainedDialogPage(Page):
    """Model Flet's managed stack, which retains closing routes until on_dismiss."""
    def __init__(self):super().__init__()
    def show_dialog(self,dialog):
        if dialog in self.dialogs:raise RuntimeError('Dialog is already opened')
        dialog.open=True;self.dialogs.append(dialog)
    def update(self):pass
    def localize(self,control):pass


class LoginLifecycleTests(unittest.IsolatedAsyncioTestCase):
    async def test_registered_credentials_survive_restart_and_consecutive_workspace_reads(self):
        with tempfile.TemporaryDirectory() as directory:
            store=WorkspaceStore(directory)
            initial=store.register('Test user','test@example.org','test-long-password')
            store.logout(initial)
            reopened=WorkspaceStore(directory)
            ui=SimpleNamespace(page=RetainedDialogPage())
            for attempt in range(2):
                token=await WorkspaceUI.call(ui,reopened.login,'TEST@example.org','test-long-password')
                self.assertEqual(await WorkspaceUI.call(ui,reopened.list_projects,token),[])
                await WorkspaceUI.call(ui,reopened.logout,token)
            dialogs=ui.page.dialogs
            self.assertGreaterEqual(len(dialogs),4)
            self.assertEqual(len({id(d) for d in dialogs}),len(dialogs))
            self.assertTrue(all(not d.open for d in dialogs))
            self.assertFalse(ui._configuration_loading.requests)
            with self.assertRaises(AccessDenied):reopened.login('test@example.org','wrong-password')

    async def test_concurrent_requests_share_open_dialog_but_new_cycle_uses_new_route(self):
        page=RetainedDialogPage();splash=ProcessingSplash(page)
        splash.show('first','First',modal=True)
        first=splash.dialog
        splash.show('second','Second',modal=True)
        self.assertIs(splash.dialog,first)
        self.assertEqual(len(page.dialogs),1)
        splash.hide('first');self.assertTrue(first.open)
        splash.hide('second');self.assertFalse(first.open)
        splash.show('third','Third',modal=True)
        self.assertIsNot(splash.dialog,first)
        self.assertTrue(splash.dialog.open)
        splash.hide('third')

    async def test_repeated_saving_uses_fresh_dialogs_before_client_dismissal(self):
        page=RetainedDialogPage();splash=ProcessingSplash(page,saving=True)
        for request in range(3):
            splash.show(request);splash.hide(request)
        self.assertEqual(len(page.dialogs),3)
        self.assertTrue(all(not dialog.open for dialog in page.dialogs))

    async def test_login_button_opens_workspace_with_real_persisted_credentials(self):
        with tempfile.TemporaryDirectory() as directory:
            store=WorkspaceStore(directory)
            registered=store.register('Test user','test@example.org','test-long-password')
            store.logout(registered)
            page=RetainedDialogPage()
            ui=WorkspaceUI(page,WorkspaceStore(directory),SimpleNamespace(),language='en',structure_viewer=SimpleNamespace())
            ui.show_login()
            controls=list(walk(page.controls))
            next(c for c in controls if isinstance(c,ft.TextField) and c.label=='E-mail').value='test@example.org'
            next(c for c in controls if isinstance(c,ft.TextField) and c.label=='Password').value='test-long-password'
            button=next(c for c in controls if isinstance(c,ft.Button) and c.content=='Sign in to workspace')
            await button.on_click(None)
            self.assertIsNotNone(ui.token)
            self.assertEqual(store.user(ui.token)['email'],'test@example.org')
            self.assertTrue(any(isinstance(c,ft.Text) and c.value=='Your workspace' for c in walk(page.controls)))
            self.assertTrue(all(not dialog.open for dialog in page.dialogs))
