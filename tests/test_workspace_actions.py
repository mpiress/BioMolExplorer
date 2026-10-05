"""Workspace deletion and clear separation of current runs from historical errors."""
import asyncio
import json
import tempfile
import unittest
from types import SimpleNamespace
from uuid import uuid4

from biomolexplorer.workspace import WorkspaceStore, AccessDenied

try:
    import flet as ft
    from biomolexplorer.ui.app import WorkspaceUI
    from biomolexplorer.ui.feedback import readable_error
except ImportError:
    WorkspaceUI = None


def controls(root):
    yield root
    for field in ('content', 'controls', 'actions', 'title', 'subtitle'):
        value = getattr(root, field, None)
        for child in value if isinstance(value, list) else [value]:
            if child is not None and not isinstance(child, (str, int, float, bool)):
                yield from controls(child)


class Page:
    def __init__(self):
        self.controls, self.dialogs = [], []
        self.width, self.height = 1440, 1000
    def update(self): pass
    def add(self, control): self.controls.append(control)
    def show_dialog(self, dialog): self.dialogs.append(dialog)
    def pop_dialog(self): return self.dialogs.pop() if self.dialogs else None


@unittest.skipIf(WorkspaceUI is None, 'Install the ui extra for Flet controls')
class WorkspaceActionsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.store = WorkspaceStore(self.temp.name)
        self.owner = self.store.register('Owner', 'owner@example.org', 'owner-password')
        self.guest = self.store.register('Guest', 'guest@example.org', 'guest-password')
        self.project = self.store.create_project(self.owner, 'Private project')
        self.pid = self.project['id']
        self.ui = object.__new__(WorkspaceUI)
        self.ui.store, self.ui.page = self.store, Page()
        self.ui.token, self.ui.current = self.owner, self.project
        self.ui.tab, self.ui.selected = 'Execuções', None
        self.ui.current_run, self.ui.run_progress = None, None
        self.ui.polling = SimpleNamespace(done=lambda: False)
        self.ui.dirty = False
        async def local_call(function, *args, **kwargs): return function(*args, **kwargs)
        self.ui.call = local_call

    def tearDown(self): self.temp.cleanup()

    def run_record(self, status, error=None, created=100):
        identifier = uuid4().hex
        stages = [{'id': uuid4().hex, 'name': 'Retrieve', 'operation': 'retrieve_compounds',
                   'status': status, 'configuration': {}, 'artifacts': [], 'error': error}]
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',
                (identifier, self.pid, self.store.user(self.owner)['id'], status,
                 json.dumps(stages), created, created, error))
        return identifier

    def test_workspace_has_delete_action_only_for_owned_projects(self):
        self.store.invite(self.owner, self.pid, 'guest@example.org', 'viewer')
        self.store.accept_invitation(self.guest, self.pid)
        self.store.create_project(self.guest, 'Guest project')
        self.ui.token = self.guest
        asyncio.run(self.ui.show_workspace())
        delete_buttons = [c for c in controls(self.ui.page.controls[0]) if isinstance(c, ft.TextButton) and c.content == 'Excluir projeto']
        self.assertEqual(len(delete_buttons), 1)
        with self.assertRaises(AccessDenied): asyncio.run(self.ui.delete_project_dialog(self.pid))

    def test_confirm_delete_removes_project_and_collaborator_access(self):
        self.store.invite(self.owner, self.pid, 'guest@example.org', 'editor')
        self.store.accept_invitation(self.guest, self.pid)
        asyncio.run(self.ui.delete_project_dialog(self.pid))
        dialog = self.ui.page.dialogs[-1]
        self.assertEqual(self.store.list_projects(self.owner)[0]['id'], self.pid)
        asyncio.run(dialog.actions[-1].on_click(None))
        self.assertEqual(self.store.list_projects(self.owner), [])
        self.assertEqual(self.store.list_projects(self.guest), [])
        with self.assertRaises(AccessDenied): self.store.project(self.guest, self.pid)

    def test_active_execution_prevents_delete_in_dialog_and_backend(self):
        self.run_record('running')
        asyncio.run(self.ui.delete_project_dialog(self.pid))
        self.assertTrue(self.ui.page.dialogs[-1].actions[-1].disabled)
        with self.assertRaisesRegex(ValueError, 'Cancele a execução'): self.store.delete_project(self.owner, self.pid)
        self.assertEqual(len(self.store.list_projects(self.owner)), 1)

    def test_old_errors_are_inside_collapsed_history_while_current_run_is_running(self):
        self.run_record('failed', 'OLD ERROR <html>server page</html>', 100)
        current_id = self.run_record('running', None, 200)
        from biomolexplorer.ui.run_progress import RunProgress
        self.ui.run_progress = RunProgress(self.ui, self.store.get_run(self.owner, current_id))
        view = asyncio.run(self.ui.runs_view(True))
        history = next(c for c in view.controls if isinstance(c, ft.ExpansionTile))
        self.assertEqual(history.title.value, 'Execuções anteriores (1)')
        self.assertFalse(history.expanded)
        self.assertNotIn('OLD ERROR', str(view.controls[1]))
        self.assertTrue(any(isinstance(c, ft.Text) and 'OLD ERROR' in (c.value or '') for c in controls(history)))
        self.assertEqual(self.ui.run_progress.run['id'], current_id)

    def test_error_page_is_not_rendered_as_a_giant_interface_message(self):
        error = 'Exception: Error getting schema with status 500 and msg <!doctype html><html>' + 'body' * 2000
        message = readable_error(error)
        self.assertIn('status 500', message)
        self.assertNotIn('<html>', message)
        self.assertLess(len(message), 500)

    def test_delete_response_cannot_remove_views_of_a_new_session(self):
        async def switched_session(function, *args, **kwargs):
            result = function(*args, **kwargs)
            self.ui.token = 'new-session'
            return result
        self.ui.call = switched_session
        asyncio.run(self.ui.delete_project_dialog(self.pid))
        self.assertEqual(self.ui.page.dialogs, [])
