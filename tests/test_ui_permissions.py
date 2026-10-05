"""Real authorization changes exercised through the Flet controller."""
import asyncio
import json
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from biomolexplorer.catalog import new_stage
from biomolexplorer.workspace import WorkspaceStore

try:
    from biomolexplorer.ui.app import WorkspaceUI
except ImportError:
    WorkspaceUI = None


class FakePage:
    def __init__(self):
        self.width,self.height=1440,1000
        self.controls=[]
        self.dialogs=[SimpleNamespace(open=True),SimpleNamespace(open=True)]
        self.closed=0
        self.messages=[]
        self.on_keyboard_event='old-editor-handler'

    def update(self): pass
    def run_task(self,*args):return SimpleNamespace(cancel=lambda:None,done=lambda:False)
    def add(self,control): self.controls.append(control)
    def pop_dialog(self):
        if not self.dialogs:return None
        self.closed+=1
        return self.dialogs.pop()
    def show_dialog(self,dialog):
        self.dialogs.append(dialog)
        if hasattr(dialog.content,'value'):self.messages.append(dialog.content.value)


@unittest.skipIf(WorkspaceUI is None,'Install the ui extra to test Flet controls')
class UIPermissionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='biomol UI permissions ')
        self.store=WorkspaceStore(self.temp.name)
        self.owner=self.store.register('Owner','owner@example.org','owner-password')
        self.guest=self.store.register('Guest','guest@example.org','guest-password')
        self.project=self.store.create_project(self.owner,'Private project')
        self.pid=self.project['id']
        self.store.invite(self.owner,self.pid,'guest@example.org','editor')
        self.store.accept_invitation(self.guest,self.pid,True,'editor')
        self.ui=object.__new__(WorkspaceUI)
        self.ui.page=FakePage()
        self.ui.store=self.store
        # These controller tests exercise real storage without an event-loop
        # self-pipe, which the execution sandbox restricts. Browser checks cover
        # the production asyncio.to_thread calls on the actual Flet server.
        async def local_call(function,*args,**kwargs):return function(*args,**kwargs)
        self.ui.call=local_call
        self.ui.service=SimpleNamespace(dock6_path=None)
        self.ui.token=self.guest
        self.ui.current=self.store.project(self.guest,self.pid)
        self.ui.tab='Pipeline'
        self.ui.dirty=True
        self.ui.base_pipeline=[]
        self.ui._save_lock=asyncio.Lock()
        self.ui.selected='local-node'
        self.ui.current_run='old-run'
        self.ui.polling=None
        self.viewer=SimpleNamespace(selection_version=0)
        self.ui.result_viewer=self.viewer
        self.ui.editing_stage=True

    def tearDown(self): self.temp.cleanup()

    def test_revoked_project_closes_private_views_and_returns_to_workspace(self):
        self.store.revoke(self.owner,self.pid,self.store.user(self.guest)['id'])
        async def action():
            await self.ui.call(self.store.read_file,self.guest,self.pid,self.store.project_dir(self.pid)/'result.csv')
        asyncio.run(self.ui.guard(action))
        self.assertIsNone(self.ui.current)
        self.assertEqual(self.ui.token,self.guest)
        self.assertIsNone(self.ui.result_viewer)
        self.assertGreater(self.viewer.selection_version,0)
        self.assertIsNone(self.ui.current_run)
        self.assertFalse(self.ui.dirty)
        self.assertEqual(self.ui.page.closed,2)
        self.assertIn('acesso não autorizado',self.ui.page.messages[-1])

    def test_expired_session_closes_dialogs_before_showing_login(self):
        with self.store.connect() as db:
            db.execute('UPDATE sessions SET expires=0 WHERE user_id=?',(self.store.user(self.guest)['id'],))
        async def action():await self.ui.call(self.store.project,self.guest,self.pid)
        asyncio.run(self.ui.guard(action))
        self.assertIsNone(self.ui.token)
        self.assertIsNone(self.ui.current)
        self.assertIsNone(self.ui.result_viewer)
        self.assertEqual(self.ui.page.closed,2)
        self.assertFalse(self.ui.editing_stage)

    def test_downgrade_discards_editor_draft_and_rebuilds_reader_controls(self):
        self.ui.current['pipeline']=[new_stage('admet')]
        self.store.invite(self.owner,self.pid,'guest@example.org','viewer')
        self.store.accept_invitation(self.guest,self.pid,True,'viewer')
        async def action():await self.ui.call(self.store.save_pipeline,self.guest,self.pid,[],0)
        asyncio.run(self.ui.guard(action))
        self.assertEqual(self.ui.current['role'],'viewer')
        self.assertEqual(self.ui.current['pipeline'],[])
        self.assertFalse(self.ui.dirty)
        self.assertFalse(self.ui.flow_editor.writable)
        self.assertIsNone(self.ui.result_viewer)
        self.assertEqual(self.ui.page.closed,2)

    def test_changed_invitation_refreshes_workspace_without_accepting_wrong_role(self):
        self.store.revoke(self.owner,self.pid,self.store.user(self.guest)['id'])
        self.ui.current=None
        self.store.invite(self.owner,self.pid,'guest@example.org','viewer')
        self.store.invite(self.owner,self.pid,'guest@example.org','editor')
        asyncio.run(self.ui.guard(lambda:self.ui.accept_invite(self.pid,True,'viewer')))
        self.assertEqual(self.store.invitations(self.guest)[0]['role'],'editor')
        self.assertEqual(self.store.list_projects(self.guest),[])
        self.assertIn('permissão do convite mudou',self.ui.page.messages[-1])

    def test_polling_does_not_query_a_run_after_logout_during_wait(self):
        async def logout_during_wait(delay):
            self.ui.token=None
            self.ui.current_run=None
        with patch('biomolexplorer.ui.app.asyncio.sleep',side_effect=logout_during_wait), \
             patch.object(self.store,'get_run') as get_run:
            asyncio.run(self.ui.poll_runs())
            get_run.assert_not_called()

    def test_project_response_from_previous_session_is_discarded(self):
        async def switched_session(function,*args):
            result=function(*args)
            self.ui.token='new-session'
            self.ui.current=None
            return result
        self.ui.call=switched_session
        asyncio.run(self.ui.open_project(self.pid))
        self.assertIsNone(self.ui.current)
        self.assertEqual(self.ui.page.controls,[])

    def test_private_visualization_response_is_not_opened_after_session_change(self):
        path=self.store.project_dir(self.pid)/'private.biomol-view.json'
        path.write_text(json.dumps({'version':1,'kind':'graph','title':'Private molecule graph',
                                    'nodes':[],'edges':[],'mcc':[]}))
        async def switched_session(function,*args):
            result=function(*args)
            self.ui.token='new-session'
            self.ui.current=None
            return result
        self.ui.call=switched_session
        asyncio.run(self.ui.preview_artifact(self.pid,path))
        self.assertEqual(len(self.ui.page.dialogs),2)
        self.assertIs(self.ui.result_viewer,self.viewer)

    def test_workspace_and_save_results_are_discarded_after_logout(self):
        original_call=self.ui.call
        async def logout_after_response(function,*args,**kwargs):
            result=await original_call(function,*args,**kwargs)
            self.ui.token=None
            self.ui.current=None
            return result
        self.ui.call=logout_after_response
        asyncio.run(self.ui.save_draft())
        self.assertIsNone(self.ui.current)
        self.assertEqual(self.ui.page.messages,[])
        self.ui.token=self.guest
        asyncio.run(self.ui.show_workspace())
        self.assertEqual(self.ui.page.controls,[])

    def test_refresh_updates_collaborator_changes_when_no_local_draft_exists(self):
        self.ui.dirty=False
        self.ui.editing_stage=False
        stage=new_stage('admet')
        self.store.save_pipeline(self.owner,self.pid,[stage],0)
        self.store.update_project(self.owner,self.pid,'Updated name','New description','#6366F1',['shared'])
        asyncio.run(self.ui.draw_project())
        self.assertEqual(self.ui.current['revision'],1)
        self.assertEqual(self.ui.current['name'],'Updated name')
        self.assertEqual(self.ui.current['pipeline'][0]['id'],stage['id'])
