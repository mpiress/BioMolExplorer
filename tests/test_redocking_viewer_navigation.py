"""Open molecular viewers in the browser gesture, retaining scene authorization."""
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch,AsyncMock
from urllib.parse import urlsplit

import test_docking_scene as scene_fixtures
from biomolexplorer.pdb_view import PREFIX
from biomolexplorer.redocking_results import RedockingResults
from biomolexplorer.workspace import AccessDenied


class ViewerNavigationTests(unittest.TestCase):
    def setUp(self):
        self.f=scene_fixtures.SceneCapabilityTests();self.f.setUp();self.addCleanup(self.f.doCleanups)
        self.pdb,self.pose=self.f.redocking()
        self.key=RedockingResults(self.f.store).simulations(self.f.token,self.f.pid,self.f.rid,self.f.sid)[0]['id']

    def reserve(self,**kwargs):
        return self.f.viewers.reserve_result(self.f.token,self.f.pid,self.f.rid,self.f.sid,'redocking',self.key,**kwargs)

    def ui(self,web=True):
        from biomolexplorer.ui.app import WorkspaceUI
        ui=object.__new__(WorkspaceUI);ui.store=self.f.store;ui.token=self.f.token
        ui.current={'id':self.f.pid};ui.structure_viewer=self.f.viewers;ui.language='pt'
        ui.page=SimpleNamespace(web=web,url='ws://localhost:9000/ws')
        return ui

    def test_reservation_and_html_do_not_load_result_files(self):
        with patch.object(self.f.viewers,'_result_spec',side_effect=AssertionError('Scene loaded during navigation')),patch.object(self.f.store,'read_file',side_effect=AssertionError('File loaded during navigation')):
            ticket=self.reserve(name='1ABC / LIG')
            status,mime,body,_=self.f.viewers.response(PREFIX+'/'+ticket)
            self.assertEqual(status,200);self.assertIn('text/html',mime)
            self.assertIn('1ABC / LIG',body.decode())
        status,_,body,_=self.f.viewers.response(PREFIX+'/'+ticket+'/scene')
        self.assertEqual(status,200);self.assertIn(b'coordinate_frame',body)

    def test_reservation_rejects_unrelated_user_and_invalid_kind(self):
        guest=self.f.store.register('Guest','guest@example.org','viewer-test-password')
        with self.assertRaises(AccessDenied):self.f.viewers.reserve_result(guest,self.f.pid,self.f.rid,self.f.sid,'redocking',self.key)
        with self.assertRaises(AccessDenied):self.f.viewers.reserve_result(self.f.token,self.f.pid,self.f.rid,self.f.sid,'other',self.key)

    def test_lazy_scene_revalidates_deleted_pose(self):
        ticket=self.reserve(language='pt');self.pose.unlink()
        self.assertEqual(self.f.viewers.response(PREFIX+'/'+ticket)[0],200)
        status,_,body,_=self.f.viewers.response(PREFIX+'/'+ticket+'/scene')
        self.assertEqual(status,403);self.assertIn('indisponível',body.decode())
        self.assertNotIn(str(self.pose),body.decode())

    def test_lazy_html_and_scene_revalidate_session(self):
        ticket=self.reserve();self.f.store.logout(self.f.token)
        for suffix in ('','/scene'):
            self.assertEqual(self.f.viewers.response(PREFIX+'/'+ticket+suffix)[0],403)

    def test_web_result_button_uses_direct_client_navigation(self):
        from biomolexplorer.ui.redocking_results import RedockingResultsTable
        ui=self.ui();r=SimpleNamespace(ui=ui,token=self.f.token,project_id=self.f.pid,run_id=self.f.rid,stage={'id':self.f.sid},valid=lambda:True,writable=False)
        simulations=RedockingResults(self.f.store).simulations(self.f.token,self.f.pid,self.f.rid,self.f.sid)
        with patch.object(self.f.viewers,'_result_spec',side_effect=AssertionError('Eager scene load')):
            table=RedockingResultsTable(r,simulations)
        button=table.table.rows[0].cells[-3].content
        self.assertIsNone(button.on_click)
        self.assertEqual(button.url.target.value,'_blank')
        self.assertTrue(button.url.url.startswith('http://localhost:9000'+PREFIX+'/'))
        self.assertNotIn(self.f.token,button.url.url)
        self.assertFalse(hasattr(ui,'_configuration_loading'))
        self.assertEqual(self.f.viewers.response(urlsplit(button.url.url).path+'/scene')[0],200)

    def test_native_result_button_keeps_callback(self):
        from biomolexplorer.ui.redocking_results import RedockingResultsTable
        ui=self.ui(web=False);r=SimpleNamespace(ui=ui,token=self.f.token,project_id=self.f.pid,run_id=self.f.rid,stage={'id':self.f.sid},valid=lambda:True,writable=False)
        table=RedockingResultsTable(r,RedockingResults(self.f.store).simulations(self.f.token,self.f.pid,self.f.rid,self.f.sid))
        button=table.table.rows[0].cells[-3].content
        self.assertIsNone(button.url);self.assertIsNotNone(button.on_click)

    def test_stale_project_does_not_reserve_a_link(self):
        ui=self.ui();ui.current={'id':'other'}
        with patch.object(self.f.viewers,'reserve_result') as reserve:
            self.assertIsNone(ui.docking_view_url(self.f.pid,self.f.rid,self.f.sid,'redocking',self.key))
            reserve.assert_not_called()


class ResidueNavigationTests(unittest.IsolatedAsyncioTestCase):
    async def test_contact_dialog_uses_the_same_direct_navigation(self):
        from biomolexplorer.ui.docking_scene import residue_dialog
        f=scene_fixtures.SceneCapabilityTests();f.setUp();self.addCleanup(f.doCleanups);f.redocking()
        service=RedockingResults(f.store);key=service.simulations(f.token,f.pid,f.rid,f.sid)[0]['id'];dialogs=[]
        page=SimpleNamespace(width=1000,height=800,on_resize=None,show_dialog=dialogs.append,update=lambda:None)
        async def call(fn,*args):return fn(*args)
        ui=SimpleNamespace(page=page,store=f.store,token=f.token,call=call)
        await residue_dialog(ui,f.pid,lambda:service.scene(f.token,f.pid,f.rid,f.sid,key),lambda:True,AsyncMock(),open_3d_url='http://localhost/molecular-viewer/test')
        button=dialogs[0].actions[0]
        self.assertIsNone(button.on_click);self.assertEqual(button.url.target.value,'_blank')
