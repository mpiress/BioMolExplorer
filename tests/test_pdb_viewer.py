"""Browser viewer authorization, package integrity and web/desktop integration."""
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from biomolexplorer.pdb_view import StructureViewers,viewer_document,RESOURCE,PREFIX
from biomolexplorer.workspace import WorkspaceStore,AccessDenied
from biomolexplorer.visualizations import MAX_VIEW_BYTES


class PDBViewerTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.store=WorkspaceStore(self.temp.name)
        self.token=self.store.register('Owner','owner@example.org','owner-password')
        self.pid=self.store.create_project(self.token,'Structures')['id']
        self.path=self.store.project_dir(self.pid)/'1CRN.pdb'
        self.path.write_bytes((Path(__file__).parent/'fixtures/dms/1CRN.pdb').read_bytes())
        self.now=100.;self.viewers=StructureViewers(self.store,clock=lambda:self.now)
    def tearDown(self):self.viewers.close();self.temp.cleanup()
    def key(self):return self.viewers.issue(self.token,self.pid,str(self.path),'en')

    def test_urls_contain_capability_and_never_login_token_or_project_path(self):
        key=self.key();url=self.viewers.url(key,True,'https://example.org:8443/')
        self.assertEqual(url,'https://example.org:8443'+PREFIX+'/'+key)
        self.assertNotIn(self.token,url);self.assertNotIn(self.pid,url)
        response=self.viewers.response(PREFIX+'/'+key)
        self.assertEqual(response[0],200);document=response[2].decode()
        self.assertIn('Drag to rotate',document);self.assertNotIn(self.token,document)
        self.assertNotIn(str(self.path),document);self.assertEqual(response[3]['Cache-Control'],'no-store')
        response=self.viewers.response(PREFIX+'/'+key+'/structure')
        self.assertEqual(response[2],self.path.read_bytes())

    def test_flet_websocket_addresses_open_as_http_pages(self):
        key=self.key()
        for page_url,expected in (
            ('ws://127.0.0.1:8550/', 'http://127.0.0.1:8550'),
            ('wss://example.org:8443/ws?session=private', 'https://example.org:8443'),
            ('http://localhost:8550/', 'http://localhost:8550')):
            self.assertEqual(self.viewers.url(key,True,page_url),expected+PREFIX+'/'+key)
        for invalid in ('file:///tmp/test', 'javascript:alert(1)', 'ws:///missing-host'):
            with self.assertRaises(ValueError):self.viewers.url(key,True,invalid)

    def test_expired_revoked_and_logged_out_access_is_denied(self):
        key=self.key();self.now+=self.viewers.TTL+1
        self.assertEqual(self.viewers.response(PREFIX+'/'+key)[0],403)
        guest=self.store.register('Reader','reader@example.org','reader-password')
        self.store.invite(self.token,self.pid,'reader@example.org','viewer');self.store.accept_invitation(guest,self.pid)
        reader_key=self.viewers.issue(guest,self.pid,str(self.path))
        self.assertEqual(self.viewers.response(PREFIX+'/'+reader_key+'/structure')[0],200)
        self.store.revoke(self.token,self.pid,self.store.user(guest)['id'])
        self.assertEqual(self.viewers.response(PREFIX+'/'+reader_key+'/structure')[0],403)
        key=self.key();self.store.logout(self.token)
        self.assertEqual(self.viewers.response(PREFIX+'/'+key+'/structure')[0],403)

    def test_capability_cannot_be_issued_for_other_projects_or_oversized_files(self):
        guest=self.store.register('Guest','guest@example.org','guest-password')
        with self.assertRaises(AccessDenied):self.viewers.issue(guest,self.pid,str(self.path))
        with self.assertRaises(AccessDenied):self.viewers.issue(self.token,self.pid,str(Path(self.temp.name)/'outside.pdb'))
        with self.path.open('r+b') as stream:stream.truncate(MAX_VIEW_BYTES+1)
        with self.assertRaisesRegex(ValueError,'32 MB'):self.key()

    def test_bundled_assets_are_pinned_and_never_expose_arbitrary_files(self):
        provenance=json.loads((RESOURCE/'third_party.json').read_text())
        for name,digest in provenance['sha256'].items():self.assertEqual(hashlib.sha256((RESOURCE/name).read_bytes()).hexdigest(),digest)
        self.assertEqual(self.viewers.response(PREFIX+'/assets/3Dmol-min.js')[0],200)
        self.assertEqual(self.viewers.response(PREFIX+'/assets/../../project.json')[0],404)
        self.assertEqual(self.viewers.response(PREFIX+'/assets/third_party.json')[0],404)
        self.assertEqual(self.viewers.response('/other')[0],404)
        document=viewer_document('</script><img src=x onerror=alert(1)>.pdb')
        self.assertNotIn('</script><img',document)
        self.assertIn('\\u003c/script',document)

    def test_desktop_server_is_loopback_and_rechecks_permissions_on_fetch(self):
        import requests
        key=self.key()
        try:url=self.viewers.url(key)
        except PermissionError:self.skipTest('Loopback sockets are restricted in this sandbox.')
        self.assertTrue(url.startswith('http://127.0.0.1:'))
        self.assertEqual(requests.get(url+'/structure',timeout=5).content,self.path.read_bytes())
        self.store.logout(self.token)
        self.assertEqual(requests.get(url+'/structure',timeout=5).status_code,403)

    def test_same_origin_web_routes_preserve_flet_and_authorize_structures(self):
        import socket
        try:
            probe=socket.socket(socket.AF_INET,socket.SOCK_STREAM);probe.close()
        except PermissionError:self.skipTest('ASGI test transport is restricted in this sandbox.')
        from fastapi.testclient import TestClient
        from biomolexplorer.ui.web_host import create_web_app
        async def session(page):pass
        app=create_web_app(session,self.store,self.viewers,'test-secret')
        key=self.key()
        with TestClient(app) as client:
            self.assertEqual(client.get('/').status_code,200)
            self.assertEqual(client.get(PREFIX+'/'+key+'/structure').content,self.path.read_bytes())
            self.assertEqual(client.get(PREFIX+'/no-such-capability').status_code,403)
            self.assertIn("frame-ancestors 'none'",client.get(PREFIX+'/'+key).headers['content-security-policy'])
            self.store.logout(self.token)
            self.assertEqual(client.get(PREFIX+'/'+key+'/structure').status_code,403)

    def test_pdb_button_opens_protected_viewer_directly(self):
        from biomolexplorer.ui.app import WorkspaceUI
        from biomolexplorer.ui.pdb_results import PDBActions
        ui=object.__new__(WorkspaceUI);ui.store=self.store;ui.token=self.token;ui.language='en'
        ui.structure_viewer=self.viewers;ui.page=SimpleNamespace(web=True,url='wss://example.org/ws')
        results=SimpleNamespace(ui=ui,project_id=self.pid,valid=lambda:True)
        actions=PDBActions(results,True).actions({'name':'1CRN.pdb','path':str(self.path)})
        self.assertTrue(actions[1].url.url.startswith('https://example.org'+PREFIX+'/'))
        self.assertEqual(actions[1].url.target.value,'_blank')
        self.assertIsNone(actions[1].on_click)
        self.assertEqual(actions[2].url,'https://www.rcsb.org/structure/1CRN')
