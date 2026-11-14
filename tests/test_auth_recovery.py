"""Password recovery and remembered devices preserve authentication boundaries."""
import asyncio
import sqlite3
import tempfile
import time
import unittest
from types import SimpleNamespace
from unittest.mock import AsyncMock,MagicMock,Mock,patch

from biomolexplorer.auth_recovery import PasswordRecoveryService,SMTPEmailSender
from biomolexplorer.workspace import AccessDenied,WorkspaceStore
from biomolexplorer.ui.app import WorkspaceUI
from test_frontend_navigation import Page


class RecoveryServiceTests(unittest.TestCase):
    def test_smtp_sender_uses_starttls_authentication_and_recipient(self):
        client=MagicMock();client.__enter__.return_value=client
        environment={'BIOMOL_SMTP_HOST':'smtp.example.org','BIOMOL_SMTP_PORT':'587',
            'BIOMOL_SMTP_USERNAME':'sender@example.org','BIOMOL_SMTP_PASSWORD':'secret',
            'BIOMOL_SMTP_FROM':'BioMolExplorer <sender@example.org>','BIOMOL_SMTP_SECURITY':'starttls'}
        with patch.dict('os.environ',environment,clear=True),patch('smtplib.SMTP',return_value=client) as smtp:
            SMTPEmailSender.from_environment()('recipient@example.org','12345678')
        smtp.assert_called_once_with('smtp.example.org',587,timeout=15)
        client.starttls.assert_called_once();client.login.assert_called_once_with('sender@example.org','secret')
        message=client.send_message.call_args.args[0]
        self.assertEqual(message['To'],'recipient@example.org');self.assertIn('12345678',message.get_content())

    def test_missing_smtp_configuration_and_delivery_failure_do_not_issue_otp(self):
        with patch.dict('os.environ',{},clear=True),self.assertRaisesRegex(RuntimeError,'não configurado'):
            SMTPEmailSender.from_environment()
        store=SimpleNamespace(account_exists=lambda email:True)
        service=PasswordRecoveryService(store,sender=Mock(side_effect=RuntimeError('delivery failed')))
        with self.assertRaisesRegex(RuntimeError,'delivery failed'):service.request('user@example.org')
        self.assertNotIn('user@example.org',service._requests)

    def test_otp_is_eight_digits_single_use_and_resets_through_store(self):
        store=SimpleNamespace(account_exists=lambda email:email=='user@example.org',reset_password=Mock())
        delivered=[];service=PasswordRecoveryService(store,sender=lambda email,code:delivered.append((email,code)))
        code=service.request('USER@example.org')
        self.assertRegex(code,r'^\d{8}$');self.assertEqual(delivered,[('user@example.org',code)])
        with self.assertRaisesRegex(ValueError,'Código incorreto'):service.validate('user@example.org','00000000' if code!='00000000' else '11111111')
        grant=service.validate('user@example.org',code)
        service.reset('user@example.org',grant,'new-password')
        store.reset_password.assert_called_once_with('user@example.org','new-password')
        with self.assertRaisesRegex(ValueError,'expirou'):service.reset('user@example.org',grant,'another-password')

    def test_invalid_email_unknown_account_and_expired_code_are_rejected(self):
        store=SimpleNamespace(account_exists=lambda email:False)
        service=PasswordRecoveryService(store,sender=lambda *args:None)
        with self.assertRaisesRegex(ValueError,'e-mail válido'):service.request('invalid')
        with self.assertRaisesRegex(ValueError,'Não existe'):service.request('nobody@example.org')
        store.account_exists=lambda email:True
        code=service.request('user@example.org')
        digest,_,attempts=service._requests['user@example.org']
        service._requests['user@example.org']=(digest,time.time()-1,attempts)
        with self.assertRaisesRegex(ValueError,'expirou'):service.validate('user@example.org',code)


class RememberedSessionTests(unittest.IsolatedAsyncioTestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.store=WorkspaceStore(self.temp.name)
        self.user_id='a'*32;self.session='session-token'
        with self.store.connect() as db:
            db.execute('INSERT INTO users VALUES (?,?,?,?,?)',
                (self.user_id,'user@example.org','User','00'*16,'unused'))
            db.execute('INSERT INTO sessions VALUES (?,?,?)',
                (self.store._digest(self.session),self.user_id,time.time()+3600))

    async def test_remembered_device_resumes_and_can_be_fully_revoked(self):
        device=self.store.remember_session(self.session)
        resumed=self.store.resume_remembered(device)
        self.assertEqual(self.store.user(resumed)['email'],'user@example.org')
        self.store.logout(resumed)
        resumed_again=self.store.resume_remembered(device)
        self.assertEqual(self.store.user(resumed_again)['id'],self.user_id)
        self.store.forget_remembered(device)
        with self.assertRaises(AccessDenied):self.store.resume_remembered(device)

    async def test_password_reset_rehashes_and_revokes_every_session_and_device(self):
        device=self.store.remember_session(self.session)
        with patch.object(self.store,'_password',return_value='new-hash'):
            self.store.reset_password('USER@example.org','new-password')
        with self.store.connect() as db:
            user=db.execute('SELECT salt,password FROM users WHERE id=?',(self.user_id,)).fetchone()
            self.assertEqual(user['password'],'new-hash');self.assertNotEqual(user['salt'],'00'*16)
            self.assertEqual(db.execute('SELECT count(*) FROM sessions WHERE user_id=?',(self.user_id,)).fetchone()[0],0)
            self.assertEqual(db.execute('SELECT count(*) FROM remembered_sessions WHERE user_id=?',(self.user_id,)).fetchone()[0],0)
        with self.assertRaises(AccessDenied):self.store.resume_remembered(device)

    async def test_login_checkbox_persists_device_and_restore_opens_workspace(self):
        page=Page();ui=WorkspaceUI(page,self.store,SimpleNamespace(),language='en',structure_viewer=SimpleNamespace())
        ui.show_login();experience=ui._login_experience;experience.remember.value=True
        experience.email.value='user@example.org';experience.password.value='password'
        ui.store.login=Mock(return_value=self.session)
        ui.show_workspace=AsyncMock()
        async def direct(function,*args,**kwargs):return function(*args,**kwargs)
        ui.call=direct
        await experience.submit.on_click(None)
        device=await ui.session_persistence.get()
        self.assertEqual(device,ui._remember_device)
        ui.token=None;ui.show_workspace.reset_mock()
        async def inline(function,*args):return function(*args)
        with patch('biomolexplorer.ui.app.asyncio.to_thread',side_effect=inline):
            self.assertTrue(await ui.restore_session())
        ui.show_workspace.assert_awaited_once_with()

    async def test_logout_without_remember_is_direct_and_remembered_logout_offers_choices(self):
        page=Page();ui=WorkspaceUI(page,self.store,SimpleNamespace(),language='en',structure_viewer=SimpleNamespace())
        async def direct(function,*args,**kwargs):return function(*args,**kwargs)
        ui.call=direct;ui.token=self.session
        async def inline(function,*args):return function(*args)
        with patch('biomolexplorer.ui.app.asyncio.to_thread',side_effect=inline):
            await ui.logout(None)
        self.assertIsNone(ui.token);self.assertFalse(page.dialogs)

        with self.store.connect() as db:
            db.execute('INSERT OR REPLACE INTO sessions VALUES (?,?,?)',
                (self.store._digest(self.session),self.user_id,time.time()+3600))
        ui.token=self.session;ui._remember_device=self.store.remember_session(self.session)
        await ui.session_persistence.set(ui._remember_device)
        await ui.logout(None)
        self.assertTrue(page.dialogs);dialog=page.dialogs[-1]
        labels=[getattr(action,'content',None) for action in dialog.actions]
        self.assertIn('Sign out only',labels);self.assertIn('Sign out and disable Remember me',labels)
        device=ui._remember_device
        with patch('biomolexplorer.ui.app.asyncio.to_thread',side_effect=inline):
            await ui._perform_logout(False)
        self.assertEqual(await ui.session_persistence.get(),device)
        self.assertIsNotNone(self.store.resume_remembered(device))

        with self.store.connect() as db:
            db.execute('INSERT OR REPLACE INTO sessions VALUES (?,?,?)',
                (self.store._digest(self.session),self.user_id,time.time()+3600))
        ui.token=self.session
        with patch('biomolexplorer.ui.app.asyncio.to_thread',side_effect=inline):
            await ui._perform_logout(True)
        self.assertIsNone(await ui.session_persistence.get());self.assertIsNone(ui._remember_device)
        with self.assertRaises(AccessDenied):self.store.resume_remembered(device)

