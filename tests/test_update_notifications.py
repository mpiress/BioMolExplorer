import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch
from biomolexplorer.updates import UpdateChecker, CHECK_INTERVAL, REMIND_INTERVAL


class UpdateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.checker = UpdateChecker(self.tmp.name, revision='1' * 40)

    def response(self, sha):
        response = MagicMock()
        response.__enter__.return_value.read.return_value = json.dumps({'sha': sha,
            'commit': {'message': 'Changes\nDetails'}}).encode()
        return response

    @patch('biomolexplorer.updates.subprocess.run', return_value=MagicMock(returncode=1))
    def test_changed_master_pins_download_and_throttles_checks(self, git):
        with patch('biomolexplorer.updates.urlopen', return_value=self.response('2' * 40)) as http:
            latest = self.checker.check()
            self.assertTrue(latest['available'])
            self.assertTrue(latest['download_url'].endswith('/' + '2' * 40 + '.zip'))
            self.checker.check()
            self.assertEqual(http.call_count, 1)
            self.checker.check(force=True)
            self.assertEqual(http.call_count, 2)

    def test_same_revision_and_network_failure(self):
        with patch('biomolexplorer.updates.urlopen', return_value=self.response('1' * 40)):
            self.assertFalse(self.checker.check()['available'])
        with patch('biomolexplorer.updates.urlopen', side_effect=OSError('offline')):
            self.assertIsNone(self.checker.check(force=True))
            self.assertEqual(self.checker.error, 'offline')

    def test_reminder_persists_per_user_and_expires(self):
        self.checker.latest = {'available': True}
        with patch('biomolexplorer.updates.time.time', return_value=100):
            self.checker.remind_later('alice')
            self.assertFalse(self.checker.visible('alice'))
            self.assertTrue(self.checker.visible('bob'))
        checker = UpdateChecker(self.tmp.name, revision='1' * 40)
        checker.latest = {'available': True}
        with patch('biomolexplorer.updates.time.time', return_value=101+REMIND_INTERVAL):
            self.assertTrue(checker.visible('alice'))

    def test_ui_visibility_does_not_wait_for_network_check_lock(self):
        import threading
        self.checker.latest={'available':True}
        result=[]
        with self.checker.lock:
            worker=threading.Thread(target=lambda:result.append(self.checker.visible('alice')))
            worker.start();worker.join(timeout=.5)
            self.assertFalse(worker.is_alive())
        self.assertEqual(result,[True])

    def test_unversioned_zip_records_first_observation(self):
        self.checker.revision = None
        with patch('biomolexplorer.updates.urlopen', return_value=self.response('2' * 40)):
            self.assertFalse(self.checker.check()['available'])
        self.assertEqual(json.loads(Path(self.tmp.name, 'updates.json').read_text())['baseline'], '2' * 40)

    @patch('biomolexplorer.updates.subprocess.run', return_value=MagicMock(returncode=0))
    def test_local_checkout_ahead_of_master_does_not_offer_downgrade(self, git):
        with patch('biomolexplorer.updates.urlopen', return_value=self.response('2' * 40)):
            self.assertFalse(self.checker.check()['available'])


class UpdateNoticeTests(unittest.IsolatedAsyncioTestCase):
    async def test_download_uses_pinned_archive_and_persists_reminder(self):
        from types import SimpleNamespace
        from unittest.mock import AsyncMock
        from biomolexplorer.ui.updates import UpdateNotice
        checker = MagicMock()
        checker.check.return_value = {'sha':'2'*40, 'available':True,
            'message':'New update', 'download_url':'https://github.com/mpiress/BioMolExplorer/archive/'+'2'*40+'.zip'}
        dialogs = []
        page = SimpleNamespace(show_dialog=dialogs.append, update=lambda:None, launch_url=AsyncMock())
        async def call(function,*args):return function(*args)
        ui = SimpleNamespace(page=page, token='token', store=SimpleNamespace(user=lambda token:{'id':'alice'}),
            call=call, tr=lambda value:value, notify=MagicMock())
        notice=UpdateNotice(ui,checker)
        notice.button()
        await notice.open(None)
        await dialogs[0].actions[1].on_click(None)
        page.launch_url.assert_awaited_once_with(checker.check.return_value['download_url'])
        checker.remind_later.assert_called_once_with('alice')
        self.assertFalse(dialogs[0].open)

    async def test_offline_check_reports_retry_without_update_dialog(self):
        from types import SimpleNamespace
        from biomolexplorer.ui.updates import UpdateNotice
        checker=MagicMock();checker.check.return_value=None
        async def call(function,*args):return function(*args)
        ui=SimpleNamespace(call=call,notify=MagicMock())
        notice=UpdateNotice(ui,checker)
        await notice.open(None)
        ui.notify.assert_called_once()
