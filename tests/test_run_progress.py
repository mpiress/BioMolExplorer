"""Real worker boundary, progress lifecycle and novice-facing run feedback."""
import asyncio
import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import AsyncMock, patch

from biomolexplorer.catalog import new_stage
from biomolexplorer.diagnostics import log_directory
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.progress import report_progress, read_progress
from biomolexplorer.workspace import WorkspaceStore

try:
    from biomolexplorer.ui.app import WorkspaceUI
    from biomolexplorer.ui.run_progress import RunProgress
except ImportError:
    RunProgress = None


class WorkerRetrievalTests(unittest.TestCase):
    def test_single_chembl220_block_uses_all_catalog_defaults_in_real_worker(self):
        original = subprocess.Popen
        fixture = Path(__file__).parent / 'fixtures' / 'chembl_worker.py'
        commands = []
        def launch(args, **kwargs):
            commands.append(args)
            return original([args[0], str(fixture), *args[3:]], **kwargs)
        with tempfile.TemporaryDirectory(prefix='biomol default retrieval ') as temp:
            store = WorkspaceStore(temp)
            token = store.register('Test', 'test@example.org', 'test-password')
            project = store.create_project(token, 'CHEMBL220 defaults')
            stage = new_stage('retrieve_compounds')
            self.assertTrue(stage['parameters']['include_pubchem'])
            self.assertEqual(stage['parameters']['search_term'], 'CHEMBL220')
            self.assertEqual(stage['templates'], {})
            self.assertEqual(stage['depends_on'], [])
            store.save_pipeline(token, project['id'], [stage], project['revision'])
            with patch('biomolexplorer.jobs.subprocess.Popen', side_effect=launch):
                service = PipelineService(store, worker_python=Path(sys.executable))
                try:
                    run = service.submit(token, project['id'], reuse_results=False)
                    deadline = time.monotonic() + 20
                    while time.monotonic() < deadline:
                        run = store.get_run(token, run['id'])
                        if run['status'] not in ('queued', 'running'):
                            break
                        time.sleep(.05)
                    self.assertEqual(run['status'], 'succeeded', run.get('error'))
                    result = run['stages'][0]
                    dataset = next(path for path in result['artifacts'] if path.endswith('/compounds/CHEMBL220/compounds.csv'))
                    import pandas as pd
                    compounds = pd.read_csv(dataset)
                    self.assertEqual(compounds['molecule_chembl_id'].tolist(), ['CHEMBL4087364', 'PUBCHEM3'])
                    self.assertEqual(compounds['source'].tolist(), ['ChEMBL', 'PubChem'])
                finally:
                    service.close()
        self.assertEqual(len(commands), 1)

    def test_provider_failure_cancels_pending_molecular_requests(self):
        from concurrent.futures import Future
        from crawlers.molecules import _monitor_queries
        failed, pending = Future(), Future()
        failed.set_exception(RuntimeError('Provider failed'))
        with self.assertRaisesRegex(RuntimeError, 'Provider failed'):
            _monitor_queries({failed: 'CHEMBL1', pending: 'CHEMBL2'}, 'Retrieval')
        self.assertTrue(pending.cancelled())

    def test_pipeline_calls_real_worker_business_function_with_ic50_and_reports_provider_errors(self):
        original = subprocess.Popen
        fixture = Path(__file__).parent / 'fixtures' / 'chembl_worker.py'
        commands = []
        def launch(args, **kwargs):
            commands.append(args)
            self.assertEqual(args[1:3], ['-m', 'biomolexplorer.worker'])
            return original([args[0], str(fixture), *args[3:]], **kwargs)
        with tempfile.TemporaryDirectory(prefix='biomol retrieval integration ') as temp:
            store = WorkspaceStore(temp)
            token = store.register('Test', 'test@example.org', 'test-password')
            project = store.create_project(token, 'CHEMBL220 IC50')
            stage = new_stage('retrieve_compounds')
            stage['parameters'].update(search_term='CHEMBL220', include_pubchem=False,
                chembl_filters={'bioactivity': {'standard_type__in': ['IC50'], 'standard_units': 'nM', 'max_value_ref': 5000}})
            store.save_pipeline(token, project['id'], [stage], project['revision'])
            with patch('biomolexplorer.jobs.subprocess.Popen', side_effect=launch):
                service = PipelineService(store, worker_python=Path(sys.executable))
                try:
                    for failure in ('', 'timeout', '500'):
                        with self.subTest(provider_failure=failure), patch.dict(os.environ, {'BIOMOL_TEST_CHEMBL_FAILURE': failure}):
                            run = service.submit(token, project['id'], reuse_results=False)
                            deadline = time.monotonic() + 20
                            while time.monotonic() < deadline:
                                run = store.get_run(token, run['id'])
                                if run['status'] not in ('queued', 'running'):
                                    break
                                time.sleep(.05)
                            self.assertEqual(run['status'], 'failed' if failure else 'succeeded', run.get('error'))
                            result = run['stages'][0]
                            self.assertIn('job_id', result)
                            self.assertGreaterEqual(result['finished_at'], result['started_at'])
                            self.assertIn('message', result['progress'])
                            log = Path(result['log_path']).read_text()
                            diagnostic = log_directory() / 'jobs' / result['job_id'] / 'backend.log'
                            self.assertIn('Worker iniciado; python=', diagnostic.read_text())
                            if failure:
                                self.assertIn('tempo limite' if failure == 'timeout' else 'HTTP 500', run['error'])
                            else:
                                dataset = next(p for p in result['artifacts'] if p.endswith('/compounds.csv'))
                                self.assertIn('CHEMBL4087364', Path(dataset).read_text())
                finally:
                    service.close()
        self.assertEqual(len(commands), 3)

    def test_progress_is_optional_atomic_and_invalid_feedback_is_ignored(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'progress.json'
            with patch.dict(os.environ, {'BIOMOL_PROGRESS_FILE': str(path)}):
                report_progress('Molecules', 2, 3)
            self.assertEqual(read_progress(path)['completed'], 2)
            path.write_text('{')
            self.assertIsNone(read_progress(path))
            path.write_text('[]')
            self.assertIsNone(read_progress(path))
            self.assertIsNone(read_progress(path.parent / 'missing'))


@unittest.skipIf(RunProgress is None, 'Install the ui extra for Flet controls')
class RunFeedbackTests(unittest.TestCase):
    def setUp(self):
        self.ui = object.__new__(WorkspaceUI)
        self.ui.token = 'session'
        self.ui.current = {'id': 'project', 'role': 'editor'}
        self.ui.page = SimpleNamespace(height=900, update=lambda: None, show_dialog=lambda d: self.dialogs.append(d))
        self.dialogs = []
        self.ui.draw_project = AsyncMock()
        self.ui.cancel_run = AsyncMock()
        self.run = {'id': 'run', 'project_id': 'project', 'status': 'queued', 'created': 100, 'updated': 100,
            'stages': [{'id': 'stage', 'name': 'Recuperar compostos', 'operation': 'retrieve_compounds', 'status': 'queued',
                        'configuration': {'enabled': True}, 'artifacts': []}], 'error': None}
        self.progress = RunProgress(self.ui, self.run)
        self.ui.run_progress = self.progress

    def test_running_stage_count_and_elapsed_are_visible_without_fake_percentage(self):
        self.run['status'] = 'running'
        self.run['stages'][0].update(status='running', started_at=110,
            progress={'message': 'Baixando moléculas', 'completed': 119, 'total': 133})
        self.progress.update(self.run, now=130)
        self.assertIn('Recuperar compostos', self.progress.message.value)
        self.assertIn('119 de 133', self.progress.message.value)
        self.assertIn('00:00:30', self.progress.timing.value)
        self.assertIn('00:00:20', self.progress.timing.value)
        self.assertIsNone(self.progress.bar.value)
        self.assertTrue(self.progress.bar.visible)

    def test_english_badge_translates_status_and_preserves_custom_stage_name(self):
        self.ui.language = 'en'
        self.run['status'] = 'running'
        self.run['stages'][0]['status'] = 'running'
        self.progress.update(self.run, now=130)
        self.assertIn('Running', self.progress.badge_text.value)
        self.assertIn('Retrieve compounds', self.progress.badge_text.value)
        self.run['stages'][0]['name'] = 'Executando meu estudo'
        self.progress.update(self.run, now=130)
        self.assertIn('Executando meu estudo', self.progress.badge_text.value)

    def test_minimizing_keeps_tracking_and_terminal_error_has_log_access(self):
        self.progress.show()
        asyncio.run(self.progress.dismiss())
        self.assertIs(self.ui.run_progress, self.progress)
        self.assertFalse(self.progress.is_open)
        self.run.update(status='failed', updated=125, error='ChEMBL indisponível (HTTP 500)')
        self.run['stages'][0].update(status='failed', log_path='/private/log')
        self.progress.update(self.run)
        self.progress.show()
        self.progress.dismissed(None)
        self.assertTrue(self.progress.is_open)
        self.assertIn('HTTP 500', self.progress.error.value)
        self.assertTrue(self.progress.log.visible)
        self.assertFalse(self.progress.bar.visible)
        self.assertFalse(self.progress.cancel.visible)
        asyncio.run(self.progress.dismiss())
        self.assertIsNone(self.ui.run_progress)

    def test_readers_cannot_cancel_and_stale_dialogs_do_not_reopen(self):
        self.ui.current['role'] = 'viewer'
        self.progress.update(self.run)
        self.assertFalse(self.progress.cancel.visible)
        self.ui.token = 'new-session'
        self.progress.show()
        asyncio.run(self.progress.request_cancel(None))
        self.assertEqual(self.dialogs, [])
        self.ui.cancel_run.assert_not_called()

    def test_waiting_feedback_names_the_block_and_offers_configuration_to_editors(self):
        self.run.update(status='awaiting_input',error='Configure o bloco “Avaliar ADMET” e escolha os arquivos de entrada.')
        self.run['stages'][0].update(status='awaiting_input',name='Avaliar ADMET')
        self.progress.update(self.run)
        self.assertIn('Avaliar ADMET',self.progress.message.value)
        self.assertTrue(self.progress.configure.visible)
        self.assertTrue(self.progress.cancel.visible)
        self.assertFalse(self.progress.bar.visible)
        self.assertFalse(self.progress.error.visible)
        self.ui.current['role']='viewer'
        self.progress.update(self.run)
        self.assertFalse(self.progress.configure.visible)

    def test_cancel_is_requested_once_and_feedback_is_immediate(self):
        asyncio.run(self.progress.request_cancel(None))
        asyncio.run(self.progress.request_cancel(None))
        self.ui.cancel_run.assert_awaited_once_with('run')
        self.assertTrue(self.progress.cancel.disabled)
        self.assertIn('Cancelamento solicitado', self.progress.hint.value)

    def test_poller_updates_terminal_popup_and_stops(self):
        self.ui.current_run = 'run'
        self.ui.tab = 'Execuções'
        self.ui.store = SimpleNamespace(get_run=lambda *args: self.run)
        async def local_call(function, *args):
            return function(*args)
        self.ui.call = local_call
        self.run.update(status='succeeded', updated=125)
        self.run['stages'][0]['status'] = 'succeeded'
        with patch('biomolexplorer.ui.app.asyncio.sleep', new=AsyncMock()):
            asyncio.run(self.ui.poll_runs())
        self.assertIsNone(self.ui.current_run)
        self.assertTrue(self.progress.is_open)
        self.assertEqual(self.progress.heading.value, 'Pipeline concluído')

    def test_reopen_waits_for_flutter_to_finish_removing_the_previous_dialog(self):
        self.progress.show()
        asyncio.run(self.progress.dismiss())
        self.progress.show()
        self.assertTrue(self.progress.reopen_requested)
        self.assertEqual(len(self.dialogs), 1)
        self.progress.dismissed(None)
        self.assertTrue(self.progress.is_open)
        self.assertEqual(len(self.dialogs), 2)

    def test_completed_run_does_not_clear_a_new_run_started_during_refresh(self):
        self.ui.current_run = 'run'
        self.ui.tab = 'Execuções'
        self.ui.store = SimpleNamespace(get_run=lambda *args: self.run)
        async def local_call(function, *args):
            return function(*args)
        self.ui.call = local_call
        self.run.update(status='succeeded', updated=125)
        async def start_next_run():
            self.ui.current_run = 'new-run'
        self.ui.draw_project = start_next_run
        waits = []
        async def wait_once(delay):
            waits.append(self.ui.current_run)
            if len(waits) == 2:
                self.assertEqual(self.ui.current_run, 'new-run')
                self.ui.token = None
        with patch('biomolexplorer.ui.app.asyncio.sleep', side_effect=wait_once):
            asyncio.run(self.ui.poll_runs())
        self.assertEqual(waits, ['run', 'new-run'])
