"""Application contracts, local process lifecycle and scientific regressions."""
import inspect
import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from importlib import import_module
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import pandas as pd
from biomolexplorer.config import AppConfig
from biomolexplorer.jobs import JobManager, TERMINAL
from biomolexplorer.operations import OPERATIONS, execute_operation, validate_operation
from biomolexplorer.paths import resolve_path
from biomolexplorer.processes import run_command
from kernel.utilities import fileHandling, fileReading
from wrappers.admet import ADMETWrapper


class ApplicationTests(unittest.TestCase):
    def test_config_and_parameter_validation(self):
        with self.assertRaises(ValueError):
            AppConfig(Path('/tmp'), cpu_workers=0)
        for params in ({'search_term': '../bad'}, {'search_term': 'A', 'pubchem_threshold': True},
                       {'search_term': 'A', 'extra': 1}):
            with self.assertRaises(ValueError):
                validate_operation('retrieve_compounds', params)
        with self.assertRaises(ValueError):
            validate_operation('unknown', {})
        validate_operation('retrieve_compounds', {'search_term': 'CHEMBL220', 'include_pubchem': True})
        for timeout in (True, '60', float('nan'), float('inf'), 0):
            with self.assertRaises(ValueError):
                AppConfig(Path('/tmp'), job_timeout=timeout)

    def test_dock6_radius_accepts_decimal_values(self):
        params = {'base_input_path': '/tmp/pdb', 'target': 'TARGET',
                  'base_selected_mols': '/tmp/mols', 'dock6_app_path': '/tmp/dock6',
                  'charge_type': 'gas', 'mol_filename': 'molecules',
                  'pdb_code': ['1ABC', 'LIG', 1, 'A'], 'base_vina_path': '/tmp/vina', 'radius': 1.4}
        validate_operation('docking_dock6', params)
        for radius in (False, 0, -1, float('nan')):
            with self.assertRaises(ValueError):
                validate_operation('docking_dock6', dict(params, radius=radius))
        with self.assertRaises(ValueError):
            validate_operation('fingerprints', {'base_input_path': '/tmp/mols', 'radius': 1.4})

    def test_registry_matches_wrapper_signatures(self):
        for name, spec in OPERATIONS.items():
            with self.subTest(operation=name):
                signature = inspect.signature(getattr(import_module(spec.module), spec.function))
                self.assertTrue(set(spec.required + spec.optional).issubset(signature.parameters))
                supplied = set(spec.required) | {'base_output_path'}
                required = {p.name for p in signature.parameters.values() if p.default is inspect.Parameter.empty}
                self.assertTrue(required.issubset(supplied), (name, required - supplied))

    def test_paths_resources_and_files(self):
        with tempfile.TemporaryDirectory(prefix='biomol with spaces ') as temp:
            root = Path(temp)
            with patch.dict(os.environ, {'BIOMOL_WORKSPACE': temp}):
                self.assertEqual(resolve_path('/datasets/compound.csv'), root / 'datasets/compound.csv')
                self.assertEqual(resolve_path(root / 'input'), root / 'input')
                self.assertTrue(resolve_path('src/scripts/crawlers/target.json').is_file())
            frame = pd.DataFrame({'value': [1, 2]})
            handler = fileHandling(input_path=str(root), output_path=str(root))
            handler.dataframe_to_csv('test', frame)
            pd.testing.assert_frame_equal(handler.csv_to_dataframe('test'), frame)
            with self.assertRaises(FileNotFoundError):
                handler.csv_to_dataframe('missing')
            reader = fileReading(str(root), 'test.csv')
            self.assertEqual(list(reader), ['value', '1', '2'])

    def test_fragment_filter_and_stable_deduplication(self):
        from kernel.filters import Molecule
        filters = Molecule()
        self.assertEqual(filters.clean_fragments(['CCO', 'CC.O', 'invalid']), ['CCO'])
        self.assertEqual(filters.remove_duplicates(['CCO', 'OCC', 'CCC']), ['CCO', 'CCC'])

    def test_redocking_service_preserves_original_inputs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            original = root / 'source/TARGET'
            original.mkdir(parents=True)
            pdb = original / '1ABC.pdb'
            pdb.write_text('original')
            def redock(**kwargs):
                working = Path(kwargs['base_input_path']) / 'TARGET/1ABC.pdb'
                self.assertNotEqual(working, pdb)
                working.unlink()
            with patch('wrappers.redocking.perform_redocking', side_effect=redock):
                execute_operation('redocking', {'target': 'TARGET', 'base_input_path': str(root / 'source'), 'pdb_codes': [['1ABC','LIG',1,'A']]}, root / 'out')
            self.assertEqual(pdb.read_text(), 'original')

    def test_admet_excludes_flagged_compounds_and_handles_empty(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            pd.DataFrame({'molecule_chembl_id': ['safe', 'flagged'], 'canonical_smiles': ['CCO', 'CCO']}).to_csv(root / 'input.csv', index=False)
            wrapper = ADMETWrapper(str(root / 'out'), str(root), 'input.csv')
            with patch.object(wrapper.evaluator, 'is_toxic', side_effect=[False, True]), patch.object(wrapper, 'generate_plot'):
                result = wrapper.run_pipeline()
            self.assertEqual(result['molecule_chembl_id'].tolist(), ['safe'])
            self.assertEqual(wrapper.excluded_count, 1)
            with patch.object(wrapper.evaluator, 'is_toxic', return_value=True), patch.object(wrapper, 'generate_plot'):
                result = wrapper.run_pipeline()
            self.assertTrue(result.empty)
            self.assertIn('BBB', result.columns)

    def test_commands_handle_spaces_failures_and_timeout(self):
        self.assertTrue(run_command([sys.executable, '-c', 'import sys; assert sys.argv[1] == "a b"', 'a b']))
        with self.assertRaisesRegex(RuntimeError, r'código 2.*invalid input'):
            run_command([sys.executable, '-c', 'import sys; sys.stderr.write("invalid input"); raise SystemExit(2)'])
        with self.assertRaises(subprocess.TimeoutExpired):
            run_command([sys.executable, '-c', 'import time; time.sleep(10)'], timeout=0.05)

    def test_fingerprints_similarity_and_graphs_can_chain_outputs(self):
        with tempfile.TemporaryDirectory() as temp, patch.dict(os.environ, {'BIOMOL_CPU_WORKERS': '1'}):
            root = Path(temp)
            source = root / 'input'
            source.mkdir()
            pd.DataFrame({'molecule_chembl_id': ['A', 'B'], 'canonical_smiles': ['CCO', 'CCCO']}).to_csv(
                source / 'compounds.csv', index=False)
            fingerprints = execute_operation('fingerprints', {'base_input_path': str(source),
                'maccs': False, 'pharmacophore': False, 'chunk_size': 1}, root / 'fp')
            self.assertEqual(len(fingerprints.artifacts), 1)
            similarity = execute_operation('similarity', {'base_input_path': str(root / 'fp'),
                'threshold': 1, 'approximate': False}, root / 'similarity')
            self.assertEqual(len(similarity.artifacts), 1)
            graphs = execute_operation('graphs', {'base_input_path': str(source),
                'similarity_path': str(Path(similarity.artifacts[0]).parent)}, root / 'graphs')
            self.assertTrue(any('maxcomp' in path for path in graphs.artifacts))

    def test_identical_fingerprints_retain_distinct_identifiers(self):
        from kernel.descriptors import MolSimilarity, similarityFunctions
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            pd.DataFrame({'molecule_chembl_id': ['A', 'B', 'C'],
                          'fingerprint': ['[1, 0, 1]'] * 3}).to_csv(root / 'morgan_test.csv', index=False)
            similarity = MolSimilarity(str(root), str(root / 'out'), approximate=False)
            similarity.perform_similarity(metric=similarityFunctions.TanimotoSimilarity)
            edges = pd.read_csv(root / 'out/Tanimoto_morgan_test.csv')
            self.assertEqual(len(edges), 6)
            self.assertFalse((edges['source'] == edges['target']).any())
            self.assertTrue((edges['value'] == 1).all())

    def test_docking_preparation_bootstraps_centers_and_keeps_each_chain(self):
        from caad.docking import Docking
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            dock = Docking(complex_input_path=str(root), output_path=str(root / 'prepared'))
            records = [('1AAA', 'LIG', '1', 'A'), ('2BBB', 'LIG', '2', 'B')]
            with patch.object(dock, 'generate_docking_script'), patch.object(dock, 'process_in_parallel'), \
                    patch.object(dock, 'prepare_on_obabel') as convert, \
                    patch.object(dock, 'calculate_ligand_centerofmass', return_value=[1, 2, 3]):
                prepared = dock.prepare_for_docking(records, 'gas', 7.4, True)
            self.assertEqual(prepared, records)
            converted = [call.args[0] for call in convert.call_args_list]
            self.assertIn('1AAA_A.dockprep.mol2', converted)
            self.assertIn('2BBB_B.dockprep.mol2', converted)
            centers = pd.read_csv(root / 'prepared/centers.csv')
            self.assertEqual(set(centers.columns), {'1AAA_LIG_1A', '2BBB_LIG_2B'})


class JobTests(unittest.TestCase):
    def test_initialization_failure_releases_workspace_lock(self):
        handles = []
        def fail(manager):
            handles.append(manager._owner)
            raise RuntimeError('Database unavailable')
        with tempfile.TemporaryDirectory() as temp:
            config = AppConfig(Path(temp))
            with patch.object(JobManager, '_connect', fail), self.assertRaises(RuntimeError):
                JobManager(config)
            self.assertTrue(handles[0].closed)
            with JobManager(config):
                pass

    def test_scheduler_failure_does_not_leave_job_queued(self):
        with tempfile.TemporaryDirectory() as temp, JobManager(AppConfig(Path(temp))) as manager:
            with patch.object(manager._executor, 'submit', side_effect=RuntimeError('Executor unavailable')):
                with self.assertRaises(RuntimeError):
                    manager.submit('admet', {'base_input_path': temp})
            job = manager.list()[0]
            self.assertEqual(job['status'], 'failed')
            self.assertIn('Could not schedule worker', job['error'])

    def test_shutdown_finds_active_jobs_outside_recent_history(self):
        with tempfile.TemporaryDirectory() as temp:
            manager = JobManager(AppConfig(Path(temp)))
            with manager._connect() as connection:
                connection.execute('INSERT INTO jobs VALUES (?,?,?,?,?,?,?,?,?,?)',
                    ('older-active', 'admet', '{}', 'queued', '0000', '0000', temp, temp, None, None))
                connection.executemany('INSERT INTO jobs VALUES (?,?,?,?,?,?,?,?,?,?)',
                    [(f'finished-{index}', 'admet', '{}', 'succeeded', '9999', '9999', temp, temp, None, None)
                     for index in range(1001)])
            manager.close()
            self.assertEqual(manager.get('older-active')['status'], 'cancelled')

    def wait_for(self, manager, job_id, predicate, timeout=20):
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            job = manager.get(job_id)
            if predicate(job):
                return job
            time.sleep(0.03)
        self.fail('Job did not reach expected state')

    def test_real_worker_success_failure_and_persistence(self):
        with tempfile.TemporaryDirectory(prefix='biomol jobs ') as temp:
            root = Path(temp)
            pd.DataFrame({'molecule_chembl_id': ['test'], 'canonical_smiles': ['CCO']}).to_csv(root / 'input.csv', index=False)
            config = AppConfig(root, job_timeout=30)
            with JobManager(config) as manager:
                job = manager.submit('admet', {'base_input_path': str(root), 'input_file': 'input.csv'})
                completed = self.wait_for(manager, job['id'], lambda j: j['status'] in TERMINAL)
                self.assertEqual(completed['status'], 'succeeded', completed['error'])
                self.assertEqual(completed['result']['details']['rows'], 1)
                self.assertTrue(any(path.endswith('_egg.png') for path in completed['result']['artifacts']))
                failed = manager.submit('admet', {'base_input_path': str(root / 'missing')})
                failed = self.wait_for(manager, failed['id'], lambda j: j['status'] in TERMINAL)
                self.assertEqual(failed['status'], 'failed')
                self.assertIn('FileNotFoundError', failed['error'])
                with self.assertRaises(RuntimeError):
                    JobManager(config)
            from biomolexplorer.diagnostics import log_directory
            central=log_directory()/'jobs'/failed['id']
            self.assertIn('FileNotFoundError',(central/'execution.log').read_text())
            self.assertIn('FileNotFoundError',(central/'backend.log').read_text())
            self.assertIn('FileNotFoundError',(central/'errors.log').read_text())
            with JobManager(config) as manager:
                self.assertEqual(manager.get(job['id'])['status'], 'succeeded')
                self.assertEqual(len(manager.list()), 2)

    def test_cancel_running_and_queued_jobs_and_timeout(self):
        original = subprocess.Popen
        def slow_worker(args, **kwargs):
            return original([sys.executable, '-c', 'import time; time.sleep(60)'], **kwargs)
        with tempfile.TemporaryDirectory() as temp, patch('biomolexplorer.jobs.subprocess.Popen', side_effect=slow_worker):
            config = AppConfig(Path(temp), job_timeout=0.4)
            with JobManager(config) as manager:
                running = manager.submit('admet', {'base_input_path': temp})
                self.wait_for(manager, running['id'], lambda j: j['status'] == 'running')
                queued = manager.submit('admet', {'base_input_path': temp})
                self.assertEqual(manager.cancel(queued['id'])['status'], 'cancelled')
                self.assertEqual(manager.cancel(running['id'])['status'], 'cancelled')
                timed = manager.submit('admet', {'base_input_path': temp})
                timed = self.wait_for(manager, timed['id'], lambda j: j['status'] in TERMINAL)
                self.assertEqual(timed['status'], 'failed')
                self.assertIn('TimeoutError', timed['error'])

    def test_restart_marks_unfinished_jobs_interrupted(self):
        with tempfile.TemporaryDirectory() as temp:
            config = AppConfig(Path(temp))
            with JobManager(config) as manager:
                with manager._connect() as connection:
                    connection.execute('INSERT INTO jobs VALUES (?,?,?,?,?,?,?,?,?,?)',
                        ('interrupted', 'admet', '{}', 'running', 'date', 'date', temp, temp, None, None))
            # close cancels outstanding records; emulate a crash after shutdown.
            import sqlite3
            with sqlite3.connect(config.state_dir / 'jobs.sqlite3') as connection:
                connection.execute("UPDATE jobs SET status='running' WHERE id='interrupted'")
            with JobManager(config) as manager:
                self.assertEqual(manager.get('interrupted')['status'], 'interrupted')


if __name__ == '__main__':
    unittest.main()
