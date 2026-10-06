"""Correlation, privacy, process-safe rotation and diagnostic reports."""
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import threading
import unittest
from unittest.mock import patch
from biomolexplorer.diagnostics import (get_logger, log_context, event, write_summary,
    diagnose_exception)
from biomolexplorer.log_report import read_events, main
from biomolexplorer.processes import run_command, ScientificToolError
from kernel.loggers import LoggerManager


class StructuredDiagnosticsTests(unittest.TestCase):
    def test_context_is_scoped_and_errors_not_duplicated(self):
        with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ, {'BIOMOL_LOG_DIR': folder}):
            logger = get_logger('context_test')
            with log_context(job_id='job-a', operation='redocking'):
                with log_context(pair='4M0E|1YL|604|A'):
                    logger.error('Failed token=private-value\nsecond line')
                logger.info('Outside pair')
            logger.info('Outside job')
            records=read_events(folder,level='INFO')
            self.assertEqual(len(records),3)
            self.assertEqual(records[0]['job_id'],'job-a')
            self.assertNotIn('pair',records[1]);self.assertNotIn('job_id',records[2])
            self.assertEqual(records[0]['message'],'Failed token=[REDACTED]\nsecond line')
            self.assertEqual((Path(folder)/'errors.log').read_text().count('second line'),1)
            self.assertNotIn('private-value',(Path(folder)/'context_test.log').read_text())

    def test_threads_have_independent_context(self):
        with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ, {'BIOMOL_LOG_DIR':folder}):
            logger=get_logger('thread_test');barrier=threading.Barrier(2)
            def emit(identity):
                with log_context(job_id=identity):
                    barrier.wait();event(logger,'thread.test',identity)
            threads=[threading.Thread(target=emit,args=(i,)) for i in ('one','two')]
            for thread in threads:thread.start()
            for thread in threads:thread.join()
            self.assertEqual({r['message']:r['job_id'] for r in read_events(folder,level='INFO')},{'one':'one','two':'two'})

    def test_legacy_science_shares_json_and_exception_cause(self):
        with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ,{'BIOMOL_LOG_DIR':folder}):
            logger=LoggerManager.get_logger('legacy-diagnostic-test','logs/docking.log')
            try:
                try:raise FileNotFoundError('receptor.pdb')
                except FileNotFoundError as exc:raise RuntimeError('Preparation failed') from exc
            except RuntimeError as exc:
                logger.exception('Preparation failure',extra=diagnose_exception(exc))
                write_summary('failed',error=exc,job_id='job-1')
            row=read_events(folder,level='ERROR')[0]
            self.assertEqual(row['error_code'],'INPUT_NOT_FOUND')
            self.assertIn('FileNotFoundError',row['exception']['traceback'])
            self.assertEqual(json.loads((Path(folder)/'diagnostic.json').read_text())['status'],'failed')

    def test_command_failure_and_timeout_have_codes_and_context(self):
        with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ,{'BIOMOL_LOG_DIR':folder}):
            with log_context(job_id='tool-job'):
                with self.assertRaises(ScientificToolError):
                    run_command([sys.executable,'-c','import sys;sys.stderr.write("bad input");sys.exit(3)'])
                with self.assertRaises(subprocess.TimeoutExpired):
                    run_command([sys.executable,'-c','import time;time.sleep(5)'],timeout=.05)
            failures=[r for r in read_events(folder,level='INFO') if r['event']=='tool.failed']
            self.assertEqual([r['error_code'] for r in failures],['TOOL_EXIT_FAILED','TOOL_TIMEOUT'])
            self.assertTrue(all(r['job_id']=='tool-job' and r['command_id'] and r['duration_ms']>=0 for r in failures))
            self.assertEqual(failures[0]['returncode'],3)

    def test_filter_report_and_rotation_tolerate_partial_records(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder);base={'timestamp':'2026-10-06T00:00:00Z','level':'ERROR','event':'tool.failed','job_id':'a','message':'bad'}
            (path/'events.jsonl.1').write_text(json.dumps(base)+'\n')
            (path/'events.jsonl').write_text(json.dumps({**base,'job_id':'b'})+'\n{"interrupted":')
            self.assertEqual([r['job_id'] for r in read_events(folder,job_id='a')],['a'])
            self.assertEqual([r['job_id'] for r in read_events(folder,limit=1)],['b'])
            with patch('builtins.print') as output:
                self.assertEqual(main(['--directory',folder,'--json']),0)
                self.assertEqual(len(json.loads(output.call_args.args[0])),2)

    def test_process_pool_shared_rotation_writes_valid_json_once(self):
        with tempfile.TemporaryDirectory() as folder:
            script='''import sys,logging
from pathlib import Path
from biomolexplorer.diagnostics import DeferredRotatingHandler,JsonFormatter,ContextFilter
logger=logging.getLogger('rotation');logger.setLevel(logging.INFO)
h=DeferredRotatingHandler(Path(sys.argv[1])/'events.jsonl',maxBytes=2500,backupCount=20,delay=True,encoding='utf-8')
h.setFormatter(JsonFormatter());h.addFilter(ContextFilter());logger.addHandler(h)
for i in range(15):logger.info('%s-%s',sys.argv[2],i)
'''
            env={**os.environ,'PYTHONPATH':str(Path(__file__).resolve().parents[1]/'src')}
            processes=[subprocess.Popen([sys.executable,'-c',script,folder,str(i)],env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE) for i in range(3)]
            for process in processes:
                stdout,stderr=process.communicate(timeout=10)
                self.assertEqual(process.returncode,0,stderr.decode())
                self.assertNotIn(b'Logging error',stderr)
            rows=[json.loads(line) for p in Path(folder).glob('events.jsonl*') if not p.name.endswith('.lock') for line in p.read_text().splitlines()]
            self.assertEqual(len(rows),45)
            self.assertEqual(len({r['message'] for r in rows}),45)

    def test_real_worker_keeps_pipeline_context_and_failure_category(self):
        from biomolexplorer.config import AppConfig
        from biomolexplorer.jobs import JobManager, TERMINAL
        with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ, {'BIOMOL_LOG_DIR': str(Path(folder)/'logs')}):
            with JobManager(AppConfig(Path(folder)/'workspace', job_timeout=20)) as manager:
                job=manager.submit('admet',{'base_input_path':str(Path(folder)/'missing')},
                                   diagnostic_context={'project_id':'project-1','run_id':'run-1','stage_id':'stage-1'})
                deadline=time.monotonic()+20
                while time.monotonic()<deadline:
                    job=manager.get(job['id'])
                    if job['status'] in TERMINAL:break
                    time.sleep(.05)
                self.assertEqual(job['status'],'failed')
                self.assertEqual(job['result']['diagnostic']['error_code'],'INPUT_NOT_FOUND')
            directory=Path(folder)/'logs/jobs'/job['id']
            summary=json.loads((directory/'diagnostic.json').read_text())
            self.assertEqual(summary['run_id'],'run-1')
            self.assertEqual(summary['error_code'],'INPUT_NOT_FOUND')
            failure=next(r for r in read_events(directory,level='ERROR') if r['event']=='worker.failed')
            self.assertEqual(failure['project_id'],'project-1')
            self.assertEqual(failure['stage_id'],'stage-1')
            self.assertEqual(failure['job_id'],job['id'])
            self.assertEqual(failure['operation'],'admet')
