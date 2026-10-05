"""Persistent local jobs, with scientific work isolated in child processes."""
import json
import os
import signal
import shutil
import sqlite3
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

from .operations import validate_operation
from .paths import SOURCE_ROOT
from .diagnostics import log_directory, get_logger
from .progress import read_progress

TERMINAL = {'succeeded', 'failed', 'cancelled', 'interrupted'}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


class JobManager:
    """One manager per workspace; suitable for a desktop/local Flet backend.

    submit/get/cancel expose JSON-compatible dictionaries. UI code polls get()
    asynchronously; it never executes scientific work on the UI thread.
    """
    def __init__(self, config):
        self.config = config
        config.state_dir.mkdir(parents=True, exist_ok=True)
        self.database = config.state_dir / 'jobs.sqlite3'
        self._lock = threading.RLock()
        self._processes = {}
        self._closed = False
        self._owner_file = config.state_dir / 'manager.lock'
        # OS lock is released on crashes, and prevents two supervisors writing
        # the same workspace. Backend currently targets the Linux scientific stack.
        import fcntl
        self._owner = self._owner_file.open('a+')
        try:
            fcntl.flock(self._owner, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            self._owner.close()
            raise RuntimeError('A JobManager is already using this workspace') from None
        try:
            with self._connect() as connection:
                connection.execute('''CREATE TABLE IF NOT EXISTS jobs (
                    id TEXT PRIMARY KEY, operation TEXT NOT NULL, parameters TEXT NOT NULL,
                    status TEXT NOT NULL, created_at TEXT NOT NULL, updated_at TEXT NOT NULL,
                    output_path TEXT NOT NULL, log_path TEXT NOT NULL, result TEXT, error TEXT)''')
                connection.execute("UPDATE jobs SET status='interrupted', error=?, updated_at=? WHERE status IN ('queued','running')",
                                   ('Previous supervisor stopped; submit a new job to resume using cached files', utc_now()))
            self._executor = ThreadPoolExecutor(max_workers=config.max_jobs, thread_name_prefix='biomol-job')
        except Exception:
            self._owner.close()
            raise

    @contextmanager
    def _connect(self):
        # Each calling thread gets its own connection, with short transactions.
        connection = sqlite3.connect(self.database, timeout=30)
        try:
            with connection:
                yield connection
        finally:
            connection.close()

    def _update(self, job_id, **values):
        values['updated_at'] = utc_now()
        with self._connect() as connection:
            connection.execute('UPDATE jobs SET ' + ', '.join(f'{name}=?' for name in values) + ' WHERE id=?',
                               (*values.values(), job_id))

    def submit(self, operation, parameters):
        validate_operation(operation, parameters)
        encoded = json.dumps(parameters, allow_nan=False)
        with self._lock:
            if self._closed:
                raise RuntimeError('JobManager is closed')
            with self._connect() as connection:
                queued = connection.execute("SELECT COUNT(*) FROM jobs WHERE status='queued'").fetchone()[0]
            if queued >= self.config.max_queued_jobs:
                raise RuntimeError('Job queue is full')
            job_id = uuid4().hex
            job_dir = self.config.state_dir / 'jobs' / job_id
            output = job_dir / 'artifacts'
            output.mkdir(parents=True)
            now = utc_now()
            with self._connect() as connection:
                connection.execute('INSERT INTO jobs VALUES (?,?,?,?,?,?,?,?,?,?)',
                    (job_id, operation, encoded, 'queued', now, now, str(output),
                     str(job_dir / 'execution.log'), None, None))
            try:
                self._executor.submit(self._run, job_id)
            except Exception as exc:
                self._update(job_id, status='failed', error=f'Could not schedule worker: {exc}')
                raise
        return self.get(job_id)

    def get(self, job_id):
        with self._connect() as connection:
            connection.row_factory = sqlite3.Row
            row = connection.execute('SELECT * FROM jobs WHERE id=?', (job_id,)).fetchone()
        if row is None:
            raise KeyError(job_id)
        job = dict(row)
        for name in ('parameters', 'result'):
            job[name] = json.loads(job[name]) if job[name] is not None else None
        job['progress'] = read_progress(Path(job['output_path']).parent / 'progress.json')
        return job

    def list(self, limit=100):
        if type(limit) is not int or not 1 <= limit <= 1000:
            raise ValueError('limit must be between 1 and 1000')
        with self._connect() as connection:
            ids = connection.execute('SELECT id FROM jobs ORDER BY created_at DESC LIMIT ?', (limit,)).fetchall()
        return [self.get(row[0]) for row in ids]

    @staticmethod
    def _stop(process):
        if process.poll() is None:
            # Stop nested docking commands/pools as well as the Python worker.
            try:
                os.killpg(process.pid, signal.SIGTERM)
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait(timeout=5)
            except ProcessLookupError:
                pass

    def cancel(self, job_id):
        with self._lock:
            job = self.get(job_id)
            if job['status'] in TERMINAL:
                return job
            self._update(job_id, status='cancelled', error='Cancelled by user')
            process = self._processes.get(job_id)
            if process is not None:
                self._stop(process)
        return self.get(job_id)

    def _run(self, job_id):
        process = None
        job = None
        try:
            with self._lock:
                job = self.get(job_id)
                if job['status'] != 'queued':
                    return
                job_dir = Path(job['output_path']).parent
                request = job_dir / 'request.json'
                response = job_dir / 'result.json'
                request.write_text(json.dumps({'operation': job['operation'], 'parameters': job['parameters'],
                                              'output_path': job['output_path']}))
                environment = dict(os.environ)
                environment['PYTHONPATH'] = str(SOURCE_ROOT) + os.pathsep + environment.get('PYTHONPATH', '')
                environment['BIOMOL_WORKSPACE'] = str(self.config.workspace)
                environment['BIOMOL_CPU_WORKERS'] = str(self.config.cpu_workers)
                environment['BIOMOL_LOG_DIR'] = str(log_directory() / 'jobs' / job_id)
                environment['MPLBACKEND'] = 'Agg'
                environment['MPLCONFIGDIR'] = str(job_dir / 'matplotlib')
                environment['BIOMOL_CACHE_DIR'] = str(self.config.state_dir / 'cache')
                environment['BIOMOL_PROGRESS_FILE'] = str(job_dir / 'progress.json')
                if self.config.resource_dir:
                    environment['BIOMOL_RESOURCE_DIR'] = str(self.config.resource_dir)
                else:
                    environment.pop('BIOMOL_RESOURCE_DIR', None)
                environment['OMP_NUM_THREADS'] = '1'
                environment['OPENBLAS_NUM_THREADS'] = '1'
                self._update(job_id, status='running')
                with Path(job['log_path']).open('w') as log:
                    process = subprocess.Popen([str(self.config.worker_python or sys.executable), '-m', 'biomolexplorer.worker',
                                                str(request), str(response)],
                        cwd=self.config.workspace, env=environment, stdout=log, stderr=log, start_new_session=True)
                self._processes[job_id] = process
            try:
                process.wait(timeout=self.config.job_timeout)
            except subprocess.TimeoutExpired:
                self._stop(process)
                raise TimeoutError(f'Job exceeded {self.config.job_timeout} seconds') from None
            with self._lock:
                if self.get(job_id)['status'] == 'cancelled':
                    return
                if process.returncode != 0:
                    error = json.loads(response.read_text()).get('error') if response.exists() else 'Worker failed; see execution.log'
                    get_logger('backend').error('Worker falhou; job=%s operation=%s error=%s diagnostics=%s',job_id,job['operation'],error,environment['BIOMOL_LOG_DIR'])
                    self._update(job_id, status='failed', error=error)
                else:
                    result = json.loads(response.read_text())
                    self._update(job_id, status='succeeded', result=json.dumps(result), error=None)
        except Exception as exc:
            get_logger('backend').exception('Falha ao supervisionar job=%s',job_id)
            with self._lock:
                if self.get(job_id)['status'] != 'cancelled':
                    self._update(job_id, status='failed', error=f'{type(exc).__name__}: {exc}')
        finally:
            if process is not None and process.poll() is None:
                self._stop(process)
            if job and Path(job['log_path']).is_file():
                try:
                    destination=log_directory()/'jobs'/job_id/'execution.log'
                    destination.parent.mkdir(parents=True,exist_ok=True)
                    shutil.copyfile(job['log_path'],destination)
                except OSError:
                    get_logger('backend').exception('Falha ao copiar log de execução; job=%s',job_id)
            with self._lock:
                self._processes.pop(job_id, None)

    def close(self):
        with self._lock:
            if self._closed:
                return
            self._closed = True
            with self._connect() as connection:
                pending = connection.execute("SELECT id FROM jobs WHERE status IN ('queued', 'running')").fetchall()
            for (job_id,) in pending:
                self.cancel(job_id)
        self._executor.shutdown(wait=True)
        self._owner.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
