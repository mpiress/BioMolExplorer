"""Large scientific jobs run to completion without implicit wall-clock deadlines."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import Mock,patch
from types import SimpleNamespace

from biomolexplorer.config import AppConfig
from biomolexplorer.jobs import JobManager,TERMINAL
from biomolexplorer.processes import run_command
from biomolexplorer.catalog import new_stage
from biomolexplorer.operations import validate_operation
from caad.graph_results import common_fragment
from biomolexplorer.docking_data import convert_structure,structure_smiles

class UnboundedExecutionTests(unittest.TestCase):
    def test_conversion_commands_have_no_deadline_and_keep_output_validation(self):
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'input.sdf';source.write_text('structure')
            destination=Path(temp)/'output.mol2'
            def convert(command,**options):
                self.assertIsNone(options.get('timeout'))
                destination.write_text('converted structure')
                return SimpleNamespace(returncode=0,stderr='')
            with patch('biomolexplorer.docking_data.subprocess.run',side_effect=convert):
                self.assertEqual(convert_structure(source,destination),destination)
            with patch('biomolexplorer.docking_data.shutil.which',return_value='obabel'),patch('biomolexplorer.docking_data.subprocess.run',return_value=SimpleNamespace(returncode=0,stdout='CCO\tinput')) as run:
                self.assertEqual(structure_smiles(source),'CCO')
                self.assertIsNone(run.call_args.kwargs.get('timeout'))
            with patch('biomolexplorer.docking_data.subprocess.run',return_value=SimpleNamespace(returncode=1,stderr='bad input')):
                with self.assertRaisesRegex(ValueError,'bad input'):convert_structure(source,destination)

    def test_command_outlives_legacy_environment_limit_and_still_reports_real_failures(self):
        with patch.dict(os.environ,{'BIOMOL_COMMAND_TIMEOUT':'0.01'}):
            self.assertTrue(run_command([sys.executable,'-c','import time;time.sleep(.1)']))
            with self.assertRaisesRegex(RuntimeError,'código 2'):
                run_command([sys.executable,'-c','raise SystemExit(2)'])

    def test_default_supervisor_waits_without_deadline(self):
        real_popen=subprocess.Popen;processes=[]
        def delayed_worker(args,**options):
            code='import sys,time,pathlib;time.sleep(.15);pathlib.Path(sys.argv[1]).write_text(\'{"artifacts":[],"details":{}}\')'
            process=real_popen([sys.executable,'-c',code,args[-1]],**options)
            process.wait=Mock(wraps=process.wait);processes.append(process)
            return process
        with tempfile.TemporaryDirectory() as temp,patch('biomolexplorer.jobs.subprocess.Popen',side_effect=delayed_worker):
            config=AppConfig(Path(temp));self.assertIsNone(config.job_timeout)
            with JobManager(config) as manager:
                job=manager.submit('admet',{'base_input_path':temp})
                deadline=time.monotonic()+5
                while job['status'] not in TERMINAL and time.monotonic()<deadline:
                    time.sleep(.02);job=manager.get(job['id'])
                self.assertEqual(job['status'],'succeeded',job['error'])
                processes[0].wait.assert_called_once_with(timeout=None)

    def test_cancellation_stops_an_unbounded_worker(self):
        real_popen=subprocess.Popen
        def long_worker(args,**options):return real_popen([sys.executable,'-c','import time;time.sleep(60)'],**options)
        with tempfile.TemporaryDirectory() as temp,patch('biomolexplorer.jobs.subprocess.Popen',side_effect=long_worker):
            with JobManager(AppConfig(Path(temp))) as manager:
                job=manager.submit('admet',{'base_input_path':temp})
                deadline=time.monotonic()+5
                while job['status']=='queued' and time.monotonic()<deadline:
                    time.sleep(.02);job=manager.get(job['id'])
                self.assertEqual(job['status'],'running')
                self.assertEqual(manager.cancel(job['id'])['status'],'cancelled')

    def test_cancellation_stops_child_that_ignores_term_after_worker_exits(self):
        import signal
        real_popen=subprocess.Popen;processes=[]
        with tempfile.TemporaryDirectory() as temp:
            ready=Path(temp)/'child.pid'
            child='import os,signal,time,pathlib;signal.signal(signal.SIGTERM,signal.SIG_IGN);pathlib.Path('+repr(str(ready))+').write_text(str(os.getpid()));time.sleep(60)'
            worker='import subprocess,sys,time;subprocess.Popen([sys.executable,"-c",'+repr(child)+']);time.sleep(60)'
            def launch(args,**options):
                process=real_popen([sys.executable,'-c',worker],**options)
                processes.append(process);return process
            try:
                with patch('biomolexplorer.jobs.subprocess.Popen',side_effect=launch),JobManager(AppConfig(Path(temp))) as manager:
                    job=manager.submit('admet',{'base_input_path':temp})
                    deadline=time.monotonic()+5
                    while not ready.exists() and time.monotonic()<deadline:time.sleep(.02)
                    self.assertTrue(ready.exists(),'Worker child did not start')
                    pid=int(ready.read_text())
                    self.assertEqual(manager.cancel(job['id'])['status'],'cancelled')
                    deadline=time.monotonic()+3
                    def running():
                        stat=Path('/proc')/str(pid)/'stat'
                        try:return stat.read_text().split()[2] not in ('Z','X')
                        except FileNotFoundError:return False
                    while running() and time.monotonic()<deadline:time.sleep(.02)
                    self.assertFalse(running(),'Scientific child survived cancellation')
            finally:
                for process in processes:
                    try:os.killpg(process.pid,signal.SIGKILL)
                    except ProcessLookupError:pass
                    process.wait(timeout=5)

    def test_graph_fragment_default_is_unlimited_and_explicit_limits_remain_valid(self):
        stage=new_stage('graphs');self.assertEqual(stage['parameters']['mcs_timeout'],0)
        validate_operation('graphs',{'mcs_timeout':0})
        validate_operation('graphs',{'mcs_timeout':1})
        result=common_fragment(['CCO','CCCO'])
        self.assertEqual(result['timeout'],0);self.assertEqual(result['status'],'complete')
