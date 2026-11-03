"""Resource admission and cross-worker leases, independent of chemistry mocks."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from biomolexplorer import preparation_resources as resources


class CapacityTests(unittest.TestCase):
    def test_cpu_reserve_memory_reserve_and_explicit_caps(self):
        cases=[(8,8,16,4),(8,2,16,1),(8,16,2,1),(8,16,3,2),(2,16,16,2)]
        with patch.dict(os.environ,{},clear=True):
            for requested,cpus,gib,expected in cases:
                with self.subTest(cpus=cpus,gib=gib),patch.object(resources,'available_cpus',return_value=cpus),patch.object(resources,'available_memory',return_value=gib*resources.GIB):
                    self.assertEqual(resources.capacity(requested),expected)
            with patch.object(resources,'available_cpus',return_value=16),patch.object(resources,'available_memory',return_value=None):
                self.assertEqual(resources.capacity(8),1)

    def test_environment_can_reduce_but_never_expand_cpu_budget(self):
        with patch.object(resources,'available_cpus',return_value=32),patch.object(resources,'available_memory',return_value=32*resources.GIB):
            for value,expected in [('1',1),('4',2)]:
                with patch.dict(os.environ,{'BIOMOL_PREPARATION_WORKERS':value}):self.assertEqual(resources.capacity(2),expected)
            for invalid in ('0','-1','invalid'):
                with patch.dict(os.environ,{'BIOMOL_PREPARATION_WORKERS':invalid}),self.assertRaises(ValueError):resources.capacity(2)

    def test_cgroup_limit_and_affinity_are_respected(self):
        data={'/sys/fs/cgroup/cpu.max':'150000 100000','/proc/meminfo':'MemAvailable: 8388608 kB',
              '/sys/fs/cgroup/memory.max':str(4*resources.GIB),'/sys/fs/cgroup/memory.current':str(3*resources.GIB)}
        with patch.object(resources,'_text',side_effect=lambda p:data.get(p,'')),patch.object(os,'sched_getaffinity',return_value=set(range(16))):
            self.assertEqual(resources.available_cpus(),1)
            self.assertEqual(resources.available_memory(),resources.GIB)
            data.update({'/proc/self/cgroup':'0::/users/study','/sys/fs/cgroup/cpu.max':'max 100000',
                         '/sys/fs/cgroup/users/cpu.max':'300000 100000',
                         '/sys/fs/cgroup/users/memory.max':str(2*resources.GIB),
                         '/sys/fs/cgroup/users/memory.current':str(3*resources.GIB//2)})
            self.assertEqual(resources.available_cpus(),3)
            self.assertEqual(resources.available_memory(),resources.GIB//2)

    def test_native_thread_limits_do_not_mutate_parent_environment(self):
        with patch.dict(os.environ,{'OMP_NUM_THREADS':'8'}):
            child=resources.single_thread_environment()
            self.assertEqual(os.environ['OMP_NUM_THREADS'],'8')
            self.assertTrue(all(child[k]=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')))


class LeaseTests(unittest.TestCase):
    def test_distinct_workers_share_two_slots_and_crash_releases_them(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            setup="from biomolexplorer import preparation_resources as r; import os,time,json; from pathlib import Path; r.tempfile.gettempdir=lambda:"+repr(folder)+"; r.capacity=lambda:2; "
            code=setup+"\nwith r.chimera_slot():\n with open("+repr(str(root/'events.jsonl'))+",'a') as f:f.write(json.dumps(['start',os.getpid(),time.monotonic()])+'\\n')\n time.sleep(.25)\n with open("+repr(str(root/'events.jsonl'))+",'a') as f:f.write(json.dumps(['end',os.getpid(),time.monotonic()])+'\\n')\n"
            children=[subprocess.Popen([sys.executable,'-c',code]) for _ in range(4)]
            try:
                for child in children:self.assertEqual(child.wait(timeout=10),0)
            finally:
                for child in children:
                    if child.poll() is None:child.kill();child.wait()
            active=maximum=0
            for kind,pid,stamp in sorted([json.loads(v) for v in (root/'events.jsonl').read_text().splitlines()],key=lambda v:v[2]):
                active+=1 if kind=='start' else -1;maximum=max(maximum,active)
            self.assertEqual(maximum,2);self.assertEqual(active,0)
            ready=root/'ready'
            hold=setup.replace('lambda:2','lambda:1')+'\nwith r.chimera_slot():\n Path('+repr(str(ready))+').write_text("ready")\n time.sleep(60)\n'
            child=subprocess.Popen([sys.executable,'-c',hold],start_new_session=True)
            try:
                import time
                deadline=time.monotonic()+5
                while not ready.exists() and time.monotonic()<deadline:time.sleep(.02)
                self.assertTrue(ready.exists())
                os.killpg(child.pid,signal.SIGTERM);child.wait(timeout=5)
                release=setup.replace('lambda:2','lambda:1')+'\nwith r.chimera_slot():pass\n'
                self.assertEqual(subprocess.run([sys.executable,'-c',release],timeout=5).returncode,0)
            finally:
                if child.poll() is None:os.killpg(child.pid,signal.SIGKILL);child.wait()
