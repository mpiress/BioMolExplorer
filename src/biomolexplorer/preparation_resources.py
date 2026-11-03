"""Conservative scheduling for independent classic-Chimera ligand preparations.

Resource policy is not a scientific setting. Linux file leases coordinate the
application's workers for the same OS user, and survive neither crash nor kill.
"""
from contextlib import contextmanager
import math
import os
from pathlib import Path
import tempfile
import time

GIB = 1024 ** 3


def _text(path):
    try:return Path(path).read_text().strip()
    except OSError:return ''


def _groups():
    root=Path('/sys/fs/cgroup')
    groups=[root]
    for line in _text('/proc/self/cgroup').splitlines():
        if line.startswith('0::'):
            relative=Path(line[3:].lstrip('/'))
            if '..' in relative.parts:continue
            path=root/relative
            while path!=root:
                groups.append(path);path=path.parent
    return groups


def available_cpus():
    try:count=len(os.sched_getaffinity(0))
    except (AttributeError,OSError):count=os.cpu_count() or 1
    for group in _groups():
        quota=_text(str(group/'cpu.max')).split()
        if len(quota)==2 and quota[0]!='max':
            try:count=min(count,max(1,math.floor(int(quota[0])/int(quota[1]))))
            except (ValueError,ZeroDivisionError):pass
    return max(1,count)


def available_memory():
    """Available bytes, respecting a cgroup-v2 limit when one is present."""
    values={}
    for line in _text('/proc/meminfo').splitlines():
        fields=line.split()
        if len(fields)>=2:
            try:values[fields[0].rstrip(':')]=int(fields[1])*1024
            except ValueError:pass
    available=values.get('MemAvailable')
    for group in _groups():
        limit,current=_text(str(group/'memory.max')),_text(str(group/'memory.current'))
        if limit and limit!='max' and current:
            try:
                free=max(0,int(limit)-int(current));available=min(available,free) if available is not None else free
            except ValueError:pass
    return available


def capacity(requested=None):
    """Reserve half the CPUs and >=1 GiB RAM; estimate 1 GiB per Chimera.

    These are admission estimates, not an OS memory limit. Unknown memory or
    a machine below the reserve runs serially rather than guessing capacity.
    """
    from .paths import worker_count
    requested=worker_count() if requested is None else requested
    override=os.environ.get('BIOMOL_PREPARATION_WORKERS')
    if override:
        try:
            value=int(override)
            if value<1:raise ValueError()
            requested=min(requested,value)
        except ValueError:raise ValueError('BIOMOL_PREPARATION_WORKERS must be a positive integer.') from None
    cpus=available_cpus();memory=available_memory()
    memory_slots=max(1,(memory-max(GIB,memory//4))//GIB) if memory is not None else 1
    return max(1,min(4,requested,max(1,cpus//2),memory_slots))


def single_thread_environment():
    environment=os.environ.copy()
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        environment[name]='1'
    return environment


@contextmanager
def chimera_slot():
    """Acquire a bounded, per-user slot immediately before launching Chimera.

    Descriptors are not inherited by native subprocesses. Worker process-group
    cancellation releases the OS lock without deleting or stealing live slots.
    """
    import fcntl
    folder=Path(tempfile.gettempdir())/('biomolexplorer-chimera-'+str(os.getuid()))
    folder.mkdir(mode=0o700,exist_ok=True)
    if folder.is_symlink() or folder.stat().st_uid!=os.getuid():
        raise OSError('Unsafe Chimera scheduling directory: '+str(folder))
    held=None
    try:
        while held is None:
            for index in range(capacity()):
                fd=os.open(folder/(str(index)+'.lock'),os.O_CREAT|os.O_RDWR|os.O_NOFOLLOW,0o600)
                stream=os.fdopen(fd,'r+')
                try:fcntl.flock(stream,fcntl.LOCK_EX|fcntl.LOCK_NB)
                except BlockingIOError:stream.close();continue
                held=stream;break
            if held is None:time.sleep(.1)
        yield
    finally:
        if held is not None:
            fcntl.flock(held,fcntl.LOCK_UN);held.close()
