"""Discover the local DOCK6 installation shared by UI validation and workers."""
import os
import shutil
import sys
from pathlib import Path


def dock6_root(configured=None, worker_python=None):
    if configured:
        return Path(configured).expanduser().resolve()
    candidates=[]
    if os.environ.get('BIOMOL_DOCK6_ROOT'):
        candidates.append(Path(os.environ['BIOMOL_DOCK6_ROOT']).expanduser())
    search_path=str(Path(worker_python or sys.executable).parent)+os.pathsep+os.environ.get('PATH','')
    executable=shutil.which('dock6',path=search_path)
    if executable:candidates.append(Path(executable).resolve().parent.parent)
    candidates.extend((Path.home()/'progs/dock6',Path.home()/'dock6',Path('/opt/dock6'),Path('/usr/local/dock6')))
    for candidate in candidates:
        root=candidate.resolve()
        if all((root/'bin'/name).is_file() and os.access(root/'bin'/name,os.X_OK)
               for name in ('dock6','grid','sphgen','sphere_selector','showbox')) and (root/'parameters/vdw_AMBER_parm99.defn').is_file():
            return root
    return None
