"""Check that the next BioMolExplorer UI start imports the current checkout.

Run with the Python that owns the UI launcher, from any directory:
    /path/to/env/bin/python scripts/check_ui_installation.py
"""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
MODULES = (
    'biomolexplorer.ui.app',
    'biomolexplorer.ui.guided',
    'biomolexplorer.ui.pose_settings',
    'biomolexplorer.ui.selected_items_table',
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    probe = ('import importlib,json; names='+repr(MODULES)+'; '
             'print(json.dumps({n:importlib.import_module(n).__file__ for n in names}))')
    with tempfile.TemporaryDirectory() as outside:
        result = subprocess.run([sys.executable,'-c',probe],cwd=outside,
                                capture_output=True,text=True,check=False)
    if result.returncode:
        print(result.stderr.strip(),file=sys.stderr)
        return 1
    paths = json.loads(result.stdout)
    failed = False
    for name,installed in paths.items():
        source = ROOT/'src'/Path(*name.split('.')).with_suffix('.py')
        equal = source.is_file() and digest(source) == digest(installed)
        print(('OK' if equal else 'DIVERGENTE'),name,installed)
        failed |= not equal
    launch = Path(sys.executable).parent/'biomolexplorer-ui'
    entry = next((e for e in importlib.metadata.distribution('biomolexplorer').entry_points
                  if e.name == 'biomolexplorer-ui'),None)
    launcher_ok = (launch.is_file() and entry is not None
                   and entry.value == 'biomolexplorer.ui.app:main'
                   and Path(launch.open(encoding='utf-8').readline().strip()[2:]).resolve()
                   == Path(sys.executable).resolve())
    print(('OK' if launcher_ok else 'DIVERGENTE'),'launcher',launch)
    print('Python:',sys.executable)
    print('Reinicie o processo da aplicação para carregar os módulos atuais.')
    return 1 if failed or not launcher_ok else 0


if __name__ == '__main__':
    raise SystemExit(main())
