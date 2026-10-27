"""Version-gated docking capabilities; unsupported methods never reach templates."""
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import re
import subprocess


def _version(value):
    match=re.search(r'(\d+)\.(\d+)(?:\.(\d+))?',value or '')
    return tuple(int(v or 0) for v in match.groups()) if match else None


@lru_cache(maxsize=8)
def dock6_version(executable):
    path=Path(executable)
    try:
        result=subprocess.run([str(path),'--version'],capture_output=True,text=True,timeout=15,check=False)
    except (OSError,subprocess.SubprocessError):return None
    match=re.search(r'DOCK\s+v(\d+\.\d+(?:\.\d+)?)',result.stdout+'\n'+result.stderr,re.I)
    return match.group(1) if match else None


def dock6_capabilities(root=None):
    from .docking_tools import dock6_root
    root=dock6_root(root)
    executable=root/'bin/dock6' if root else None
    version=dock6_version(str(executable)) if executable else None
    parsed=_version(version)
    return dict(backend='DOCK6',version=version,available=bool(executable),standard_docking=bool(executable),
        covalent_attach_and_grow=bool(parsed and parsed>=(6,12,0)),
        covalent_requirement='DOCK >= 6.12')


def vina_capabilities():
    try:
        import vina
        version=getattr(vina,'__version__',None)
    except ImportError:version=None
    return dict(backend='AutoDock Vina',version=version,available=version is not None,
        standard_docking=version is not None,covalent_docking=False,
        reason='No validated native covalent-docking protocol is enabled for Vina 1.2.5.')


@dataclass(frozen=True)
class StandardDocking:
    capabilities: dict

    def ensure_available(self):
        if not self.capabilities.get('standard_docking'):raise ValueError('Standard docking backend is unavailable.')


@dataclass(frozen=True)
class CovalentDocking:
    capabilities: dict

    def ensure_available(self):
        if not self.capabilities.get('covalent_attach_and_grow'):
            raise ValueError('Covalent Attach-and-Grow requires DOCK >= 6.12; the installed backend must not receive conformer_search_type covalent.')

    def conformer_search_type(self):
        self.ensure_available()
        return 'covalent'
