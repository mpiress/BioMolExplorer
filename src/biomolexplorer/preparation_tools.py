"""Compatibility API restricted to Chimera bundled charge utilities."""
from dataclasses import dataclass
import os
from pathlib import Path
import shutil
from .processes import ScientificToolError,run_command
@dataclass(frozen=True)
class AmberTool:
    name:str
    executable:Path
    environment:str
    prefix:Path

def resolve_amber_tool(name,prefix='',environment=''):
    chimera=shutil.which('chimera')
    if not chimera:raise ScientificToolError('UCSF Chimera is required for molecular preparation.','TOOL_NOT_FOUND','Make classic UCSF Chimera available on PATH.')
    installation=Path(chimera).resolve().parent.parent
    if name not in ('antechamber','sqm','am1bcc','atomtype','bondtype','parmchk2'):
        raise ScientificToolError('The retired standalone preparation protocol is unavailable: '+name,'UNSUPPORTED_PREPARATION_PROTOCOL','Use Chimera preparation for docking, or native DOCK6 6.11 Amber Score. Existing results are preserved.')
    candidates=sorted((installation/'bin').glob('amber*/bin/'+name),reverse=True)
    executable=next((p for p in candidates if p.is_file() and os.access(p,os.X_OK)),None)
    if executable is None:raise ScientificToolError('Chimera bundled charge utility is unavailable: '+name,'TOOL_NOT_FOUND','Check the classic Chimera installation.')
    return AmberTool(name,executable.resolve(),'UCSF Chimera bundled utilities',executable.resolve().parent.parent)

def amber_tool(name,prefix='',search_path=None,environment=''):
    return resolve_amber_tool(name,prefix,environment).executable

def amber_command(name,args=(),*,prefix='',environment=''):
    resolved=resolve_amber_tool(name,prefix,environment)
    return [str(resolved.executable),*[str(v) for v in args]],resolved

def run_amber_tool(name,args=(),cwd=None,check=True,timeout=None,*,env=None,prefix='',environment=''):
    command,resolved=amber_command(name,args,prefix=prefix,environment=environment)
    child=dict(os.environ if env is None else env);child['AMBERHOME']=str(resolved.prefix)
    library=resolved.prefix.parent.parent/'lib'
    child['LD_LIBRARY_PATH']=str(library)+os.pathsep+child.get('LD_LIBRARY_PATH','')
    child['PATH']=str(resolved.prefix/'bin')+os.pathsep+child.get('PATH','')
    return run_command(command,cwd=cwd,check=check,timeout=timeout,env=child,tool_name=name,tool_environment=resolved.environment)
