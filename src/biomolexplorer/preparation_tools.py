"""Central resolution and execution of AmberTools programs.

The BioMolExplorer interpreter stays in its main environment. Amber programs
are invoked through ``conda run`` when the configured installation is a Conda
environment, so no interactive activation is required.
"""
from dataclasses import dataclass
from functools import lru_cache
import json
import os
from pathlib import Path
import shutil
import subprocess

from .processes import ScientificToolError, run_command

DEFAULT_AMBERTOOLS_ENV='BioMolExplorer-AmberTools'
AMBERTOOLS_ENV_VAR='BIOMOL_AMBERTOOLS_ENV'
AMBERTOOLS_PATH_VAR='BIOMOL_AMBERTOOLS_PATH'


@dataclass(frozen=True)
class AmberTool:
    name: str
    executable: Path
    environment: str
    conda: Path | None=None
    prefix: Path | None=None


def conda_executable():
    """Find Conda without assuming that the launching shell activated it."""
    configured=os.environ.get('CONDA_EXE','').strip()
    candidates=[configured,shutil.which('conda')]
    import sys
    current=Path(sys.executable).resolve()
    candidates.append(current.parent/'conda')
    candidates.extend(parent/'bin'/'conda' for parent in current.parents)
    for value in candidates:
        if value:
            path=Path(value).expanduser()
            if path.is_file() and os.access(path,os.X_OK):return path.resolve()
    return None


@lru_cache(maxsize=8)
def _conda_environments(conda):
    try:
        result=subprocess.run([str(conda),'env','list','--json'],capture_output=True,
            text=True,check=False,timeout=30)
        if result.returncode:return []
        return [Path(value).expanduser().resolve() for value in json.loads(result.stdout).get('envs',[])]
    except (OSError,ValueError,subprocess.TimeoutExpired):
        return []


def _environment_prefix(conda,name):
    candidate=Path(name).expanduser()
    if (candidate/'conda-meta').is_dir():return candidate.resolve()
    return next((prefix for prefix in _conda_environments(conda) if prefix.name==name),None)


def _explicit_tool(name,prefix):
    if not prefix:return None
    root=Path(prefix).expanduser().resolve()
    candidates=(root,root/name,root/'bin'/name) if root.is_file() else (root/name,root/'bin'/name)
    executable=next((path for path in candidates if path.is_file() and os.access(path,os.X_OK)),None)
    if executable is None:
        raise ScientificToolError(
            f'Explicit AmberTools installation was found, but "{name}" is unavailable: {root}',
            'TOOL_NOT_FOUND','Set ambertools_path or BIOMOL_AMBERTOOLS_PATH to an installation containing bin/'+name+'.')
    return AmberTool(name,executable.resolve(),'explicit:'+str(root),prefix=executable.resolve().parent.parent)


def resolve_amber_tool(name,prefix='',environment=''):
    """Resolve explicit path -> configured Conda environment -> current PATH.

    Chimera-bundled Amber executables are intentionally not searched.  The
    legacy Chimera backend is opt-in at its call site and must never become an
    implicit dependency of the native preparation path.
    """
    explicit=_explicit_tool(name,prefix or os.environ.get(AMBERTOOLS_PATH_VAR,''))
    if explicit:return explicit
    env_name=(environment or os.environ.get(AMBERTOOLS_ENV_VAR,'') or DEFAULT_AMBERTOOLS_ENV).strip()
    env_was_explicit=bool(environment or os.environ.get(AMBERTOOLS_ENV_VAR,''))
    conda=conda_executable()
    if conda:
        env_prefix=_environment_prefix(conda,env_name)
        if env_prefix:
            executable=env_prefix/'bin'/name
            if not executable.is_file() or not os.access(executable,os.X_OK):
                raise ScientificToolError(
                    f'AmberTools environment was found, but "{name}" is unavailable.',
                    'TOOL_NOT_FOUND',f'Install {name} in Conda environment "{env_name}" or select another AmberTools environment.')
            return AmberTool(name,executable.resolve(),env_name,conda.resolve(),env_prefix)
        if env_was_explicit:
            raise ScientificToolError(
                f'AmberTools Conda environment "{env_name}" was not found.',
                'AMBER_ENV_NOT_FOUND','Create the environment from requirements-ambertools.yml or configure BIOMOL_AMBERTOOLS_ENV.')
    elif env_was_explicit:
        raise ScientificToolError(
            f'Conda could not be located while resolving AmberTools environment "{env_name}".',
            'CONDA_NOT_FOUND','Set CONDA_EXE to the Conda executable or make conda available on PATH.')
    executable=shutil.which(name)
    if executable:return AmberTool(name,Path(executable).resolve(),'current-environment')
    if not conda:
        raise ScientificToolError(
            f'Conda could not be located and AmberTools tool "{name}" is not on PATH.',
            'CONDA_NOT_FOUND','Set CONDA_EXE, install the environment from requirements-ambertools.yml, or configure BIOMOL_AMBERTOOLS_PATH.')
    raise ScientificToolError(
        f'AmberTools Conda environment "{env_name}" was not found and "{name}" is unavailable in the current environment.',
        'AMBER_ENV_NOT_FOUND','Create the environment from requirements-ambertools.yml or configure BIOMOL_AMBERTOOLS_ENV.')


def amber_tool(name,prefix='',search_path=None,environment=''):
    """Compatibility resolver returning the real executable for identities."""
    try:return resolve_amber_tool(name,prefix,environment).executable
    except ScientificToolError:
        if search_path and not (prefix or environment or os.environ.get(AMBERTOOLS_ENV_VAR)):
            executable=shutil.which(name,path=search_path)
            if executable:return Path(executable).resolve()
        raise


def amber_command(name,args=(),*,prefix='',environment=''):
    resolved=resolve_amber_tool(name,prefix,environment)
    arguments=[str(value) for value in args]
    if resolved.conda:
        selector=['-n',resolved.environment] if not Path(resolved.environment).is_absolute() else ['-p',str(resolved.prefix)]
        # Use the already validated absolute executable. An unrelated program
        # with the same name (notably DOCK6's bundled tleap launcher) may be
        # earlier on the parent process PATH even after ``conda run``.
        command=[str(resolved.conda),'run','--no-capture-output',*selector,str(resolved.executable),*arguments]
    else:command=[str(resolved.executable),*arguments]
    return command,resolved


def run_amber_tool(name,args=(),cwd=None,check=True,timeout=None,*,env=None,prefix='',environment=''):
    """Run an AmberTools executable with environment-aware diagnostics."""
    command,resolved=amber_command(name,args,prefix=prefix,environment=environment)
    try:
        return run_command(command,cwd=cwd,check=check,timeout=timeout,env=env,
            tool_name=name,tool_environment=resolved.environment)
    except ScientificToolError as error:
        if error.error_code in ('TOOL_EXIT_FAILED','TOOL_REPORTED_ERROR'):
            error.action=(f'AmberTools tool={name}; environment={resolved.environment}; '
                f'working_directory={Path(cwd or ".").resolve()}. '+error.action)
        raise
