"""Run a Vina batch in a crash-isolated child with per-ligand checkpoints."""
import argparse
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import time

from .storage import write_json
from .calculation_cache import CalculationCache, identity
from .stage_cache import file_digest


def _options(configuration):
    values={}
    try:text=Path(configuration).read_text()
    except OSError:return values
    for line in text.splitlines():
        line=line.split('#',1)[0].strip()
        if not line:continue
        name,separator,value=line.partition('=')
        if separator:values[name.strip()]=value.strip()
    return values


def _output(configuration,directory):
    value=_options(configuration).get('out')
    return (Path(directory)/value).resolve() if value else None


def _identity(configuration):
    options=_options(configuration);name=Path(options.get('ligand') or Path(configuration).stem).name
    return name.removesuffix('.lig.pdbqt').removesuffix('.pdbqt')


def validate_pose(role, path):
    from .ligand_cache import validate_prepared_ligand
    from .docking_data import score
    validate_prepared_ligand(path, 'pdbqt')
    score(path, 'vina')


def calculation_request(configuration, directory, fingerprints=None):
    """Replace transport paths with content hashes; preserve scientific options."""
    import shutil
    options = _options(configuration)
    request = dict(options)
    fingerprints={} if fingerprints is None else fingerprints
    def fingerprint(path):
        path=Path(path);stat=path.stat()
        stamp=(stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns)
        previous=fingerprints.get(path)
        if previous is None or previous[0]!=stamp:
            previous=(stamp,file_digest(path));fingerprints[path]=previous
        return previous[1]
    for name in ('receptor', 'ligand', 'flex'):
        if name in request:
            request[name] = fingerprint(Path(directory) / request[name])
    # Custom map-prefix configurations need a separate map manifest.
    if 'maps' in request or not {'receptor', 'ligand', 'out'} <= options.keys():
        return None
    request.pop('out', None); request.pop('log', None)
    try:request['vina_version'] = importlib.metadata.version('vina')
    except importlib.metadata.PackageNotFoundError:request['vina_version'] = None
    executable = shutil.which('vina')
    request['vina_executable'] = fingerprint(executable) if executable else None
    try:
        specification=importlib.util.find_spec('vina')
        root=Path(specification.origin).parent if specification and specification.origin else None
        request['vina_extension']={p.name:fingerprint(p) for p in sorted(root.glob('*.so'))} if root else {}
    except (ImportError,ValueError,AttributeError):request['vina_extension']={}
    request['adapter'] = identity([fingerprint(Path(__file__)),
                                  fingerprint(Path(__file__).with_name('vina_session.py'))])
    return request


def execute(request_path,state_path):
    """Child entry point. A native abort leaves the active config checkpointed."""
    request=json.loads(Path(request_path).read_text());directory=Path(request['directory'])
    state=dict(active=None,completed=[],failed=[])
    from .ligand_cache import validate_prepared_ligand
    from .vina_session import VinaSessions
    sessions=VinaSessions()
    cache=CalculationCache('vina',directory)
    fingerprints={}
    for name in request['configurations']:
        configuration=directory/name;state['active']=name;write_json(state,state_path)
        started=time.monotonic()
        try:
            options=_options(configuration)
            if options.get('ligand'):
                validate_prepared_ligand(directory/options['ligand'],'pdbqt')
            output=_output(configuration,directory)
            recipe=calculation_request(configuration,directory,fingerprints)
            key=identity(recipe) if recipe else identity(options)
            with cache.locked(key):
                reused=recipe is not None and cache.restore(key,{'pose':output},validate_pose) is not None
                if not reused:
                    used_api=sessions.run(configuration,directory)
                    if not used_api:
                        result=subprocess.run(['vina','--config',name],cwd=directory,capture_output=True,text=True)
                        if result.returncode:
                            raise RuntimeError('Vina CLI terminou com código '+str(result.returncode)+': '+(result.stderr or result.stdout)[-2000:])
                    if output is None or not output.is_file() or not output.stat().st_size:
                        raise RuntimeError('Vina não produziu a pose configurada.')
                    validate_pose('pose',output)
                    if recipe is not None:
                        cache.remember(key,{'pose':output},dict(recipe=recipe,compound_identifier=_identity(configuration),
                            receptor_identifier=Path(options['receptor']).name.removesuffix('.pdbqt').removesuffix('.dockprep')),validate_pose)
            state['completed'].append(dict(configuration=name,compound_identifier=_identity(configuration),
                                           cache_key=key,reused=bool(reused),
                                           duration_ms=round((time.monotonic()-started)*1000)))
        except Exception as error:
            output=_output(configuration,directory)
            if output:output.unlink(missing_ok=True)
            state['failed'].append(dict(configuration=name,compound_identifier=_identity(configuration),
                                        error=type(error).__name__+': '+str(error),
                                        duration_ms=round((time.monotonic()-started)*1000)))
        state['active']=None;write_json(state,state_path)
    return state


def run_resilient(configurations,directory,python=None):
    """Restart after an uncatchable native failure and skip only its ligand."""
    from .diagnostics import event,get_logger
    from .progress import report_progress
    directory=Path(directory);names=[Path(path).name for path in configurations]
    pending=list(names);completed=[];failed=[];attempts=[];python=str(python or sys.executable)
    request=directory/'.vina-batch-request.json';state_path=directory/'.vina-batch-state.json'
    logger=get_logger('scientific_tools')
    cache=CalculationCache('vina',directory);fingerprints={}
    # Resolve indexed hits before starting a native worker. A full hit needs no
    # Vina import, map generation or subprocess at all.
    unresolved=[]
    for name in pending:
        try:
            recipe=calculation_request(directory/name,directory,fingerprints)
            key=identity(recipe) if recipe else None
            if recipe is not None and cache.restore(key,{'pose':_output(directory/name,directory)},validate_pose) is not None:
                completed.append(dict(configuration=name,compound_identifier=_identity(directory/name),
                                      cache_key=key,reused=True,duration_ms=0))
                continue
        except (OSError,ValueError,KeyError):pass
        unresolved.append(name)
    pending=unresolved
    try:version=importlib.metadata.version('vina')
    except importlib.metadata.PackageNotFoundError:version='CLI or unknown'
    event(logger,'vina.batch.started','Lote Vina iniciado em processo isolado.',total=len(names),vina_version=version)
    attempt=0
    while pending:
        attempt+=1;write_json(dict(directory=str(directory.resolve()),configurations=pending),request)
        state_path.unlink(missing_ok=True);started=time.monotonic()
        result=subprocess.run([python,'-m','biomolexplorer.vina_batch','--request',str(request),'--state',str(state_path)],
                              cwd=directory,capture_output=True,text=True)
        try:state=json.loads(state_path.read_text())
        except (OSError,ValueError,TypeError):state=dict(active=None,completed=[],failed=[])
        done={entry['configuration'] for entry in state.get('completed',[])}
        known_completed={entry['configuration'] for entry in completed}
        known_failed={entry['configuration'] for entry in failed}
        completed.extend(entry for entry in state.get('completed',[]) if entry['configuration'] not in known_completed)
        failed.extend(entry for entry in state.get('failed',[]) if entry['configuration'] not in known_failed)
        active=state.get('active')
        attempts.append(dict(attempt=attempt,returncode=result.returncode,duration_ms=round((time.monotonic()-started)*1000),
                             active=active,stderr=(result.stderr or '')[-4000:]))
        if result.returncode and active and active not in {entry['configuration'] for entry in failed}:
            output=_output(directory/active,directory)
            if output:output.unlink(missing_ok=True)
            failed.append(dict(configuration=active,compound_identifier=_identity(directory/active),
                               error='Falha nativa isolada do Vina (código '+str(result.returncode)+'): '+(result.stderr or '')[-1000:]))
            event(logger,'vina.ligand.crashed','Falha nativa do Vina isolada; os demais compostos continuarão.',
                  configuration=active,compound_identifier=_identity(directory/active))
        accounted=done|{entry['configuration'] for entry in failed}
        new_pending=[name for name in pending if name not in accounted]
        if not active and new_pending==pending:
            failed.extend(dict(configuration=name,compound_identifier=_identity(directory/name),
                               error='O processo Vina terminou antes de identificar o composto ativo.') for name in pending)
            pending=[]
        else:pending=new_pending
        report_progress('Vina: '+str(len(completed))+' concluídos; '+str(len(failed))+' falhas; '+str(len(pending))+' pendentes.',
                        len(completed)+len(failed),len(names))
    experiments={}
    for name in names:
        try:experiments[name]=dict(options=_options(directory/name),recipe=calculation_request(directory/name,directory,fingerprints))
        except (OSError,ValueError):experiments[name]=dict(options=_options(directory/name))
    summary=dict(version=1,vina_version=version,configurations=len(names),completed=completed,failed=failed,
                 attempts=attempts,experiments=experiments)
    write_json(summary,directory/'vina_execution.json')
    request.unlink(missing_ok=True);state_path.unlink(missing_ok=True)
    event(logger,'vina.batch.completed','Lote Vina concluído.',completed=len(completed),failed=len(failed),attempts=attempt)
    return summary


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--request',required=True);parser.add_argument('--state',required=True)
    args=parser.parse_args();execute(args.request,args.state)


if __name__=='__main__':main()
