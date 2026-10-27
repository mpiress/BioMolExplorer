"""Content-based reuse of project stage results, independent of canvas positions."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
from functools import lru_cache

from .catalog import PATH_FIELDS, template_names
from .paths import SOURCE_ROOT
from .templates import RESOURCE_ROOT
from .bindings import dependencies


def file_digest(path):
    with Path(path).open('rb') as stream:
        digest = hashlib.sha256()
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
        return digest.hexdigest()


def artifact_manifest(artifacts):
    return {str(Path(path).resolve()): file_digest(path) for path in artifacts}


def manifest_matches(manifest):
    if not manifest:
        return False
    try:
        return all(file_digest(path) == digest for path, digest in manifest.items())
    except OSError:
        return False


@lru_cache(maxsize=256)
def _runtime_digest(path, signature):
    return file_digest(path)


def runtime_manifest(operation, parameters, worker_python=None):
    """Native software changes also invalidate whole-stage hits."""
    if operation not in ('redocking','prepare_structures','docking_vina','docking_dock6'):
        return {}
    binary=Path(worker_python or sys.executable).resolve()
    search=str(binary.parent)+os.pathsep+os.environ.get('PATH','')
    # Docking and redocking prepare exclusively with Chimera, including saved
    # auto/python options. Hash the tools actually executed, not external AmberTools.
    names=['obabel','chimera']
    if operation=='docking_vina' or operation=='redocking' and parameters.get('redocking_engine','vina')=='vina':
        names.append('vina')
    files={name:shutil.which(name,path=search) for name in names}
    if (operation=='docking_dock6' or operation=='redocking' and parameters.get('redocking_engine','vina')=='dock6') and parameters.get('dock6_app_path'):
        root=Path(parameters['dock6_app_path'])
        for folder in ('bin','parameters'):
            for path in sorted((root/folder).glob('*')):
                if path.is_file():files[str(path.relative_to(root))]=str(path)
    # Conda records identify the worker's Python extension builds without
    # importing scientific libraries on the UI/scheduler thread.
    for path in sorted((binary.parent.parent/'conda-meta').glob('*.json')):
        if path.name.startswith(('vina-','rdkit-','openbabel-','numpy-','pymol-open-source-')):
            files['package:'+path.name]=str(path)
    result={}
    for name,path in files.items():
        if path is None:result[name]=None;continue
        stat=Path(path).stat()
        result[name]=_runtime_digest(path,(stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns))
    return result


def implementation_digest():
    digest = hashlib.sha256()
    for name in ('artifact_naming.py','artifacts.py'):
        digest.update((SOURCE_ROOT/'biomolexplorer'/name).read_bytes())
    for name in ('structure_selection.py','amber_preparation.py','native_preparation.py'):
        digest.update((SOURCE_ROOT/'biomolexplorer'/name).read_bytes())
    digest.update((SOURCE_ROOT/'biomolexplorer/preparation_failures.py').read_bytes())
    for name in ('chimera_charge_repair.py','chimera_mol2.py','dock6_redocking.py','redocking_config.py','docking_preparation.py','ligand_cache.py','docking_data.py','dock6_validation.py','pdb_ligands.py','redocking_preparation.py','redocking_cache.py','structure_cache.py','dock6_cache.py','molecular_core.py','molecular_adapters.py','preparation_tools.py','molecular_parsing.py'):
        digest.update((SOURCE_ROOT/'biomolexplorer'/name).read_bytes())
    for name in ('md_config.py','md_engine.py','md_analysis.py','md_inputs.py','md_cache.py','result_visibility.py'):
        digest.update((SOURCE_ROOT/'biomolexplorer'/name).read_bytes())
    for folder in ('caad', 'crawlers', 'kernel', 'wrappers'):
        for path in sorted((SOURCE_ROOT / folder).rglob('*.py')):
            digest.update(str(path.relative_to(SOURCE_ROOT)).encode())
            digest.update(path.read_bytes())
    for name in ('operations.py','pipeline.py','input_validation.py','molecule_quality.py','fingerprint_selection.py','graph_inputs.py','docking_inputs.py','visualizations.py','admet_config.py','admet_predictors.py','admet_ai_runner.py','egg_plot.py','vina_session.py','vina_batch.py','calculation_cache.py'):
        digest.update((SOURCE_ROOT / 'biomolexplorer' / name).read_bytes())
    return digest.hexdigest()


def path_manifest(path):
    path = Path(path)
    if path.is_file():
        return {path.name: file_digest(path)}
    ignored = {'.biomolexplorer', '.calculation-cache', '.curation', '__pycache__', 'cache', 'logs', '.history', '.exports', '.input-cache', 'project.json'}
    return {p.relative_to(path).as_posix(): file_digest(p) for p in sorted(path.rglob('*'))
            if p.is_file() and not (set(p.relative_to(path).parts) & ignored)}


def stage_key(stage, parameters, inputs, software):
    config = {'version': 1, 'operation': stage['operation'], 'parameters': dict(parameters),
              'bindings': stage.get('bindings',{}), 'provided_results':stage.get('provided_results'),
              'inputs': dict(inputs), 'software': software, 'templates': {},
              'input_processing':stage.get('input_processing','merge')}
    if stage['operation']=='redocking':config['selection_scope']='requested_pairs_v1'
    if stage['operation'] in ('redocking','prepare_structures','docking_vina','docking_dock6'):
        # Effective, materialized inputs carry exact selections. Auto versus
        # the equivalent confirmed selector must not force a redundant worker.
        from .bindings import sources
        config['bindings']={field:sorted({(ref.get('stage') or '',ref.get('asset') or '') for ref in sources(group)})
                            for field,group in stage.get('bindings',{}).items()}
    if stage.get('provided_results'):
        config['bindings']={}
    if stage['operation']=='graphs' and parameters.get('graph_inputs'):
        # Source run paths must not change an otherwise identical experiment.
        config['parameters']['graph_inputs']=[dict(entry,
            file={'name':Path(entry['file']).name,'sha256':file_digest(entry['file'])},
            compound_files=[{'name':Path(p).name,'sha256':file_digest(p)} for p in entry['compound_files']])
            for entry in parameters['graph_inputs']]
    if stage['operation']=='dynamics_dock6':
        from .md_inputs import molecular_identity
        config['parameters']['md_inputs']=[dict(id=p['id'],**{role:dict(identity=p[role]['identity'],sha256=molecular_identity(p[role]))
            for role in ('receptor','ligand')}) for p in parameters.get('md_inputs',[])]
        config['parameters'].pop('md_seed_path',None)
    # The installation directory is a machine setting. Scientific options and
    # immutable input/output hashes still determine whether an experiment fits.
    if stage['operation'] in ('docking_dock6','dynamics_dock6'):
        config['parameters'].pop('dock6_app_path',None)
    # Bound paths change with run IDs; input contents, selectors and inferred
    # filenames determine the scientific request instead of those output paths.
    for field in stage.get('bindings', {}):
        config['parameters'].pop(field, None)
    for field in PATH_FIELDS:
        if parameters.get(field) and field not in stage.get('bindings', {}):
            config['inputs']['path:' + field] = path_manifest(parameters[field])
            config['parameters'].pop(field,None)
    for name in template_names(stage['operation']):
        config['templates'][name] = stage.get('templates', {}).get(name, (RESOURCE_ROOT / name).read_text())
    # Extra valid overlays can influence preparation protocols indirectly.
    config['templates'].update(stage.get('templates', {}))
    encoded = json.dumps(config, sort_keys=True, ensure_ascii=False, allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def input_manifests(stage, results):
    inputs = {}
    for identifier in sorted(dependencies(stage)):
        files = []
        for path in results[identifier]:
            parts = Path(path).parts
            # Keep the output's internal path while discarding its source run ID.
            start = len(parts) - 1 - list(reversed(parts)).index('artifacts') if 'artifacts' in parts else len(parts) - 2
            files.append(('/'.join(parts[start + 1:]), file_digest(path)))
        inputs['stage:' + identifier] = sorted(files)
    return inputs
