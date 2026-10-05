"""Content-based reuse of project stage results, independent of canvas positions."""
import hashlib
import json
from pathlib import Path

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


def implementation_digest():
    digest = hashlib.sha256()
    for folder in ('caad', 'crawlers', 'kernel', 'wrappers'):
        for path in sorted((SOURCE_ROOT / folder).rglob('*.py')):
            digest.update(str(path.relative_to(SOURCE_ROOT)).encode())
            digest.update(path.read_bytes())
    for name in ('operations.py','pipeline.py','input_validation.py','molecule_quality.py','fingerprint_selection.py','graph_inputs.py','docking_inputs.py','visualizations.py'):
        digest.update((SOURCE_ROOT / 'biomolexplorer' / name).read_bytes())
    return digest.hexdigest()


def path_manifest(path):
    path = Path(path)
    if path.is_file():
        return {path.name: file_digest(path)}
    ignored = {'.biomolexplorer', '.curation', '__pycache__', 'cache', 'logs', '.history', '.exports', '.input-cache', 'project.json'}
    return {p.relative_to(path).as_posix(): file_digest(p) for p in sorted(path.rglob('*'))
            if p.is_file() and not (set(p.relative_to(path).parts) & ignored)}


def stage_key(stage, parameters, inputs, software):
    config = {'version': 1, 'operation': stage['operation'], 'parameters': dict(parameters),
              'bindings': stage.get('bindings',{}), 'provided_results':stage.get('provided_results'),
              'inputs': dict(inputs), 'software': software, 'templates': {},
              'input_processing':stage.get('input_processing','merge')}
    if stage.get('provided_results'):
        config['bindings']={}
    if stage['operation']=='graphs' and parameters.get('graph_inputs'):
        # Source run paths must not change an otherwise identical experiment.
        config['parameters']['graph_inputs']=[dict(entry,
            file={'name':Path(entry['file']).name,'sha256':file_digest(entry['file'])},
            compound_files=[{'name':Path(p).name,'sha256':file_digest(p)} for p in entry['compound_files']])
            for entry in parameters['graph_inputs']]
    # The installation directory is a machine setting. Scientific options and
    # immutable input/output hashes still determine whether an experiment fits.
    if stage['operation']=='docking_dock6':
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
