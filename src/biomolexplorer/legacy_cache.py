"""Read-only verification of the audited version-1 cache contract."""
import hashlib
import json
from pathlib import Path
from .stage_cache import file_digest,path_manifest
from .catalog import PATH_FIELDS,template_names
from .templates import RESOURCE_ROOT


def compatible(operation,provenance,parameters,current_software=None):
    """Adopt only an explicitly audited source build with unchanged science."""
    from .stage_cache import operation_digest,runtime_signature
    import sys
    path=RESOURCE_ROOT/'legacy_cache_compatibility.json'
    if not path.is_file() or not provenance:return False
    policy=json.loads(path.read_text())
    if operation=='admet' and parameters.get('prediction_backend')=='admet_ai':return False
    if provenance.split(':',1)[0] not in policy.get('source_builds',[]) or policy.get('operations',{}).get(operation)!=operation_digest(operation):return False
    runtime=runtime_signature(operation,current_software.split(':',1)[-1] if current_software else sys.executable,parameters)
    approved=policy.get('runtime',{})
    return runtime.get('python')==approved.get('python') and all(
        value==approved.get('packages',{}).get(name) for name,value in runtime.get('packages',{}).items()) and all(
        value==policy.get('tools',{}).get(name) for name,value in runtime.get('tools',{}).items())

def legacy_stage_key(stage, parameters, inputs, software):
    config = {'version': 1, 'operation': stage['operation'], 'parameters': dict(parameters),
              'bindings': stage.get('bindings',{}), 'provided_results':stage.get('provided_results'),
              'inputs': dict(inputs), 'software': software, 'templates': {},
              'input_processing':stage.get('input_processing','merge')}
    if stage.get('provided_results'):
        config['bindings']={}
    elif stage['operation']=='admet':
        # Schema 3 omitted independently selected rules/alerts. These outputs
        # must be regenerated once, never adopted as a complete schema 4 run.
        from .admet_outputs import SCHEMA_VERSION
        config['admet_output_schema']=SCHEMA_VERSION
    elif stage['operation']=='redocking':
        from .ligand_mol2 import PREPARATION_SCHEMA_VERSION
        config['redocking_preparation_schema']=PREPARATION_SCHEMA_VERSION
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
