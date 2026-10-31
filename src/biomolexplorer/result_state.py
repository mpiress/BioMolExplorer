"""Scientific request equivalence, independent of presentation and run paths.

Only declared independent axes may grow monotonically. An arbitrary list of
scientific parameters is never assumed to represent independent work.
"""
import copy
import hashlib
import json
from functools import lru_cache
from pathlib import Path
from .bindings import sources

FINGERPRINTS=('morgan','maccs','pharmacophore')


def preparation(config,charge='gas'):
    from .redocking_config import DEFAULTS
    config=config or {}
    value=dict(cofactors=sorted(set(config.get('cofactors',[]))),**{
        role:{**DEFAULTS,'charge_type':charge if role=='ligand' else 'gas',**config.get(role,{})}
        for role in ('receptor','ligand')})
    if config.get('ligand_chain') is not None:value['ligand_chain']=config['ligand_chain']
    return value


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def request(stage):
    keys=('operation','parameters','provided_results','bindings','templates','depends_on','input_processing')
    payload=json.dumps({key:stage[key] for key in keys if key in stage},sort_keys=True,allow_nan=False)
    return copy.deepcopy(_request(payload))


@lru_cache(maxsize=512)
def _request(payload):
    from .catalog import new_stage,operation_fields,template_names
    from .templates import RESOURCE_ROOT
    from .pipeline import configured_inputs
    stage=configured_inputs(json.loads(payload))
    op=stage['operation']
    parameters=(copy.deepcopy(new_stage(op)['parameters']) if op=='import_results' else
                {field['name']:copy.deepcopy(field['default']) for field in operation_fields(op)
                 if field['default'] is not None and not field['path']})
    parameters.update(stage.get('parameters',{}))
    provided=stage.get('provided_results')
    if provided:parameters=copy.deepcopy(provided)
    parameters.pop('verbose',None)
    if not provided and op=='admet':
        from .admet_config import resolve_config,PARAMETERS
        effective=resolve_config(parameters);effective.pop('admet_profile',None)
        effective['alert_catalogs']=sorted(set(effective['alert_catalogs']))
        parameters={key:value for key,value in parameters.items() if key not in PARAMETERS}
        parameters.update(effective)
    if not provided and op in ('redocking','prepare_structures','docking_vina','docking_dock6'):
        parameters['preparation_options']=preparation(parameters.get('preparation_options'),parameters.get('charge_type','gas'))
    # Filenames selected in a binding have one canonical representation.
    for field in stage.get('bindings',{}):parameters.pop(field,None)
    if stage.get('bindings'):
        for field in ('files','filename','input_file','dataset_origin'):
            parameters.pop(field,None)
    if parameters.get('graph_inputs'):
        parameters['graph_inputs']=sorted([
            {key:([Path(p).name for p in value] if key=='compound_files' else Path(value).name if key=='file' else value)
             for key,value in entry.items() if key!='label'} for entry in parameters['graph_inputs']],key=digest)
    axes={}
    if op=='import_results' or provided:
        axes['assets']=sorted(set(parameters.pop('asset_ids',[])))
        types=parameters.pop('asset_types',{})
        axes['asset_types']={key:types.get(key,parameters.get('kind')) for key in axes['assets']}
    if not provided and op=='fingerprints':
        axes['actions']=[key for key in FINGERPRINTS if parameters.pop(key, key=='morgan')]
        if 'morgan' not in axes['actions']:
            parameters.pop('radius',None);parameters.pop('morgan_n_bits',None)
        parameters.pop('chunk_size',None) # Scheduling, not fingerprint science.
    if not provided and op in ('redocking','prepare_structures'):
        records=parameters.pop('pdb_codes',[]) or []
        selected=[Path(r.get('selector','')).name for group in stage.get('bindings',{}).values() for r in sources(group) if r.get('selector','auto')!='auto']
        pdbs={name.split('.',1)[0].split('_',1)[0] for name in selected if name.endswith(('.pdb','.pdbqt','.mol2'))}
        if pdbs:records=[row for row in records if row[0] in pdbs]
        from .redocking_config import pair_key
        axes['pairs']=sorted({pair_key(row) for row in records})
        parameters.pop('preparation_pairs',None)
    if op=='redocking':
        engine=parameters.get('redocking_engine','vina')
        unused=('dock6_app_path','density','radius','distance','conformer_search_type','max_orientations','num_scored_conformers','dock6_seed') if engine=='vina' else ('sizeof_box','exhaustiveness','num_modes','dock6_app_path')
        for field in unused:parameters.pop(field,None)
    if op=='docking_dock6':parameters.pop('dock6_app_path',None)
    templates={}
    for name in template_names(op) if not provided else []:
        if op=='redocking' and name.startswith('vina/' if parameters.get('redocking_engine','vina')=='dock6' else 'dock6/'):continue
        templates[name]=stage.get('templates',{}).get(name,(RESOURCE_ROOT/name).read_text())
    bindings={} if provided else stage.get('bindings',{})
    inputs={field:sorted((dict(r,selector=r.get('selector','auto')) if 'stage' in r else r for r in sources(group)),key=digest)
            for field,group in bindings.items()}
    return dict(operation=op,provided=bool(provided),parameters=parameters,axes=axes,templates=templates,
                inputs=inputs,dependencies=sorted(stage.get('depends_on',[])),
                input_processing=stage.get('input_processing','merge'))


def signature(stage):return digest(request(stage))


def removed_axes(previous,current):
    old,new=request(previous)['axes'],request(current)['axes']
    return any(set(values)-set(new.get(axis,[])) for axis,values in old.items() if axis!='asset_types')


def fingerprint_units(batch):
    """Expose identifiable outputs of old combined jobs as immutable units."""
    config=batch.get('configuration',{})
    if config.get('operation')!='fingerprints':return [batch]
    actions=[key for key in FINGERPRINTS if config['parameters'].get(key,True)]
    if len(actions)<2:return [batch]
    groups={key:[p for p in batch.get('artifacts',[]) if Path(p).name.startswith(key+'_')] for key in actions}
    if not all(groups.values()):return [batch] # Do not guess an unindexed legacy mapping.
    shared=[p for p in batch['artifacts'] if not any(Path(p).name.startswith(key+'_') for key in FINGERPRINTS)]
    units=[]
    for action in actions:
        unit=copy.deepcopy(batch);unit['configuration']['parameters'].update({key:key==action for key in FINGERPRINTS})
        unit['artifacts']=groups[action]+shared;allowed=set(unit['artifacts'])
        unit['artifact_manifest']={p:sha for p,sha in unit.get('artifact_manifest',{}).items() if p in allowed}
        unit['artifact_records']=[r for r in unit.get('artifact_records',[]) if r['path'] in allowed]
        unit['label']=batch.get('label',config.get('name',''))+' · '+action
        units.append(unit)
    return units


@lru_cache(maxsize=16384)
def _file_digest(path,stamp):
    from .stage_cache import file_digest
    return file_digest(path)


def input_digest(path):
    path=Path(path);stat=path.stat()
    return _file_digest(str(path),(stat.st_dev,stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns))


def inputs_changed(entry):
    """Compare persisted content hashes; old unversioned inputs are not guessed."""
    manifest=entry.get('input_manifest',{})
    try:
        for field,paths in entry.get('input_files',{}).items():
            expected=manifest.get('selected:'+field)
            if isinstance(expected,list) and sorted(input_digest(p) for p in paths)!=sorted(value[1] for value in expected):return True
    except OSError:return True
    return False


def incompatible(previous,current,*,include_axes=True):
    """Can this earlier input/action unit still be requested unchanged?"""
    old,new=request(previous),request(current)
    # A batch for one algorithm must be compared with that same algorithm,
    # rather than with options of newly enabled algorithms.
    if old['operation']==new['operation']=='fingerprints':
        selected=old['axes'].get('actions',[])
        if 'morgan' not in selected:
            for field in ('radius','morgan_n_bits'):
                old['parameters'].pop(field,None);new['parameters'].pop(field,None)
    for key in ('operation','provided','parameters','templates','input_processing'):
        if old[key]!=new[key]:return True
    if not include_axes:return False
    for axis,values in old['axes'].items():
        if axis=='asset_types':
            if any(new['axes'].get(axis,{}).get(key)!=value for key,value in values.items() if key in new['axes'].get('assets',[])):return True
        elif not set(values)<=set(new['axes'].get(axis,[])):return True
    if not set(old['dependencies'])<=set(new['dependencies']):return True
    if previous.get('operation') in ('redocking','prepare_structures'):
        from .redocking_config import effective_preparation,pair_key
        before,after=effective_preparation(previous['parameters']),effective_preparation(current['parameters'])
        selected={pair_key(row) for row in current['parameters'].get('pdb_codes',[]) or []}
        for row in previous['parameters'].get('pdb_codes',[]) or []:
            key=pair_key(row)
            if key in selected and preparation(before.get(key),previous['parameters'].get('charge_type','gas'))!=preparation(after.get(key),current['parameters'].get('charge_type','gas')):return True
    return False
