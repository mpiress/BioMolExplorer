"""Verified, immutable completed results, separate from current-version job cache.

Choosing reuse means retaining a completed experiment with its original producer
version. It never certifies that an old experiment used the current algorithm.
"""
import copy
import json
from collections import Counter
from pathlib import Path

from .bindings import sources,asset_references,dependencies
from .result_state import request,signature,input_digest,digest
from .artifacts import validate_reference,locator


def _scoped(root,path):
    resolved=Path(path).resolve()
    if not resolved.is_relative_to(root):raise ValueError('O caminho deve pertencer ao projeto.')
    return resolved


def artifact_integrity(store,pid,entry,*,include_batches=False,root=None):
    """Verify output bytes without treating an unfinished stage as completed."""
    paths=entry.get('artifacts',[]);manifest=dict(entry.get('artifact_manifest',{}))
    if include_batches:
        for batch in entry.get('batches',[]):manifest.update(batch.get('artifact_manifest',{}))
    if not paths:return False
    if not include_batches and set(manifest)!={str(Path(p).resolve()) for p in paths}:return False
    try:
        root=store.project_dir(pid).resolve() if root is None else root
        return all(input_digest(_scoped(root,p))==manifest.get(str(Path(p).resolve())) for p in paths)
    except (OSError,ValueError):return False


def intact(store,pid,entry,root=None):
    if entry.get('status')!='succeeded' or any(entry.get(k) for k in
            ('outputs_removed','dependency_obsolete','partial_results','redocking_failures','docking_failures','md_failures')):
        return False
    return artifact_integrity(store,pid,entry,root=root)


def _input_hashes(entry):
    """Read selected-input hashes of a stage, nested batch aggregate or unit."""
    values={}
    def visit(manifest):
        if not isinstance(manifest,dict):return
        for name,data in manifest.items():
            if name.startswith('selected:') and isinstance(data,list):
                values.setdefault(name[9:],Counter()).update(sha for _,sha in data)
            elif isinstance(data,dict):visit(data)
    visit(entry.get('input_manifest',{}))
    return values


def _protocol(stage):
    params=stage.get('parameters',{})
    # These explicit protocol choices must not disappear in request defaults.
    return {key:params.get(key) for key in ('preparation_backend','amber_score_type')}, {
        key:(params.get(key) or {}).get('preparation_backend') for key in ('preparation_options',)}


def equivalent(stage,previous):
    return request(stage)==request(previous) and _protocol(stage)==_protocol(previous)


def verify(store,pid,stage,entry,snapshot,results):
    """Prove request, exact files, original inputs and current source equivalence."""
    try:
        if not equivalent(stage,entry.get('configuration',{})):return False
        root=store.project_dir(pid).resolve()
        if not intact(store,pid,entry,root):return False
    except (ValueError,KeyError,TypeError,AttributeError):return False
    prior={item['id']:item for item in snapshot}
    expected={p:sha for item in snapshot for p,sha in item.get('artifact_manifest',{}).items()}
    selected_hashes=_input_hashes(entry)
    try:
        # Materialized original inputs must still match their persisted hashes.
        # Old redocking snapshots lack input_manifest, but have producer manifests.
        def original_inputs(unit):
            if unit.get('batches'):
                if not all(original_inputs(batch) for batch in unit['batches']):return False
                indexed={str(Path(p).resolve()) for batch in unit['batches'] for paths in batch.get('input_files',{}).values() for p in paths}
                return all(str(Path(p).resolve()) in indexed for paths in unit.get('input_files',{}).values() for p in paths)
            hashes=_input_hashes(unit)
            for field,paths in unit.get('input_files',{}).items():
                actual=Counter(input_digest(_scoped(root,p)) for p in paths)
                if field in hashes:
                    if actual!=hashes[field]:return False
                elif any(expected.get(str(Path(p).resolve()))!=input_digest(_scoped(root,p)) for p in paths):return False
            return True
        if not original_inputs(entry):return False
        asset_ids=asset_references(stage)
        for asset in asset_ids:
            with store.connect() as db:
                row=db.execute('SELECT path FROM assets WHERE project_id=? AND id=?',(pid,asset)).fetchone()
            if row is None:return False
            asset_path=str(_scoped(root,row[0]).resolve())
            sha=entry.get('source_manifest',{}).get(asset)
            if sha is not None:
                if sha!=input_digest(asset_path):return False
            elif asset_path not in {str(Path(p).resolve()) for paths in entry.get('input_files',{}).values() for p in paths}:return False
        represented=set()
        for group in stage.get('bindings',{}).values():
            for ref in sources(group):
                owner=ref.get('stage')
                if not owner:continue
                represented.add(owner)
                old=prior.get(owner,{})
                old_files=old.get('artifacts',[]);new_files=results.get(owner,[])
                for path in old_files:_scoped(root,path)
                if not old_files or not new_files:return False
                if ref.get('selector','auto')!='auto':
                    previous_file=validate_reference(ref,old_files);current_file=validate_reference(ref,new_files)
                    sha=old.get('artifact_manifest',{}).get(str(Path(previous_file).resolve()))
                    if not sha or input_digest(_scoped(root,current_file))!=sha:return False
                else:
                    original=Counter((locator(owner,p),old.get('artifact_manifest',{}).get(str(Path(p).resolve()))) for p in old_files)
                    current=Counter((locator(owner,p),input_digest(_scoped(root,p))) for p in new_files)
                    if any(sha is None for _,sha in original) or original!=current:return False
        # Direct paths need original per-file evidence, not only a path string.
        from .catalog import PATH_FIELDS
        for field in PATH_FIELDS-stage.get('bindings',{}).keys():
            value=stage.get('parameters',{}).get(field)
            if not value:continue
            direct=_scoped(root,value)
            from .stage_cache import path_manifest
            actual=Counter(path_manifest(direct).values())
            recorded=Counter(input_digest(_scoped(root,p)) for p in entry.get('input_files',{}).get(field,[])) or selected_hashes.get(field)
            if recorded is None or actual!=recorded:return False
        for owner in dependencies(stage)-represented:
            old=prior.get(owner,{})
            original={str(Path(p).resolve()):sha for p,sha in old.get('artifact_manifest',{}).items()}
            if not original or Counter(original.values())!=Counter(input_digest(_scoped(root,p)) for p in results.get(owner,[])):return False
    except (OSError,ValueError,KeyError,TypeError,AttributeError):return False
    return True


def adopted(entry,run_id):
    result=copy.deepcopy(entry)
    producer=entry.get('producer_software') or entry.get('cache_software')
    result.update(reused_from_run=entry.get('reused_from_run') or run_id,reuse_kind='published_result',
        producer_software=producer,publication_run_id=run_id,
        reuse_provenance=dict(policy='verified_completed_request_v1',producer_version=entry.get('processing_version') or entry.get('adopted_processing_version'),
            producer_software=producer,source_run_id=run_id,request_signature=signature(entry['configuration'])))
    return result


def histories(store,pid,*,include_active=False):
    clause='' if include_active else " AND status NOT IN ('queued','running','awaiting_input')"
    with store.connect() as db:
        return [(row[0],json.loads(row[1])) for row in db.execute(
            'SELECT id,stages FROM runs WHERE project_id=?'+clause+' ORDER BY created DESC',(pid,))]


def catalog(store,pid,pipeline,*,include_checkpoints=False):
    """Current block results, optionally collecting verified completed checkpoints."""
    by_id={s['id']:s for s in pipeline};pending={s['id']:dependencies(s) for s in pipeline};order=[]
    while pending:
        ready=[identifier for identifier,owners in pending.items() if owners<=set(order)]
        if not ready:break
        for identifier in ready:order.append(identifier);pending.pop(identifier)
    history=histories(store,pid,include_active=include_checkpoints)
    results={};found={};projected_history={}
    for identifier in order:
        stage=by_id[identifier]
        if not stage.get('enabled',True):continue
        for run_id,snapshot in history:
            entry=next((item for item in snapshot if item['id']==identifier),None)
            if entry and retained_configuration(store,pid,stage,entry,snapshot,results) is not None:
                found[identifier]=adopted(entry,run_id);results[identifier]=list(entry['artifacts']);break
        if include_checkpoints and identifier not in found:
            from .result_visibility import project_run,retained
            from .result_state import incompatible
            parts=[];seen_units=set()
            for run_id,snapshot in history:
                if run_id not in projected_history:
                    projected_history[run_id]=project_run(dict(status='interrupted',stages=snapshot),pipeline)
                projected=projected_history[run_id]
                entry=next((s for s in projected['stages'] if s['id']==identifier),None)
                if not entry or entry.get('outputs_removed') or entry.get('dependency_obsolete'):continue
                if entry.get('configuration',{}).get('input_processing','merge')!=stage.get('input_processing','merge'):continue
                valid=[]
                for unit in entry.get('batches') or [entry]:
                    previous=unit.get('configuration',entry.get('configuration',{}))
                    current=copy.deepcopy(stage)
                    if entry.get('batches'):current['input_processing']=previous.get('input_processing','merge')
                    if not previous or incompatible(previous,current) or not retained(current,previous):continue
                    unit_identity=digest([signature(previous),{field:sorted(counts.items()) for field,counts in _input_hashes(unit).items()}])
                    if unit_identity in seen_units:continue
                    if verify(store,pid,previous,unit,snapshot,results):
                        valid.append(unit);seen_units.add(unit_identity)
                paths=list(dict.fromkeys(p for unit in valid for p in unit.get('artifacts',[])))
                if not paths:continue
                part=adopted(entry,run_id);part['coverage_pending']=True
                part['reuse_kind']='verified_checkpoint'
                provenance=part['reuse_provenance']
                provenance['policy']='verified_completed_units_v1'
                provenance['source_request_signature']=provenance.pop('request_signature')
                provenance['units']=[dict(request_signature=signature(unit['configuration']),processing_version=unit.get('processing_version'),
                    producer_software=unit.get('producer_software'),cache_key=unit.get('cache_key')) for unit in valid]
                if entry.get('batches'):part['batches']=valid
                part['_checkpoint_subset']=set(paths)!=set(entry.get('artifacts',[]))
                part['artifacts']=paths
                part['artifact_manifest']={p:sha for unit in valid for p,sha in unit.get('artifact_manifest',{}).items()}
                parts.append(part)
            if parts:
                combined=copy.deepcopy(parts[0])
                combined['artifacts']=list(dict.fromkeys(p for part in parts for p in part['artifacts']))
                combined['artifact_manifest']={p:sha for part in parts for p,sha in part['artifact_manifest'].items()}
                if len(parts)>1:combined['_result_sources']=parts
                found[identifier]=combined;results[identifier]=combined['artifacts']
    return found


def compatible_configuration(stage,entry,snapshot,results):
    """Restore a prior confirmed selector only when its entire source is unchanged."""
    from .pipeline import PipelineService
    previous=entry.get('configuration',{})
    if previous.get('operation')!=stage.get('operation'):return None
    try:configuration=PipelineService._reuse_configuration(stage,previous)
    except (ValueError,TypeError,KeyError,AttributeError):return None
    if configuration is None:return None
    prior={item['id']:item for item in snapshot}
    try:
        for field,group in stage.get('bindings',{}).items():
            refs=sources(group)
            if configuration.get('bindings',{}).get(field)==group:continue
            if not all(ref.get('stage') and ref.get('selector','auto')=='auto' for ref in refs):return None
            for ref in refs:
                owner=ref['stage'];original=prior.get(owner,{})
                expected=Counter((locator(owner,p),original.get('artifact_manifest',{}).get(str(Path(p).resolve()))) for p in original.get('artifacts',[]))
                actual=Counter((locator(owner,p),input_digest(p)) for p in results.get(owner,[]))
                if not expected or any(sha is None for _,sha in expected) or expected!=actual:return None
    except (OSError,ValueError,KeyError,TypeError):return None
    return configuration


def retained_configuration(store,pid,stage,entry,snapshot,results):
    if entry.get('status')!='succeeded' or not entry.get('configuration'):return None
    if verify(store,pid,stage,entry,snapshot,results):return copy.deepcopy(stage)
    configuration=compatible_configuration(stage,entry,snapshot,results)
    if configuration is not None and verify(store,pid,configuration,entry,snapshot,results):return configuration
    return None
