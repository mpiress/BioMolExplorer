"""Project current selections onto saved results without destroying cache/history."""
import copy
from pathlib import Path
from .bindings import sources
from .artifact_choices import matches_selector


def covered(current,previous):
    """A current automatic reference includes every old file from that origin."""
    if (current.get('stage'),current.get('asset'))!=(previous.get('stage'),previous.get('asset')):return False
    if current.get('compound_id') and current.get('compound_id')!=previous.get('compound_id'):return False
    selector=current.get('selector','auto');old=previous.get('selector','auto')
    if selector=='auto':return True
    if old=='auto':return False
    if current.get('artifact_id') and previous.get('artifact_id'):return current['artifact_id']==previous['artifact_id']
    return matches_selector(old,selector) or matches_selector(selector,old)


def retained(current,previous):
    if current is None or not current.get('enabled',True):return False
    for field,group in previous.get('bindings',{}).items():
        refs=sources(current.get('bindings',{}).get(field,{}))
        if any(not any(covered(now,old) for now in refs) for old in sources(group)):return False
    for field in ('asset_ids','files','pdb_codes','pdb_code'):
        before=previous['parameters'].get(field);after=current['parameters'].get(field)
        if before and after is not None:
            # Records have nested lists; compare their stable textual identity.
            if not set(map(str,before))<=set(map(str,after)):return False
    for field in ('input_file','filename','mol_filename','compound_id','base_input_path','base_selected_mols','base_vina_path','base_dock6_path','similarity_path'):
        before=previous['parameters'].get(field);after=current['parameters'].get(field)
        if before and after and before!=after:return False
    return True


def project_run(run,pipeline):
    """Hide removed inputs and their derived outputs; archives remain immutable."""
    if run['status'] in ('queued','running','awaiting_input'):return run
    original_files={item['id']:set(item.get('artifacts',[])) for item in run['stages']}
    run=copy.deepcopy(run);current={s['id']:s for s in pipeline};removed=set()
    for item in run['stages']:
        stage=current.get(item['id']);previous=item.get('configuration')
        if not previous:continue
        # Curation has already rewritten the authoritative artifact list and
        # tables transactionally. Comparing its old selection hides surviving
        # redocking pairs after removal of a failed sibling.
        if item.get('outputs_removed') or item.get('dependency_obsolete'):continue
        if stage is None and not (item.get('project_stage') or item.get('selection_removed')):continue
        old_files=set(item.get('artifacts',[]))
        if item.get('batches'):
            item['batches']=[b for b in item['batches'] if retained(stage,b.get('configuration',previous))]
            item['artifacts']=list(dict.fromkeys(p for b in item['batches'] for p in b.get('artifacts',[])))
        elif not retained(stage,previous):item['artifacts']=[]
        if stage and item['operation']=='dynamics_dock6':
            params=stage['parameters'];excluded=set(params.get('excluded_pairs') or []);selected=params.get('selected_pairs')
            def active(path):
                parts=Path(path).parts
                if 'pairs' not in parts:return True
                key=parts[parts.index('pairs')+1]
                return key not in excluded and (selected is None or key in selected)
            item['artifacts']=[p for p in item.get('artifacts',[]) if active(p)]
        removed.update(old_files-set(item.get('artifacts',[])))
    # Propagate exact removed file provenance through individual batches. A
    # merged computation is hidden as a unit until its new result is available.
    changed=True
    while changed:
        changed=False
        for item in run['stages']:
            old_files=set(item.get('artifacts',[]))
            def affected(entry):return bool(removed.intersection(p for paths in entry.get('input_files',{}).values() for p in paths))
            if item.get('batches'):
                item['batches']=[b for b in item['batches'] if not affected(b)]
                item['artifacts']=list(dict.fromkeys(p for b in item['batches'] for p in b.get('artifacts',[])))
            elif affected(item):item['artifacts']=[]
            missing=old_files-set(item.get('artifacts',[]))
            if missing:removed.update(missing);changed=True
    for item in run['stages']:
        item['artifact_records']=[r for r in item.get('artifact_records',[]) if r['path'] in item.get('artifacts',[])]
        if not item.get('artifacts') and item.get('status')=='succeeded' and original_files.get(item['id'],set()).intersection(removed):item['results_stale']=True
    return run


def consolidated_run(runs,pipeline,available,*,valid_result=None):
    """Read-only latest execution plus verified published results, once per block.

    Each displayed result retains its real producer run ID for authorization,
    downloads, molecular scenes and provenance. No history row is rewritten.
    """
    if not runs:return None
    latest=copy.deepcopy(runs[0]);current={s['id']:s for s in latest['stages']}
    published={s['id']:s for s in available};entries=[]
    for config in pipeline:
        identifier=config['id'];execution=current.get(identifier);saved=published.get(identifier)
        execution_intact=not execution or not execution.get('artifacts') or valid_result is None or valid_result(execution)
        if execution and execution_intact and execution.get('artifacts') and (execution.get('status')=='succeeded' or not saved or set(saved.get('artifacts',[]))<=set(execution['artifacts'])) and not execution.get('results_stale') and not execution.get('dependency_obsolete'):
            entry=copy.deepcopy(execution);producer=latest['id']
        elif saved:
            entry=copy.deepcopy(saved);producer=saved['publication_run_id']
            entry['reused']=True
        elif execution:
            entry=copy.deepcopy(execution);producer=latest['id']
        else:
            entry=dict(id=identifier,name=config['name'],operation=config['operation'],configuration=config,
                status='skipped',artifacts=[],error=None)
            producer=latest['id']
        if execution and not execution_intact:
            entry['_integrity_notice']='Os arquivos de resultado foram alterados ou estão ausentes. Consulte o histórico ou execute novamente.'
            if not saved:
                entry.update(artifacts=[],artifact_manifest={},results_stale=True,dependency_obsolete=True)
        entry['name']=config['name'];entry['_result_run_id']=producer
        if execution and producer!=latest['id']:
            entry['_execution_status']=execution['status']
        entries.append(entry)
    latest['stages']=entries
    return latest
