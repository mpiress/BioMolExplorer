"""Plan and atomically reconcile exclusive results after confirmed graph edits."""
import copy
import json
from pathlib import Path
import shutil
from uuid import uuid4
from collections import defaultdict,deque

from .bindings import sources,pack
from .clear_results import unregistered_original
from .result_state import incompatible,signature,digest,inputs_changed,input_digest,removed_axes,fingerprint_units,preparation


def refs(stage):
    bindings={} if stage.get('provided_results') else stage.get('bindings',{})
    references=[r for group in bindings.values() for r in sources(group)]
    imported=stage.get('provided_results') or (stage.get('parameters',{}) if stage.get('operation')=='import_results' else {})
    return references+[{'asset':asset} for asset in imported.get('asset_ids',[])]


def retained(reference,current):
    for other in refs(current):
        if (reference.get('stage'),reference.get('asset'))!=(other.get('stage'),other.get('asset')):continue
        if reference.get('compound_id')!=other.get('compound_id'):continue
        if reference.get('selector','auto')=='auto':return True
        if other.get('selector','auto')=='auto':return True
        if reference.get('artifact_id') and reference.get('artifact_id')==other.get('artifact_id'):return True
        if reference.get('selector','auto')==other.get('selector','auto'):return True
    return False


def pair_removals(entry,current):
    """Keep verified files of surviving pairs, including shared receptors.

    Aggregate tables become obsolete. The original immutable pair index remains
    provenance; seed readers filter it to the surviving current pairs.
    """
    from .stage_cache import file_digest
    from .redocking_config import pair_key,effective_preparation
    previous=entry.get('configuration',{});before=effective_preparation(previous.get('parameters',{}))
    after=effective_preparation(current.get('parameters',{}));selected={pair_key(r) for r in current['parameters'].get('pdb_codes',[]) or []}
    keep=set();indexes=set();found=False
    for filename in entry.get('artifacts',[]):
        path=Path(filename)
        if path.name!='redocking_cache.json':continue
        try:
            if file_digest(path)!=entry.get('artifact_manifest',{}).get(str(path.resolve())):continue
            value=json.loads(path.read_text())
            if value.get('version')!=1:continue
            found=True
            for pair in value['entries']:
                key=pair_key(pair['record'])
                if key not in selected or preparation(before.get(key),previous.get('parameters',{}).get('charge_type','gas'))!=preparation(after.get(key),current['parameters'].get('charge_type','gas')):continue
                paths={str((path.parent/name).resolve()) for name in pair['files']}
                # Index paths must remain within their original artifact tree.
                if any(not Path(p).is_relative_to(path.parent.resolve()) for p in paths):continue
                keep.update(paths);indexes.add(str(path))
        except (OSError,ValueError,KeyError,TypeError):continue
    lost=set(entry.get('artifacts',[]))-(keep|indexes) if found else set(entry.get('artifacts',[]))
    entry['partial_results']=bool(keep)
    entry['artifacts']=[p for p in entry.get('artifacts',[]) if p not in lost]
    entry['artifact_manifest']={p:sha for p,sha in entry.get('artifact_manifest',{}).items() if p not in lost}
    entry['artifact_records']=[r for r in entry.get('artifact_records',[]) if r['path'] not in lost]
    entry.pop('cache_key',None)
    return lost


def plan(db,pid,pipeline,project_root=None):
    pipeline=copy.deepcopy(pipeline)
    current={s['id']:s for s in pipeline};runs=[];candidates=set();protected=set();affected=0;lost_owners={}
    changed_blocks={};active_conflicts=set();signatures={};selection_owners=set()
    def discard(item,lost):
        nonlocal affected
        if not lost:return
        affected+=1;candidates.update(lost)
        changed_blocks[item['id']]=current.get(item['id'],{}).get('name',item.get('name',item['id']))
        for path in lost:lost_owners.setdefault(path,set()).add(item['id'])
        item['artifacts']=[p for p in item.get('artifacts',[]) if p not in lost]
        item['artifact_manifest']={p:d for p,d in item.get('artifact_manifest',{}).items() if p not in lost}
        item['artifact_records']=[r for r in item.get('artifact_records',[]) if r['path'] not in lost]
        item['outputs_removed']=not item['artifacts'];item['reused']=False
        item['dependency_obsolete']=True;item.pop('cache_key',None)
    active_ids=set()
    for row in db.execute('SELECT id,stages,status FROM runs WHERE project_id=?',(pid,)):
        items=json.loads(row[1]);runs.append((row[0],items))
        if row[2] in ('queued','running','awaiting_input'):active_ids.add(row[0])
        for item in items:
            if item.get('outputs_removed') and not item.get('artifacts'):continue
            config=item.get('configuration',{});new=current.get(item['id'])
            if not config.get('operation'):continue
            provided=unregistered_original(item)
            if provided:
                protected.update(item.get('artifacts',[]))
                if new is not None:continue
            parameter_change=bool(new is not None and incompatible(config,new,include_axes=False))
            changed_inputs=inputs_changed(item) or any(inputs_changed(batch) for batch in item.get('batches',[]))
            changed_assets=set()
            for asset,sha in item.get('source_manifest',{}).items():
                row_asset=db.execute('SELECT path FROM assets WHERE project_id=? AND id=?',(pid,asset)).fetchone()
                try:valid=row_asset is not None and input_digest(row_asset[0])==sha
                except OSError:valid=False
                if not valid:changed_assets.add(asset)
            removed=bool(new is None or incompatible(config,new) or changed_inputs or changed_assets or any(not retained(r,new) for r in refs(config)))
            if not removed:continue
            if new is None or removed_axes(config,new) or any(not retained(r,new) for r in refs(config)):
                selection_owners.add(item['id'])
            if item.get('artifacts'):
                signatures[item['id']]=None if new is None else signature(new)
            lost=set()
            if new is not None:
                old_import=config.get('provided_results') or (config.get('parameters',{}) if config.get('operation')=='import_results' else {})
                new_import=new.get('provided_results') or (new.get('parameters',{}) if new.get('operation')=='import_results' else {})
                removed_assets=(set(old_import.get('asset_ids',[]))-set(new_import.get('asset_ids',[])))|changed_assets
                names=set()
                for asset in removed_assets:
                    source=db.execute('SELECT path FROM assets WHERE project_id=? AND id=?',(pid,asset)).fetchone()
                    if source:
                        filename=Path(source[0]);names.update((filename.name,filename.with_suffix('.biomol-view.json').name))
                # Import validation guarantees distinct source filenames. These
                # are normalized copies; library originals remain protected.
                lost.update(p for p in item.get('artifacts',[]) if Path(p).name in names)
            if new is None:
                lost.update(item.get('artifacts',[]));item['batches']=[]
                if project_root:
                    directory=Path(project_root)/'runs'/row[0]/item['id']
                    for name in ('artifacts','inputs','redocking-inputs'):
                        for folder in directory.rglob(name) if directory.is_dir() else []:
                            if folder.is_dir() and not folder.is_symlink():
                                lost.update(str(p.resolve()) for p in folder.rglob('*') if p.is_file() and p.resolve().is_relative_to(Path(project_root).resolve()))
            elif item.get('batches'):
                batches=[]
                units=[unit for batch in item['batches'] for unit in fingerprint_units(batch)]
                for batch in units:
                    old_config=batch.get('configuration',config);oldrefs=refs(old_config)
                    unit=copy.deepcopy(new)
                    same_mode=config.get('input_processing','merge')==new.get('input_processing','merge')
                    if same_mode:unit['input_processing']=old_config.get('input_processing','merge')
                    invalid=not same_mode or incompatible(old_config,unit) or inputs_changed(batch) or oldrefs and any(not retained(r,new) for r in oldrefs)
                    if invalid:
                        partial=(config['operation']=='redocking' and same_mode and not incompatible(old_config,unit,include_axes=False) and not inputs_changed(batch) and all(retained(r,new) for r in oldrefs))
                        lost.update(pair_removals(batch,unit) if partial else batch.get('artifacts',[]))
                        if partial and batch.get('artifacts'):batches.append(batch)
                    else:batches.append(batch)
                item['batches']=batches
                lost.difference_update(p for batch in batches for p in batch.get('artifacts',[]))
                if any(b.get('partial_results') for b in batches):item['partial_results']=True
            elif config['operation']=='redocking' and not parameter_change and not changed_inputs and all(retained(r,new) for r in refs(config)):
                lost.update(pair_removals(item,new))
            elif parameter_change or changed_inputs or refs(config) and not any(retained(r,new) for r in refs(config)) or config['operation'] not in ('import_results','redocking','prepare_structures') and not config.get('provided_results') and incompatible(config,new):
                lost.update(item.get('artifacts',[]))
            if row[0] in active_ids and lost:active_conflicts.add(item['id'])
            discard(item,lost)
            # A partial native index can seed valid pairs, never a whole stage.
            item['dependency_obsolete']=True;item.pop('cache_key',None)
    # Follow actual file lineage, rather than deleting every descendant of a
    # changed block. Unaffected batches and direct independent file users survive.
    from .bindings import dependencies
    consumers=defaultdict(list);legacy_consumers=defaultdict(list)
    for run_id,items in runs:
        for item in items:
            if unregistered_original(item):continue
            for unit in item.get('batches') or [item]:
                config=unit.get('configuration',item.get('configuration',{}));inputs=unit.get('input_files',{})
                declared=dependencies(config);represented={r.get('stage') for r in refs(config) if r.get('stage')}
                for paths in inputs.values():
                    for path in paths:
                        for owner in declared:consumers[owner,path].append((run_id,item,unit))
                for owner in declared-represented if any(inputs.values()) else declared:
                    legacy_consumers[owner].append((run_id,item,unit))
    events=deque((owner,path) for path,owners in lost_owners.items() for owner in owners)
    seen_events=set(events);seen_owners=set();invalid_units=set()
    while events:
        owner,path=events.popleft();users=list(consumers.get((owner,path),[]))
        if owner not in seen_owners:
            users.extend(legacy_consumers.get(owner,[]));seen_owners.add(owner)
        for run_id,item,unit in users:
            if id(unit) in invalid_units:continue
            invalid_units.add(id(unit));lost=set(unit.get('artifacts',[]))
            if not lost:continue
            if run_id in active_ids:active_conflicts.add(item['id'])
            if owner in selection_owners:selection_owners.add(item['id'])
            discard(item,lost)
            for filename in lost:
                event=(item['id'],filename)
                if event not in seen_events:seen_events.add(event);events.append(event)
    for _,items in runs:
        for item in items:
            if item.get('batches'):item['batches']=[b for b in item['batches'] if id(b) not in invalid_units]
    from .artifact_choices import matches_selector
    from .artifacts import slot_id,locator
    def removed_reference(reference):
        if reference.get('selector','auto')=='auto':return False
        owner=reference.get('stage')
        if owner not in selection_owners:return False
        return any(owner in owners and (reference.get('artifact_id')==slot_id(owner,locator(owner,path)) or
                   matches_selector(path,reference.get('selector',''))) for path,owners in lost_owners.items())
    for stage in pipeline:
        for field,group in list(stage.get('bindings',{}).items()):
            old=sources(group);kept=[r for r in old if not removed_reference(r)]
            if old==kept:continue
            if not kept:kept=[dict(stage=r['stage'],selector='auto') for r in old if r.get('stage') in current]
            value=pack(kept)
            if value:stage['bindings'][field]=value
            else:stage['bindings'].pop(field,None)
            for key in ('files','filename','input_file','mol_filename','graph_inputs','dataset_origin'):stage.get('parameters',{}).pop(key,None)
    for _,items in runs:
        for item in items:
            if not item.get('dependency_obsolete'):
                protected.update(item.get('artifacts',[]))
                for paths in item.get('input_files',{}).values():protected.update(paths)
            for batch in item.get('batches',[]):
                protected.update(batch.get('artifacts',[]))
                for paths in batch.get('input_files',{}).values():protected.update(paths)
    protected.update(row[0] for row in db.execute('SELECT path FROM assets WHERE project_id=?',(pid,)))
    blocks=[dict(id=key,name=value) for key,value in changed_blocks.items()]
    confirmation=digest(dict(pipeline=[signature(s) for s in pipeline],removed=sorted(candidates),blocks=blocks))
    return dict(runs=runs,pipeline=pipeline,paths=sorted(candidates-protected),affected=affected,blocks=blocks,
                active_conflicts=sorted(active_conflicts),confirmation=confirmation,signatures=signatures,
                shared_preserved=len(candidates&protected))


def apply(store,db,pid,proposal,moved,journal):
    root=store.project_dir(pid,db).resolve();backup=root/'.history'/('.sync-'+uuid4().hex)
    selected=[]
    for value in proposal['paths']:
        path=Path(value).resolve()
        if not path.is_relative_to(root) or not path.is_file():continue
        selected.append((path,backup/str(len(selected))))
    journal.prepare(selected)
    for path,saved in selected:
        backup.mkdir(parents=True,exist_ok=True)
        journal.move(path,saved);moved.append((path,saved))
    journal.flush()
    for rid,items in proposal['runs']:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(items),rid))


def finish(moved,rollback=False):
    for path,saved in reversed(moved):
        if rollback and saved.exists():path.parent.mkdir(parents=True,exist_ok=True);saved.replace(path)
    for parent in {saved.parent for _,saved in moved}:shutil.rmtree(parent,ignore_errors=True)
