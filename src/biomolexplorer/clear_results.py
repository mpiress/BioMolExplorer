"""Audited project result cleanup, keeping source assets and pipeline edges."""
import copy
import json
from pathlib import Path
import shutil
import time
from uuid import uuid4

from .bindings import sources,pack
from .workspace import AccessDenied


def reset_selections(pipeline):
    pipeline=copy.deepcopy(pipeline);count=0
    for stage in pipeline:
        bindings=stage.get('bindings',{})
        for field,group in list(bindings.items()):
            refs=[]
            for ref in sources(group):
                if 'stage' in ref:refs.append(dict(stage=ref['stage'],selector='auto'))
                elif 'asset' in ref:refs.append(dict(asset=ref['asset']))
            updated=pack(refs)
            if updated!=group:count+=1
            if updated is None:bindings.pop(field)
            else:bindings[field]=updated
        params=stage.get('parameters',{})
        if stage['operation'] in ('redocking','prepare_structures','docking_vina','docking_dock6'):
            for key in ('pdb_codes','pdb_code','preparation_pairs'):params.pop(key,None)
        if not bindings:continue
        for field in bindings:params.pop(field,None)
        for key in ('input_file','files','filename','mol_filename','compound_id','graph_inputs','dataset_origin'):
            params.pop(key,None)
        if 'selection_mode' in params:params['selection_mode']='all'
        for key in ('input_expansions','input_files','artifact_records'):stage.pop(key,None)
    return pipeline,count


def clear_project_results(store,token,pid):
    """Commit file removal and selection reset together, restoring both on failure."""
    moved=[];backup=None;files_removed=0;stages_cleared=0
    try:
        with store.change(token,pid,'editor','clear_all_results','Resultados e seleções de entrada limpos') as db:
            if db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')",(pid,)).fetchone():
                raise ValueError('Aguarde ou cancele a execução antes de limpar os resultados.')
            root=store.project_dir(pid,db).resolve()
            def scoped(value):
                path=Path(value)
                path=(root/path).resolve() if not path.is_absolute() else path.resolve()
                return path if path.is_relative_to(root) and path!=root else None
            runs=[(row[0],json.loads(row[1])) for row in db.execute('SELECT id,stages FROM runs WHERE project_id=?',(pid,))]
            protected={p for row in db.execute('SELECT path FROM assets WHERE project_id=?',(pid,)) if (p:=scoped(row[0]))}
            def provided(item):return item.get('operation')=='import_results' or item.get('provided') or bool(item.get('configuration',{}).get('provided_results'))
            for _,stages in runs:
                for item in stages:
                    if provided(item):protected.update(p for value in item.get('artifacts',[]) if (p:=scoped(value)))
            candidates=set()
            for rid,stages in runs:
                for item in stages:
                    if provided(item):continue
                    stages_cleared+=1
                    candidates.update(p for value in item.get('artifacts',[]) if (p:=scoped(value)))
                    stage_root=root/'runs'/rid/item['id']
                    if stage_root.is_dir():
                        for name in ('artifacts','inputs','redocking-inputs'):
                            candidates.update(p.resolve() for p in stage_root.rglob(name) if p.is_dir() and not p.is_symlink())
                    item.update(artifacts=[],artifact_manifest={},artifact_records=[],input_files={},batches=[],
                        outputs_removed=True,curated_at=time.time(),reused=False,redocking_failures=0,redocking_reused=0,excluded_records=0)
                    for key in ('cache_key','cache_software','reused_from_run','current_batch','progress'):item.pop(key,None)
            if (root/'.input-cache').is_dir():candidates.add(root/'.input-cache')
            owned=set()
            for path in candidates:
                if not path.is_relative_to(root) or path in protected:continue
                if path.is_dir() and any(p.is_relative_to(path) for p in protected):
                    owned.update(p for p in path.rglob('*') if p.is_file() and p.resolve().is_relative_to(root) and p.resolve() not in protected)
                else:owned.add(path)
            # Move parent directories once; never move the same bytes twice.
            selected=[]
            for path in sorted(owned,key=lambda p:(len(p.parts),str(p))):
                if path.exists() and not any(path.is_relative_to(parent) for parent in selected if parent.is_dir()):selected.append(path)
            backup=root/'.history'/('.clear-'+uuid4().hex);backup.mkdir(parents=True)
            for index,path in enumerate(selected):
                files_removed+=sum(p.is_file() for p in path.rglob('*')) if path.is_dir() else 1
                destination=backup/str(index);path.replace(destination);moved.append((path,destination))
            pipeline=json.loads(db.execute('SELECT pipeline FROM projects WHERE id=?',(pid,)).fetchone()[0])
            pipeline,links_reset=reset_selections(pipeline)
            # Run timestamps describe execution, not the later cleanup duration.
            for rid,stages in runs:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),rid))
            db.execute('UPDATE projects SET pipeline=?,revision=revision+1,updated=? WHERE id=?',(json.dumps(pipeline),time.time(),pid))
    except BaseException:
        for path,saved in reversed(moved):
            if saved.exists():path.parent.mkdir(parents=True,exist_ok=True);saved.replace(path)
        if backup:shutil.rmtree(backup,ignore_errors=True)
        raise
    if backup:shutil.rmtree(backup,ignore_errors=True)
    return dict(files_removed=files_removed,stages_cleared=stages_cleared,links_reset=links_reset)
