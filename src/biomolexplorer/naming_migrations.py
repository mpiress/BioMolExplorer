"""Conservative metadata-only repairs for already published scientific results."""
import copy
import json
import time
from pathlib import Path


def repair_artifact_lineage_names(store, project_id=None):
    """Restore complete semantic lineage in existing analytical public names.

    Artifact bytes, physical paths, identities and scientific results are not
    changed.  A recovery document is written before publication manifests and
    run metadata are updated.
    """
    from .artifacts import index_run
    from .storage import write_json
    with store.connect() as db:
        values=() if project_id is None else (project_id,)
        where='' if project_id is None else ' WHERE project_id=?'
        rows=[dict(row) for row in db.execute('SELECT id,project_id,stages FROM runs'+where+' ORDER BY created',values)]
        assets={pid:[dict(item) for item in db.execute('SELECT id,name,path FROM assets WHERE project_id=?',(pid,))]
            for pid in {row['project_id'] for row in rows}}
    changes=[];manifest_original={};updated_runs={}
    for row in rows:
        stages=json.loads(row['stages'])
        if not any(stage.get('artifact_records') for stage in stages):continue
        candidate=copy.deepcopy(stages)
        try:index_run(candidate,row['id'],assets.get(row['project_id'],[]))
        except (OSError,ValueError,KeyError,TypeError):continue
        recalculated={(stage['id'],record['path']):record for stage in candidate
            for record in stage.get('artifact_records',[])}
        changed=False
        for stage in stages:
            for record in stage.get('artifact_records',[]):
                replacement=recalculated.get((stage['id'],record['path']))
                if not replacement or replacement.get('public_name')==record.get('public_name'):continue
                old={key:copy.deepcopy(record.get(key)) for key in ('public_name','display_name','naming')}
                for key in ('public_name','display_name','naming'):record[key]=copy.deepcopy(replacement.get(key))
                path=Path(record['path']);manifest=path.parent/'.artifact-publication.json'
                if manifest.is_file():
                    manifest_original.setdefault(str(manifest),json.loads(manifest.read_text(encoding='utf-8')))
                changes.append({'project_id':row['project_id'],'run_id':row['id'],'stage_id':stage['id'],
                    'occurrence_id':record.get('occurrence_id'),'path':record['path'],'old':old,
                    'new':{key:copy.deepcopy(record.get(key)) for key in ('public_name','display_name','naming')}})
                changed=True
        if changed:updated_runs[row['id']]=(row['project_id'],stages)
    if not changes:return {'changed':0,'runs':0,'backup':None}
    stamp=time.strftime('%Y%m%d-%H%M%S');projects=sorted({item['project_id'] for item in changes})
    backups=[]
    for pid in projects:
        folder=store.project_dir(pid)/'.history'/'naming-migrations';folder.mkdir(parents=True,exist_ok=True)
        backup=folder/f'artifact-lineage-{stamp}.json'
        write_json({'schema':'biomolexplorer.naming-repair','version':1,'created':time.time(),
            'changes':[item for item in changes if item['project_id']==pid],
            'publication_manifests':{path:value for path,value in manifest_original.items()
                if Path(path).is_relative_to(store.project_dir(pid))}},backup)
        backups.append(str(backup))
    touched=[]
    try:
        for manifest_name,original in manifest_original.items():
            manifest=Path(manifest_name);updated=copy.deepcopy(original)
            for item in changes:
                if Path(item['path']).parent!=manifest.parent:continue
                entry=updated.setdefault(Path(item['path']).name,{})
                entry.update(public_name=item['new']['public_name'],naming=item['new']['naming'])
            write_json(updated,manifest);touched.append(manifest)
        with store.connect() as db:
            db.execute('BEGIN IMMEDIATE')
            for run_id,(_,stages) in updated_runs.items():
                db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages,ensure_ascii=False),run_id))
            for item in changes:
                if not item.get('occurrence_id'):continue
                row=db.execute('SELECT metadata FROM artifact_instances WHERE project_id=? AND indexed_run_id=? AND occurrence_id=?',
                    (item['project_id'],item['run_id'],item['occurrence_id'])).fetchone()
                if row is None:continue
                metadata=json.loads(row[0]);metadata.update(public_name=item['new']['public_name'],naming=item['new']['naming'])
                db.execute('UPDATE artifact_instances SET metadata=? WHERE project_id=? AND indexed_run_id=? AND occurrence_id=?',
                    (json.dumps(metadata,sort_keys=True),item['project_id'],item['run_id'],item['occurrence_id']))
    except BaseException:
        for manifest in touched:write_json(manifest_original[str(manifest)],manifest)
        raise
    return {'changed':len(changes),'runs':len(updated_runs),'backup':backups}


def repair_pubchem_lineage_names(store, project_id=None):
    """Backward-compatible entry point for the broader lineage repair."""
    return repair_artifact_lineage_names(store,project_id)
