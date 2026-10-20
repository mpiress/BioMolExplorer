"""Audited failed-row deletion follows exact input provenance, never filenames alone."""
import csv
import io
import json
from pathlib import Path
import time
from .storage import write_text
from .stage_cache import artifact_manifest


def remove_failure(service,token,pid,rid,sid,key):
    store=service.store;backups={};removed=set();changed_files=set()
    def remember(path):
        if path not in backups:backups[path]=path.read_bytes() if path.is_file() else None
    def discard(path):
        remember(path)
        path.unlink(missing_ok=True);removed.add(str(path))
    def rewrite(path,predicate):
        path=store.scoped_path(pid,path)
        if not path.is_file():return
        with path.open(encoding='utf-8-sig',newline='') as stream:
            reader=csv.DictReader(stream);fields=reader.fieldnames;rows=list(reader)
        kept=[r for r in rows if not predicate(r)]
        if len(kept)==len(rows):return
        remember(path);output=io.StringIO(newline='');writer=csv.DictWriter(output,fieldnames=fields)
        writer.writeheader();writer.writerows(kept);write_text(output.getvalue(),path);changed_files.add(str(path))
    try:
        with store.change(token,pid,'editor','remove_failed_redocking','Complexo com falha removido do redocking e da origem') as db:
            if db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')",(pid,)).fetchone():
                raise ValueError('Aguarde a execução terminar antes de remover o complexo.')
            simulation=service.simulation(token,pid,rid,sid,key)
            if simulation['status']!='preparation_failed':raise ValueError('Somente linhas com falha de preparação podem ser excluídas por esta ação.')
            stage=service.results.stage(token,pid,rid,sid)
            metadata=store.scoped_path(pid,simulation['metadata']);code=simulation['pdb']
            batch=next((b for b in stage.get('batches',[]) if str(metadata) in b.get('artifacts',[])),stage)
            paths=list(batch.get('input_files',{}).get('base_input_path',[]))
            for record in stage.get('artifact_records',[]):
                if record['path']==str(metadata):
                    paths.extend(r['path'] for r in record.get('inputs',{}).get('base_input_path',[]) if r.get('path'))
            raw_name=code+'.pdb'
            origins={store.scoped_path(pid,p) for p in paths if Path(p).name==raw_name}
            snapshot=store.scoped_path(pid,metadata.parent/raw_name)
            if snapshot.is_file():origins.add(snapshot)
            config=batch.get('configuration',stage.get('configuration',{})).get('parameters',{})
            if not paths and config.get('base_input_path'):
                candidate=store.scoped_path(pid,Path(config['base_input_path'])/config.get('target','MeuAlvo').replace(' ','')/raw_name)
                if candidate.is_file():origins.add(candidate)
            if not origins:raise ValueError('A origem do complexo não está disponível para exclusão segura.')
            runs=[(r['id'],json.loads(r['stages'])) for r in db.execute('SELECT id,stages FROM runs WHERE project_id=?',(pid,))]
            origin_strings={str(p) for p in origins};producer_ids=set()
            for _,stages in runs:
                for item in stages:
                    if origin_strings.intersection(item.get('artifacts',[])) and item['operation']=='retrieve_structures':producer_ids.add(item['id'])
            # Delete only registered source/snapshot files, plus verified copies
            # in the disposable input cache. Unrelated collections are untouched.
            from .project_state import digest as file_digest
            if snapshot.is_file() and any(p.is_file() and file_digest(p)!=file_digest(snapshot) for p in origins):
                raise ValueError('O complexo de origem foi alterado. Reabra os resultados antes de excluir.')
            digests={file_digest(p) for p in origins if p.is_file()}
            cache=store.project_dir(pid)/'.input-cache'
            if cache.is_dir():
                for p in cache.rglob(raw_name):
                    p=store.scoped_path(pid,p)
                    if p.is_file() and file_digest(p) in digests:origins.add(p)
            for p in origins:discard(p)
            identity=(code,simulation['ligand'],simulation['residue'],simulation['chain'])
            def same_row(r):return tuple(r.get(k,'') for k in ('PDB_CODE','LIGAND','RESNUM','CHAIN'))==identity
            rewrite(metadata,same_row)
            for source in origins:
                if source.parent!=metadata.parent:
                    for name in ('pdb_codes.csv','ligand_candidates.csv'):
                        rewrite(source.parent/name,lambda r:r.get('PDB_CODE','').upper()==code.upper())
            for item in stage.get('artifacts',[]):
                path=store.scoped_path(pid,item)
                if path.name=='redocking_failures.csv' and path.is_relative_to(service.root(metadata)):rewrite(path,same_row)
            asset_ids={r['id'] for r in db.execute('SELECT id,path FROM assets WHERE project_id=?',(pid,)) if str(store.scoped_path(pid,r['path'])) in removed}
            for asset in asset_ids:db.execute('DELETE FROM assets WHERE id=? AND project_id=?',(asset,pid))
            project=db.execute('SELECT pipeline FROM projects WHERE id=?',(pid,)).fetchone()
            pipeline=json.loads(project['pipeline'])
            for item in pipeline:
                params=item['parameters']
                if item['id'] in producer_ids:
                    exclusions=params.get('excluded_pdb_ids') or []
                    from .retrieval import identifiers
                    exclusions=identifiers(exclusions,'pdb') if exclusions else []
                    params['excluded_pdb_ids']=list(dict.fromkeys(exclusions+[code]))
                # Every pair from this physical PDB is unavailable for another
                # execution; preserve completed historical simulation outputs.
                if item['id']==sid and params.get('pdb_codes'):
                    params['pdb_codes']=[r for r in params['pdb_codes'] if str(r[0]).upper()!=code.upper()]
                    keys={'|'.join(str(v) for v in r[:4]) for r in params['pdb_codes']}
                    params['preparation_pairs']={k:v for k,v in (params.get('preparation_pairs') or {}).items() if k in keys}
                if asset_ids and item['operation']=='import_results':
                    params['asset_ids']=[a for a in params.get('asset_ids',[]) if a not in asset_ids]
                    params['asset_types']={a:k for a,k in params.get('asset_types',{}).items() if a not in asset_ids}
            db.execute('UPDATE projects SET pipeline=?,revision=revision+1,updated=? WHERE id=?',(json.dumps(pipeline),time.time(),pid))
            for run_id,stages in runs:
                changed=False
                for item in stages:
                    affected=(removed|changed_files).intersection(item.get('artifacts',[]))
                    if affected:
                        item['artifacts']=[p for p in item['artifacts'] if p not in removed]
                        item['artifact_records']=[r for r in item.get('artifact_records',[]) if r['path'] not in removed]
                        item['artifact_manifest']=artifact_manifest(item['artifacts'])
                        item['outputs_removed']=True;item['curated_at']=time.time();changed=True
                        if item['id']==sid:item['redocking_failures']=max(0,item.get('redocking_failures',1)-1)
                        for b in item.get('batches',[]):b['artifacts']=[p for p in b.get('artifacts',[]) if p not in removed]
                if changed:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),run_id))
    except Exception:
        for path,data in backups.items():
            if data is None:path.unlink(missing_ok=True)
            else:path.write_bytes(data)
        raise
    return dict(removed_pdb=code,removed_files=len(removed))
