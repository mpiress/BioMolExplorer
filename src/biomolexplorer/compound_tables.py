"""Authorized, paginated compound tables and atomic, audited dataset curation."""
import csv
import hashlib
import io
import json
import os
import shutil
import time
from pathlib import Path
from tempfile import NamedTemporaryFile
from uuid import uuid4

from .stage_cache import artifact_manifest
from .workspace import AccessDenied


class CompoundTables:
    def __init__(self,store):
        self.store=store

    def _stage(self,token,project_id,run_id,stage_id):
        run=self.store.get_run(token,run_id)
        if run['project_id']!=project_id:
            raise AccessDenied('Resultado de outro projeto.')
        stage=next((s for s in run['stages'] if s['id']==stage_id and s['status']=='succeeded'),None)
        if stage is None:
            raise ValueError('A tabela estará disponível após concluir a recuperação.')
        return run,stage

    def tables(self,token,project_id,run_id,stage_id):
        _,stage=self._stage(token,project_id,run_id,stage_id)
        result=[]
        for filename in stage.get('artifacts',[]):
            path=self.store.scoped_path(project_id,filename)
            name=path.name
            is_summary=path.suffix.lower()=='.csv' and path.stem.upper().endswith(('_FULL','_MOLS','_SIMS'))
            is_integrated=name=='compounds.csv' and path.parent.parent.name=='compounds'
            is_provided=bool(stage.get('configuration',{}).get('provided_results')) and path.suffix.lower()=='.csv'
            if (is_summary or is_integrated or is_provided) and path.is_file():
                result.append({'path':str(path),'name':name,'integrated':is_integrated})
        result.sort(key=lambda t:(not t['integrated'],t['name']))
        return result

    def _path(self,token,project_id,run_id,stage_id,filename):
        choices=self.tables(token,project_id,run_id,stage_id)
        path=self.store.scoped_path(project_id,filename)
        if str(path) not in {t['path'] for t in choices}:
            raise AccessDenied('Tabela não autorizada para esta etapa.')
        return path

    @staticmethod
    def _digest(stream):
        digest=hashlib.sha256()
        for chunk in iter(lambda:stream.read(1024*1024),b''):
            digest.update(chunk)
        stream.seek(0)
        return digest.hexdigest()

    def page(self,token,project_id,run_id,stage_id,filename,offset=0,limit=25,query=''):
        if type(offset) is not int or offset<0 or type(limit) is not int or not 1<=limit<=100:
            raise ValueError('Página de tabela inválida.')
        if not isinstance(query,str) or len(query)>200:
            raise ValueError('Use até 200 caracteres na busca.')
        path=self._path(token,project_id,run_id,stage_id,filename)
        rows=[];matched=total=0
        with path.open('rb') as source:
            version=self._digest(source)
            reader=csv.DictReader(io.TextIOWrapper(source,encoding='utf-8-sig',newline=''))
            if not {'molecule_chembl_id','canonical_smiles'}.issubset(reader.fieldnames or []):
                raise ValueError('Esta tabela não contém código do composto e SMILES.')
            for index,row in enumerate(reader):
                total+=1
                identifier=row['molecule_chembl_id'] or ''
                smiles=row['canonical_smiles'] or ''
                if query.lower() not in (identifier+' '+smiles).lower():
                    continue
                if offset<=matched<offset+limit:
                    rows.append({'index':index,'id':identifier,'smiles':smiles})
                matched+=1
        project=self.store.project(token,project_id)
        with self.store.connect() as db:
            active=db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')",(project_id,)).fetchone()
        return {'rows':rows,'total':total,'matched':matched,'version':version,'offset':offset,'limit':limit,'name':path.name,
                'can_edit':project['role']!='viewer' and active is None}

    def remove(self,token,project_id,run_id,stage_id,filename,index,expected_version):
        self.store.project(token,project_id,'editor')
        if type(index) is not int or index<0 or not isinstance(expected_version,str):
            raise ValueError('Seleção inválida.')
        path=self._path(token,project_id,run_id,stage_id,filename)
        temporary=None;backup=None
        with self.store.connect() as db:
            db.execute('BEGIN IMMEDIATE')
            user=self.store.user(token)
            self.store._require_user(user['id'],project_id,'editor')
            from .project_state import snapshot,record
            before=snapshot(self.store,db,project_id)
            if db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')",(project_id,)).fetchone():
                raise ValueError('Aguarde a execução terminar antes de remover compostos.')
            try:
                with path.open('rb') as source:
                    if self._digest(source)!=expected_version:
                        raise ValueError('A tabela foi alterada. Atualize os dados antes de remover o composto.')
                    reader=csv.DictReader(io.TextIOWrapper(source,encoding='utf-8-sig',newline=''))
                    if not {'molecule_chembl_id','canonical_smiles'}.issubset(reader.fieldnames or []):
                        raise ValueError('Esta tabela não contém código do composto e SMILES.')
                    removed=None
                    with NamedTemporaryFile(mode='w',encoding='utf-8',newline='',dir=path.parent,delete=False,suffix='.tmp') as out:
                        temporary=Path(out.name)
                        writer=csv.DictWriter(out,fieldnames=reader.fieldnames)
                        writer.writeheader()
                        for row_index,row in enumerate(reader):
                            if row_index==index:
                                removed=row
                            else:
                                writer.writerow(row)
                        out.flush();os.fsync(out.fileno())
                    if removed is None:
                        raise ValueError('O composto não está mais nesta tabela.')
                # Retain the original bytes and an audit trail outside scientific
                # artifacts; failures restore the original CSV.
                backup_dir=self.store.project_dir(project_id)/'.curation'
                backup_dir.mkdir(mode=0o700,exist_ok=True)
                edit_id=uuid4().hex
                backup=backup_dir/(edit_id+'.csv')
                shutil.copyfile(path,backup)
                temporary.replace(path)
                db.execute('CREATE TABLE IF NOT EXISTS compound_edits(id TEXT PRIMARY KEY, project_id TEXT, user_id TEXT, run_id TEXT, stage_id TEXT, path TEXT, compound_id TEXT, backup TEXT, created REAL)')
                db.execute('INSERT INTO compound_edits VALUES (?,?,?,?,?,?,?,?,?)',
                    (edit_id,project_id,user['id'],run_id,stage_id,str(path),removed.get('molecule_chembl_id'),str(backup),time.time()))
                # Reused runs reference the same files. Refresh their manifests so
                # retrieval stays reusable; downstream input hashes invalidate.
                for row in db.execute('SELECT id,stages FROM runs WHERE project_id=?',(project_id,)).fetchall():
                    stages=json.loads(row[1]);changed=False
                    for item in stages:
                        if str(path) in item.get('artifacts',[]):
                            item['artifact_manifest']=artifact_manifest(item['artifacts'])
                            item['curated_at']=time.time();item['curated_by']=user['id'];changed=True
                    if changed:
                        db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),row[0]))
                record(self.store,db,project_id,user['id'],'curation',f'Composto {removed.get("molecule_chembl_id")} removido de {path.name}',before)
                db.commit()
                return {'id':removed.get('molecule_chembl_id'),'table':path.name}
            except Exception:
                if backup is not None and backup.exists():
                    with NamedTemporaryFile(dir=path.parent,delete=False,suffix='.tmp') as restored:
                        restore_path=Path(restored.name)
                    try:
                        shutil.copyfile(backup,restore_path)
                        restore_path.replace(path)
                    finally:
                        restore_path.unlink(missing_ok=True)
                raise
            finally:
                if temporary is not None:
                    temporary.unlink(missing_ok=True)
