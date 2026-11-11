"""Authorized, paginated compound tables and atomic, audited dataset curation."""
import csv
import hashlib
import io
import json
import os
import re
import shutil
import time
from pathlib import Path
from tempfile import NamedTemporaryFile
from uuid import uuid4
from contextlib import contextmanager
from types import SimpleNamespace
from itertools import chain
from .table_sort import ordered,is_numeric

from .stage_cache import artifact_manifest
from .workspace import AccessDenied


def compound_links(row):
    """Build provider links only from recognized compound identifiers."""
    identifier = str(row.get('molecule_chembl_id') or '').strip().upper()
    links = []
    if re.fullmatch(r'CHEMBL[1-9][0-9]*', identifier):
        links.append({'provider': 'ChEMBL', 'url': 'https://www.ebi.ac.uk/chembl/explore/compound/' + identifier})
    pubchem = re.fullmatch(r'PUBCHEM([1-9][0-9]*)', identifier)
    cid = pubchem.group(1) if pubchem else str(row.get('PubChem_CID') or '').strip()
    # CSVs produced by pandas can serialize an optional numeric CID as 2244.0.
    match = re.fullmatch(r'([1-9][0-9]*)(?:\.0+)?', cid)
    if match:
        links.append({'provider': 'PubChem', 'url': 'https://pubchem.ncbi.nlm.nih.gov/compound/' + match.group(1)})
    return links


class CompoundTables:
    requires_compounds=True
    available_statuses=('succeeded',)
    def __init__(self,store):
        self.store=store
        self._snapshot=None
        self._view=None
        self._stage_snapshot=None
        self._table_access=None

    @contextmanager
    def _reader(self,path):
        stat=path.stat();stamp=(str(path),stat.st_dev,stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns)
        cached=self._snapshot
        if cached and cached[0]==stamp:
            yield cached[1],SimpleNamespace(fieldnames=cached[2],rows=cached[3])
            return
        with path.open('rb') as source:
            version=self._digest(source)
            reader=csv.DictReader(io.TextIOWrapper(source,encoding='utf-8-sig',newline=''))
            if stat.st_size<=8*1024*1024:
                rows=[]
                for row in reader:
                    rows.append(row)
                    if len(rows)>50000:
                        self._snapshot=None
                        yield version,SimpleNamespace(fieldnames=reader.fieldnames,rows=chain(rows,reader))
                        return
                fields=list(reader.fieldnames or [])
                self._snapshot=(stamp,version,fields,rows)
                yield version,SimpleNamespace(fieldnames=fields,rows=rows)
            else:
                self._snapshot=None
                yield version,SimpleNamespace(fieldnames=reader.fieldnames,rows=reader)

    def _stage(self,token,project_id,run_id,stage_id):
        project=self.store.project(token,project_id)
        with self.store.connect() as db:
            revision=db.execute('''SELECT r.project_id,r.updated,r.status,coalesce(v.revision,0) AS result_revision
                FROM runs r LEFT JOIN result_read_versions v ON v.run_id=r.id WHERE r.id=?''',(run_id,)).fetchone()
        if revision is None or revision['project_id']!=project_id:
            raise AccessDenied('Resultado de outro projeto.')
        key=(project_id,run_id,stage_id,revision['updated'],revision['status'],revision['result_revision'],project.get('revision'))
        if self._stage_snapshot and self._stage_snapshot[0]==key:
            return self._stage_snapshot[1:]
        run=getattr(self.store,'get_result_run',self.store.get_run)(token,run_id)
        if run['project_id']!=project_id:
            raise AccessDenied('Resultado de outro projeto.')
        stage=next((s for s in run['stages'] if s['id']==stage_id and s['status'] in self.available_statuses),None)
        if stage is None:
            raise ValueError('A tabela estará disponível após concluir a recuperação.')
        self._stage_snapshot=(key,run,stage)
        return run,stage

    def tables(self,token,project_id,run_id,stage_id):
        run,stage=self._stage(token,project_id,run_id,stage_id)
        result=[]
        for filename in stage.get('artifacts',[]):
            path=self.store.scoped_path(project_id,filename)
            name=path.name
            is_summary=path.suffix.lower()=='.csv' and path.stem.upper().endswith(('_FULL','_MOLS','_SIMS'))
            is_integrated=name=='compounds.csv' and path.parent.parent.name=='compounds'
            from .artifact_naming import PATTERN
            match=PATTERN.match(path.stem)
            is_summary=is_summary or bool(match and match[3] in ('mol','sim','merged','compounds'))
            is_provided=bool(stage.get('configuration',{}).get('provided_results')) and path.suffix.lower()=='.csv'
            if (is_summary or is_integrated or is_provided) and path.is_file():
                result.append({'path':str(path),'name':name,'integrated':is_integrated})
        result.sort(key=lambda t:(not t['integrated'],t['name']))
        from .artifacts import annotate_files
        result=annotate_files(stage,result,run_id)
        from .result_origins import ResultOrigins
        origins=ResultOrigins(run['stages'],run_id)
        for item in result:
            info=origins.describe(item['path'],stage_id,run_id)
            item.update(public_name=info['public_name'],display_name=info['public_name'],label=info['public_name'],naming=info['naming'])
        return result

    def _path(self,token,project_id,run_id,stage_id,filename):
        path=self.store.scoped_path(project_id,filename)
        _,stage=self._stage(token,project_id,run_id,stage_id)
        from .result_index import file_stamp
        try:stamp=file_stamp(path)
        except OSError:raise AccessDenied('Tabela não autorizada para esta etapa.') from None
        scope=[]
        for value in stage.get('artifacts',[]):
            if not str(value).lower().endswith('.csv'):continue
            candidate=self.store.scoped_path(project_id,value)
            try:current=stamp if candidate==path else file_stamp(candidate)
            except OSError:current=None
            scope.append((str(candidate),current))
        key=(self._stage_snapshot[0],tuple(scope))
        cached=self._table_access
        if cached and cached[0]==key and cached[2].get(str(path))==stamp:return path
        choices=self.tables(token,project_id,run_id,stage_id)
        allowed={table['path'] for table in choices}
        if str(path) not in allowed:
            raise AccessDenied('Tabela não autorizada para esta etapa.')
        stamps=dict(cached[2]) if cached and cached[0]==key else {}
        stamps[str(path)]=stamp
        self._table_access=(key,allowed,stamps)
        return path

    @staticmethod
    def _digest(stream):
        digest=hashlib.sha256()
        for chunk in iter(lambda:stream.read(1024*1024),b''):
            digest.update(chunk)
        stream.seek(0)
        return digest.hexdigest()

    def page(self,token,project_id,run_id,stage_id,filename,offset=0,limit=25,query='',sort_by=None,ascending=True):
        if type(offset) is not int or offset<0 or type(limit) is not int or not 1<=limit<=100:
            raise ValueError('Página de tabela inválida.')
        if not isinstance(query,str) or len(query)>200:
            raise ValueError('Use até 200 caracteres na busca.')
        if sort_by is not None and (not isinstance(sort_by,str) or type(ascending) is not bool):
            raise ValueError('Ordenação de tabela inválida.')
        path=self._path(token,project_id,run_id,stage_id,filename)
        from .result_index import MIN_INDEX_BYTES,ResultIndex
        if path.stat().st_size>=MIN_INDEX_BYTES:
            indexed=ResultIndex(self.store.project_dir(project_id),path,self._digest).page(
                offset,limit,query,sort_by,ascending,self.requires_compounds)
            rows=[]
            for index,row in indexed.pop('visible'):
                identifier=row.get('molecule_chembl_id',row.get(indexed['columns'][0],'')) or ''
                rows.append(dict(index=index,id=identifier,smiles=row.get('canonical_smiles') or '',links=compound_links(row),properties=row))
            project=self.store.project(token,project_id)
            with self.store.connect() as db:
                active=db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')",(project_id,)).fetchone()
            return dict(indexed,rows=rows,offset=offset,limit=limit,name=path.name,can_edit=project['role']!='viewer' and active is None)
        rows=[];matched=total=0;selected=[]
        with self._reader(path) as (version,reader):
            if self.requires_compounds and not {'molecule_chembl_id','canonical_smiles'}.issubset(reader.fieldnames or []):
                raise ValueError('Esta tabela não contém código do composto e SMILES.')
            fields=list(reader.fieldnames or [])
            if 'vina' in fields and 'dock6' in fields and 'min-max' in fields and 'normalized_score' not in fields:
                fields=fields+['normalized_score']
            if sort_by is not None and sort_by not in fields:
                raise ValueError('Coluna de ordenação inválida.')
            # Keep one bounded view of the cached CSV. Authorization and the file
            # stamp are checked on every request, including cache hits.
            key=(version,query.lower(),sort_by,ascending,self.requires_compounds)
            cacheable=isinstance(reader.rows,list)
            view=self._view
            if cacheable and view is not None and view[0]==key:
                total,selected=view[1:]
                matched=len(selected)
            else:
                needle=query.lower()
                for index,original in enumerate(reader.rows):
                    row=original
                    if 'normalized_score' in fields and 'normalized_score' not in row:
                        row=dict(row,normalized_score=row.get('min-max',''))
                    total+=1
                    if needle:
                        identifier=row.get('molecule_chembl_id',row.get(fields[0],'')) or ''
                        smiles=row.get('canonical_smiles') or ''
                        searchable=identifier+' '+smiles if self.requires_compounds else ' '.join(str(v or '') for v in row.values())
                        if needle not in searchable.lower():continue
                    if cacheable or sort_by is not None or offset<=matched<offset+limit:
                        selected.append((index,row))
                    matched+=1
                if sort_by is not None:
                    numeric=sort_by in ('score','vina','dock6','normalized_score','z-score','min-max') or (
                        sort_by not in ('molecule_chembl_id','canonical_smiles','receptor_id') and is_numeric(row.get(sort_by) for _,row in selected))
                    selected=ordered(selected,lambda item:item[1].get(sort_by),ascending,numeric)
                self._view=(key,total,selected) if cacheable else None
            visible=selected[offset:offset+limit] if cacheable or sort_by is not None else selected
            for index,original in visible:
                row=dict(original)
                identifier=row.get('molecule_chembl_id',row.get(fields[0],'')) or ''
                smiles=row.get('canonical_smiles') or ''
                rows.append({'index':index,'id':identifier,'smiles':smiles,'links':compound_links(row),'properties':row})
        project=self.store.project(token,project_id)
        with self.store.connect() as db:
            active=db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')",(project_id,)).fetchone()
        return {'rows':rows,'total':total,'matched':matched,'version':version,'offset':offset,'limit':limit,'name':path.name,
                'can_edit':project['role']!='viewer' and active is None,'columns':fields}

    def _record(self,path,index,version,*,root):
        if type(index) is not int or index<0:raise ValueError('Seleção de linha inválida.')
        from .result_index import MIN_INDEX_BYTES,ResultIndex
        if path.stat().st_size>=MIN_INDEX_BYTES:
            return ResultIndex(root,path,self._digest).row(index,version)
        with self._reader(path) as (current,reader):
            if current!=version:raise ValueError('Atualize a tabela antes de visualizar a conformação.')
            if isinstance(reader.rows,list):
                if index>=len(reader.rows):raise IndexError('A linha não está mais nesta tabela.')
                return dict(reader.rows[index])
            from itertools import islice
            row=next(islice(reader.rows,index,index+1),None)
            if row is None:raise IndexError('A linha não está mais nesta tabela.')
            return dict(row)

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
