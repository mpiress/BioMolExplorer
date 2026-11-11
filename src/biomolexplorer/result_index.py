"""Disposable, versioned SQLite indexes for large scientific CSV tables.

Original CSVs remain authoritative. Row ordinals and SHA-256 versions retain
the existing curation/preview contracts; no chemical engine runs on indexing.
"""
from contextlib import contextmanager
import csv
import hashlib
import io
import json
from pathlib import Path
import sqlite3
import time

from .table_sort import value_key

MIN_INDEX_BYTES=512*1024
SCHEMA_VERSION=1


def file_stamp(path):
    stat=Path(path).stat()
    return (stat.st_dev,stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns)


class ResultIndex:
    def __init__(self,root,path,digest):
        self.source=Path(path)
        key=hashlib.sha256(str(self.source.resolve()).encode()).hexdigest()
        self.root=Path(root)/'.biomolexplorer/cache/result-index'
        self.root.mkdir(parents=True,exist_ok=True)
        self.database=self.root/(key+'.sqlite3')
        self.digest=digest

    @contextmanager
    def _open(self):
        import fcntl
        with self.database.with_suffix('.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            for attempt in range(2):
                db=None
                try:
                    db=sqlite3.connect(self.database,timeout=30)
                    db.execute('PRAGMA temp_store=FILE')
                    db.execute('PRAGMA cache_size=-4096')
                    db.executescript('''
                        CREATE TABLE IF NOT EXISTS metadata(id INTEGER PRIMARY KEY,value TEXT NOT NULL);
                        CREATE TABLE IF NOT EXISTS rows(row_idx INTEGER PRIMARY KEY,payload TEXT NOT NULL,
                            compound_search TEXT NOT NULL,all_search TEXT NOT NULL);
                        CREATE TABLE IF NOT EXISTS sort_values(column_idx INTEGER NOT NULL,row_idx INTEGER NOT NULL,
                            text_value TEXT,number_value REAL,PRIMARY KEY(column_idx,row_idx));
                        CREATE TABLE IF NOT EXISTS sort_columns(column_idx INTEGER PRIMARY KEY,numeric INTEGER NOT NULL);
                        CREATE TABLE IF NOT EXISTS views(key TEXT PRIMARY KEY,matched INTEGER NOT NULL,used REAL NOT NULL);
                        CREATE TABLE IF NOT EXISTS view_rows(key TEXT NOT NULL,position INTEGER NOT NULL,row_idx INTEGER NOT NULL,
                            PRIMARY KEY(key,position));
                    ''')
                    self._ensure(db)
                    break
                except (sqlite3.DatabaseError,json.JSONDecodeError):
                    if db:db.close()
                    if attempt:raise
                    for suffix in ('','-wal','-shm','-journal'):
                        Path(str(self.database)+suffix).unlink(missing_ok=True)
                except BaseException:
                    if db:db.close()
                    raise
            try:
                with db:yield db
            finally:db.close()

    def _ensure(self,db):
        stamp=file_stamp(self.source)
        found=db.execute('SELECT value FROM metadata WHERE id=1').fetchone()
        if found:
            metadata=json.loads(found[0])
            if metadata.get('schema')==SCHEMA_VERSION and metadata.get('stamp')==list(stamp):
                self.metadata=metadata
                return
        with db:
            for table in ('metadata','rows','sort_values','sort_columns','views','view_rows'):
                db.execute('DELETE FROM '+table)
            with self.source.open('rb') as stream:
                version=self.digest(stream)
                reader=csv.DictReader(io.TextIOWrapper(stream,encoding='utf-8-sig',newline=''))
                fields=list(reader.fieldnames or [])
                alias={'vina','dock6','min-max'}<=set(fields) and 'normalized_score' not in fields
                if alias:fields.append('normalized_score')
                batch=[];total=0
                for index,original in enumerate(reader):
                    row=dict(original)
                    if alias:row['normalized_score']=row.get('min-max','')
                    identifier=row.get('molecule_chembl_id',row.get(fields[0],'')) or ''
                    smiles=row.get('canonical_smiles') or ''
                    batch.append((index,json.dumps(row,ensure_ascii=False),
                        (identifier+' '+smiles).lower(),' '.join(str(value or '') for value in row.values()).lower()))
                    total+=1
                    if len(batch)==512:
                        db.executemany('INSERT INTO rows VALUES (?,?,?,?)',batch);batch.clear()
                if batch:db.executemany('INSERT INTO rows VALUES (?,?,?,?)',batch)
            if file_stamp(self.source)!=stamp:
                raise ValueError('A tabela foi alterada durante a indexação. Atualize os resultados.')
            self.metadata=dict(schema=SCHEMA_VERSION,stamp=list(stamp),version=version,fields=fields,total=total)
            db.execute('INSERT INTO metadata VALUES (1,?)',(json.dumps(self.metadata),))

    def _sort_column(self,db,column):
        index=self.metadata['fields'].index(column)
        previous=db.execute('SELECT numeric FROM sort_columns WHERE column_idx=?',(index,)).fetchone()
        if previous:return index,bool(previous[0])
        nonempty=numbers=0;batch=[]
        for ordinal,payload in db.execute('SELECT row_idx,payload FROM rows ORDER BY row_idx'):
            value=json.loads(payload).get(column);text=value_key(value);number=value_key(value,True)
            nonempty+=text is not None;numbers+=number is not None
            batch.append((index,ordinal,text,number))
            if len(batch)==512:
                db.executemany('INSERT INTO sort_values VALUES (?,?,?,?)',batch);batch.clear()
        if batch:db.executemany('INSERT INTO sort_values VALUES (?,?,?,?)',batch)
        numeric=column in ('score','vina','dock6','normalized_score','z-score','min-max') or (
            column not in ('molecule_chembl_id','canonical_smiles','receptor_id') and nonempty>0 and numbers==nonempty)
        db.execute('INSERT INTO sort_columns VALUES (?,?)',(index,int(numeric)))
        return index,numeric

    def page(self,offset,limit,query,sort_by,ascending,requires_compounds):
        with self._open() as db:
            fields=self.metadata['fields']
            if requires_compounds and not {'molecule_chembl_id','canonical_smiles'}<=set(fields):
                raise ValueError('Esta tabela não contém código do composto e SMILES.')
            if sort_by is not None and sort_by not in fields:raise ValueError('Coluna de ordenação inválida.')
            if not query and sort_by is None:
                visible=db.execute('SELECT row_idx,payload FROM rows ORDER BY row_idx LIMIT ? OFFSET ?',(limit,offset)).fetchall()
                matched=self.metadata['total']
            else:
                key=json.dumps([query.lower(),sort_by,ascending,requires_compounds],ensure_ascii=False)
                cached=db.execute('SELECT matched,used FROM views WHERE key=?',(key,)).fetchone()
                if cached is None:
                    params=[];join='';order='r.row_idx ASC'
                    if sort_by is not None:
                        column,numeric=self._sort_column(db,sort_by)
                        join=' JOIN sort_values s ON s.row_idx=r.row_idx AND s.column_idx=?'
                        params.append(column)
                        value='s.number_value' if numeric else 's.text_value'
                        order=f'({value} IS NULL) ASC,{value} '+('ASC' if ascending else 'DESC')+',r.row_idx ASC'
                    search='r.compound_search' if requires_compounds else 'r.all_search'
                    where=' WHERE instr('+search+',?)>0' if query else ''
                    if query:params.append(query.lower())
                    db.execute('INSERT INTO view_rows SELECT ?,ROW_NUMBER() OVER (ORDER BY '+order+')-1,r.row_idx FROM rows r'+join+where,
                        [key,*params])
                    matched=db.execute('SELECT count(*) FROM view_rows WHERE key=?',(key,)).fetchone()[0]
                    db.execute('INSERT INTO views VALUES (?,?,?)',(key,matched,time.time()))
                    obsolete=db.execute('SELECT key FROM views WHERE key<>? ORDER BY used DESC LIMIT -1 OFFSET 7',(key,)).fetchall()
                    for old, in obsolete:
                        db.execute('DELETE FROM view_rows WHERE key=?',(old,));db.execute('DELETE FROM views WHERE key=?',(old,))
                else:matched=cached[0]
                if cached is not None and cached[1]<time.time()-60:
                    db.execute('UPDATE views SET used=? WHERE key=?',(time.time(),key))
                visible=db.execute('''SELECT r.row_idx,r.payload FROM view_rows v JOIN rows r ON r.row_idx=v.row_idx
                    WHERE v.key=? AND v.position>=? AND v.position<? ORDER BY v.position''',(key,offset,offset+limit)).fetchall()
            result=dict(version=self.metadata['version'],total=self.metadata['total'],matched=matched,
                columns=list(fields),visible=[(ordinal,json.loads(payload)) for ordinal,payload in visible])
            if file_stamp(self.source)!=tuple(self.metadata['stamp']):raise ValueError('A tabela foi alterada. Atualize os resultados.')
            return result

    def row(self,index,version):
        with self._open() as db:
            if self.metadata['version']!=version:raise ValueError('Atualize a tabela antes de visualizar a conformação.')
            found=db.execute('SELECT payload FROM rows WHERE row_idx=?',(index,)).fetchone()
            if found is None:raise IndexError('A linha não está mais nesta tabela.')
            if file_stamp(self.source)!=tuple(self.metadata['stamp']):raise ValueError('A tabela foi alterada. Atualize os resultados.')
            return json.loads(found[0])
