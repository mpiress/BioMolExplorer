"""Compact, revisioned read models; authoritative execution JSON is untouched."""
import codecs
import json
import sqlite3
import threading
from weakref import WeakValueDictionary
import zlib

VERSION=4
_GATES=WeakValueDictionary()
_GATES_LOCK=threading.Lock()
STABLE_FIELDS={'artifact_records','artifact_manifest','artifacts','input_manifest','input_files','artifact_context',
               '_review_input_groups','_review_naming_groups','execution_manifest','processing_runtime','batches'}


class _Gate:
    def __init__(self):self.lock=threading.RLock()
RECORD_FIELDS=('artifact_id','block_id','output','original_name','record_count','naming',
               'occurrence_id','run_id','path','context','output_port','artifact_type',
               '_review_input_group','_review_naming_group')


def create_schema(db):
    db.executescript('''
        CREATE TABLE IF NOT EXISTS result_read_versions(run_id TEXT PRIMARY KEY,revision INTEGER NOT NULL);
        CREATE TABLE IF NOT EXISTS result_read_catalog(run_id TEXT PRIMARY KEY,revision INTEGER NOT NULL,
            schema_version INTEGER NOT NULL,stages TEXT NOT NULL,overlay TEXT NOT NULL DEFAULT '{}');
        CREATE TABLE IF NOT EXISTS result_read_stages(run_id TEXT NOT NULL,stage_id TEXT NOT NULL,
            revision INTEGER NOT NULL,schema_version INTEGER NOT NULL,position INTEGER NOT NULL,
            stage BLOB NOT NULL,proof BLOB NOT NULL,PRIMARY KEY(run_id,stage_id));
        CREATE TRIGGER IF NOT EXISTS result_stage_delete AFTER DELETE ON runs BEGIN
            DELETE FROM result_read_stages WHERE run_id=OLD.id;
        END;
        CREATE TRIGGER IF NOT EXISTS result_read_insert AFTER INSERT ON runs BEGIN
            INSERT OR REPLACE INTO result_read_versions VALUES(NEW.id,1);
            DELETE FROM result_read_catalog WHERE run_id=NEW.id;
        END;
        DROP TRIGGER IF EXISTS result_read_update;
        CREATE TRIGGER result_read_update AFTER UPDATE OF stages ON runs BEGIN
            INSERT INTO result_read_versions VALUES(NEW.id,1)
                ON CONFLICT(run_id) DO UPDATE SET revision=revision+1;
        END;
        CREATE TRIGGER IF NOT EXISTS result_read_delete AFTER DELETE ON runs BEGIN
            DELETE FROM result_read_versions WHERE run_id=OLD.id;
            DELETE FROM result_read_catalog WHERE run_id=OLD.id;
        END;
    ''')
    if 'overlay' not in {row[1] for row in db.execute('PRAGMA table_info(result_read_catalog)')}:
        db.execute("ALTER TABLE result_read_catalog ADD COLUMN overlay TEXT NOT NULL DEFAULT '{}'")


class RecordPacker:
    def __init__(self):
        self.groups={};self.previous=None;self.previous_key=None
        self.naming_groups={};self.previous_naming=None;self.naming_key=None

    def add(self,record):
        item={key:record[key] for key in RECORD_FIELDS if key in record}
        naming=record.get('naming') or {}
        parents=naming.get('source_artifact_ids')
        if parents:
            if parents!=self.previous_naming:
                self.naming_key=next((key for key,value in self.naming_groups.items() if value==parents),None)
                if self.naming_key is None:
                    self.naming_key=str(len(self.naming_groups));self.naming_groups[self.naming_key]=parents
                self.previous_naming=parents
            item['naming']={key:value for key,value in naming.items() if key!='source_artifact_ids'}
            item['_review_naming_group']=self.naming_key
        inputs=record.get('inputs') or {}
        if inputs:
            if inputs!=self.previous:
                self.previous_key=next((key for key,value in self.groups.items() if value==inputs),None)
                if self.previous_key is None:
                    self.previous_key=str(len(self.groups));self.groups[self.previous_key]=inputs
                self.previous=inputs
            item['_review_input_group']=self.previous_key
        return item


def compact_stage(original):
    stage={key:value for key,value in original.items() if key not in
        ('artifact_records','batches','execution_manifest','processing_runtime')}
    packer=RecordPacker()
    records=[packer.add(record) for record in original.get('artifact_records',[])]
    if 'artifact_records' in original:stage['artifact_records']=records
    if packer.groups:stage['_review_input_groups']=packer.groups
    if packer.naming_groups:stage['_review_naming_groups']=packer.naming_groups
    if 'batches' in original:stage['batches']=[compact_stage(batch) for batch in original['batches']]
    return stage


class ProjectionReader:
    """Decode legacy run TEXT incrementally, one descriptor at a time."""
    def __init__(self,blob):
        self.blob=blob;self.buffer='';self.position=0;self.eof=False
        self.decoder=codecs.getincrementaldecoder('utf-8')()
        self.json=json.JSONDecoder()

    def fill(self):
        if self.position:
            self.buffer=self.buffer[self.position:];self.position=0
        chunk=self.blob.read(65536)
        self.eof=not chunk
        self.buffer+=self.decoder.decode(chunk,final=self.eof)

    def peek(self):
        while True:
            while self.position<len(self.buffer) and self.buffer[self.position].isspace():self.position+=1
            if self.position<len(self.buffer):return self.buffer[self.position]
            if self.eof:raise ValueError('Metadados de execução truncados.')
            self.fill()

    def consume(self,character):
        if self.peek()!=character:raise ValueError('Metadados de execução inválidos.')
        self.position+=1

    def value(self):
        self.peek()
        while True:
            try:
                value,end=self.json.raw_decode(self.buffer,self.position)
                if end==len(self.buffer) and not self.eof:
                    self.fill();continue
                self.position=end
                return value
            except json.JSONDecodeError:
                if self.eof:raise
                self.fill()

    def array(self,reader):
        self.consume('[');items=[]
        if self.peek()==']':self.consume(']');return items
        while True:
            items.append(reader())
            if self.peek()==']':self.consume(']');return items
            self.consume(',')

    def stage(self):
        self.consume('{');stage={};packer=RecordPacker()
        if self.peek()=='}':self.consume('}');return stage
        while True:
            key=self.value();self.consume(':')
            if key=='artifact_records':
                stage[key]=self.array(lambda:packer.add(self.value()))
            elif key=='batches':stage[key]=self.array(self.stage)
            elif key in ('execution_manifest','processing_runtime'):self.value()
            else:stage[key]=self.value()
            if self.peek()=='}':
                self.consume('}')
                if packer.groups:stage['_review_input_groups']=packer.groups
                if packer.naming_groups:stage['_review_naming_groups']=packer.naming_groups
                return stage
            self.consume(',')


def store_snapshot(db,run_id,stages,*,unchanged=False):
    revision=db.execute('SELECT revision FROM result_read_versions WHERE run_id=?',(run_id,)).fetchone()
    current=revision[0] if revision else 0
    cached=db.execute('SELECT revision,schema_version FROM result_read_catalog WHERE run_id=?',(run_id,)).fetchone()
    if unchanged and cached and cached[0]==current-1 and cached[1]==VERSION:
        overlay={stage['id']:_status_overlay(stage) for stage in stages}
        db.execute('UPDATE result_read_catalog SET revision=?,overlay=? WHERE run_id=?',(current,json.dumps(overlay),run_id))
        db.execute('UPDATE result_read_stages SET revision=? WHERE run_id=?',(current,run_id))
    else:_save(db,run_id,current,[compact_stage(stage) for stage in stages])


def _status_overlay(stage):
    value={key:item for key,item in stage.items() if key not in STABLE_FIELDS}
    if 'batches' in stage:value['batches']=[_status_overlay(batch) for batch in stage['batches']]
    return value


def _apply_overlay(stage,overlay):
    if not overlay:return stage
    for key in list(stage):
        if key not in STABLE_FIELDS and key not in overlay:stage.pop(key)
    for key,value in overlay.items():
        if key!='batches':stage[key]=value
    for batch,update in zip(stage.get('batches',[]),overlay.get('batches',[])):_apply_overlay(batch,update)
    return stage


def _encode(value):
    return zlib.compress(json.dumps(value,separators=(',',':')).encode(),1)


def proof_stage(stage):
    result={key:value for key,value in stage.items() if key not in
        ('artifact_records','_review_input_groups','_review_naming_groups','batches')}
    if 'batches' in stage:result['batches']=[proof_stage(batch) for batch in stage['batches']]
    return result


def _save(db,run_id,revision,stages):
    db.execute('''INSERT OR REPLACE INTO result_read_catalog(run_id,revision,schema_version,stages,overlay)
        VALUES (?,?,?,?,?)''',(run_id,revision,VERSION,_encode(stages),'{}'))
    db.execute('DELETE FROM result_read_stages WHERE run_id=?',(run_id,))
    db.executemany('INSERT INTO result_read_stages VALUES (?,?,?,?,?,?,?)',[
        (run_id,stage['id'],revision,VERSION,index,_encode(stage),_encode(proof_stage(stage))) for index,stage in enumerate(stages)])


def read_proofs(store,run_id):
    with store.connect() as db:
        rows=db.execute('''SELECT c.proof,c.stage_id,m.overlay FROM result_read_stages c JOIN runs r ON r.id=c.run_id
            JOIN result_read_catalog m ON m.run_id=c.run_id
            LEFT JOIN result_read_versions v ON v.run_id=r.id WHERE c.run_id=? AND c.revision=coalesce(v.revision,0)
            AND c.schema_version=? AND m.revision=c.revision AND m.schema_version=c.schema_version ORDER BY c.position''',(run_id,VERSION)).fetchall()
    if rows:
        try:
            overlay=json.loads(rows[0][2])
            return [_apply_overlay(json.loads(zlib.decompress(row[0])),overlay.get(row[1],{})) for row in rows]
        except (zlib.error,json.JSONDecodeError,TypeError):rows=[]
    run=read_run(store,run_id)
    return [proof_stage(stage) for stage in run['stages']] if run else None


def read_stage(store,run_id,stage_id):
    with store.connect() as db:
        row=db.execute('''SELECT c.stage,m.overlay FROM result_read_stages c JOIN runs r ON r.id=c.run_id
            JOIN result_read_catalog m ON m.run_id=c.run_id
            LEFT JOIN result_read_versions v ON v.run_id=r.id WHERE c.run_id=? AND c.stage_id=?
            AND c.revision=coalesce(v.revision,0) AND c.schema_version=?
            AND m.revision=c.revision AND m.schema_version=c.schema_version''',(run_id,stage_id,VERSION)).fetchone()
    if row:
        try:return _apply_overlay(json.loads(zlib.decompress(row[0])),json.loads(row[1]).get(stage_id,{}))
        except (zlib.error,json.JSONDecodeError,TypeError):row=None
    run=read_run(store,run_id)
    return next((stage for stage in run['stages'] if stage['id']==stage_id),None) if run else None


def read_run(store,run_id):
    with store.connect() as db:
        warm=db.execute('''SELECT 1 FROM result_read_catalog c LEFT JOIN result_read_versions v ON v.run_id=c.run_id
            WHERE c.run_id=? AND c.revision=coalesce(v.revision,0) AND c.schema_version=?''',(run_id,VERSION)).fetchone()
    if warm:return _read_run(store,run_id)
    key=(str(store.database),run_id)
    with _GATES_LOCK:
        gate=_GATES.get(key)
        if gate is None:gate=_Gate();_GATES[key]=gate
    with gate.lock:return _read_run(store,run_id)


def _read_run(store,run_id):
    with store.connect() as db:
        db.execute('BEGIN')
        row=db.execute('''SELECT r.rowid AS sqlite_row_id,r.id,r.project_id,r.user_id,r.status,r.created,r.updated,r.error,
            coalesce(v.revision,0) AS source_revision,c.revision AS catalog_revision,c.schema_version,c.stages AS review_stages,c.overlay
            FROM runs r LEFT JOIN result_read_versions v ON v.run_id=r.id
            LEFT JOIN result_read_catalog c ON c.run_id=r.id WHERE r.id=?''',(run_id,)).fetchone()
        if row is None:return None
        result={key:row[key] for key in ('id','project_id','user_id','status','created','updated','error')}
        cached=row['schema_version']==VERSION and row['catalog_revision']==row['source_revision']
        if cached:
            try:
                overlay=json.loads(row['overlay'])
                stages=[_apply_overlay(stage,overlay.get(stage['id'],{})) for stage in json.loads(zlib.decompress(row['review_stages']))]
            except (zlib.error,json.JSONDecodeError,TypeError):cached=False
        if not cached:
            with db.blobopen('runs','stages',row['sqlite_row_id'],readonly=True) as blob:
                stages=_project_blob(blob)
        revision=row['source_revision']
    if not cached:
        try:
            with store.connect() as db:
                # A read helper may be called inside an audited deletion that
                # already owns the write transaction on another connection.
                db.execute('PRAGMA busy_timeout=0')
                db.execute('BEGIN IMMEDIATE')
                current=db.execute('SELECT revision FROM result_read_versions WHERE run_id=?',(run_id,)).fetchone()
                if (current[0] if current else 0)==revision and db.execute('SELECT 1 FROM runs WHERE id=?',(run_id,)).fetchone():
                    _save(db,run_id,revision,stages)
        except sqlite3.OperationalError as error:
            if error.sqlite_errorcode not in (sqlite3.SQLITE_BUSY,sqlite3.SQLITE_LOCKED):raise
    result['stages']=stages
    return result


def _project_blob(blob):
    reader=ProjectionReader(blob)
    return reader.array(reader.stage)
