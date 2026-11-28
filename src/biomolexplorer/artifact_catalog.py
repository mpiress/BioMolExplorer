"""Content objects, operation identity and dependency edges, as catalog rows.

``artifact_instances`` answers which files a run published.  It cannot answer
which operation produced an occurrence, what that operation consumed, or which
logical artifacts share one physical content, because the dependencies were
stored as a JSON list inside each row.  The four identities of #6 are separated
here:

* **content** — the SHA-256 of the bytes, in ``content_objects``;
* **logical artifact** — ``artifact_id`` and ``occurrence_id``, unchanged;
* **scientific entity** — ``artifact_entities``, unchanged;
* **operation** — one execution of one block, in ``operations``, with its
  inputs and outputs as edge rows.

The derivation reads the artifact records of a run, which carry the producing
run, the software identity and the exact input lineage.  Everything it cannot
establish is recorded in ``catalog_backfill_issues`` rather than guessed, and
every query here answers from indexes alone: none of them reads a directory.
"""
import hashlib
import json
import time
from pathlib import Path
from .catalog_schema import CATALOG_VERSION

SIGNATURE_VERSION=1

MEDIA_TYPES={'.pdb':'chemical/x-pdb','.ent':'chemical/x-pdb','.pdbqt':'chemical/x-pdbqt',
    '.cif':'chemical/x-cif','.mmcif':'chemical/x-cif','.mol2':'chemical/x-mol2',
    '.sdf':'chemical/x-mdl-sdfile','.mol':'chemical/x-mdl-molfile','.smi':'chemical/x-daylight-smiles',
    '.dms':'chemical/x-dms','.csv':'text/csv','.tsv':'text/tab-separated-values',
    '.json':'application/json','.txt':'text/plain','.log':'text/plain','.md':'text/markdown',
    '.html':'text/html','.png':'image/png','.svg':'image/svg+xml','.zip':'application/zip','.gz':'application/gzip'}

STATES={'succeeded':'completed','failed':'failed','interrupted':'interrupted',
        'cancelled':'cancelled','skipped':'skipped','awaiting_input':'awaiting_input'}


def media_type(path):
    """Byte format of the stored file, independent of its scientific role."""
    return MEDIA_TYPES.get(Path(path).suffix.lower(),'application/octet-stream')


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,default=str).encode()).hexdigest()


def operation_id(project_id,stage_id,run_id):
    """Deterministic identity of one execution of one block."""
    return digest(['operation',project_id,stage_id,run_id])


def _parameters_digest(stage):
    config=stage.get('configuration',stage)
    return digest(['parameters',config.get('parameters',{})])


def _software_digest(group):
    """Identity of the software that produced a group of records.

    ``None`` when the run recorded none, which is the case for every run made
    before the implementation digest existed.  A digest of absence would be a
    signature that cannot distinguish two implementations.
    """
    values=sorted({json.dumps([record.get('processing_version'),record.get('producer_software'),
        record.get('runtime_id')],sort_keys=True,default=str) for record in group})
    if all(json.loads(value)==[None,None,None] for value in values):return None
    return digest(['software',SIGNATURE_VERSION,values])


def _input_digests(group):
    """Content digests of the inputs, as the run recorded them."""
    found=set()
    for record in group:
        for entries in (record.get('dependency_manifest') or {}).values():
            for entry in entries or []:
                if isinstance(entry,(list,tuple)) and len(entry)==2:found.add(entry[1])
    return sorted(found)


def signature(operation_type,parameters_digest,software_digest,input_digests):
    """Provisional operation signature; ADR-004 fixes the definitive model."""
    return digest(['signature',SIGNATURE_VERSION,operation_type,parameters_digest,software_digest,input_digests])


def _edges(group):
    """Input edges of one operation, in a deterministic order."""
    ports={}
    for record in group:
        for port,rows in (record.get('inputs') or {}).items():
            bucket=ports.setdefault(port,{})
            for row in rows or []:
                occurrence=row.get('occurrence_id');asset=row.get('asset_id')
                role='artifact' if occurrence else 'asset' if asset else 'external'
                bucket[(role,str(occurrence or asset or row.get('path','')))]=(occurrence,asset,role)
    return {port:[bucket[key] for key in sorted(bucket)] for port,bucket in ports.items()}


def _issue(db,project_id,run_id,stage_id,reason,details):
    db.execute('INSERT OR REPLACE INTO catalog_backfill_issues VALUES (?,?,?,?,?,?)',
        (project_id,run_id,stage_id,reason,json.dumps(details,sort_keys=True,default=str),time.time()))


# An operation survives while an occurrence indexed by another run still
# declares it as its producer.  Clearing the results of one run must not erase
# the lineage of a run that reused them, and must not depend on the order in
# which the runs of a project are re-indexed.
# ``CROSS JOIN`` pins the join order: starting from the operations of this run
# keeps the work proportional to what the run published, instead of to every
# artifact the project holds.
_CLAIMED='''SELECT e.operation_id FROM operations o
    CROSS JOIN operation_outputs e ON e.operation_id=o.operation_id
    CROSS JOIN artifact_instances i ON i.project_id=? AND i.occurrence_id=e.occurrence_id AND i.indexed_run_id<>?
    WHERE o.project_id=? AND o.run_id=?'''


def forget_run(db,project_id,run_id):
    """Drop the operations a run produced, before its index is rebuilt."""
    keys=(project_id,run_id,project_id,run_id)
    for table in ('operation_inputs','operation_outputs'):
        db.execute(f'''DELETE FROM {table} WHERE operation_id IN
            (SELECT operation_id FROM operations WHERE project_id=? AND run_id=?
             EXCEPT {_CLAIMED})''',(project_id,run_id,*keys))
    db.execute(f'''DELETE FROM operations WHERE project_id=? AND run_id=?
        AND operation_id NOT IN ({_CLAIMED})''',(project_id,run_id,*keys))
    db.execute('DELETE FROM catalog_backfill_issues WHERE project_id=? AND run_id=?',(project_id,run_id))


def index_operations(db,project_id,run_id,stages,sizes=None):
    """Record content objects, operations and edges for one indexed run."""
    sizes=sizes or {}
    now=time.time()
    for stage in stages:
        stage_id=stage.get('id')
        if not stage_id:continue
        records=stage.get('artifact_records') or []
        if not records:
            if stage.get('artifacts'):
                # A published file with no descriptor has no occurrence, so it
                # has no edge either.  The gap is recorded under its stage.
                _issue(db,project_id,run_id,stage_id,'missing_artifact_records',
                    {'artifacts':len(stage.get('artifacts') or [])})
            continue
        groups={}
        for record in records:
            if not record.get('occurrence_id'):
                _issue(db,project_id,run_id,stage_id,'missing_occurrence',{'path':record.get('path')})
                continue
            groups.setdefault(record.get('run_id') or run_id,[]).append(record)
        operation_type=stage.get('operation') or stage.get('configuration',{}).get('operation') or 'unknown'
        parameters=_parameters_digest(stage)
        for producer,group in groups.items():
            identity=operation_id(project_id,stage_id,producer)
            own=producer==run_id
            software=_software_digest(group)
            if software is None:
                _issue(db,project_id,run_id,stage_id,'missing_software_identity',
                    {'producer_run':producer,'records':len(group)})
            mark=signature(operation_type,parameters,software,_input_digests(group)) if software else None
            state=STATES.get(stage.get('status'),'unknown') if own else 'completed'
            row=(identity,project_id,operation_type,stage_id,
                 next((record.get('block_id') for record in group if record.get('block_id')),None),
                 producer,mark,SIGNATURE_VERSION,parameters,software,state,
                 stage.get('started_at') if own else None,stage.get('finished_at') if own else None,now)
            # The run that produced an operation is the authority on it.  A run
            # that only reuses its result records what it knows without
            # overwriting the producer's own row.
            db.execute(('INSERT OR REPLACE' if own else 'INSERT OR IGNORE')+
                ' INTO operations VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)',row)
            keyword='OR REPLACE' if own else 'OR IGNORE'
            for record in group:
                db.execute(f'INSERT {keyword} INTO operation_outputs VALUES (?,?,?)',
                    (identity,record['occurrence_id'],record.get('output_port')))
                content=record.get('content_sha256')
                if content:
                    db.execute('INSERT OR IGNORE INTO content_objects VALUES (?,?,?,?,?)',
                        (content,sizes.get(record['path']),media_type(record['path']),None,now))
                    if sizes.get(record['path']) is not None:
                        db.execute('UPDATE content_objects SET size=? WHERE sha256=? AND size IS NULL',
                            (sizes[record['path']],content))
            # The producing run is the authority on what its operation consumed,
            # so its edges are replaced as a set: a positional overwrite alone
            # would leave behind the tail of a longer earlier registration.
            if own:db.execute('DELETE FROM operation_inputs WHERE operation_id=?',(identity,))
            for port,rows in _edges(group).items():
                for position,(occurrence,asset,role) in enumerate(rows):
                    db.execute(f'INSERT {keyword} INTO operation_inputs VALUES (?,?,?,?,?,?)',
                        (identity,port,position,occurrence,asset,role))
    db.execute('INSERT OR REPLACE INTO catalog_indexed_runs VALUES (?,?,?,?)',
        (project_id,run_id,CATALOG_VERSION,now))


def backfill_operations(db):
    """Derive the operations of runs indexed before the catalog recorded them.

    Idempotent and resumable: a run is processed once per catalog version, and
    an interrupted pass leaves the runs it had finished marked.
    """
    rows=db.execute('''SELECT r.id,r.project_id,r.stages FROM runs r
        WHERE NOT EXISTS (SELECT 1 FROM catalog_indexed_runs c
            WHERE c.project_id=r.project_id AND c.run_id=r.id AND c.version>=?)''',(CATALOG_VERSION,)).fetchall()
    for row in rows:
        try:
            stages=json.loads(row['stages'])
        except (TypeError,ValueError):
            _issue(db,row['project_id'],row['id'],'','unreadable_run',{})
            continue
        index_operations(db,row['project_id'],row['id'],stages)
    return len(rows)


def _instances(db,project_id,occurrences):
    """Artifact rows for occurrences, one row each, from the index alone."""
    occurrences=[value for value in dict.fromkeys(occurrences) if value]
    if not occurrences:return []
    found={}
    for start in range(0,len(occurrences),400):
        chunk=occurrences[start:start+400]
        placeholders=','.join('?' for _ in chunk)
        for row in db.execute(f'''SELECT * FROM artifact_instances
            WHERE project_id=? AND occurrence_id IN ({placeholders}) ORDER BY created''',(project_id,*chunk)):
            found.setdefault(row['occurrence_id'],dict(row))
    return [found[value] for value in occurrences if value in found]


def artifacts_by_entity(db,project_id,entity_type,entity_id,qualifiers=None):
    """Artifacts published for one scientific entity (REQ-CAT-003)."""
    rows=db.execute('''SELECT occurrence_id,qualifiers FROM artifact_entities
        WHERE project_id=? AND entity_type=? AND entity_id=?''',(project_id,entity_type,entity_id)).fetchall()
    if qualifiers:
        rows=[row for row in rows if all(str(json.loads(row['qualifiers'] or '{}').get(key,''))==str(value)
                                         for key,value in qualifiers.items())]
    return _instances(db,project_id,[row['occurrence_id'] for row in rows])


def artifacts_by_content(db,project_id,sha256):
    """Every logical artifact that references one physical content object."""
    return [dict(row) for row in db.execute('''SELECT * FROM artifact_instances
        WHERE project_id=? AND content_sha256=? ORDER BY created,occurrence_id''',(project_id,sha256))]


def artifacts_by_producer(db,project_id,stage_id,run_id=None,port=None):
    """Artifacts a block published, optionally in one run and on one port."""
    where=['project_id=?','stage_id=?'];values=[project_id,stage_id]
    if run_id:where.append('indexed_run_id=?');values.append(run_id)
    if port:where.append('output_port=?');values.append(port)
    return [dict(row) for row in db.execute('SELECT * FROM artifact_instances WHERE '+
        ' AND '.join(where)+' ORDER BY created,occurrence_id',values)]


def content_object(db,sha256):
    row=db.execute('SELECT * FROM content_objects WHERE sha256=?',(sha256,)).fetchone()
    return dict(row) if row else None


def producing_operation(db,project_id,occurrence_id):
    """The operation that published an occurrence."""
    row=db.execute('''SELECT o.* FROM operation_outputs e
        JOIN operations o ON o.operation_id=e.operation_id
        WHERE e.occurrence_id=? AND o.project_id=?''',(occurrence_id,project_id)).fetchone()
    return dict(row) if row else None


def operation_edges(db,operation_id):
    """Recorded inputs and outputs of one operation."""
    return {'inputs':[dict(row) for row in db.execute('''SELECT * FROM operation_inputs
                WHERE operation_id=? ORDER BY port,position''',(operation_id,))],
            'outputs':[dict(row) for row in db.execute('''SELECT * FROM operation_outputs
                WHERE operation_id=? ORDER BY occurrence_id''',(operation_id,))]}


def operations_of_run(db,project_id,run_id):
    return [dict(row) for row in db.execute('''SELECT * FROM operations
        WHERE project_id=? AND run_id=? ORDER BY started,stage_id''',(project_id,run_id))]


def _walk(db,project_id,occurrence_id,depth,query):
    frontier=[occurrence_id];seen={occurrence_id};order=[];level=0
    while frontier and (depth is None or level<depth):
        found=[]
        for value in frontier:
            for row in db.execute(query,(value,project_id)):
                if row[0] and row[0] not in seen:
                    seen.add(row[0]);found.append(row[0]);order.append(row[0])
        frontier=found;level+=1
    return _instances(db,project_id,order)


_UPSTREAM='''SELECT i.occurrence_id FROM operation_outputs e
    JOIN operations o ON o.operation_id=e.operation_id
    JOIN operation_inputs i ON i.operation_id=e.operation_id
    WHERE e.occurrence_id=? AND o.project_id=? AND i.occurrence_id IS NOT NULL'''

_DOWNSTREAM='''SELECT e.occurrence_id FROM operation_inputs i
    JOIN operations o ON o.operation_id=i.operation_id
    JOIN operation_outputs e ON e.operation_id=i.operation_id
    WHERE i.occurrence_id=? AND o.project_id=?'''


def upstream(db,project_id,occurrence_id,depth=None):
    """Occurrences an artifact was computed from, nearest first."""
    return _walk(db,project_id,occurrence_id,depth,_UPSTREAM)


def downstream(db,project_id,occurrence_id,depth=None):
    """Occurrences computed from an artifact, nearest first."""
    return _walk(db,project_id,occurrence_id,depth,_DOWNSTREAM)
