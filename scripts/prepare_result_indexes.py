"""Build disposable read indexes for an existing local project, without rerunning it."""
import argparse
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import sqlite3
import sys
import time

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from biomolexplorer.result_catalog import create_schema,read_run,read_stage
from biomolexplorer.result_index import ResultIndex,MIN_INDEX_BYTES
from biomolexplorer.docking_data import result_tables
from biomolexplorer.input_validation import columns


def digest(stream):
    value=hashlib.sha256()
    for chunk in iter(lambda:stream.read(1048576),b''):value.update(chunk)
    stream.seek(0)
    return value.hexdigest()


class LocalStore:
    def __init__(self,database):self.database=database

    @contextmanager
    def connect(self):
        connection=sqlite3.connect(self.database,timeout=30);connection.row_factory=sqlite3.Row
        try:
            with connection:yield connection
        finally:connection.close()


def source_digest(db,run_id):
    row=db.execute('SELECT rowid FROM runs WHERE id=?',(run_id,)).fetchone()
    with db.blobopen('runs','stages',row[0],readonly=True) as blob:return digest(blob)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--database',required=True,type=Path)
    parser.add_argument('--project-id',required=True)
    parser.add_argument('--run-id')
    parser.add_argument('--limit',type=int,default=8)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    if not 1<=args.limit<=100:raise ValueError('Use a limit between 1 and 100 tables')
    if not args.database.is_file():raise FileNotFoundError(args.database)
    store=LocalStore(args.database.resolve());started=time.monotonic()
    with store.connect() as db:
        run=db.execute('SELECT id,status,updated FROM runs WHERE project_id=? '+
            ('AND id=?' if args.run_id else 'ORDER BY updated DESC LIMIT 1'),
            (args.project_id,args.run_id) if args.run_id else (args.project_id,)).fetchone()
        if run is None:raise ValueError('Project run not found')
        if run['status'] in ('queued','running','awaiting_input'):raise ValueError('Finish or cancel the run before prewarming its indexes')
        location=db.execute('SELECT path FROM project_locations WHERE project_id=?',(args.project_id,)).fetchone()
        root=Path(location[0]) if location else args.database.parent/'projects'/args.project_id
        root=root.resolve();before=source_digest(db,run['id']);create_schema(db)
    review=read_run(store,run['id']);report=dict(run_id=run['id'],tables=[])
    priority={'docking_vina':0,'docking_dock6':1,'consensus':2,'admet':3}
    for stage in sorted(review['stages'],key=lambda item:priority.get(item.get('operation'),10)):
        paths=[Path(value) for value in stage.get('artifacts',[]) if str(value).lower().endswith('.csv')]
        if stage.get('operation') in ('docking_vina','docking_dock6'):paths=result_tables(paths)
        for path in paths:
            resolved=path.resolve()
            if not resolved.is_relative_to(root):raise ValueError('Source is outside the project')
            if not resolved.is_file() or resolved.stat().st_size<MIN_INDEX_BYTES:continue
            fields=columns(resolved)
            if not {'molecule_chembl_id','canonical_smiles'}<=fields:continue
            sort='score' if 'score' in fields else 'normalized_score' if 'normalized_score' in fields else None
            before_table=time.monotonic()
            result=ResultIndex(root,resolved,digest).page(0,25,'',sort,True,True)
            expected=stage.get('artifact_manifest',{}).get(str(resolved))
            if expected and expected!=result['version']:raise ValueError('Source checksum differs from its execution manifest')
            report['tables'].append(dict(path=str(resolved),rows=result['total'],seconds=time.monotonic()-before_table))
            if len(report['tables'])>=args.limit:break
        if len(report['tables'])>=args.limit:break
    before_warm=time.monotonic();read_stage(store,run['id'],review['stages'][0]['id'])
    report['cached_stage_seconds']=time.monotonic()-before_warm
    with store.connect() as db:
        after=source_digest(db,run['id'])
        report['authoritative_stages_unchanged']=before==after
        report['catalog_bytes']=db.execute('SELECT length(stages) FROM result_read_catalog WHERE run_id=?',(run['id'],)).fetchone()[0]
    if before!=after:raise ValueError('Run changed concurrently; repeat the preparation')
    report['seconds']=time.monotonic()-started
    args.report.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
