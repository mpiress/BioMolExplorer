"""Read-only retrieval baseline for PHASE-00 (#7): discovery, catalog and result reads.

Every case runs ``--repeats`` times in two modes:

* ``cold``: the files involved are evicted from the OS page cache with
  ``posix_fadvise(DONTNEED)`` before each repetition. Directory metadata
  (dentries/inodes) cannot be evicted without root, so traversal cases are
  reported as ``warm`` only.
* ``warm``: one untimed priming call, then the timed repetitions.

The workspace database is opened with ``mode=ro``. Projection readers are only
called when their cached read model is valid, so they never write. Disposable
CSV indexes are built under ``--scratch``, never inside a project. No engine runs.
"""
import argparse
import csv
import hashlib
import json
import os
import random
import sqlite3
import statistics
import sys
import tempfile
import time
import zlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from biomolexplorer import result_catalog  # noqa: E402
from biomolexplorer.result_index import ResultIndex  # noqa: E402

csv.field_size_limit(sys.maxsize)


class ReadOnlyStore:
    """The subset of WorkspaceStore used by result_catalog readers, read-only."""
    def __init__(self, database):
        self.database = Path(database).resolve()

    def connect(self):
        db = sqlite3.connect('file:' + str(self.database) + '?mode=ro', uri=True, timeout=30)
        db.row_factory = sqlite3.Row
        return db


def evict(paths):
    for path in paths:
        try:
            fd = os.open(path, os.O_RDONLY)
        except OSError:
            continue
        try:
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
        finally:
            os.close(fd)


def stats(values):
    return {'n': len(values), 'mean': statistics.fmean(values), 'stdev': statistics.stdev(values) if len(values) > 1 else 0.0,
            'min': min(values), 'median': statistics.median(values), 'max': max(values)}


def measure(function, repeats, cold_paths=None):
    result = {}
    outputs = set()
    modes = (('cold', True), ('warm', False)) if cold_paths else (('warm', False),)
    for mode, cold in modes:
        if not cold:
            outputs.add(repr(function()))
        timings = []
        for _ in range(repeats):
            if cold:
                evict(cold_paths)
            started = time.perf_counter()
            output = function()
            timings.append(time.perf_counter() - started)
            outputs.add(repr(output))
        result[mode] = stats(timings)
    result['consistent_output'] = len(outputs) == 1
    result['output'] = next(iter(outputs)) if len(outputs) == 1 else sorted(outputs)[:3]
    return result


def database_files(database):
    base = Path(database)
    return [str(base)] + [str(base) + suffix for suffix in ('-wal', '-shm') if Path(str(base) + suffix).exists()]


def plan(db, sql, args):
    return [row[3] for row in db.execute('EXPLAIN QUERY PLAN ' + sql, args)]


def digest(stream):
    value = hashlib.sha256()
    for chunk in iter(lambda: stream.read(1 << 20), b''):
        value.update(chunk)
    stream.seek(0)
    return value.hexdigest()


def project_cases(label, root, project_id, database, repeats, scratch, seed):
    store = ReadOnlyStore(database)
    db_files = database_files(database)
    cases = {}
    runs_dir = Path(root) / 'runs'

    def walk():
        total = 0
        for _, _, names in os.walk(runs_dir):
            total += len(names)
        return total

    def rglob_csv():
        return sum(1 for _ in runs_dir.rglob('*.csv'))

    cases['discovery_os_walk_runs'] = measure(walk, repeats)
    cases['discovery_rglob_csv'] = measure(rglob_csv, repeats)

    with store.connect() as db:
        run = db.execute('''SELECT id,rowid,length(CAST(stages AS BLOB)) AS bytes FROM runs WHERE project_id=? AND status='succeeded'
                            ORDER BY created DESC LIMIT 1''', (project_id,)).fetchone()
        if run is None:
            return {'skipped': 'no succeeded run'}
        stage = db.execute('''SELECT stage_id,count(*) AS n FROM artifact_instances WHERE project_id=? AND indexed_run_id=?
                              GROUP BY stage_id ORDER BY n DESC LIMIT 1''', (project_id, run['id'])).fetchone()
        rng = random.Random(seed)
        entities = db.execute('SELECT entity_type,entity_id FROM artifact_entities WHERE project_id=? ORDER BY rowid LIMIT 20000',
                              (project_id,)).fetchall()
        entities = rng.sample(entities, min(50, len(entities)))
        contents = [r[0] for r in db.execute('''SELECT DISTINCT content_sha256 FROM artifact_instances WHERE project_id=?
                    AND content_sha256 IS NOT NULL ORDER BY rowid LIMIT 20000''', (project_id,))]
        contents = rng.sample(contents, min(50, len(contents)))
        warm_projection = db.execute('''SELECT 1 FROM result_read_catalog c LEFT JOIN result_read_versions v ON v.run_id=c.run_id
            WHERE c.run_id=? AND c.revision=coalesce(v.revision,0) AND c.schema_version=?''',
                                     (run['id'], result_catalog.VERSION)).fetchone() is not None
        stage_sql = 'SELECT occurrence_id,storage_path,content_sha256 FROM artifact_instances WHERE project_id=? AND stage_id=? AND indexed_run_id=?'
        entity_sql = 'SELECT occurrence_id,artifact_id FROM artifact_entities WHERE project_id=? AND entity_type=? AND entity_id=?'
        content_sql = 'SELECT occurrence_id,storage_path FROM artifact_instances WHERE project_id=? AND content_sha256=?'
        plans = {'stage_artifacts': plan(db, stage_sql, (project_id, stage['stage_id'] if stage else '', run['id'])),
                 'entity_lookup': plan(db, entity_sql, (project_id, 'compound', 'x')),
                 'content_lookup': plan(db, content_sql, (project_id, 'x'))}
    cases['_context'] = {'run_id': run['id'], 'run_json_bytes': run['bytes'], 'largest_stage': dict(stage) if stage else None,
                         'sampled_entities': len(entities), 'sampled_contents': len(contents),
                         'projection_cached': warm_projection, 'query_plans': plans}

    def stage_artifacts():
        with store.connect() as db:
            return len(db.execute(stage_sql, (project_id, stage['stage_id'], run['id'])).fetchall())

    def entity_lookup():
        with store.connect() as db:
            return sum(len(db.execute(entity_sql, (project_id, kind, ident)).fetchall()) for kind, ident in entities)

    def content_lookup():
        with store.connect() as db:
            return sum(len(db.execute(content_sql, (project_id, key)).fetchall()) for key in contents)

    def run_json():
        with store.connect() as db:
            return len(json.loads(db.execute('SELECT stages FROM runs WHERE id=?', (run['id'],)).fetchone()[0]))

    def run_stream_projection():
        with store.connect() as db:
            with db.blobopen('runs', 'stages', run['rowid'], readonly=True) as blob:
                return len(result_catalog._project_blob(blob))

    if stage:
        cases['indexed_stage_artifacts'] = measure(stage_artifacts, repeats, db_files)
    if entities:
        cases['indexed_entity_lookup_x50'] = measure(entity_lookup, repeats, db_files)
    if contents:
        cases['indexed_content_lookup_x50'] = measure(content_lookup, repeats, db_files)
    cases['run_full_json_load'] = measure(run_json, repeats, db_files)
    cases['run_streamed_projection'] = measure(run_stream_projection, repeats, db_files)
    if warm_projection:
        cases['run_cached_projection_read_run'] = measure(lambda: len(result_catalog.read_run(store, run['id'])['stages']), repeats, db_files)
        cases['run_cached_proofs'] = measure(lambda: len(result_catalog.read_proofs(store, run['id'])), repeats, db_files)

    with store.connect() as db:
        registered = {Path(r[0]) for r in db.execute('''SELECT DISTINCT storage_path FROM artifact_instances
            WHERE project_id=? AND indexed_run_id=? AND storage_path LIKE '%.csv' ''', (project_id, run['id']))}
    tables = sorted((p for p in registered if p.is_file()), key=lambda p: -p.stat().st_size)[:3]
    cases['_tables'] = [{'path': str(p.relative_to(root)), 'bytes': p.stat().st_size} for p in tables]
    for number, table in enumerate(tables, 1):
        def full_read(table=table):
            with table.open(encoding='utf-8-sig', newline='') as stream:
                return sum(1 for _ in csv.DictReader(stream))
        cases[f'csv{number}_full_read'] = measure(full_read, repeats, [str(table)])
        with tempfile.TemporaryDirectory(dir=scratch, prefix='result-index-') as folder:
            index = ResultIndex(Path(folder), table, digest)
            evict([str(table)])
            started = time.perf_counter()
            first = index.page(0, 25, '', None, True, False)
            build = time.perf_counter() - started
            offset = first['total'] // 2
            pages = measure(lambda index=index: [i for i, _ in index.page(offset, 25, '', None, True, False)['visible']], repeats,
                            [str(index.database)])
            cases[f'csv{number}_indexed'] = {'cold_build_and_first_page_seconds': build, 'rows': first['total'],
                                             'index_bytes': index.database.stat().st_size, 'page_offset': offset,
                                             'page_25_rows': pages}

    snaps = sorted((Path(root) / '.history' / 'snapshots').glob('*.json'), key=lambda p: p.stat().st_size)
    snaps = [p for p in snaps if p.stat().st_size <= 1536 << 20]
    if snaps:
        middle = snaps[len(snaps) // 2]
        cases['_snapshot'] = {'file': middle.name, 'bytes': middle.stat().st_size}
        cases['history_snapshot_parse_median'] = measure(lambda: len(json.loads(middle.read_bytes())), max(2, repeats // 2), [str(middle)])
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--database', type=Path, required=True)
    parser.add_argument('--project', action='append', default=[], metavar='LABEL=DIR')
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--seed', type=int, default=6)
    parser.add_argument('--scratch', type=Path, help='Directory for disposable CSV indexes (default: next to --output).')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    scratch = args.scratch or args.output.parent
    scratch.mkdir(parents=True, exist_ok=True)
    with ReadOnlyStore(args.database).connect() as db:
        locations = {Path(r[1]).resolve(): r[0] for r in db.execute('SELECT project_id,path FROM project_locations')}
    report = {'generated_unix': time.time(), 'repeats': args.repeats, 'seed': args.seed, 'read_only': True,
              'cold_method': 'posix_fadvise(DONTNEED) on the involved files before each repetition; directory metadata is not evicted',
              'projects': {}}
    started = time.monotonic()
    for item in args.project:
        label, folder = item.split('=', 1)
        folder = Path(folder).resolve()
        project_id = locations.get(folder) or folder.name
        print(f'[{label}] {folder}', flush=True)
        report['projects'][label] = {'project_id': project_id, **project_cases(label, folder, project_id, args.database,
                                                                               args.repeats, scratch, args.seed)}
    report['elapsed_seconds'] = time.monotonic() - started
    args.output.write_text(json.dumps(report, indent=1, default=str) + '\n')
    print(f'Wrote {args.output} in {report["elapsed_seconds"]:.1f} s')


if __name__ == '__main__':
    main()
