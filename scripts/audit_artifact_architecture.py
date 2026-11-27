"""Read-only inventory of artifact storage, catalog and reuse for PHASE-00 (#7).

The workspace database is opened with ``mode=ro`` and project folders are only
walked and read. Nothing is written except the JSON report given by ``--output``.
No scientific engine, application service or migration is started.

Example::

    python scripts/audit_artifact_architecture.py \
        --database ~/.local/share/biomolexplorer/workspace.sqlite3 \
        --project large=/path/to/project --hash-duplicates --output audit.json
"""
import argparse
import collections
import hashlib
import json
import os
import platform
import shutil
import sqlite3
import statistics
import subprocess
import sys
import time
from pathlib import Path

RECEPTOR_SUFFIXES = ('.dockprep.pdbqt', '.dockprep.mol2', '.rec.pdb', '.rec.mol2', '.rec.pdbqt', '.dms', '.sph')
HISTORY_AREAS = ('blobs', 'run-stages', 'snapshots', '.removed', '.result-transactions')


def connect(database):
    db = sqlite3.connect('file:' + str(Path(database).resolve()) + '?mode=ro', uri=True)
    db.row_factory = sqlite3.Row
    return db


def sha256(path):
    value = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(4 << 20), b''):
            value.update(block)
    return value.hexdigest()


def describe(values):
    values = list(values)
    if not values:
        return {'n': 0}
    return {'n': len(values), 'sum': sum(values), 'min': min(values), 'median': statistics.median(values),
            'mean': statistics.fmean(values), 'max': max(values)}


def environment(projects, output):
    tools = {}
    for name in ('vina', 'obabel', 'chimera', 'dock6', 'antechamber', 'sander', 'tleap'):
        tools[name] = shutil.which(name)
    mounts = {}
    for label, root in projects.items():
        try:
            mounts[label] = subprocess.run(['findmnt', '-no', 'FSTYPE,SOURCE', '-T', str(root)], capture_output=True,
                                           text=True, check=True).stdout.split()
        except (OSError, subprocess.CalledProcessError):
            mounts[label] = None
    return {'platform': platform.platform(), 'python': sys.version.split()[0], 'cpu_count': os.cpu_count(),
            'memory_bytes': os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES'), 'tools_on_path': tools,
            'project_filesystems': mounts, 'reflink_probe': reflink_probe(output.parent)}


def reflink_probe(directory):
    """Try FICLONE in the report directory (never in a project)."""
    import fcntl
    import tempfile
    try:
        with tempfile.TemporaryDirectory(dir=directory) as folder:
            source, target = Path(folder) / 'a', Path(folder) / 'b'
            source.write_bytes(b'x' * 4096)
            with source.open('rb') as src, target.open('wb') as dst:
                fcntl.ioctl(dst.fileno(), 0x40049409, src.fileno())
            return {'directory': str(directory), 'supported': True}
    except OSError as error:
        return {'directory': str(directory), 'supported': False, 'error': error.strerror}


def database_inventory(db, path, table_sizes):
    report = {'path': str(path), 'bytes': Path(path).stat().st_size, 'tables': {}}
    for (name,) in db.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"):
        report['tables'][name] = {'rows': db.execute(f'SELECT count(*) FROM "{name}"').fetchone()[0]}
    for pragma in ('page_size', 'page_count', 'freelist_count', 'journal_mode'):
        report[pragma] = db.execute(f'PRAGMA {pragma}').fetchone()[0]
    report['free_bytes'] = report['page_size'] * report['freelist_count']
    report['indexes'] = [dict(row) for row in db.execute(
        "SELECT name,tbl_name AS table_name,sql FROM sqlite_master WHERE type='index' ORDER BY tbl_name,name")]
    report['triggers'] = [row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='trigger'")]
    if table_sizes:
        started = time.monotonic()
        try:
            report['object_bytes'] = {row[0]: row[1] for row in db.execute(
                'SELECT name,sum(pgsize) FROM dbstat GROUP BY name ORDER BY 2 DESC')}
        except sqlite3.OperationalError as error:
            report['object_bytes'] = {'error': str(error)}
        report['dbstat_seconds'] = time.monotonic() - started
    projects = []
    for row in db.execute('''SELECT p.id,p.name,p.deleted,p.archived,l.path,
            (SELECT count(*) FROM runs r WHERE r.project_id=p.id) AS runs,
            (SELECT coalesce(sum(length(CAST(stages AS BLOB))),0) FROM runs r WHERE r.project_id=p.id) AS run_json_bytes,
            (SELECT count(*) FROM project_history h WHERE h.project_id=p.id) AS history_events,
            (SELECT count(*) FROM artifact_instances a WHERE a.project_id=p.id) AS artifact_instances
            FROM projects p LEFT JOIN project_locations l ON l.project_id=p.id'''):
        projects.append(dict(row))
    report['projects'] = projects
    return report


def catalog_inventory(db, project_id):
    def one(sql, *args):
        return db.execute(sql, (project_id, *args)).fetchone()
    report = {}
    report['instances'] = one('SELECT count(*) FROM artifact_instances WHERE project_id=?')[0]
    report['distinct_artifact_ids'] = one('SELECT count(DISTINCT artifact_id) FROM artifact_instances WHERE project_id=?')[0]
    report['distinct_occurrences'] = one('SELECT count(DISTINCT occurrence_id) FROM artifact_instances WHERE project_id=?')[0]
    report['distinct_content'] = one('SELECT count(DISTINCT content_sha256) FROM artifact_instances WHERE project_id=?')[0]
    report['without_content_hash'] = one('SELECT count(*) FROM artifact_instances WHERE project_id=? AND content_sha256 IS NULL')[0]
    report['distinct_storage_paths'] = one('SELECT count(DISTINCT storage_path) FROM artifact_instances WHERE project_id=?')[0]
    report['indexed_runs'] = one('SELECT count(DISTINCT indexed_run_id) FROM artifact_instances WHERE project_id=?')[0]
    report['rows_reindexed_from_other_runs'] = one(
        'SELECT count(*) FROM artifact_instances WHERE project_id=? AND producer_run_id IS NOT NULL AND producer_run_id<>indexed_run_id')[0]
    row = one('''SELECT count(*) AS content_with_many_paths,coalesce(sum(paths-1),0) AS extra_paths FROM (
        SELECT content_sha256,count(DISTINCT storage_path) AS paths FROM artifact_instances
        WHERE project_id=? AND content_sha256 IS NOT NULL GROUP BY content_sha256 HAVING paths>1)''')
    report['content_with_multiple_storage_paths'] = dict(row)
    report['json_bytes'] = dict(one('''SELECT sum(length(metadata)) AS metadata,sum(length(dependencies)) AS dependencies,
        sum(length(provenance)) AS provenance,sum(length(storage_path)) AS storage_path FROM artifact_instances WHERE project_id=?'''))
    report['by_type'] = {r[0]: r[1] for r in db.execute(
        'SELECT artifact_type,count(*) FROM artifact_instances WHERE project_id=? GROUP BY 1 ORDER BY 2 DESC', (project_id,))}
    report['by_state'] = {r[0]: r[1] for r in db.execute(
        'SELECT state,count(*) FROM artifact_instances WHERE project_id=? GROUP BY 1', (project_id,))}
    report['entities_by_type'] = {r[0]: r[1] for r in db.execute(
        'SELECT entity_type,count(*) FROM artifact_entities WHERE project_id=? GROUP BY 1 ORDER BY 2 DESC', (project_id,))}
    report['entity_index_states'] = {r[0]: r[1] for r in db.execute(
        'SELECT state,count(*) FROM artifact_entity_indexes WHERE project_id=? GROUP BY 1', (project_id,))}
    report['execution_manifests'] = dict(one(
        'SELECT count(*) AS rows,coalesce(sum(length(manifest)),0) AS bytes FROM execution_manifests WHERE project_id=?'))
    report['migration_issues'] = one('SELECT count(*) FROM artifact_migration_issues WHERE project_id=?')[0]
    return report


def run_inventory(db, project_id):
    """Per-stage execution, reuse and descriptor-size facts from runs.stages."""
    runs, stages = [], []
    for row in db.execute('SELECT id,status,created,updated,length(CAST(stages AS BLOB)) AS bytes FROM runs WHERE project_id=? ORDER BY created',
                          (project_id,)).fetchall():
        raw = db.execute('SELECT stages FROM runs WHERE id=?', (row['id'],)).fetchone()[0]
        started = time.monotonic()
        items = json.loads(raw)
        parse = time.monotonic() - started
        del raw
        record_bytes = 0
        for stage in items:
            size = len(json.dumps(stage.get('artifact_records', []), separators=(',', ':')))
            record_bytes += size
            begin, end = stage.get('started_at'), stage.get('finished_at')
            stages.append({'run_id': row['id'], 'stage_id': stage.get('id'), 'operation': stage.get('operation'),
                           'status': stage.get('status'), 'seconds': (end - begin) if begin and end else None,
                           'reused': bool(stage.get('reused')), 'reuse_kind': stage.get('reuse_kind'),
                           'reused_from_run': stage.get('reused_from_run'), 'cache_key': stage.get('cache_key'),
                           'request_signature': (stage.get('reuse_provenance') or {}).get('request_signature'),
                           'artifacts': len(stage.get('artifacts', [])), 'artifact_records': len(stage.get('artifact_records', [])),
                           'artifact_record_bytes': size, 'batches': len(stage.get('batches', []))})
        runs.append({'run_id': row['id'], 'status': row['status'], 'wall_seconds': row['updated'] - row['created'],
                     'json_bytes': row['bytes'], 'artifact_record_bytes': record_bytes, 'stages': len(items),
                     'json_parse_seconds': parse})
        del items
    return runs, stages, reuse_summary(stages)


def reuse_summary(stages):
    """Historical reuse, from recorded stage states only (no inference of timings)."""
    finished = [s for s in stages if s['status'] == 'succeeded']
    computed = [s for s in finished if not s['reused']]
    reused = [s for s in finished if s['reused']]
    seen, redundant = {}, []
    for stage in finished:
        key = stage['cache_key']
        if not key:
            continue
        if key in seen and not stage['reused']:
            redundant.append(stage)
        seen.setdefault(key, stage)
    by_operation = collections.defaultdict(lambda: {'computed': 0, 'reused': 0, 'computed_seconds': 0.0,
                                                    'redundant_recomputations': 0, 'redundant_seconds': 0.0})
    for stage in finished:
        item = by_operation[stage['operation']]
        if stage['reused']:
            item['reused'] += 1
        else:
            item['computed'] += 1
            item['computed_seconds'] += stage['seconds'] or 0.0
    for stage in redundant:
        by_operation[stage['operation']]['redundant_recomputations'] += 1
        by_operation[stage['operation']]['redundant_seconds'] += stage['seconds'] or 0.0
    eligible = len(reused) + len(redundant)
    return {'succeeded_stages': len(finished), 'computed': len(computed), 'reused': len(reused),
            'stage_hit_rate_percent': 100 * len(reused) / len(finished) if finished else None,
            'redundant_recomputations': len(redundant),
            'redundant_seconds': sum(s['seconds'] or 0.0 for s in redundant),
            'historical_avoided_computation_ratio_percent': 100 * len(reused) / eligible if eligible else None,
            'by_operation': dict(by_operation),
            'note': 'A stage is redundant when its cache_key already succeeded earlier and it was computed again. '
                    'Recomputations requested by the user (reuse declined) are included; they are not necessarily defects.'}


def area_of(relative):
    parts = relative.parts
    if not parts:
        return '(root)'
    if parts[0] == '.history' and len(parts) > 1:
        return '.history/' + (parts[1] if parts[1] in HISTORY_AREAS else 'other')
    if parts[0] == 'runs' and '.biomolexplorer' in parts and 'jobs' in parts:
        return 'runs (job folders)'
    if parts[0] == '.input-cache' and len(parts) > 1 and parts[1] == '.calculation-cache':
        return '.input-cache/.calculation-cache'
    return parts[0] if len(parts) > 1 else '(root files)'


def extension(name):
    lower = name.lower()
    for suffix in ('.csv.gz', '.biomol-view.json', '.dockprep.pdbqt', '.dockprep.mol2'):
        if lower.endswith(suffix):
            return suffix
    return Path(lower).suffix or '(none)'


def project_inventory(root, hash_duplicates, snapshot_sample):
    root = Path(root)
    areas = collections.defaultdict(lambda: {'files': 0, 'bytes': 0, 'directories': 0})
    extensions = collections.defaultdict(lambda: {'files': 0, 'bytes': 0})
    by_size = collections.defaultdict(list)
    receptors = collections.defaultdict(list)
    started = time.monotonic()
    for folder, subfolders, names in os.walk(root, followlinks=False):
        relative = Path(folder).relative_to(root)
        areas[area_of(relative / 'x')]['directories'] += 1
        for name in names:
            path = Path(folder) / name
            try:
                info = path.lstat()
            except OSError:
                continue
            if not path.is_file() or path.is_symlink():
                continue
            area = area_of((relative / name))
            areas[area]['files'] += 1
            areas[area]['bytes'] += info.st_size
            if area.startswith('runs'):
                ext = extension(name)
                extensions[ext]['files'] += 1
                extensions[ext]['bytes'] += info.st_size
            if info.st_size and info.st_nlink == 1:
                by_size[info.st_size].append((path, area))
            if name.lower().endswith(RECEPTOR_SUFFIXES):
                receptors[name].append((path, info.st_size))
    walk_seconds = time.monotonic() - started
    report = {'root': str(root), 'walk_seconds_warm_metadata': walk_seconds,
              'total_bytes': sum(a['bytes'] for a in areas.values()), 'total_files': sum(a['files'] for a in areas.values()),
              'areas': dict(sorted(areas.items(), key=lambda item: -item[1]['bytes'])),
              'run_extensions': dict(sorted(extensions.items(), key=lambda item: -item[1]['bytes'])[:25])}
    if hash_duplicates:
        report['duplicates'] = duplicates(root, by_size)
        report['receptors'] = receptor_copies(root, receptors)
    report['snapshots'] = snapshots(root / '.history' / 'snapshots', snapshot_sample)
    return report


def duplicates(root, by_size):
    started = time.monotonic()
    hashed_bytes = 0
    groups = collections.defaultdict(list)
    for size, entries in by_size.items():
        if len(entries) < 2:
            continue
        for path, area in entries:
            try:
                groups[(size, sha256(path))].append((path, area))
                hashed_bytes += size
            except OSError:
                continue
    redundant = 0
    by_areas = collections.Counter()
    by_ext = collections.Counter()
    top = []
    for (size, key), entries in groups.items():
        if len(entries) < 2:
            continue
        extra = size * (len(entries) - 1)
        redundant += extra
        by_areas[' + '.join(sorted({area for _, area in entries}))] += extra
        by_ext[extension(entries[0][0].name)] += extra
        top.append((extra, size, key, entries))
    top.sort(key=lambda item: -item[0])
    return {'hashed_bytes': hashed_bytes, 'hash_seconds': time.monotonic() - started,
            'duplicate_groups': sum(1 for entries in groups.values() if len(entries) > 1),
            'redundant_bytes': redundant,
            'redundant_bytes_by_areas': dict(by_areas.most_common(15)),
            'redundant_bytes_by_extension': dict(by_ext.most_common(15)),
            'largest_groups': [{'redundant_bytes': extra, 'size': size, 'copies': len(entries), 'sha256': key,
                                'examples': [str(path.relative_to(root)) for path, _ in entries[:4]]}
                               for extra, size, key, entries in top[:15]],
            'note': 'Files with st_nlink>1 are excluded (already shared). Zero-byte files are ignored.'}


def receptor_copies(root, receptors):
    rows = []
    for name, entries in receptors.items():
        hashes = collections.Counter()
        for path, _ in entries:
            try:
                hashes[sha256(path)] += 1
            except OSError:
                pass
        rows.append({'name': name, 'copies': len(entries), 'distinct_contents': len(hashes),
                     'bytes': sum(size for _, size in entries)})
    rows.sort(key=lambda row: -row['copies'])
    return {'distinct_names': len(rows), 'files': sum(r['copies'] for r in rows),
            'distinct_contents': sum(r['distinct_contents'] for r in rows),
            'bytes': sum(r['bytes'] for r in rows), 'most_copied': rows[:20]}


def snapshots(folder, sample, max_parse_bytes=1536 << 20):
    if not folder.is_dir():
        return None
    files = sorted(folder.glob('*.json'), key=lambda p: p.stat().st_size)
    sizes = [p.stat().st_size for p in files]
    report = {'files': len(files), 'size': describe(sizes)}
    references = compact = 0
    for path in files:
        with path.open('rb') as stream:
            head = stream.read(4096)
        if path.stat().st_size < 4096 and b'_snapshot_reference' in head:
            references += 1
        if head.startswith(b'{"_compact_runs": true'):
            compact += 1
    report['reference_only_files'] = references
    report['compact_runs_files'] = compact
    report['compact_runs_bytes'] = sum(p.stat().st_size for p in files if p.open('rb').read(24).startswith(b'{"_compact_runs": true'))
    report['larger_than_1GiB'] = sum(1 for size in sizes if size > 1 << 30)
    parseable = [p for p in files if p.stat().st_size <= max_parse_bytes]
    chosen = parseable[-sample:] if sample else []
    report['sample_limit_bytes'] = max_parse_bytes
    report['sample'] = []
    for path in chosen:
        started = time.monotonic()
        state = json.loads(path.read_text(encoding='utf-8'))
        parse = time.monotonic() - started
        parts = {key: len(json.dumps(value, ensure_ascii=False)) for key, value in state.items()}
        report['sample'].append({'file': path.name, 'bytes': path.stat().st_size, 'parse_seconds': parse,
                                 'compact_runs': bool(state.get('_compact_runs')), 'files_listed': len(state.get('files', {})),
                                 'runs': len(state.get('runs', [])), 'bytes_by_key': dict(sorted(parts.items(), key=lambda i: -i[1]))})
        del state
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--database', type=Path, required=True)
    parser.add_argument('--project', action='append', default=[], metavar='LABEL=DIR',
                        help='Project folder to inventory; the project id is resolved from project_locations or the folder name.')
    parser.add_argument('--hash-duplicates', action='store_true', help='Hash size-colliding files to measure duplicated bytes.')
    parser.add_argument('--snapshot-sample', type=int, default=3, help='Parse the N largest history snapshots.')
    parser.add_argument('--table-sizes', action='store_true', help='Measure per-table bytes with dbstat (slow on large databases).')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    projects = dict(item.split('=', 1) for item in args.project)
    started = time.monotonic()
    db = connect(args.database)
    report = {'generated_unix': time.time(), 'read_only': True, 'environment': environment(projects, args.output),
              'database': database_inventory(db, args.database, args.table_sizes), 'projects': {}}
    locations = {Path(row['path']).resolve(): row['project_id'] for row in db.execute('SELECT project_id,path FROM project_locations')}
    for label, folder in projects.items():
        folder = Path(folder).resolve()
        project_id = locations.get(folder) or folder.name
        print(f'[{label}] {folder} ({project_id})', flush=True)
        item = {'project_id': project_id}
        item['catalog'] = catalog_inventory(db, project_id)
        item['runs'], item['stages'], item['reuse'] = run_inventory(db, project_id)
        item['storage'] = project_inventory(folder, args.hash_duplicates, args.snapshot_sample)
        report['projects'][label] = item
    report['elapsed_seconds'] = time.monotonic() - started
    args.output.write_text(json.dumps(report, indent=1, default=str) + '\n')
    print(f'Wrote {args.output} in {report["elapsed_seconds"]:.1f} s')


if __name__ == '__main__':
    main()
