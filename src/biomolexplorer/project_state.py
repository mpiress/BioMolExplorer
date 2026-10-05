"""Portable project snapshots, audit history and content-addressed file versions."""
import hashlib
import json
import os
import shutil
import time
from pathlib import Path
from uuid import uuid4

EXCLUDED = {'.history', '.exports', '.input-cache', 'project.json'}
_DIGESTS = {}


def translate(value, old, new):
    if isinstance(value, str):
        return new + value[len(old):] if value == old or value.startswith(old + '/') else value
    if isinstance(value, list):
        return [translate(v, old, new) for v in value]
    if isinstance(value, dict):
        return {translate(k, old, new): translate(v, old, new) for k, v in value.items()}
    return value


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def snapshot(store, db, project_id):
    root = store.project_dir(project_id,db)
    project = dict(db.execute('SELECT * FROM projects WHERE id=?', (project_id,)).fetchone())
    project['pipeline'] = json.loads(project['pipeline'])
    project['tags'] = json.loads(project['tags'])
    assets = [dict(r) for r in db.execute('SELECT * FROM assets WHERE project_id=?', (project_id,))]
    runs = [dict(r) for r in db.execute("SELECT * FROM runs WHERE project_id=? AND status NOT IN ('queued','running','awaiting_input')", (project_id,))]
    for run in runs:
        run['stages'] = json.loads(run['stages'])
    members = [dict(r) for r in db.execute('SELECT m.*,u.email FROM members m JOIN users u ON u.id=m.user_id WHERE project_id=?', (project_id,))]
    active = {r[0] for r in db.execute("SELECT id FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')", (project_id,))}
    registered_assets = {a['id'] for a in assets}
    files, directories = {}, []
    blobs = root / '.history' / 'blobs'
    blobs.mkdir(parents=True, exist_ok=True, mode=0o700)
    for folder, subfolders, filenames in os.walk(root, followlinks=False):
        relative = Path(folder).relative_to(root)
        subfolders[:] = [name for name in subfolders if name not in EXCLUDED and not (Path(folder) / name).is_symlink()
                         and not (relative == Path('runs') and name in active)
                         and not (relative == Path('assets') and name not in registered_assets)]
        if relative.parts:
            directories.append(relative.as_posix())
        for name in filenames:
            if name in EXCLUDED:
                continue
            path = Path(folder) / name
            if path.is_symlink() or not path.is_file():
                raise ValueError('O projeto contém um arquivo simbólico. Use arquivos regulares para versionar e exportar.')
            stat = path.stat()
            signature=(stat.st_dev,stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns)
            cached=_DIGESTS.get(str(path))
            key=cached[1] if cached and cached[0]==signature else digest(path)
            _DIGESTS[str(path)]=(signature,key)
            blob = blobs / key
            if not blob.exists():
                temporary = blobs / (key + '.' + uuid4().hex)
                shutil.copyfile(path, temporary)
                temporary.replace(blob)
            files[path.relative_to(root).as_posix()] = {'sha256': key, 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns}
    edits=[dict(r) for r in db.execute('SELECT * FROM compound_edits WHERE project_id=?',(project_id,))] if db.execute("SELECT 1 FROM sqlite_master WHERE name='compound_edits'").fetchone() else []
    state = translate({'project': project, 'assets': assets, 'runs': runs, 'members': members, 'edits':edits}, str(root), '@project')
    state.update(files=files, directories=directories)
    return state


def mirror(store, db, project_id):
    row = dict(db.execute('SELECT * FROM projects WHERE id=?', (project_id,)).fetchone())
    for key in ('tags', 'pipeline'):
        row[key] = json.loads(row[key])
    root = store.project_dir(project_id,db)
    value = translate(row, str(root), '@project')
    value.pop('owner', None)
    temporary = root / ('.project-' + uuid4().hex)
    temporary.write_text(json.dumps({'format': 'biomolexplorer-project', 'version': 1, **value}, ensure_ascii=False, indent=2), encoding='utf-8')
    temporary.replace(root / 'project.json')


def record(store, db, project_id, user_id, action, summary, before=None):
    root = store.project_dir(project_id,db)
    state = snapshot(store, db, project_id)
    folder = root / '.history' / 'snapshots'
    folder.mkdir(parents=True, exist_ok=True)
    event_id = uuid4().hex
    for suffix, value in (('after', state), ('before', before)):
        if value is not None:
            (folder / f'{event_id}-{suffix}.json').write_text(json.dumps(value, ensure_ascii=False), encoding='utf-8')
    actor = db.execute('SELECT name,email FROM users WHERE id=?', (user_id,)).fetchone()
    db.execute('INSERT INTO project_history VALUES (?,?,?,?,?,?,?,?,?)',
               (event_id, project_id, user_id, actor['name'] if actor else 'Sistema', actor['email'] if actor else '',
                action, summary, time.time(), int(before is not None)))
    db.execute('UPDATE projects SET updated=? WHERE id=?', (time.time(), project_id))
    mirror(store, db, project_id)
    return event_id


def read_snapshot(store, project_id, event_id, side='before'):
    path = store.project_dir(project_id) / '.history' / 'snapshots' / f'{event_id}-{side}.json'
    return json.loads(path.read_text(encoding='utf-8'))


def restore_database(store, db, project_id, state, owner=None):
    root = store.project_dir(project_id)
    state = translate(state, '@project', str(root))
    project = state['project']
    db.execute('UPDATE projects SET name=?,description=?,color=?,tags=?,pipeline=?,archived=?,deleted=?,revision=revision+1,updated=? WHERE id=?',
               (project['name'], project['description'], project['color'], json.dumps(project['tags']),
                json.dumps(project['pipeline']), project['archived'], 0, time.time(), project_id))
    for table in ('assets', 'runs', 'members'):
        db.execute(f'DELETE FROM {table} WHERE project_id=?', (project_id,))
    for row in state['assets']:
        db.execute('INSERT INTO assets VALUES (?,?,?,?,?,?,?)', (row['id'], project_id, row['name'], row['path'], row['size'], row['kind'], row['created']))
    current_owner = owner or db.execute('SELECT owner FROM projects WHERE id=?', (project_id,)).fetchone()[0]
    for row in state['runs']:
        user_id = row['user_id'] if db.execute('SELECT 1 FROM users WHERE id=?', (row['user_id'],)).fetchone() else current_owner
        db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)', (row['id'], project_id, user_id, row['status'], json.dumps(row['stages']), row['created'], row['updated'], row['error']))
    for row in state['members']:
        user = db.execute('SELECT id FROM users WHERE email=?', (row['email'],)).fetchone()
        if user and user[0] != current_owner:
            db.execute('INSERT INTO members VALUES (?,?,?,?)', (project_id, user[0], row['role'], row['accepted'] if owner is None else 0))
    if state.get('edits') or db.execute("SELECT 1 FROM sqlite_master WHERE name='compound_edits'").fetchone():
        db.execute('CREATE TABLE IF NOT EXISTS compound_edits(id TEXT PRIMARY KEY, project_id TEXT, user_id TEXT, run_id TEXT, stage_id TEXT, path TEXT, compound_id TEXT, backup TEXT, created REAL)')
        db.execute('DELETE FROM compound_edits WHERE project_id=?',(project_id,))
        for row in state.get('edits',[]):
            db.execute('INSERT INTO compound_edits VALUES (?,?,?,?,?,?,?,?,?)',(row['id'],project_id,current_owner if owner else row['user_id'],row['run_id'],row['stage_id'],row['path'],row['compound_id'],row['backup'],row['created']))


def rollback(store, token, project_id, event_id):
    actor = store.user(token)
    store.project(token, project_id, 'owner')
    root = store.project_dir(project_id)
    with store.connect() as db:
        db.execute('BEGIN IMMEDIATE')
        store._require_user(actor['id'], project_id, 'owner')
        event = db.execute('SELECT * FROM project_history WHERE id=? AND project_id=? AND has_before=1', (event_id, project_id)).fetchone()
        if not event:
            raise ValueError('Esta alteração não possui uma versão anterior disponível.')
        if db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')", (project_id,)).fetchone():
            raise ValueError('Aguarde ou cancele a execução antes de restaurar uma versão.')
        before = snapshot(store, db, project_id)
        target = read_snapshot(store, project_id, event_id)
        recovery = root / '.history' / ('restore-' + uuid4().hex)
        staged, backup = recovery / 'staged', recovery / 'backup'
        staged.mkdir(parents=True)
        backup.mkdir()
        moved = []
        try:
            for directory in target['directories']:
                (staged / directory).mkdir(parents=True, exist_ok=True)
            for name, info in target['files'].items():
                path = Path(name)
                if path.is_absolute() or '..' in path.parts or path.parts[0] in EXCLUDED:
                    raise ValueError('Caminho inválido na versão do projeto.')
                blob = root / '.history' / 'blobs' / info['sha256']
                if digest(blob) != info['sha256']:
                    raise ValueError('Um arquivo do histórico está incompleto ou foi alterado.')
                destination = staged / path
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(blob, destination)
                os.utime(destination, ns=(info['mtime_ns'], info['mtime_ns']))
            for path in list(root.iterdir()):
                if path.name not in EXCLUDED:
                    path.replace(backup / path.name)
                    moved.append(path.name)
            for path in list(staged.iterdir()):
                path.replace(root / path.name)
            restore_database(store, db, project_id, target)
            record(store, db, project_id, actor['id'], 'rollback', f'Restauração para antes de: {event["summary"]}', before)
            db.commit()
        except BaseException:
            for path in list(root.iterdir()):
                if path.name not in EXCLUDED and path.name not in moved:
                    shutil.rmtree(path) if path.is_dir() else path.unlink()
            for path in backup.iterdir():
                destination = root / path.name
                if destination.exists():
                    shutil.rmtree(destination) if destination.is_dir() else destination.unlink()
                path.replace(destination)
            db.rollback()
            mirror(store,db,project_id)
            raise
        finally:
            shutil.rmtree(recovery, ignore_errors=True)
    return store.project(token, project_id)
