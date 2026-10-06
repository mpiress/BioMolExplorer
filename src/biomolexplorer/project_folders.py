"""Named project folders and permanent, resumable filesystem removal."""
import hashlib
import json
import re
import shutil
import sqlite3
import time
import sys
from pathlib import Path
from uuid import uuid4

from .workspace import AccessDenied, COLORS


def folder_name(name):
    if not isinstance(name,str):
        raise ValueError('Preencha um nome de projeto válido antes de escolher a pasta, sem separadores de caminho.')
    value = name.strip()
    if (not value or len(value) > 120 or value in ('.', '..') or value.startswith('.')
            or value.endswith(('.', ' ')) or any(ord(c) < 32 or c in '/\\<>:"|?*' for c in value)
            or re.fullmatch(r'(CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])(?:\..*)?', value, re.I)):
        raise ValueError('Preencha um nome de projeto válido antes de escolher a pasta, sem separadores de caminho.')
    return value


def check_path(store, db, folder, user_id, replacing=False, deleting=False):
    if folder.resolve() != folder or folder.is_symlink():
        raise ValueError('Escolha uma pasta sem links simbólicos ou segmentos .. no caminho.')
    if (Path(__file__).resolve().is_relative_to(folder) or Path(sys.executable).resolve().is_relative_to(folder)
            or folder == Path.home().resolve() or folder == Path(folder.anchor)
            or store.database.is_relative_to(folder) or folder.is_relative_to(store.staging)):
        raise ValueError('Escolha uma pasta exclusiva para os arquivos do projeto.')
    if folder.exists() and not folder.is_dir():
        raise ValueError('O caminho do projeto já existe e não é uma pasta.')
    existing = None
    for row in db.execute('SELECT id,owner FROM projects'):
        other = store.project_dir(row['id'], db)
        if folder == other and replacing:
            if row['owner'] != user_id:
                raise AccessDenied('Somente o proprietário pode substituir este projeto.')
            if db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')", (row['id'],)).fetchone():
                raise ValueError('Cancele a execução antes de excluir o projeto.' if deleting else
                                 'Conclua ou cancele a execução antes de substituir o projeto.')
            existing = row['id']
        elif folder.is_relative_to(other) or other.is_relative_to(folder):
            raise ValueError('A pasta não pode conter nem pertencer a outro projeto.')
    return existing


def plan(store, token, name, parent, db=None):
    user = store.user(token)
    component = folder_name(name)
    parent = Path(parent).expanduser().absolute()
    if not parent.is_dir():
        raise ValueError('Selecione uma pasta principal existente.')
    folder = parent / component
    if db is None:
        with store.connect() as connection:
            return plan(store, token, name, str(parent), connection)
    existing = check_path(store, db, folder, user['id'], replacing=True)
    # Bind confirmation to the exact directory, contents and current registered project.
    entries = []
    if folder.exists():
        for path in [folder, *sorted(folder.rglob('*'))]:
            info = path.lstat()
            entries.append((str(path.relative_to(folder)), info.st_dev, info.st_ino,
                            info.st_size, info.st_mtime_ns, info.st_ctime_ns))
    project = db.execute('SELECT updated,revision FROM projects WHERE id=?', (existing,)).fetchone() if existing else None
    signature = hashlib.sha256(json.dumps([str(folder), existing, list(project) if project else None, entries]).encode()).hexdigest()
    return {'directory': str(folder), 'parent': str(parent), 'name': component,
            'project_id': existing, 'replace': bool(existing or entries), 'signature': signature}


def purge_rows(store, db, project_id):
    for table in ('members', 'assets', 'uploads', 'runs', 'project_history', 'project_locations'):
        db.execute(f'DELETE FROM {table} WHERE project_id=?', (project_id,))
    if db.execute("SELECT 1 FROM sqlite_master WHERE name='compound_edits'").fetchone():
        db.execute('DELETE FROM compound_edits WHERE project_id=?', (project_id,))
    db.execute('DELETE FROM projects WHERE id=?', (project_id,))


def park(store, db, folder, user_id, project_id):
    parked = folder.parent / ('.biomol-delete-' + uuid4().hex)
    tickets=[r['id'] for r in db.execute('SELECT id FROM uploads WHERE project_id=?', (project_id,))]
    db.execute('INSERT INTO pending_deletions VALUES (?,?,?,?,?,?)',
               (project_id, user_id, str(folder), str(parked), time.time(), json.dumps(tickets)))
    if folder.exists():folder.rename(parked)
    return parked


def finish(store, project_id=None, token=None):
    user = store.user(token) if token else None
    with store.connect() as db:
        rows = db.execute('SELECT * FROM pending_deletions' + (' WHERE id=?' if project_id else ''),
                          (project_id,) if project_id else ()).fetchall()
    for row in rows:
        if user and row['owner'] != user['id']:
            raise AccessDenied('Somente o proprietário pode concluir a exclusão.')
        path = Path(row['parked'])
        if not path.name.startswith('.biomol-delete-') or path.parent != Path(row['original']).parent or path.resolve() != path:
            raise ValueError('Caminho de exclusão pendente inválido.')
        for ticket in json.loads(row['uploads']):
            if not re.fullmatch('[a-f0-9]{32}', ticket):
                raise ValueError('Identificador de upload pendente inválido.')
            (store.staging / ticket).unlink(missing_ok=True)
        if path.exists():
            shutil.rmtree(path)  # unlink internal symlinks; never follows their targets
        with store.connect() as db:
            db.execute('DELETE FROM pending_deletions WHERE id=?', (row['id'],))


def create(store, token, name, parent, description='', color=COLORS[0], tags=None, confirmation=None):
    user = store.user(token)
    store._validate_metadata(name, color, tags or [])
    project_id, now = uuid4().hex, time.time()
    parked = None
    created = False
    folder = None
    try:
        with store.connect() as db:
            db.execute('BEGIN IMMEDIATE')
            current = plan(store, token, name, parent, db)
            folder = Path(current['directory'])
            if current['replace'] and confirmation != current:
                raise ValueError('A pasta já existe ou mudou. Confirme novamente sua substituição permanente.')
            old_id = current['project_id'] or uuid4().hex
            parked = park(store, db, folder, user['id'], old_id)
            folder.mkdir(mode=0o700)
            created = True
            if current['project_id']:
                purge_rows(store, db, current['project_id'])
            db.execute('INSERT INTO projects VALUES (?,?,?,?,?,?,?,0,0,0,?,?)',
                       (project_id, user['id'], name.strip(), description[:4000], color, json.dumps(tags or []), '[]', now, now))
            db.execute('INSERT INTO project_locations VALUES (?,?)', (project_id, str(folder)))
            from .project_state import record
            record(store, db, project_id, user['id'], 'create', 'Projeto criado')
    except Exception:
        # Restore the original first: cleanup failures must not strand its data.
        failed = None
        if created and folder.exists():
            failed = folder.parent / ('.biomol-delete-' + uuid4().hex)
            folder.rename(failed)
        if parked and parked.exists():
            parked.rename(folder)
        if failed:
            try:
                cleanup_id = uuid4().hex
                with store.connect() as db:
                    db.execute('INSERT INTO pending_deletions VALUES (?,?,?,?,?,?)',
                               (cleanup_id, user['id'], str(folder), str(failed), time.time(), '[]'))
                finish(store, cleanup_id)
            except (OSError, sqlite3.Error):
                import logging
                logging.getLogger(__name__).exception('Limpeza pendente da criação incompleta: %s', failed)
        raise
    if parked:
        try:
            finish(store, old_id, token)
        except OSError as error:
            raise ValueError('O novo projeto foi criado, mas a limpeza dos dados anteriores não terminou. '
                             'Verifique as permissões de '+str(parked)+'. A limpeza será retomada ao reiniciar.') from error
    return store.project(token, project_id)


def delete(store, token, project_id, expected_directory=None):
    user = store.user(token)
    parked = None
    folder = None
    try:
        with store.connect() as db:
            db.execute('BEGIN IMMEDIATE')
            row = db.execute('SELECT * FROM projects WHERE id=?', (project_id,)).fetchone()
            if row is None:
                pending = db.execute('SELECT owner FROM pending_deletions WHERE id=?', (project_id,)).fetchone()
                if not pending or pending['owner'] != user['id']:
                    raise AccessDenied('Projeto não encontrado ou acesso não autorizado.')
            else:
                if row['owner'] != user['id']:
                    raise AccessDenied('Somente o proprietário pode excluir o projeto.')
                folder = store.project_dir(project_id, db)
                if expected_directory is not None and str(folder) != expected_directory:
                    raise ValueError('A pasta do projeto mudou. Reabra a confirmação de exclusão.')
                check_path(store, db, folder, user['id'], replacing=True, deleting=True)
                parked = park(store, db, folder, user['id'], project_id)
                purge_rows(store, db, project_id)
    except Exception:
        if parked and parked.exists():
            parked.rename(folder)
        raise
    try:
        finish(store, project_id, token)
    except OSError as error:
        raise ValueError('Projeto removido do workspace, mas a remoção física não terminou. '
                         'Verifique as permissões da pasta e tente excluir novamente. '
                         'A limpeza pendente também será retomada ao reiniciar.') from error
