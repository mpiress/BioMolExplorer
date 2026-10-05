"""Authenticated project storage; authorization is enforced on every public call."""
import hashlib
import hmac
import json
import re
import secrets
import shutil
import sqlite3
import time
from contextlib import contextmanager
from pathlib import Path
from uuid import uuid4

ROLES = {'viewer': 0, 'editor': 1, 'owner': 2}
COLORS = ['#14B8A6', '#6366F1', '#F59E0B', '#F43F5E', '#3B82F6', '#8B5CF6']


class AccessDenied(PermissionError):
    pass


class WorkspaceStore:
    def __init__(self, root, max_upload_bytes=200 * 1024 * 1024):
        self.root = Path(root).expanduser().resolve()
        self.root.mkdir(parents=True, exist_ok=True, mode=0o700)
        self.projects_root = self.root / 'projects'
        self.staging = self.root / 'staging'
        self.projects_root.mkdir(exist_ok=True, mode=0o700)
        self.staging.mkdir(exist_ok=True, mode=0o700)
        self.database = self.root / 'workspace.sqlite3'
        self.max_upload_bytes = max_upload_bytes
        with self.connect() as db:
            db.executescript('''
                CREATE TABLE IF NOT EXISTS users(id TEXT PRIMARY KEY, email TEXT UNIQUE NOT NULL,
                    name TEXT NOT NULL, salt TEXT NOT NULL, password TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS sessions(token TEXT PRIMARY KEY, user_id TEXT REFERENCES users(id), expires REAL NOT NULL);
                CREATE TABLE IF NOT EXISTS projects(id TEXT PRIMARY KEY, owner TEXT REFERENCES users(id), name TEXT NOT NULL,
                    description TEXT NOT NULL, color TEXT NOT NULL, tags TEXT NOT NULL, pipeline TEXT NOT NULL,
                    revision INTEGER NOT NULL DEFAULT 0, archived INTEGER NOT NULL DEFAULT 0, deleted INTEGER NOT NULL DEFAULT 0,
                    created REAL NOT NULL, updated REAL NOT NULL);
                CREATE TABLE IF NOT EXISTS members(project_id TEXT REFERENCES projects(id), user_id TEXT REFERENCES users(id),
                    role TEXT NOT NULL, accepted INTEGER NOT NULL DEFAULT 0, PRIMARY KEY(project_id,user_id));
                CREATE TABLE IF NOT EXISTS assets(id TEXT PRIMARY KEY, project_id TEXT REFERENCES projects(id),
                    name TEXT NOT NULL, path TEXT NOT NULL, size INTEGER NOT NULL, kind TEXT NOT NULL, created REAL NOT NULL);
                CREATE TABLE IF NOT EXISTS uploads(id TEXT PRIMARY KEY, project_id TEXT REFERENCES projects(id),
                    user_id TEXT REFERENCES users(id), name TEXT NOT NULL, kind TEXT NOT NULL, expires REAL NOT NULL);
                CREATE TABLE IF NOT EXISTS runs(id TEXT PRIMARY KEY, project_id TEXT REFERENCES projects(id),
                    user_id TEXT REFERENCES users(id), status TEXT NOT NULL, stages TEXT NOT NULL,
                    created REAL NOT NULL, updated REAL NOT NULL, error TEXT);
                CREATE INDEX IF NOT EXISTS runs_project ON runs(project_id,created);
                CREATE INDEX IF NOT EXISTS assets_project ON assets(project_id);
                CREATE TABLE IF NOT EXISTS project_locations(project_id TEXT PRIMARY KEY REFERENCES projects(id), path TEXT UNIQUE NOT NULL);
                CREATE TABLE IF NOT EXISTS project_history(id TEXT PRIMARY KEY, project_id TEXT REFERENCES projects(id),
                    user_id TEXT, actor_name TEXT NOT NULL, actor_email TEXT NOT NULL, action TEXT NOT NULL,
                    summary TEXT NOT NULL, created REAL NOT NULL, has_before INTEGER NOT NULL);
                CREATE INDEX IF NOT EXISTS history_project ON project_history(project_id,created);
            ''')

    @contextmanager
    def connect(self):
        db = sqlite3.connect(self.database, timeout=30)
        db.row_factory = sqlite3.Row
        db.execute('PRAGMA foreign_keys=ON')
        try:
            with db:
                yield db
        finally:
            db.close()

    @staticmethod
    def _password(password, salt):
        return hashlib.scrypt(password.encode(), salt=bytes.fromhex(salt), n=16384, r=8, p=1).hex()

    def register(self, name, email, password):
        name, email = name.strip(), email.strip().lower()
        if not name or len(name) > 100 or not re.fullmatch(r'[^\s@]+@[^\s@]+\.[^\s@]+', email):
            raise ValueError('Informe seu nome e um e-mail válido.')
        if not isinstance(password, str) or not 10 <= len(password) <= 1024:
            raise ValueError('A senha deve ter entre 10 e 1024 caracteres.')
        salt, user_id = secrets.token_hex(16), uuid4().hex
        hashed = self._password(password, salt)
        try:
            with self.connect() as db:
                db.execute('INSERT INTO users VALUES (?,?,?,?,?)', (user_id, email, name, salt, hashed))
        except sqlite3.IntegrityError:
            raise ValueError('Este e-mail já possui uma conta.') from None
        return self.login(email, password)

    def login(self, email, password):
        if not isinstance(password, str) or len(password) > 1024:
            raise AccessDenied('E-mail ou senha inválidos.')
        with self.connect() as db:
            user = db.execute('SELECT * FROM users WHERE email=?', (email.strip().lower(),)).fetchone()
        # Derive even for unknown accounts to avoid a cheap account timing oracle.
        derived = self._password(password, user['salt'] if user else '00' * 16)
        if not user or not hmac.compare_digest(derived, user['password']):
            raise AccessDenied('E-mail ou senha inválidos.')
        token = secrets.token_urlsafe(32)
        with self.connect() as db:
            db.execute('DELETE FROM sessions WHERE expires<?', (time.time(),))
            db.execute('INSERT INTO sessions VALUES (?,?,?)', (self._digest(token), user['id'], time.time() + 8 * 3600))
        return token

    @staticmethod
    def _digest(token):
        if not isinstance(token,str) or not token:
            raise AccessDenied('Entre na aplicação para continuar.')
        return hashlib.sha256(token.encode()).hexdigest()

    def user(self, token):
        with self.connect() as db:
            user = db.execute('''SELECT u.id,u.name,u.email FROM users u JOIN sessions s ON u.id=s.user_id
                WHERE s.token=? AND s.expires>?''', (self._digest(token), time.time())).fetchone()
        if user is None:
            raise AccessDenied('Sua sessão expirou. Entre novamente.')
        return dict(user)

    def logout(self, token):
        with self.connect() as db:
            db.execute('DELETE FROM sessions WHERE token=?', (self._digest(token),))

    def _require_user(self, user_id, project_id, minimum='viewer'):
        with self.connect() as db:
            project = db.execute('SELECT * FROM projects WHERE id=? AND deleted=0', (project_id,)).fetchone()
            member = db.execute('SELECT role FROM members WHERE project_id=? AND user_id=? AND accepted=1',
                                (project_id, user_id)).fetchone()
        role = 'owner' if project and project['owner'] == user_id else member['role'] if member else None
        if project is None or role is None or ROLES[role] < ROLES[minimum]:
            raise AccessDenied('Projeto não encontrado ou acesso não autorizado.')
        result = dict(project)
        from .graph_contract import normalize_legacy_graphs
        result.update(role=role, tags=json.loads(result['tags']), pipeline=normalize_legacy_graphs(json.loads(result['pipeline'])))
        result['directory'] = str(self.project_dir(project_id))
        return result

    def project(self, token, project_id, minimum='viewer'):
        return self._require_user(self.user(token)['id'], project_id, minimum)

    def project_dir(self, project_id, db=None):
        if not re.fullmatch('[a-f0-9]{32}', project_id):
            raise ValueError('Identificador de projeto inválido.')
        if db is None:
            with self.connect() as connection:
                row = connection.execute('SELECT path FROM project_locations WHERE project_id=?', (project_id,)).fetchone()
        else:
            row=db.execute('SELECT path FROM project_locations WHERE project_id=?',(project_id,)).fetchone()
        return Path(row['path']) if row else self.projects_root / project_id

    def create_project(self, token, name, description='', color=COLORS[0], tags=None, directory=None):
        user = self.user(token)
        self._validate_metadata(name, color, tags or [])
        project_id, now = uuid4().hex, time.time()
        if directory is not None and not str(directory).strip():
            raise ValueError('Informe a pasta do projeto.')
        folder = Path(directory).expanduser().absolute() if directory is not None else self.projects_root / project_id
        if folder.resolve() != folder:
            raise ValueError('Escolha uma pasta sem links simbólicos ou segmentos .. no caminho.')
        if folder.exists() and (not folder.is_dir() or any(folder.iterdir())):
            raise ValueError('Escolha uma pasta nova ou vazia para evitar sobrescrever arquivos existentes.')
        with self.connect() as db:
            db.execute('BEGIN IMMEDIATE')
            for row in db.execute('SELECT id FROM projects'):
                other = self.project_dir(row['id'])
                if folder.is_relative_to(other) or other.is_relative_to(folder):
                    raise ValueError('A pasta não pode conter nem pertencer a outro projeto.')
            if self.database.is_relative_to(folder) or folder.is_relative_to(self.staging):
                raise ValueError('Escolha uma pasta exclusiva para os arquivos do projeto.')
            folder.mkdir(parents=True, exist_ok=True, mode=0o700)
            db.execute('INSERT INTO projects VALUES (?,?,?,?,?,?,?,0,0,0,?,?)',
                (project_id, user['id'], name.strip(), description[:4000], color, json.dumps(tags or []), '[]', now, now))
            db.execute('INSERT INTO project_locations VALUES (?,?)', (project_id, str(folder)))
            from .project_state import record
            record(self, db, project_id, user['id'], 'create', 'Projeto criado')
        return self.project(token, project_id)

    @contextmanager
    def change(self, token, project_id, role, action, summary):
        from .project_state import snapshot, record
        user = self.user(token)
        with self.connect() as db:
            db.execute('BEGIN IMMEDIATE')
            self._require_user(user['id'], project_id, role)
            before = snapshot(self, db, project_id)
            yield db
            record(self, db, project_id, user['id'], action, summary, before)

    def history(self, token, project_id):
        self.project(token, project_id)
        with self.connect() as db:
            return [dict(row) for row in db.execute('SELECT * FROM project_history WHERE project_id=? ORDER BY created DESC', (project_id,))]

    def rollback(self, token, project_id, event_id):
        from .project_state import rollback
        return rollback(self, token, project_id, event_id)

    def export_project(self, token, project_id):
        from .project_archive import export_project
        return export_project(self, token, project_id)

    def import_project(self, token, archive, directory):
        from .project_archive import import_project
        return import_project(self, token, archive, directory)

    @staticmethod
    def _validate_metadata(name, color, tags):
        if not isinstance(name, str) or not name.strip() or len(name) > 120:
            raise ValueError('Dê ao projeto um nome de até 120 caracteres.')
        if color not in COLORS or not isinstance(tags, list) or len(tags) > 20 or not all(isinstance(tag, str) and len(tag) <= 40 for tag in tags):
            raise ValueError('Cor ou tags inválidas.')

    def list_projects(self, token, archived=False):
        user = self.user(token)
        with self.connect() as db:
            ids = db.execute('''SELECT DISTINCT p.id FROM projects p LEFT JOIN members m ON m.project_id=p.id
                WHERE p.deleted=0 AND p.archived=? AND (p.owner=? OR (m.user_id=? AND m.accepted=1)) ORDER BY p.updated DESC''',
                (int(archived), user['id'], user['id'])).fetchall()
        return [self._require_user(user['id'], row['id']) for row in ids]

    def update_project(self, token, project_id, name, description, color, tags, archived=False):
        self.project(token, project_id, 'editor')
        self._validate_metadata(name, color, tags)
        with self.change(token, project_id, 'editor', 'metadata', 'Nome, descrição, tags ou cor alterados') as db:
            db.execute('UPDATE projects SET name=?,description=?,color=?,tags=?,archived=?,updated=? WHERE id=?',
                       (name.strip(), description[:4000], color, json.dumps(tags), int(archived), time.time(), project_id))

    def delete_project(self, token, project_id):
        self.project(token, project_id, 'owner')
        with self.change(token, project_id, 'owner', 'delete', 'Projeto excluído do workspace') as db:
            active = db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')", (project_id,)).fetchone()
            if active:
                raise ValueError('Cancele a execução antes de excluir o projeto.')
            db.execute('UPDATE projects SET deleted=1,updated=? WHERE id=?', (time.time(), project_id))
        # Tombstone rather than deleting large scientific datasets during a UI request.

    def save_pipeline(self, token, project_id, stages, expected_revision, base_pipeline=None):
        self.project(token, project_id, 'editor')
        from .pipeline import validate_pipeline
        validate_pipeline(stages)
        with self.change(token, project_id, 'editor', 'pipeline', 'Pipeline e configurações dos blocos alterados') as db:
            current = db.execute('SELECT pipeline,revision FROM projects WHERE id=?', (project_id,)).fetchone()
            if current['revision'] != expected_revision:
                if base_pipeline is None:
                    raise ValueError('O projeto foi alterado por outro usuário. Reabra-o antes de salvar.')
                from .project_merge import merge_pipeline
                from .graph_contract import normalize_legacy_graphs
                stages = merge_pipeline(base_pipeline, stages, normalize_legacy_graphs(json.loads(current['pipeline'])))
                validate_pipeline(stages)
            db.execute('UPDATE projects SET pipeline=?,revision=revision+1,updated=? WHERE id=?',
                (json.dumps(stages, allow_nan=False), time.time(), project_id))
        return self.project(token, project_id)

    def invite(self, token, project_id, email, role='viewer'):
        project = self.project(token, project_id, 'owner')
        if role not in ('viewer', 'editor'):
            raise ValueError('Escolha leitor ou editor.')
        with self.change(token, project_id, 'owner', 'invite', f'Convite para {email.strip().lower()} como {role}') as db:
            recipient = db.execute('SELECT id FROM users WHERE email=?', (email.strip().lower(),)).fetchone()
            if not recipient:
                raise ValueError('O convidado precisa criar uma conta antes de receber o convite.')
            if recipient['id'] == project['owner']:
                raise ValueError('Você já é o proprietário.')
            existing = db.execute('SELECT role,accepted FROM members WHERE project_id=? AND user_id=?',
                                  (project_id,recipient['id'])).fetchone()
            if existing and existing['accepted'] and existing['role'] == role:
                raise ValueError('Este usuário já tem acesso ao projeto com essa permissão.')
            db.execute('INSERT INTO members VALUES (?,?,?,0) ON CONFLICT(project_id,user_id) DO UPDATE SET role=excluded.role,accepted=0',
                       (project_id, recipient['id'], role))

    def invitations(self, token):
        user = self.user(token)
        with self.connect() as db:
            return [dict(row) for row in db.execute('''SELECT p.id,p.name,u.name AS owner_name,m.role FROM members m
                JOIN projects p ON p.id=m.project_id JOIN users u ON u.id=p.owner
                WHERE m.user_id=? AND m.accepted=0 AND p.deleted=0''', (user['id'],))]

    def accept_invitation(self, token, project_id, accept=True, expected_role=None):
        user = self.user(token)
        if type(accept) is not bool or expected_role not in (None,'viewer','editor'):
            raise ValueError('Resposta ou permissão do convite inválida.')
        with self.connect() as db:
            db.execute('BEGIN IMMEDIATE')
            pending = db.execute('''SELECT m.role FROM members m JOIN projects p ON p.id=m.project_id
                WHERE m.project_id=? AND m.user_id=? AND m.accepted=0 AND p.deleted=0''',
                (project_id,user['id'])).fetchone()
            if not pending:
                raise AccessDenied('Convite inexistente, revogado ou já respondido.')
            if expected_role is not None and pending['role'] != expected_role:
                raise ValueError('A permissão do convite mudou. Confira o convite atualizado antes de responder.')
            if accept:
                db.execute('UPDATE members SET accepted=1 WHERE project_id=? AND user_id=? AND accepted=0', (project_id, user['id']))
            else:
                db.execute('DELETE FROM members WHERE project_id=? AND user_id=? AND accepted=0', (project_id, user['id']))
            from .project_state import snapshot, record
            # Acceptance is also an attributable project change. Reconstruct its previous membership.
            before = snapshot(self, db, project_id)
            before['members'] = [m for m in before['members'] if m['user_id'] != user['id']] + [
                {'project_id':project_id,'user_id':user['id'],'email':user['email'],'role':pending['role'],'accepted':0}]
            record(self, db, project_id, user['id'], 'invitation', 'Convite aceito' if accept else 'Convite recusado', before)

    def members(self, token, project_id):
        self.project(token, project_id, 'owner')
        with self.connect() as db:
            return [dict(row) for row in db.execute('''SELECT u.id,u.name,u.email,m.role,m.accepted FROM members m
                JOIN users u ON u.id=m.user_id WHERE m.project_id=?''', (project_id,))]

    def revoke(self, token, project_id, user_id):
        self.project(token, project_id, 'owner')
        with self.change(token, project_id, 'owner', 'revoke', 'Acesso de colaborador revogado') as db:
            db.execute('DELETE FROM members WHERE project_id=? AND user_id=?', (project_id, user_id))

    def prepare_upload(self, token, project_id, name, kind):
        user = self.user(token)
        self._require_user(user['id'], project_id, 'editor')
        if Path(name).name != name or not name or '\\' in name or len(name) > 180 or any(c in name for c in '\r\n\x00'):
            raise ValueError('Nome de arquivo inválido.')
        if kind not in ('compounds', 'structures', 'prepared_structures', 'fingerprints', 'similarity', 'vina', 'dock6', 'scores', 'visualization', 'other'):
            raise ValueError('Tipo de arquivo inválido.')
        ticket = uuid4().hex
        with self.connect() as db:
            db.execute('INSERT INTO uploads VALUES (?,?,?,?,?,?)',
                       (ticket, project_id, user['id'], name, kind, time.time() + 600))
        return ticket

    def finish_upload(self, token, ticket):
        user = self.user(token)
        with self.connect() as db:
            row = db.execute('SELECT * FROM uploads WHERE id=? AND user_id=? AND expires>?',
                             (ticket, user['id'], time.time())).fetchone()
        if row is None:
            raise AccessDenied('Envio expirado ou não autorizado.')
        self._require_user(user['id'], row['project_id'], 'editor')
        source = self.staging / ticket
        if source.is_symlink() or not source.is_file() or source.stat().st_size > self.max_upload_bytes:
            source.unlink(missing_ok=True)
            raise ValueError('Arquivo inexistente ou maior que o limite de envio.')
        asset_id = uuid4().hex
        from .input_validation import validate_file
        # Validate with the original filename: the upload ticket has no extension.
        folder = self.project_dir(row['project_id']) / 'assets' / asset_id
        folder.mkdir(parents=True, mode=0o700)
        target = folder / row['name']
        shutil.move(source, target)
        try:
            validate_file(target, row['kind'])
            with self.change(token, row['project_id'], 'editor', 'upload', f'Arquivo fornecido: {row["name"]}') as db:
                db.execute('INSERT INTO assets VALUES (?,?,?,?,?,?,?)',
                           (asset_id, row['project_id'], row['name'], str(target), target.stat().st_size, row['kind'], time.time()))
                db.execute('DELETE FROM uploads WHERE id=?', (ticket,))
        except BaseException:
            shutil.rmtree(folder, ignore_errors=True)
            raise
        return asset_id

    def import_local_file(self, token, project_id, path, kind='other'):
        # Native desktop picker supplies this path; web clients use signed uploads.
        source = Path(path)
        if not source.is_file() or source.stat().st_size > self.max_upload_bytes:
            raise ValueError('Arquivo inexistente ou maior que o limite de envio.')
        ticket = self.prepare_upload(token, project_id, source.name, kind)
        shutil.copyfile(source, self.staging / ticket)
        return self.finish_upload(token, ticket)

    def assets(self, token, project_id):
        self.project(token, project_id)
        with self.connect() as db:
            return [dict(row) for row in db.execute('SELECT id,name,size,kind,created FROM assets WHERE project_id=? ORDER BY created DESC', (project_id,))]

    def asset_path(self, token, project_id, asset_id):
        self.project(token, project_id)
        with self.connect() as db:
            row = db.execute('SELECT path FROM assets WHERE id=? AND project_id=?', (asset_id, project_id)).fetchone()
        if row is None:
            raise AccessDenied('Arquivo não autorizado.')
        return self.scoped_path(project_id, row['path'])

    def scoped_path(self, project_id, path):
        resolved = Path(path).resolve()
        if not resolved.is_relative_to(self.project_dir(project_id).resolve()):
            raise AccessDenied('O caminho deve pertencer ao projeto.')
        return resolved

    def read_file(self, token, project_id, path, max_bytes=None):
        self.project(token, project_id)
        path = self.scoped_path(project_id, path)
        if not path.is_file():
            raise FileNotFoundError('Arquivo indisponível.')
        with path.open('rb') as stream:
            return stream.read(max_bytes) if max_bytes else stream.read()

    def list_runs(self, token, project_id):
        self.project(token, project_id)
        with self.connect() as db:
            return [self._run_row(row) for row in db.execute('SELECT * FROM runs WHERE project_id=? ORDER BY created DESC LIMIT 100', (project_id,))]

    @staticmethod
    def _run_row(row):
        result = dict(row)
        result['stages'] = json.loads(result['stages'])
        return result

    def get_run(self, token, run_id):
        with self.connect() as db:
            row = db.execute('SELECT * FROM runs WHERE id=?', (run_id,)).fetchone()
        if row is None:
            raise AccessDenied('Execução não encontrada.')
        self.project(token, row['project_id'])
        return self._run_row(row)
