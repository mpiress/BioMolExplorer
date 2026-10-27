"""Project-scoped relational index of verified scientific checkpoints.

Objects are content addressed, shared across experiments and never hard-linked
to mutable output files. Publication is transactional and occurs after validation.
"""
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import tempfile

from .stage_cache import file_digest


def identity(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def copy_artifact(source, destination):
    """Use copy-on-write where supported, otherwise an independent copy."""
    import fcntl
    source, destination = Path(source), Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as stream:
        temporary = Path(stream.name)
    try:
        try:
            with source.open('rb') as src, temporary.open('wb') as dst:
                fcntl.ioctl(dst.fileno(), 0x40049409, src.fileno())  # FICLONE
        except OSError:
            shutil.copyfile(source, temporary)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


class CalculationCache:
    def __init__(self, namespace, destination):
        preparation = os.environ.get('BIOMOL_PREPARATION_CACHE')
        self.root = Path(preparation).parent / 'calculations' if preparation else Path(destination) / '.calculation-cache'
        self.root.mkdir(parents=True, exist_ok=True)
        self.namespace = namespace
        self.reuse = os.environ.get('BIOMOL_REUSE_PREPARATIONS', '1') != '0'
        self.database = self.root / 'index.sqlite3'
        with self.connect() as db:
            db.executescript('''
                CREATE TABLE IF NOT EXISTS experiments (
                    namespace TEXT NOT NULL, key TEXT NOT NULL, metadata TEXT NOT NULL,
                    PRIMARY KEY(namespace,key));
                CREATE TABLE IF NOT EXISTS objects (
                    sha256 TEXT PRIMARY KEY, size INTEGER NOT NULL);
                CREATE TABLE IF NOT EXISTS artifacts (
                    namespace TEXT NOT NULL, experiment_key TEXT NOT NULL,
                    role TEXT NOT NULL, sha256 TEXT NOT NULL REFERENCES objects(sha256),
                    PRIMARY KEY(namespace,experiment_key,role),
                    FOREIGN KEY(namespace,experiment_key) REFERENCES experiments(namespace,key));
                CREATE INDEX IF NOT EXISTS artifacts_object ON artifacts(sha256);
                CREATE TABLE IF NOT EXISTS subjects (
                    namespace TEXT NOT NULL, experiment_key TEXT NOT NULL,
                    role TEXT NOT NULL, identifier TEXT NOT NULL, input_sha256 TEXT,
                    PRIMARY KEY(namespace,experiment_key,role),
                    FOREIGN KEY(namespace,experiment_key) REFERENCES experiments(namespace,key));
                CREATE INDEX IF NOT EXISTS subjects_identifier ON subjects(role,identifier);
                CREATE INDEX IF NOT EXISTS subjects_content ON subjects(input_sha256);
            ''')

    @contextmanager
    def connect(self):
        db = sqlite3.connect(self.database, timeout=60)
        db.execute('PRAGMA foreign_keys=ON')
        try:
            with db:
                yield db
        finally:
            db.close()

    @contextmanager
    def locked(self, key):
        import fcntl
        # Independent scientific layers acquire nested locks (structure ->
        # charges/topology). Sharing stripes across namespaces can deadlock
        # even within one worker when two keys land on the same inode.
        locks = self.root / 'locks' / identity(self.namespace)[:16]
        locks.mkdir(parents=True,exist_ok=True)
        stripe=str(int(identity([self.namespace,key])[:2],16) % 64)+'.lock'
        with (locks / stripe).open('a') as stream:
            fcntl.flock(stream, fcntl.LOCK_EX)
            yield

    def object_path(self, digest):
        return self.root / 'objects' / digest[:2] / digest

    def metadata(self, key):
        if not self.reuse:
            return None
        with self.connect() as db:
            row = db.execute('SELECT metadata FROM experiments WHERE namespace=? AND key=?',
                             (self.namespace, key)).fetchone()
        try:
            metadata=json.loads(row[0]) if row else None
            return metadata if isinstance(metadata,dict) else None
        except (TypeError,ValueError):return None

    def restore(self, key, targets, validate=None):
        if not self.reuse:
            return None
        with self.connect() as db:
            db.execute('BEGIN')
            record = db.execute('SELECT metadata FROM experiments WHERE namespace=? AND key=?',
                                (self.namespace, key)).fetchone()
            rows = db.execute('SELECT role,sha256 FROM artifacts WHERE namespace=? AND experiment_key=?',
                              (self.namespace, key)).fetchall()
        if not record or set(targets) != {row[0] for row in rows}:
            return None
        try:
            metadata=json.loads(record[0])
            if not isinstance(metadata,dict):return None
            # Validate every object before materializing any output.
            for role, digest in rows:
                source = self.object_path(digest)
                if file_digest(source) != digest:
                    return None
                if validate:
                    validate(role, source)
            for role, digest in rows:
                copy_artifact(self.object_path(digest), targets[role])
            from .diagnostics import event,get_logger
            event(get_logger('scientific_tools'),'calculation.reused','Checkpoint científico validado e reutilizado.',
                  tool=self.namespace,configuration=key,
                  compound_identifier=metadata.get('compound_identifier') or metadata.get('compound'))
            return metadata
        except (OSError, ValueError, KeyError):
            return None

    def remember(self, key, files, metadata=None, validate=None):
        from .diagnostics import current_context
        metadata=dict(metadata or {})
        metadata.setdefault('origin',current_context())
        entries = []
        for role, source in files.items():
            source = Path(source)
            if not source.is_file() or not source.stat().st_size:
                raise ValueError('Cache output is absent or empty: ' + str(source))
            if validate:
                validate(role, source)
            digest = file_digest(source); target = self.object_path(digest)
            if not target.is_file() or file_digest(target) != digest:
                copy_artifact(source, target)
            entries.append((role, digest, source.stat().st_size))
        with self.connect() as db:
            db.execute('INSERT OR REPLACE INTO experiments VALUES (?,?,?)',
                       (self.namespace, key, json.dumps(metadata, allow_nan=False)))
            db.execute('DELETE FROM artifacts WHERE namespace=? AND experiment_key=?', (self.namespace, key))
            db.execute('DELETE FROM subjects WHERE namespace=? AND experiment_key=?', (self.namespace, key))
            compound=metadata.get('compound_identifier') or metadata.get('compound')
            recipe=metadata.get('recipe',{})
            record=metadata.get('record')
            subjects={}
            if compound:subjects['ligand']=(str(compound),recipe.get('ligand'))
            if metadata.get('receptor_identifier'):
                subjects['receptor']=(str(metadata['receptor_identifier']),recipe.get('receptor'))
            if record:
                subjects['receptor']=(str(record[0])+'_'+str(record[3]),recipe.get('receptor'))
                subjects['reference']=('|'.join(map(str,record[:4])),None)
            for role in ('receptor','ligand'):
                if role in recipe and role not in subjects:
                    subjects[role]=(str(recipe[role]),str(recipe[role]))
            for role,(identifier,digest) in subjects.items():
                db.execute('INSERT INTO subjects VALUES (?,?,?,?,?)',(self.namespace,key,role,identifier,digest))
            for role, digest, size in entries:
                db.execute('INSERT OR IGNORE INTO objects VALUES (?,?)', (digest, size))
                db.execute('INSERT INTO artifacts VALUES (?,?,?,?)', (self.namespace, key, role, digest))
