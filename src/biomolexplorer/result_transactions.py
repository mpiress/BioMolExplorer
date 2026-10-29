"""Recover file moves after process termination at either side of a DB commit."""
import json
import os
from pathlib import Path
import shutil
import sqlite3
from tempfile import NamedTemporaryFile
from uuid import uuid4


def sync(path):
    if os.name!='posix':return
    fd=os.open(path,os.O_DIRECTORY)
    try:os.fsync(fd)
    finally:os.close(fd)


class ResultTransaction:
    def __init__(self,store,db,pid):
        self.store,self.pid=store,pid
        self.root=store.project_dir(pid,db).resolve();self.id=uuid4().hex;self.moves=[]
        self.planned=set()
        self.journal=self.root/'.history/.result-transactions'/(self.id+'.json')
        db.execute('INSERT INTO result_file_transactions VALUES (?,?)',(self.id,pid))

    def prepare(self,moves):
        for source,destination in moves:
            source,destination=Path(source).resolve(),Path(destination).resolve()
            if not source.is_relative_to(self.root) or source==self.root or source.is_relative_to(self.root/'.history') or not destination.is_relative_to(self.root/'.history'):
                raise ValueError('Arquivo fora do projeto.')
            pair=(str(source),str(destination))
            if pair not in self.planned:self.moves.append(pair);self.planned.add(pair)
        self.journal.parent.mkdir(parents=True,exist_ok=True)
        with NamedTemporaryFile(dir=self.journal.parent,mode='w',encoding='utf-8',delete=False) as stream:
            temp=Path(stream.name);json.dump(dict(id=self.id,project_id=self.pid,moves=self.moves),stream)
            stream.flush();os.fsync(stream.fileno())
        temp.replace(self.journal);sync(self.journal.parent);sync(self.journal.parent.parent);sync(self.root)

    def move(self,source,destination):
        source,destination=Path(source).resolve(),Path(destination).resolve()
        if (str(source),str(destination)) not in self.planned:self.prepare([(source,destination)])
        source.replace(destination)

    def flush(self):
        for parent in {Path(path).parent for pair in self.moves for path in pair}:
            if parent.is_dir():sync(parent)

    def finish(self):
        # Retain the journal if filesystem cleanup needs another attempt.
        if any(Path(saved).exists() for _,saved in self.moves):return
        try:
            self.journal.unlink(missing_ok=True)
            if self.journal.parent.is_dir():sync(self.journal.parent)
            with self.store.connect() as db:db.execute('DELETE FROM result_file_transactions WHERE id=?',(self.id,))
        except (OSError,sqlite3.Error):
            from .diagnostics import get_logger
            get_logger('backend').exception('Result cleanup journal retained for recovery')


def recover(store):
    """The write lock distinguishes a crashed transaction from a live writer."""
    with store.connect() as db:
        db.execute('BEGIN IMMEDIATE')
        for row in db.execute('SELECT id FROM projects WHERE deleted=0').fetchall():
            pid=row[0];root=store.project_dir(pid,db).resolve();folder=root/'.history/.result-transactions'
            for journal in folder.glob('*.json') if folder.is_dir() else []:
                value=json.loads(journal.read_text())
                if value.get('id')!=journal.stem or value.get('project_id')!=pid:continue
                committed=db.execute('SELECT 1 FROM result_file_transactions WHERE id=? AND project_id=?',(journal.stem,pid)).fetchone()
                moves=[]
                for original,saved in value.get('moves',[]):
                    original,saved=Path(original).resolve(),Path(saved).resolve()
                    if not original.is_relative_to(root) or original==root or original.is_relative_to(root/'.history') or not saved.is_relative_to(root/'.history') or saved==root/'.history':
                        raise ValueError('Arquivo fora do projeto.')
                    moves.append((original,saved))
                for original,saved in reversed(moves):
                    if not saved.exists():continue
                    if committed:
                        if saved.is_dir():shutil.rmtree(saved)
                        else:saved.unlink()
                    else:
                        if original.exists():raise ValueError('Conflito ao recuperar arquivos de resultados.')
                        original.parent.mkdir(parents=True,exist_ok=True);saved.replace(original);sync(original.parent)
                # project.json is a portable mirror; SQLite remains authoritative
                # if termination interrupted its publication around the commit.
                from .project_state import mirror
                mirror(store,db,pid)
                journal.unlink();sync(folder)
                db.execute('DELETE FROM result_file_transactions WHERE id=?',(journal.stem,))
            # A committed metadata-only cleanup has no file journal.
            pending={p.stem for p in folder.glob('*.json')} if folder.is_dir() else set()
            for marker in db.execute('SELECT id FROM result_file_transactions WHERE project_id=?',(pid,)).fetchall():
                if marker[0] not in pending:db.execute('DELETE FROM result_file_transactions WHERE id=?',(marker[0],))
