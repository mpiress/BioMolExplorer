"""Verified project archives contain experiments and versions, never credentials."""
import json
import os
import shutil
import stat
import zipfile
from pathlib import Path
from uuid import uuid4

from .project_state import digest, record, restore_database, snapshot

MAX_BYTES = 4 * 1024 ** 3
MAX_FILES = 100000


def export_project(store, token, project_id):
    store.project(token, project_id)
    root = store.project_dir(project_id)
    folder = root / '.exports'
    folder.mkdir(exist_ok=True)
    target = folder / f'project-{project_id}-{uuid4().hex[:8]}.bme.zip'
    with store.connect() as db:
        db.execute('BEGIN IMMEDIATE')
        store.project(token, project_id)
        if db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')", (project_id,)).fetchone():
            raise ValueError('Aguarde a conclusão ou cancele a execução antes de exportar o projeto.')
        state = snapshot(store, db, project_id)
        history = [dict(row) for row in db.execute('SELECT * FROM project_history WHERE project_id=? ORDER BY created', (project_id,))]
        files = {name: info for name, info in state['files'].items()}
        for path in (root / '.history').rglob('*'):
            if path.is_file() and not path.is_symlink() and not any(p.startswith('restore-') for p in path.relative_to(root).parts):
                files[path.relative_to(root).as_posix()] = {'sha256': digest(path), 'size': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns}
        manifest = {'format': 'biomolexplorer-project', 'version': 1, 'state': state, 'history': history, 'files': files}
        # Compress immutable version contents after releasing SQLite's write
        # lock, so exporting a large project does not stop other running jobs.
        originals={name:root/'.history'/'blobs'/info['sha256'] for name,info in state['files'].items()}
        originals.update({name:root/name for name in files if name not in originals})
    try:
        with zipfile.ZipFile(target, 'w', zipfile.ZIP_DEFLATED) as archive:
            archive.writestr('manifest.json', json.dumps(manifest, ensure_ascii=False))
            for name in files:
                archive.write(originals[name], 'data/' + name)
    except BaseException:
        target.unlink(missing_ok=True)
        raise
    return target


def _relative(name):
    path = Path(name)
    if not name or '\\' in name or '\x00' in name or path.is_absolute() or '..' in path.parts or path.as_posix() != name:
        raise ValueError('O projeto exportado contém um caminho inválido.')
    return path


def _validate_state(state):
    from .pipeline import validate_pipeline
    from .workspace import WorkspaceStore
    from .graph_contract import normalize_legacy_graphs
    WorkspaceStore._validate_metadata(state['project']['name'],state['project']['color'],state['project']['tags'])
    state['project']['pipeline']=normalize_legacy_graphs(state['project']['pipeline'])
    validate_pipeline(state['project']['pipeline'])
    import re
    for table in ('assets','runs'):
        identifiers=[r['id'] for r in state[table]]
        if len(identifiers)!=len(set(identifiers)) or any(not isinstance(i,str) or not re.fullmatch('[a-f0-9]{32}',i) for i in identifiers):
            raise ValueError('Identificadores inválidos no projeto exportado.')
    if any(m['role'] not in ('editor','viewer') or m['accepted'] not in (0,1) for m in state['members']):
        raise ValueError('Permissão de compartilhamento inválida no projeto exportado.')
    if len(state['files'])>MAX_FILES or len(state['directories'])>MAX_FILES:
        raise ValueError('A versão contém mais arquivos ou pastas do que o limite permitido.')
    for name, info in state['files'].items():
        path = _relative(name)
        if path.parts[0] in ('.history', '.exports', '.input-cache', 'project.json'):
            raise ValueError('Arquivo reservado no estado do projeto.')
        if not isinstance(info['sha256'], str) or len(info['sha256']) != 64 or any(c not in '0123456789abcdef' for c in info['sha256']):
            raise ValueError('Hash inválido no projeto exportado.')
    for name in state['directories']:
        path=_relative(name)
        if path.parts[0] in ('.history','.exports','.input-cache','project.json'):
            raise ValueError('Pasta reservada no estado do projeto.')
    for asset in state['assets']:
        if not asset['path'].startswith('@project/'):
            raise ValueError('O arquivo fornecido precisa pertencer ao projeto exportado.')
        _relative(asset['path'][9:])
    for run in state['runs']:
        if run['status'] in ('queued', 'running', 'awaiting_input'):
            raise ValueError('O arquivo exportado contém uma execução ainda ativa.')
        for stage in run['stages']:
            for path in stage.get('artifacts', []):
                if not path.startswith('@project/'):
                    raise ValueError('Resultado externo ao projeto exportado.')
                _relative(path[9:])


def import_project(store, token, archive_path, directory):
    user = store.user(token)
    project = None
    try:
        with zipfile.ZipFile(archive_path) as archive:
            infos = archive.infolist()
            names = [info.filename for info in infos]
            if len(names) != len(set(names)) or len(names) > MAX_FILES or sum(i.file_size for i in infos) > MAX_BYTES:
                raise ValueError('Arquivo exportado duplicado ou acima do limite de 4 GB / 100 mil arquivos.')
            for info in infos:
                _relative(info.filename)
                if stat.S_ISLNK(info.external_attr >> 16) or info.flag_bits & 1:
                    raise ValueError('Arquivos simbólicos ou criptografados não são aceitos.')
            manifest_info = archive.getinfo('manifest.json')
            if manifest_info.file_size > 64 * 1024 * 1024:
                raise ValueError('Manifesto acima do limite de 64 MB.')
            manifest = json.loads(archive.read('manifest.json'))
            if manifest.get('format') != 'biomolexplorer-project' or manifest.get('version') != 1:
                raise ValueError('Selecione um projeto exportado pela BioMolExplorer (.bme.zip).')
            state = manifest['state']
            _validate_state(state)
            files = manifest['files']
            if set(names) != {'manifest.json'} | {'data/' + name for name in files}:
                raise ValueError('A lista de arquivos não corresponde ao manifesto do projeto.')
            for name in files:
                path = _relative(name)
                if path.parts[0] in ('.exports', '.input-cache', 'project.json'):
                    raise ValueError('O projeto exportado contém uma pasta reservada.')
                if name.startswith('.history/blobs/') and (len(path.parts)!=3 or path.name!=files[name]['sha256']):
                    raise ValueError('Um conteúdo versionado não corresponde ao seu hash.')
            def check_version(version):
                _validate_state(version)
                for name,info in version['files'].items():
                    blob=files.get('.history/blobs/'+info['sha256'])
                    if blob is None or blob['sha256']!=info['sha256'] or blob['size']!=info['size']:
                        raise ValueError('Um arquivo necessário para restaurar uma versão está incompleto.')
            check_version(state)
            if any(name not in files or files[name]['sha256']!=info['sha256'] or files[name]['size']!=info['size'] for name,info in state['files'].items()):
                raise ValueError('Os resultados não correspondem ao estado atual do projeto.')
            # Validate all bytes and history before creating a destination project.
            import hashlib
            for name, info in files.items():
                hasher, size = hashlib.sha256(), 0
                with archive.open('data/' + name) as stream:
                    for block in iter(lambda: stream.read(1024 * 1024), b''):
                        hasher.update(block)
                        size += len(block)
                if hasher.hexdigest() != info['sha256'] or size != info['size']:
                    raise ValueError(f'O arquivo {name} foi alterado ou está incompleto.')
                if name.startswith('.history/snapshots/'):
                    check_version(json.loads(archive.read('data/' + name)))
            project = store.create_project(token, state['project']['name'], state['project']['description'],
                                           state['project']['color'], state['project']['tags'], directory=directory)
            root = store.project_dir(project['id'])
            with store.connect() as db:
                db.execute('BEGIN IMMEDIATE')
                before = snapshot(store, db, project['id'])
                mapping = {}
                for table in ('assets', 'runs', 'edits'):
                    for row in state[table]:
                        sql_table='compound_edits' if table=='edits' else table
                        if db.execute('SELECT 1 FROM sqlite_master WHERE name=?',(sql_table,)).fetchone() and db.execute(f'SELECT 1 FROM {sql_table} WHERE id=?', (row['id'],)).fetchone():
                            mapping[row['id']] = uuid4().hex
                for event in manifest['history']:
                    mapping[event['id']] = uuid4().hex

                def remap(value):
                    if isinstance(value, str):
                        if value in mapping:
                            return mapping[value]
                        value='/'.join(mapping.get(part,part) for part in value.split('/'))
                        for old, new in mapping.items():
                            value = value.replace('/' + old + '/', '/' + new + '/').replace('/' + old + '-', '/' + new + '-')
                        return value
                    if isinstance(value, list):
                        return [remap(v) for v in value]
                    if isinstance(value, dict):
                        return {remap(k): remap(v) for k, v in value.items()}
                    return value

                state = remap(state)
                for name, info in files.items():
                    destination = root / remap(name)
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    if name.startswith('.history/snapshots/'):
                        version=remap(json.loads(archive.read('data/' + name)))
                        for member in version['members']:member['accepted']=0
                        for run in version['runs']:
                            for item in run['stages']:item.pop('cache_key',None)
                        destination.write_text(json.dumps(version), encoding='utf-8')
                    else:
                        with archive.open('data/' + name) as source, destination.open('wb') as out:
                            shutil.copyfileobj(source, out)
                    os.utime(destination, ns=(info['mtime_ns'], info['mtime_ns']))
                for name in state['directories']:
                    (root / name).mkdir(parents=True, exist_ok=True)
                for run in state['runs']:
                    for stage in run['stages']:
                        stage.pop('cache_key', None)  # Adopt using rebased configuration and verified input/output hashes.
                restore_database(store, db, project['id'], state, owner=user['id'])
                for event in manifest['history']:
                    db.execute('INSERT INTO project_history VALUES (?,?,?,?,?,?,?,?,?)',
                               (mapping[event['id']], project['id'], None, event['actor_name'], event['actor_email'],
                                event['action'], event['summary'], event['created'], event['has_before']))
                record(store, db, project['id'], user['id'], 'import', 'Projeto importado com arquivos, resultados e histórico', before)
            return store.project(token, project['id'])
    except (zipfile.BadZipFile, KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError('Arquivo de projeto inválido ou incompleto. Use um .bme.zip exportado pela plataforma.') from exc
    finally:
        # A failed import must not leave a partial experiment in the workspace.
        if project:
            with store.connect() as db:
                imported = db.execute("SELECT 1 FROM project_history WHERE project_id=? AND action='import'", (project['id'],)).fetchone()
                if not imported:
                    for table in ('project_history', 'members', 'assets', 'runs', 'uploads', 'project_locations', 'projects'):
                        column = 'id' if table == 'projects' else 'project_id'
                        db.execute(f'DELETE FROM {table} WHERE {column}=?', (project['id'],))
                    shutil.rmtree(Path(project['directory']), ignore_errors=True)
