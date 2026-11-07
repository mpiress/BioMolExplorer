"""Check master without modifying the running installation or research data."""
import json
import re
import subprocess
import sys
import threading
import time
import zipfile
from pathlib import Path
from tempfile import NamedTemporaryFile
from urllib.request import Request, urlopen
from urllib.parse import urlparse

REPOSITORY = 'https://github.com/mpiress/BioMolExplorer'
API = 'https://api.github.com/repos/mpiress/BioMolExplorer/commits/master'
CHECK_INTERVAL = 3600
REMIND_INTERVAL = 24 * 3600
MAX_ARCHIVE_BYTES = 256 * 1024 * 1024
MAX_EXPANDED_BYTES = 1024 * 1024 * 1024


def trusted_download_url(url, sha):
    """Accept only the immutable GitHub archive generated for the checked SHA."""
    parsed = urlparse(str(url))
    expected = f'/mpiress/BioMolExplorer/archive/{sha}.zip'
    return parsed.scheme == 'https' and parsed.netloc == 'github.com' and parsed.path == expected \
        and not parsed.params and not parsed.query and not parsed.fragment


def prepare_install_archive(path, destination, sha):
    """Validate GitHub's tree and repack only files required by Python packaging."""
    prefix=f'BioMolExplorer-{sha}/'
    roots=('pyproject.toml','README.md','LICENSE','LICENSE.txt')
    with zipfile.ZipFile(path) as archive, zipfile.ZipFile(destination,'w',zipfile.ZIP_DEFLATED) as output:
        members=archive.infolist()
        names={item.filename for item in members}
        if not members or sum(item.file_size for item in members)>MAX_EXPANDED_BYTES:
            raise ValueError('Pacote de atualização vazio ou excessivamente grande.')
        if prefix+'pyproject.toml' not in names or prefix+'src/biomolexplorer/__init__.py' not in names:
            raise ValueError('O pacote baixado não contém uma distribuição válida do BioMolExplorer.')
        for item in members:
            name=item.filename
            if not name.startswith(prefix) or name.startswith('/') or '..' in Path(name).parts:
                raise ValueError('O pacote de atualização contém caminhos inseguros.')
            relative=name[len(prefix):]
            selected=relative.startswith('src/') or relative in roots
            if not selected or name.endswith('/'):continue
            if (item.external_attr >> 16) & 0o170000 == 0o120000:
                raise ValueError('O pacote de atualização contém links não permitidos.')
            output.writestr(item,archive.read(item))
    return destination


def validate_archive(path, sha):
    """Compatibility probe that validates through the exact installable subset."""
    with NamedTemporaryFile(suffix='.zip') as temporary:
        prepare_install_archive(path,temporary.name,sha)
    return True


def installed_revision():
    root = Path(__file__).resolve().parents[2]
    archive_revision = re.fullmatch(r'BioMolExplorer-([a-f0-9]{40})', root.name)
    if archive_revision:
        return archive_revision[1]
    if (root / '.git').exists():
        try:
            result = subprocess.run(['git', '-C', str(root), 'rev-parse', 'HEAD'],
                capture_output=True, text=True, timeout=5, check=True)
            return result.stdout.strip()
        except (OSError, subprocess.SubprocessError):
            pass
    revision = Path(__file__).parent / 'resources' / 'revision.txt'
    if revision.is_file():
        value = revision.read_text().strip()
        if re.fullmatch('[a-f0-9]{40}', value): return value
    return None


class UpdateChecker:
    def __init__(self, root, revision=None):
        self.state_file = Path(root) / 'updates.json'
        self.revision = revision or installed_revision()
        self.lock = threading.Lock()
        self.checked_at = None
        self.latest = None
        self.error = None

    def _state(self):
        try:
            state = json.loads(self.state_file.read_text())
            return state if isinstance(state, dict) else {}
        except (OSError, ValueError): return {}

    def _save(self, state):
        with NamedTemporaryFile(mode='w', encoding='utf-8', dir=self.state_file.parent,
                prefix='.updates-', suffix='.tmp', delete=False) as stream:
            json.dump(state, stream)
            temporary = Path(stream.name)
        try:
            temporary.replace(self.state_file)
        finally:
            temporary.unlink(missing_ok=True)

    def check(self, force=False):
        with self.lock:
            now = time.time()
            if not force and self.checked_at is not None and now-self.checked_at < CHECK_INTERVAL:
                return self.latest
            self.checked_at = now
            try:
                request = Request(API, headers={'Accept': 'application/vnd.github+json',
                    'User-Agent': 'BioMolExplorer-update-checker', 'X-GitHub-Api-Version': '2022-11-28'})
                with urlopen(request, timeout=10) as response:
                    data = json.loads(response.read(1024 * 1024))
                sha = data['sha']
                if not isinstance(sha, str) or not re.fullmatch('[a-f0-9]{40}', sha):
                    raise ValueError('Revisão inválida retornada pelo GitHub.')
                state = self._state()
                baseline = self.revision or state.get('baseline')
                if not baseline:
                    # Source ZIPs have no Git metadata: track subsequent master changes.
                    state['baseline'] = sha
                    self._save(state)
                    baseline = sha
                available = sha != baseline
                root = Path(__file__).resolve().parents[2]
                if available and (root / '.git').exists():
                    ancestor = subprocess.run(['git', '-C', str(root), 'merge-base', '--is-ancestor', sha, baseline],
                        capture_output=True, timeout=5)
                    if ancestor.returncode == 0: available = False
                self.latest = {'sha': sha, 'available': available,
                    'message': str(data.get('commit', {}).get('message', '')).split('\n')[0][:240],
                    'download_url': f'{REPOSITORY}/archive/{sha}.zip'}
                if not trusted_download_url(self.latest['download_url'], sha):
                    raise ValueError('URL de atualização não confiável.')
                self.error = None
            except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as exc:
                self.error = str(exc)
                # A failed request cannot assert that the installation is current.
                self.latest = None
            return self.latest

    def visible(self, user='installation'):
        # Rendering must never wait for the network worker's check lock.
        latest = self.latest
        until = self._state().get('reminders', {}).get(user, 0)
        return bool(latest and latest['available'] and time.time() >= until)

    def remind_later(self, user='installation'):
        with self.lock:
            state = self._state()
            state.setdefault('reminders', {})[user] = time.time() + REMIND_INTERVAL
            self._save(state)

    def install(self, latest=None):
        """Download, verify and atomically delegate replacement to this Python environment."""
        with self.lock:
            latest=latest or self.latest
            if not latest or not latest.get('available'):
                raise ValueError('Nenhuma atualização válida está disponível para instalação.')
            sha=latest.get('sha');url=latest.get('download_url')
            if not isinstance(sha,str) or not re.fullmatch('[a-f0-9]{40}',sha) or not trusted_download_url(url,sha):
                raise ValueError('URL de atualização não confiável.')
            directory=self.state_file.parent/'updates';directory.mkdir(parents=True,exist_ok=True)
            archive_path=directory/f'BioMolExplorer-{sha}.zip'
            install_path=directory/f'BioMolExplorer-{sha}-install.zip'
            temporary=None
            try:
                request=Request(url,headers={'User-Agent':'BioMolExplorer-updater'})
                with urlopen(request,timeout=30) as response, NamedTemporaryFile(dir=directory,
                        prefix='.download-',suffix='.zip',delete=False) as stream:
                    temporary=Path(stream.name);total=0
                    while True:
                        chunk=response.read(1024*1024)
                        if not chunk:break
                        total+=len(chunk)
                        if total>MAX_ARCHIVE_BYTES:raise ValueError('Pacote de atualização excessivamente grande.')
                        stream.write(chunk)
                prepare_install_archive(temporary,install_path,sha)
                temporary.replace(archive_path);temporary=None
                command=[sys.executable,'-m','pip','install','--disable-pip-version-check','--no-deps','--upgrade',str(install_path)]
                result=subprocess.run(command,capture_output=True,text=True,timeout=900,check=False)
                log=directory/f'install-{sha}.log'
                log.write_text('$ '+' '.join(command)+'\n\n'+(result.stdout or '')+'\n'+(result.stderr or ''),encoding='utf-8')
                if result.returncode:
                    raise RuntimeError('A instalação terminou com erro. Consulte '+str(log))
                state=self._state();state.update(baseline=sha,installed=sha,installed_at=time.time(),install_log=str(log))
                self._save(state);self.revision=sha
                self.latest=dict(latest,available=False)
                return {'sha':sha,'archive':str(install_path),'source_archive':str(archive_path),
                    'log':str(log),'restart_required':True}
            finally:
                if temporary is not None:temporary.unlink(missing_ok=True)
