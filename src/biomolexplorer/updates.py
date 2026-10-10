"""Check master without modifying the running installation or research data."""
import json
import re
import subprocess
import threading
import time
from pathlib import Path
from tempfile import NamedTemporaryFile
from urllib.request import Request, urlopen

REPOSITORY = 'https://github.com/mpiress/BioMolExplorer'
API = 'https://api.github.com/repos/mpiress/BioMolExplorer/commits/master'
CHECK_INTERVAL = 3600
REMIND_INTERVAL = 24 * 3600


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
