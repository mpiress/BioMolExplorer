"""Optional, atomic worker feedback; scientific functions also work without a UI."""
import json
import os
import threading
import time
from pathlib import Path

_LOCK = threading.Lock()


def report_progress(message, completed=None, total=None):
    filename = os.environ.get('BIOMOL_PROGRESS_FILE')
    if not filename:
        return
    event = {'message': str(message)[:1500], 'updated_at': time.time()}
    if completed is not None and total is not None:
        event.update(completed=int(completed), total=int(total))
    # Monitoring failures must never invalidate scientific output.
    try:
        with _LOCK:
            path = Path(filename)
            temporary = path.with_suffix('.tmp')
            temporary.write_text(json.dumps(event, ensure_ascii=False), encoding='utf-8')
            temporary.replace(path)
    except OSError:
        from .diagnostics import get_logger
        get_logger('backend').warning('Não foi possível atualizar o acompanhamento do worker.', exc_info=True)


def read_progress(path):
    try:
        with Path(path).open(encoding='utf-8') as stream:
            event = json.loads(stream.read(8192))
        if not isinstance(event, dict) or not isinstance(event.get('message'), str):
            return None
        return event
    except (OSError, ValueError):
        return None
