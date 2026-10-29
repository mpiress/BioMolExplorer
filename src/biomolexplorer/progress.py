"""Optional, atomic worker feedback; scientific functions also work without a UI."""
import json
import os
import threading
import time
from pathlib import Path

_LOCK = threading.Lock()


def report_progress(message=None, completed=None, total=None, *, key=None, arguments=None):
    filename = os.environ.get('BIOMOL_PROGRESS_FILE')
    if not filename:
        return
    # Monitoring failures must never invalidate scientific output.
    try:
        from .i18n import RuntimeMessages
        messages = RuntimeMessages('pt')
        if key is None and isinstance(message, str) and message in messages.catalog:key = message
        if key is not None:message = messages(key, **(arguments or {}))
        event = {'message': str(message)[:1500], 'updated_at': time.time()}
        if key is not None:event.update(message_key=key,message_args=arguments or {})
        if completed is not None and total is not None:
            event.update(completed=int(completed), total=int(total))
        payload = json.dumps(event, ensure_ascii=False)
        if len(payload)>8000:
            event.pop('message_key',None);event.pop('message_args',None)
            payload=json.dumps(event,ensure_ascii=False)
        with _LOCK:
            path = Path(filename)
            temporary = path.with_suffix('.tmp')
            temporary.write_text(payload, encoding='utf-8')
            temporary.replace(path)
    except (OSError,KeyError,TypeError,ValueError):
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
