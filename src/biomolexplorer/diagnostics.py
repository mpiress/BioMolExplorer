"""Contextual text and JSONL diagnostics, shared by UI, supervisors and science."""
import json
import logging
import os
import re
import sys
import threading
from contextlib import contextmanager
from contextvars import ContextVar
from datetime import datetime, timezone
from logging.handlers import RotatingFileHandler
from pathlib import Path
from .paths import SOURCE_ROOT

_LOCK = threading.RLock()
_CONTEXT = ContextVar('biomol_log_context', default={})
_CONTEXT_KEYS = ('project_id', 'run_id', 'stage_id', 'job_id', 'operation', 'pair', 'command_id')
_SECRET = re.compile(r'(?i)(\b(?:password|passwd|token|secret|api[_-]?key|authorization)\b\s*[=:]\s*)(?:Bearer\s+)?([^\s,;]+)')


def redact(value):
    """Defensive masking; callers must still avoid logging credentials or parameters."""
    return _SECRET.sub(r'\1[REDACTED]', str(value))


def log_directory():
    default = SOURCE_ROOT.parent/'logs' if (SOURCE_ROOT.parent/'pyproject.toml').exists() else Path.cwd()/'logs'
    return Path(os.environ.get('BIOMOL_LOG_DIR', default)).expanduser().resolve()


def current_context():
    context = {}
    try:
        inherited = json.loads(os.environ.get('BIOMOL_LOG_CONTEXT', '{}'))
        if isinstance(inherited, dict):
            context = {k: inherited[k] for k in _CONTEXT_KEYS if k in inherited}
    except (TypeError, ValueError):
        pass
    context.update(_CONTEXT.get())
    return context


@contextmanager
def log_context(**values):
    """Scoped, thread/task-local correlation; inherited explicitly by worker processes."""
    token = _CONTEXT.set({**_CONTEXT.get(), **{k: v for k, v in values.items() if k in _CONTEXT_KEYS and v is not None}})
    try:
        yield
    finally:
        _CONTEXT.reset(token)


class ContextFilter(logging.Filter):
    def filter(self, record):
        context = current_context()
        for key in _CONTEXT_KEYS:
            if not hasattr(record, key):
                setattr(record, key, context.get(key))
        if record.exc_info and record.exc_info[0] is not None and not hasattr(record, 'error_code'):
            for key, value in diagnose_exception(record.exc_info[1]).items():
                setattr(record, key, value)
        if not hasattr(record, 'event'):
            record.event = 'exception' if record.exc_info else 'message'
        return True


class TextFormatter(logging.Formatter):
    def format(self, record):
        timestamp = datetime.fromtimestamp(record.created, timezone.utc).isoformat(timespec='milliseconds')
        context = ' '.join(f'{k}={getattr(record,k)}' for k in _CONTEXT_KEYS if getattr(record,k,None) is not None)
        location = f'{record.pathname}:{record.lineno}'
        message = record.getMessage()
        details = []
        for key in ('error_code', 'action', 'duration_ms', 'returncode', 'tool', 'cwd', 'diagnostic_path', 'configuration'):
            if hasattr(record, key): details.append(f'{key}={getattr(record,key)}')
        result = f'{timestamp} {record.levelname} {record.name} event={record.event} [{context or "application"}] {message}'
        if details: result += '\n  ' + ' | '.join(details)
        if record.levelno >= logging.WARNING: result += f'\n  source={location}'
        if record.exc_info: result += '\n' + self.formatException(record.exc_info)
        return redact(result)


class JsonFormatter(logging.Formatter):
    def format(self, record):
        result = dict(schema_version=1, timestamp=datetime.fromtimestamp(record.created, timezone.utc).isoformat(timespec='milliseconds'),
                      level=record.levelname, logger=record.name, event=record.event,
                      message=redact(record.getMessage()), pid=record.process, thread=record.threadName,
                      source=dict(file=record.pathname, line=record.lineno, function=record.funcName))
        for key in (*_CONTEXT_KEYS, 'error_code', 'action', 'duration_ms', 'returncode', 'tool', 'cwd', 'diagnostic_path', 'artifacts', 'channel', 'configuration'):
            value = getattr(record, key, None)
            if value is not None: result[key] = redact(value) if isinstance(value, str) else value
        if record.exc_info and record.exc_info[0] is not None:
            result['exception'] = dict(type=record.exc_info[0].__name__, message=redact(record.exc_info[1]),
                                       traceback=redact(self.formatException(record.exc_info)))
        return json.dumps(result, ensure_ascii=False, default=str)


class DeferredRotatingHandler(RotatingFileHandler):
    def emit(self, record):
        # Workers use process pools. Serialize writes/rotation across processes,
        # and reopen after another process has rotated the shared file.
        import fcntl
        try:
            path = Path(self.baseFilename)
            path.parent.mkdir(parents=True, exist_ok=True)
            with Path(str(path)+'.lock').open('a') as lock:
                fcntl.flock(lock, fcntl.LOCK_EX)
                if self.stream is not None:
                    try:
                        stale = os.fstat(self.stream.fileno()).st_ino != path.stat().st_ino
                    except FileNotFoundError:
                        stale = True
                    if stale:
                        self.stream.close()
                        self.stream = None
                super().emit(record)
        except OSError:
            # A diagnostic write must not invalidate scientific output.
            pass

    def _open(self):
        Path(self.baseFilename).parent.mkdir(parents=True, exist_ok=True)
        return super()._open()


# Files have bounded retention; created on first event, never on import.
def _file_handler(path, level, structured=False):
    handler = DeferredRotatingHandler(path, maxBytes=5*1024*1024, backupCount=3, delay=True, encoding='utf-8')
    handler.setLevel(level)
    handler.setFormatter(JsonFormatter() if structured else TextFormatter())
    handler.addFilter(ContextFilter())
    handler._biomol_diagnostic = True
    return handler


_SHARED = {}
_CONSOLE = None


def get_logger(component, *, name=None, level=logging.INFO):
    logger = logging.getLogger(name or 'biomolexplorer.'+component)
    logger.setLevel(level)
    logger.propagate = False
    directory = log_directory()
    paths = {directory/(Path(component).name+'.log'): (logging.DEBUG, False),
             directory/'errors.log': (logging.ERROR, False), directory/'events.jsonl': (logging.DEBUG, True)}
    with _LOCK:
        for old in list(logger.handlers):
            if getattr(old, '_biomol_diagnostic', False) and isinstance(old, logging.FileHandler) and Path(old.baseFilename) not in paths:
                logger.removeHandler(old)
                if old in logging.getLogger().handlers: logging.getLogger().removeHandler(old)
                # Shared handlers can still belong to other loggers; delayed open permits reuse.
                old.close()
        for path, (threshold, structured) in paths.items():
            if not any(isinstance(h, logging.FileHandler) and h.baseFilename == str(path) for h in logger.handlers):
                key = (str(path), threshold, structured)
                handler = _SHARED.get(key)
                if handler is None:
                    handler = _SHARED[key] = _file_handler(path, threshold, structured)
                logger.addHandler(handler)
        if _CONSOLE is not None and _CONSOLE not in logger.handlers: logger.addHandler(_CONSOLE)
    return logger


def event(logger, name, message, *, level=logging.INFO, **fields):
    logger.log(level, message, extra={'event': name, **fields}, stacklevel=2)


def diagnose_exception(exc):
    """Stable categories and explicit next steps, without guessing scientific validity."""
    visited = set()
    while getattr(exc, '__cause__', None) is not None and not getattr(exc, 'error_code', None) and id(exc) not in visited:
        visited.add(id(exc))
        exc = exc.__cause__
    if getattr(exc, 'error_code', None):
        return {'error_code': exc.error_code, 'action': exc.action}
    if isinstance(exc, (TimeoutError,)):
        return {'error_code': 'EXECUTION_TIMEOUT', 'action': 'Review the timeout and the last running command before retrying.'}
    if isinstance(exc, FileNotFoundError):
        return {'error_code': 'INPUT_NOT_FOUND', 'action': 'Verify the reported file path and the selected input files.'}
    if isinstance(exc, PermissionError):
        return {'error_code': 'ACCESS_DENIED', 'action': 'Verify project permissions and filesystem access.'}
    if isinstance(exc, ValueError):
        return {'error_code': 'VALIDATION_FAILED', 'action': 'Correct the input or configuration described in the error before retrying.'}
    return {'error_code': 'UNEXPECTED_ERROR', 'action': 'Inspect the original exception and source location in errors.log.'}


def write_summary(status, *, error=None, directory=None, **fields):
    """Atomic single-job diagnostic summary. Logging failures cannot change a job result."""
    payload = {'schema_version': 1, 'updated_at': datetime.now(timezone.utc).isoformat(), 'status': status, **current_context(), **fields}
    if error is not None:
        payload.update(diagnose_exception(error), error=redact(f'{type(error).__name__}: {error}'))
    try:
        path = (Path(directory) if directory is not None else log_directory())/'diagnostic.json'
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(f'.{os.getpid()}.tmp')
        temporary.write_text(json.dumps(payload, ensure_ascii=False, default=str), encoding='utf-8')
        temporary.replace(path)
    except OSError:
        pass


def configure_logging(component='backend', *, console=False):
    global _CONSOLE
    with _LOCK:
        if console and _CONSOLE is None:
            _CONSOLE = logging.StreamHandler()
            _CONSOLE.setFormatter(TextFormatter())
            _CONSOLE.addFilter(ContextFilter())
            _CONSOLE._biomol_diagnostic = True
    logger = get_logger(component)
    root = logging.getLogger()
    for handler in list(root.handlers):
        if getattr(handler, '_biomol_diagnostic', False) and handler not in logger.handlers:
            root.removeHandler(handler)
    for handler in logger.handlers:
        if handler not in root.handlers: root.addHandler(handler)
    root.setLevel(logging.INFO)
    logging.captureWarnings(True)
    previous = sys.excepthook
    if not getattr(previous, '_biomol_hook', False):
        def exception_hook(kind, value, tb):
            get_logger(component).error('Unhandled exception', exc_info=(kind, value, tb),
                                        extra={'event': 'application.failed', **diagnose_exception(value)})
            previous(kind, value, tb)
        exception_hook._biomol_hook = True
        sys.excepthook = exception_hook
    prior_thread = threading.excepthook
    if not getattr(prior_thread, '_biomol_hook', False):
        def thread_hook(args):
            get_logger(component).error('Unhandled thread exception', exc_info=(args.exc_type, args.exc_value, args.exc_traceback),
                                        extra={'event': 'thread.failed', **diagnose_exception(args.exc_value)})
            prior_thread(args)
        thread_hook._biomol_hook = True
        threading.excepthook = thread_hook
    return logger
