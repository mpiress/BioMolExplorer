"""Central diagnostic files; project execution logs remain private and downloadable."""
import logging
import os
import sys
import threading
from logging.handlers import RotatingFileHandler
from pathlib import Path
from .paths import SOURCE_ROOT

_LOCK=threading.RLock()


def log_directory():
    default=SOURCE_ROOT.parent/'logs' if (SOURCE_ROOT.parent/'pyproject.toml').exists() else Path.cwd()/'logs'
    return Path(os.environ.get('BIOMOL_LOG_DIR',default)).expanduser().resolve()


class DeferredRotatingHandler(RotatingFileHandler):
    def _open(self):
        Path(self.baseFilename).parent.mkdir(parents=True,exist_ok=True)
        return super()._open()


def get_logger(component):
    logger=logging.getLogger('biomolexplorer.'+component)
    logger.setLevel(logging.INFO)
    logger.propagate=False
    path=log_directory()/(component+'.log')
    errors=log_directory()/'errors.log'
    with _LOCK:
        for old in list(logger.handlers):
            if getattr(old,'_biomol_diagnostic',False) and old.baseFilename not in (str(path),str(errors)):
                logger.removeHandler(old)
                if old in logging.getLogger().handlers: logging.getLogger().removeHandler(old)
                old.close()
        if not any(isinstance(h,logging.FileHandler) and h.baseFilename==str(path) for h in logger.handlers):
            handler=DeferredRotatingHandler(path,maxBytes=5*1024*1024,backupCount=3,delay=True,encoding='utf-8')
            handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(name)s [pid=%(process)d thread=%(threadName)s] %(message)s'))
            handler._biomol_diagnostic=True
            logger.addHandler(handler)
        if not any(isinstance(h,logging.FileHandler) and h.baseFilename==str(errors) for h in logger.handlers):
            error_handler=next((h for h in logging.getLogger().handlers if isinstance(h,logging.FileHandler) and h.baseFilename==str(errors)),None)
            if error_handler is None:
                error_handler=DeferredRotatingHandler(errors,maxBytes=5*1024*1024,backupCount=3,delay=True,encoding='utf-8')
                error_handler.setLevel(logging.ERROR)
                error_handler._biomol_diagnostic=True
                error_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(name)s [pid=%(process)d thread=%(threadName)s] %(message)s'))
            logger.addHandler(error_handler)
    return logger


def configure_logging(component='backend'):
    logger=get_logger(component)
    root=logging.getLogger()
    # Capture Flet/framework errors and unhandled exceptions as well as our handlers.
    for handler in logger.handlers:
        if handler not in root.handlers: root.addHandler(handler)
    errors=log_directory()/'errors.log'
    if not any(isinstance(h,logging.FileHandler) and h.baseFilename==str(errors) for h in root.handlers):
        handler=DeferredRotatingHandler(errors,maxBytes=5*1024*1024,backupCount=3,delay=True,encoding='utf-8')
        handler.setLevel(logging.ERROR)
        handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(name)s [pid=%(process)d] %(message)s'))
        root.addHandler(handler)
    root.setLevel(logging.INFO)
    logging.captureWarnings(True)
    previous=sys.excepthook
    if not getattr(previous,'_biomol_hook',False):
        def exception_hook(kind,value,tb):
            logger.error('Exceção não tratada',exc_info=(kind,value,tb))
            previous(kind,value,tb)
        exception_hook._biomol_hook=True
        sys.excepthook=exception_hook
    prior_thread=threading.excepthook
    if not getattr(prior_thread,'_biomol_hook',False):
        def thread_hook(args):
            logger.error('Exceção não tratada na thread %s',args.thread.name,exc_info=(args.exc_type,args.exc_value,args.exc_traceback))
            prior_thread(args)
        thread_hook._biomol_hook=True
        threading.excepthook=thread_hook
    return logger
