"""Logging without filesystem writes or warning suppression during import."""
import logging
import threading
from pathlib import Path


class LoggerManager:
    _lock = threading.RLock()

    @classmethod
    def get_logger(cls, name, log_file=None, level=logging.INFO):
        logger = logging.getLogger(name)
        logger.setLevel(level)
        if log_file:
            from biomolexplorer.diagnostics import log_directory
            path = log_directory() / Path(log_file).name
            with cls._lock:
                if not any(isinstance(handler, logging.FileHandler) and handler.baseFilename == str(path)
                           for handler in logger.handlers):
                    from biomolexplorer.diagnostics import DeferredRotatingHandler
                    handler = DeferredRotatingHandler(path,maxBytes=5*1024*1024,backupCount=3,delay=True,encoding='utf-8')
                    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(name)s %(message)s'))
                    logger.addHandler(handler)
        return logger
