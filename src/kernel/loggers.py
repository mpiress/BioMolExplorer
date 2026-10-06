"""Compatibility adapter for the shared application diagnostic system."""
import logging
from pathlib import Path


class LoggerManager:
    @classmethod
    def get_logger(cls, name, log_file=None, level=logging.INFO):
        from biomolexplorer.diagnostics import get_logger
        component = Path(log_file).stem if log_file else name
        return get_logger(component, name=name, level=level)
