import logging
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from biomolexplorer.diagnostics import get_logger,log_directory
from kernel.loggers import LoggerManager

class DiagnosticTests(unittest.TestCase):
    def test_caught_exception_has_traceback_in_central_folder(self):
        with tempfile.TemporaryDirectory() as folder,patch.dict(os.environ,{'BIOMOL_LOG_DIR':folder}):
            logger=get_logger('frontend')
            try: raise ValueError('Falha de teste')
            except ValueError: logger.exception('Configuração do bloco')
            for handler in logger.handlers: handler.flush()
            contents=(Path(folder)/'frontend.log').read_text()
            self.assertIn('Traceback',contents);self.assertIn('Falha de teste',contents)
            self.assertIn('Falha de teste',(Path(folder)/'errors.log').read_text())
            self.assertEqual(contents.count('Configuração do bloco'),1)
            self.assertEqual(log_directory(),Path(folder))
            for handler in list(logger.handlers): handler.close();logger.removeHandler(handler)
    def test_scientific_module_uses_same_root(self):
        with tempfile.TemporaryDirectory() as folder,patch.dict(os.environ,{'BIOMOL_LOG_DIR':folder}):
            logger=LoggerManager.get_logger('biomol-test-science','logs/loaders.log')
            logger.error('API indisponível')
            for handler in logger.handlers: handler.flush()
            self.assertIn('API indisponível',(Path(folder)/'loaders.log').read_text())
            for handler in list(logger.handlers):handler.close();logger.removeHandler(handler)
