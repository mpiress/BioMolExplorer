import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from biomolexplorer.preparation_tools import (
    amber_command, conda_executable, resolve_amber_tool, run_amber_tool,
)
from biomolexplorer.processes import ScientificToolError


class AmberToolsExecutionTests(unittest.TestCase):
    def executable(self,root,name,text='#!/bin/sh\nexit 0\n'):
        path=Path(root)/'bin'/name;path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(text);path.chmod(0o755)
        return path

    def test_conda_executable_prefers_conda_exe(self):
        with tempfile.TemporaryDirectory() as folder:
            executable=self.executable(folder,'conda')
            with patch.dict(os.environ,{'CONDA_EXE':str(executable)},clear=True):
                self.assertEqual(conda_executable(),executable.resolve())

    def test_constructs_named_environment_command_without_shell(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);conda=self.executable(root,'conda')
            environment=root/'envs'/'amber env';self.executable(environment,'antechamber')
            with patch('biomolexplorer.preparation_tools.conda_executable',return_value=conda), \
                 patch('biomolexplorer.preparation_tools._conda_environments',return_value=[environment]):
                command,resolved=amber_command('antechamber',['-i','input with spaces.mol2'],environment='amber env')
            self.assertEqual(command[:5],[str(conda.resolve()),'run','--no-capture-output','-n','amber env'])
            self.assertEqual(command[5],str((environment/'bin'/'antechamber').resolve()))
            self.assertEqual(command[-2:],['-i','input with spaces.mol2'])
            self.assertEqual(resolved.environment,'amber env')

    def test_missing_conda_has_clear_error(self):
        with patch('biomolexplorer.preparation_tools.conda_executable',return_value=None), \
             patch('biomolexplorer.preparation_tools.shutil.which',return_value=None):
            with self.assertRaisesRegex(ScientificToolError,'Conda could not be located') as raised:
                resolve_amber_tool('tleap',environment='isolated')
        self.assertEqual(raised.exception.error_code,'CONDA_NOT_FOUND')

    def test_missing_configured_environment_has_clear_error(self):
        with tempfile.TemporaryDirectory() as folder:
            conda=self.executable(folder,'conda')
            with patch('biomolexplorer.preparation_tools.conda_executable',return_value=conda), \
                 patch('biomolexplorer.preparation_tools._conda_environments',return_value=[]):
                with self.assertRaisesRegex(ScientificToolError,'"missing" was not found') as raised:
                    resolve_amber_tool('tleap',environment='missing')
        self.assertEqual(raised.exception.error_code,'AMBER_ENV_NOT_FOUND')

    def test_missing_tool_in_existing_environment_has_clear_error(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);conda=self.executable(root,'conda');environment=root/'envs'/'amber'
            (environment/'conda-meta').mkdir(parents=True)
            with patch('biomolexplorer.preparation_tools.conda_executable',return_value=conda), \
                 patch('biomolexplorer.preparation_tools._conda_environments',return_value=[environment]):
                with self.assertRaisesRegex(ScientificToolError,'"cpptraj" is unavailable'):
                    resolve_amber_tool('cpptraj',environment='amber')

    def test_explicit_execution_preserves_working_directory_with_spaces(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);installation=root/'Amber Tools';work=root/'working directory';work.mkdir()
            self.executable(installation,'probe','#!/bin/sh\npwd > observed.txt\n')
            self.assertTrue(run_amber_tool('probe',cwd=work,prefix=installation))
            self.assertEqual((work/'observed.txt').read_text().strip(),str(work.resolve()))

    def test_failed_external_execution_reports_failure(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);self.executable(root,'broken','#!/bin/sh\necho relevant-error >&2\nexit 7\n')
            with self.assertRaisesRegex(ScientificToolError,'relevant-error') as raised:
                run_amber_tool('broken',cwd=root,prefix=root)
        self.assertEqual(raised.exception.error_code,'TOOL_EXIT_FAILED')
        self.assertIn('working_directory=',raised.exception.action)


if __name__=='__main__':unittest.main()
