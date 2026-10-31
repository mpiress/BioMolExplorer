"""Preparation must never discover standalone AmberTools/Conda tools."""
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from biomolexplorer.preparation_tools import resolve_amber_tool,run_amber_tool
from biomolexplorer.processes import ScientificToolError


class BundledPreparationTests(unittest.TestCase):
    def installation(self,folder):
        root=Path(folder)/'Chimera installation';binary=root/'bin/chimera'
        binary.parent.mkdir(parents=True);binary.write_text('#!/bin/sh\n');binary.chmod(0o755)
        tool=root/'bin/amber18/bin/antechamber';tool.parent.mkdir(parents=True)
        tool.write_text('#!/bin/sh\n');tool.chmod(0o755)
        return root,binary,tool

    def test_external_prefix_and_environment_cannot_override_chimera(self):
        with tempfile.TemporaryDirectory() as folder:
            root,binary,tool=self.installation(folder)
            with patch('biomolexplorer.preparation_tools.shutil.which',return_value=str(binary)),patch.dict(os.environ,AMBERHOME='/external',BIOMOL_AMBERTOOLS_ENV='external'):
                result=resolve_amber_tool('antechamber','/external','external')
            self.assertEqual(result.executable,tool.resolve())
            self.assertEqual(result.prefix,root/'bin/amber18')

    def test_no_conda_or_path_fallback_when_chimera_is_missing(self):
        with patch('biomolexplorer.preparation_tools.shutil.which',return_value=None):
            with self.assertRaisesRegex(ScientificToolError,'Chimera'):resolve_amber_tool('antechamber')

    def test_external_md_and_topology_protocols_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root,binary,tool=self.installation(folder)
            with patch('biomolexplorer.preparation_tools.shutil.which',return_value=str(binary)):
                for name in ('sander','cpptraj','tleap'):
                    with self.subTest(name=name),self.assertRaisesRegex(ScientificToolError,'retired'):resolve_amber_tool(name)

    def test_execution_sets_only_child_environment_and_uses_absolute_program(self):
        with tempfile.TemporaryDirectory() as folder:
            root,binary,tool=self.installation(folder)
            original=dict(os.environ)
            with patch('biomolexplorer.preparation_tools.shutil.which',return_value=str(binary)),patch('biomolexplorer.preparation_tools.run_command',return_value=True) as run:
                self.assertTrue(run_amber_tool('antechamber',['-i','input with spaces.mol2'],cwd=folder))
            self.assertEqual(run.call_args.args[0],[str(tool),'-i','input with spaces.mol2'])
            self.assertEqual(run.call_args.kwargs['env']['AMBERHOME'],str(root/'bin/amber18'))
            self.assertEqual(run.call_args.kwargs['env']['LD_LIBRARY_PATH'].split(os.pathsep)[0],str(root/'lib'))
            self.assertEqual(dict(os.environ),original)
