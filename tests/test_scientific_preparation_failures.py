"""Classic Chimera's zero exit status must not hide failed preparation commands."""
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch
from biomolexplorer.processes import run_command
from caad.docking import Docking

ATOM = 'ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n'

class ScientificPreparationFailures(unittest.TestCase):
    def docking(self, root):
        dock=object.__new__(Docking)
        dock.outputpath=str(root);dock.logger=MagicMock()
        return dock

    def test_chimera_zero_exit_status_with_command_error_fails(self):
        result=subprocess.CompletedProcess(['chimera'],0,'', 'Error while sourcing prepare.com, line 0:\nNo such file or directory')
        with patch('biomolexplorer.processes.subprocess.run',return_value=result):
            with self.assertRaisesRegex(RuntimeError,'Error while sourcing'):
                run_command(['chimera','--nogui','--silent','prepare.com'])

    def test_normal_chimera_output_is_logged_and_succeeds(self):
        result=subprocess.CompletedProcess(['chimera'],0,'Wrote receptor.pdb','')
        with patch('biomolexplorer.processes.subprocess.run',return_value=result), self.assertLogs('biomolexplorer.scientific_tools') as log:
            self.assertTrue(run_command(['chimera','--nogui','--silent','prepare.com']))
        self.assertIn('Wrote receptor.pdb','\n'.join(log.output))

    def test_missing_chimera_output_keeps_script_and_blocks_next_step(self):
        with tempfile.TemporaryDirectory(prefix='prepared structures ') as folder:
            root=Path(folder);dock=self.docking(root);dock.perform_subprocess=MagicMock(return_value=True)
            script=root/'prepare.com';script.write_text(f'write format pdb #0 {root}/ligand.pdb\n')
            with self.assertRaisesRegex(ValueError,'Arquivo preparado ausente ou vazio'):dock.prepare_on_chimera(script.name)
            self.assertTrue(script.exists())
            (root/'ligand.pdb').write_text(ATOM)
            self.assertTrue(dock.prepare_on_chimera(script.name))
            self.assertFalse(script.exists())

    def test_empty_openbabel_output_fails_even_with_successful_exit(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);dock=self.docking(root)
            (root/'ligand.pdb').write_text(ATOM)
            def convert(*args): (root/'ligand.pdbqt').write_text('')
            dock.perform_subprocess=MagicMock(side_effect=convert)
            with self.assertRaisesRegex(ValueError,'Arquivo preparado ausente ou vazio'):
                dock.prepare_on_obabel('ligand.pdb','ligand.pdbqt',input_format='pdb',output_format='pdbqt')
            (root/'ligand.pdb').unlink();dock.perform_subprocess.reset_mock()
            with self.assertRaisesRegex(ValueError,'Arquivo preparado ausente ou vazio'):
                dock.prepare_on_obabel('ligand.pdb','ligand.pdbqt',input_format='pdb',output_format='pdbqt')
            dock.perform_subprocess.assert_not_called()

    def test_files_without_atoms_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            for name,content,kind in (('empty.pdb','END\n','pdb'),('empty.mol2','@<TRIPOS>ATOM\n@<TRIPOS>BOND\n','mol2')):
                path=root/name;path.write_text(content)
                with self.assertRaisesRegex(ValueError,'Arquivo preparado sem átomos'):
                    Docking.validate_prepared_file(path,kind)
