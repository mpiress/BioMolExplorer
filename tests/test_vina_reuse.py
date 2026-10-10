"""Receptor reuse, scientific configuration and shared maps regression checks."""
import sys
import importlib.util
import math
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from biomolexplorer.vina_session import VinaSessions


class VinaReuseTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root/'rec.pdbqt').write_text('receptor contents')
        self.factory = MagicMock()
        self.module = patch.dict(sys.modules, {'vina': SimpleNamespace(Vina=self.factory)})
        self.module.start()
        self.addCleanup(self.module.stop)

    def config(self, compound='M1', center=0, size=12, extra=''):
        path = self.root/(compound+'.vina')
        path.write_text(f'receptor = rec.pdbqt\nligand = {compound}.pdbqt\nout = {compound}.out.pdbqt\n'
                        f'center_x = {center}\ncenter_y = 0\ncenter_z = 0\n'
                        f'size_x = {size}\nsize_y = 12\nsize_z = 12\n'
                        'seed = 73\nverbosity = 0\nexhaustiveness = 2\nnum_modes = 3\n'+extra)
        return path

    def test_many_compounds_use_one_receptor_and_all_atom_type_maps(self):
        sessions = VinaSessions()
        for index in range(15):
            self.assertTrue(sessions.run(self.config('M'+str(index)), self.root))
        self.factory.assert_called_once_with(sf_name='vina', cpu=0, seed=73, verbosity=0)
        engine = self.factory.return_value
        engine.set_receptor.assert_called_once_with(str(self.root/'rec.pdbqt'))
        engine.compute_vina_maps.assert_called_once_with(center=[0.,0.,0.], box_size=[12.,12.,12.], spacing=.375)
        self.assertEqual(engine.set_ligand_from_file.call_count, 15)
        self.assertEqual(engine.write_poses.call_count, 15)
        calls = engine.method_calls
        self.assertLess(next(i for i,c in enumerate(calls) if c[0]=='compute_vina_maps'),
                        next(i for i,c in enumerate(calls) if c[0]=='set_ligand_from_file'))
        engine.dock.assert_called_with(exhaustiveness=2,n_poses=3,min_rmsd=1.,max_evals=0)

    def test_changed_receptor_center_size_and_scoring_create_new_maps(self):
        sessions = VinaSessions()
        sessions.run(self.config(),self.root)
        sessions.run(self.config(center=3),self.root)
        sessions.run(self.config(center=3,size=15),self.root)
        (self.root/'rec.pdbqt').write_text('changed contents')
        sessions.run(self.config(center=3,size=15),self.root)
        sessions.run(self.config(center=3,size=15,extra='scoring = vinardo\n'),self.root)
        self.assertEqual(self.factory.call_count,5)

    def test_unsupported_custom_options_and_missing_api_keep_cli(self):
        sessions=VinaSessions()
        self.assertFalse(sessions.run(self.config(extra='flex = flexible.pdbqt\n'),self.root))
        self.factory.assert_not_called()
        with patch.dict(sys.modules,{'vina':None}):
            self.assertFalse(sessions.run(self.config(),self.root))

    def test_engine_failure_is_propagated_without_retry_or_double_docking(self):
        self.factory.return_value.dock.side_effect=RuntimeError('bad ligand')
        with self.assertRaisesRegex(RuntimeError,'bad ligand'):
            VinaSessions().run(self.config(),self.root)
        self.factory.return_value.dock.assert_called_once()

    def test_redocking_receptor_reused_for_every_batch_without_preparation(self):
        from wrappers.docking import prepare_docking_receptor
        prepared=self.root/'Target/Prepared';prepared.mkdir(parents=True)
        receptor=prepared/'1ABC_A.dockprep.pdbqt';receptor.write_text('unchanged receptor')
        (prepared/'centers.csv').write_text('1ABC_LIG_1A\n1\n2\n3\n')
        (prepared.parent/'pdb_codes.csv').write_text('PDB_CODE,LIGAND,RESNUM,CHAIN\n1ABC,LIG,1,A\n')
        before=receptor.stat().st_mtime_ns
        with patch('wrappers.redocking.prepare_structures') as prepare:
            for index in range(8):
                result=prepare_docking_receptor(str(self.root),'Target',self.root/str(index),None,7.4,None,True)
                self.assertEqual(result,str(self.root))
            prepare.assert_not_called()
        self.assertEqual(receptor.read_text(),'unchanged receptor')
        self.assertEqual(receptor.stat().st_mtime_ns,before)
        receptor.unlink()
        with patch('wrappers.redocking.prepare_structures') as prepare:
            with self.assertRaisesRegex(ValueError,'dockprep.pdbqt'):
                prepare_docking_receptor(str(self.root),'Target',self.root/'out',None,7.4,None,True)
            prepare.assert_not_called()

    def test_raw_receptor_conversion_happens_once_per_chain_not_reference_ligand(self):
        from caad.docking import Docking
        obj=object.__new__(Docking)
        obj.outputpath=str(self.root);obj.complexpath=str(self.root);obj.logger=MagicMock()
        obj.generate_docking_script=MagicMock();obj.process_in_parallel=MagicMock()
        obj.prepare_on_obabel=MagicMock()
        records=[['1ABC','LIG',1,'A'],['1ABC','LIG',2,'A'],['1ABC','LIG',3,'B']]
        obj.prepare_for_docking(records,'gas',7.4,False)
        receptor_calls=[c for c in obj.prepare_on_obabel.call_args_list if '.dockprep.' in c.args[0]]
        self.assertEqual(len(receptor_calls),2)
        self.assertEqual(len(obj.generate_docking_script.call_args_list),7)

    def test_docking_loop_uses_shared_session_and_keeps_reference_pose_names(self):
        from caad.docking import DockVina
        receptor=self.root/'1ABC_A.dockprep.pdbqt';receptor.write_text('prepared receptor')
        ligands=self.root/'ligands';ligands.mkdir()
        output=self.root/'poses';output.mkdir()
        for code in ('M1','M2','M3'):
            (ligands/(code+'.lig.pdbqt')).write_text('prepared ligand')
        engine=DockVina(ligand_input_path=str(ligands),receptor_input_path=str(self.root),
                        output_path=str(output),pdb_codes=[['1ABC','LIG',1,'A'],['1ABC','LIG',2,'A']],
                        centerofmasspath=str(self.root),sizeof_box=[12]*3,exhaustiveness=2,num_modes=3)
        engine.retrieve_centerofmass_dataset=MagicMock(return_value=[0,0,0])
        engine.perform_subprocess=MagicMock()
        engine.docking(str(self.root))
        engine.perform_subprocess.assert_not_called()
        self.factory.assert_called_once()
        writes=[Path(c.args[0]).name for c in self.factory.return_value.write_poses.call_args_list]
        self.assertEqual(set(writes),{f'1ABC_LIG_{residue}A_M{mol}.lig.pdbqt' for residue in (1,2) for mol in (1,2,3)})
        self.assertEqual(receptor.read_text(),'prepared receptor')


@unittest.skipUnless(importlib.util.find_spec('vina'), 'Vina Python API is not installed')
class RealVinaReuseTests(unittest.TestCase):
    def test_real_engine_docks_different_atom_types_with_one_session(self):
        from biomolexplorer.docking_data import score
        with tempfile.TemporaryDirectory(prefix='bme vina reuse ') as folder:
            root=Path(folder)
            atom='ATOM      1  C   ALA A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\n'
            (root/'rec.pdbqt').write_text(atom)
            sessions=VinaSessions();first=None
            first_pose=None
            for kind in ('OA','NA','OA'):
                (root/(kind+'.pdbqt')).write_text('ROOT\n'+atom+
                    f'ATOM      2  X   UNL A   1       1.400   0.000   0.000  1.00  0.00     0.000 {kind}\nENDROOT\nTORSDOF 0\n')
                config=root/(kind+'.vina')
                config.write_text(f'receptor = rec.pdbqt\nligand = {kind}.pdbqt\nout = {kind}.out.pdbqt\n'
                    'center_x = 0\ncenter_y = 0\ncenter_z = 0\nsize_x = 6\nsize_y = 6\nsize_z = 6\n'
                    'cpu = 1\nseed = 73\nverbosity = 0\nexhaustiveness = 1\nnum_modes = 1\n')
                self.assertTrue(sessions.run(config,root))
                if first is None:first=sessions.session
                self.assertIs(sessions.session,first)
                self.assertTrue(math.isfinite(score(root/(kind+'.out.pdbqt'),'vina')))
                if kind=='OA':
                    pose=(root/(kind+'.out.pdbqt')).read_text()
                    if first_pose is None:first_pose=pose
                    else:self.assertEqual(pose,first_pose)
            self.assertEqual((root/'rec.pdbqt').read_text(),atom)


if __name__=='__main__':
    unittest.main()
