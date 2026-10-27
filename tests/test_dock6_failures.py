"""DOCK6 isolates native crashes and rejects unusable stage outputs."""
import csv
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from caad.docking import Dock6
from biomolexplorer.docking_data import write_csv
from biomolexplorer.dock6_validation import compound_code,validate_mol2
from biomolexplorer.preparation_failures import PreparationFailures
from biomolexplorer.processes import ScientificToolError

MOL2='@<TRIPOS>MOLECULE\nligand\n1 0 1 0 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 C 0 0 0 C.3 1 LIG 0.0\n@<TRIPOS>BOND\n@<TRIPOS>SUBSTRUCTURE\n1 LIG 1\n'
SCORED='########## Grid_Score: -5.0\n'+MOL2


class Dock6FailureTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        # External commands are mocked below; the cache still fingerprints the
        # installation directories before collecting per-compound failures.
        (self.root/'tools/bin').mkdir(parents=True)
        (self.root/'tools/parameters').mkdir()
        ligands=self.root/'input';ligands.mkdir()
        receptor=self.root/'rec';receptor.mkdir()
        (receptor/'1ABC_A.dockprep.mol2').write_text(MOL2)
        self.codes=('FIRST','BAD','LAST','A.B')
        for code in self.codes:(ligands/(code+'.lig.mol2')).write_text(MOL2)
        self.dataset=self.root/'compounds.csv'
        write_csv(self.dataset,[{'molecule_chembl_id':code,'canonical_smiles':'C'} for code in self.codes])
        self.report=PreparationFailures()
        self.engine=Dock6(dock6_path=str(self.root/'tools'),ligand_input_path=str(ligands),
                          receptor_input_path=str(receptor),base_output_path=str(self.root/'output'),
                          pdb_code='1ABC_A',binding_site_center=[1,2,3])
        self.engine.set_compound_report(self.dataset,self.report)
        self.engine.nprocess=2
        self.addCleanup(lambda:shutil.rmtree(self.engine.runtime_path,ignore_errors=True))
        self.calls=[]

    def external(self,command,cwd=None,**kwargs):
        config=Path(command[command.index('-i')+1])
        values=dict(line.split(None,1) for line in config.read_text().splitlines() if line.strip())
        prefix=Path(values['ligand_outfile_prefix']);code=compound_code(values['ligand_atom_file'])
        self.calls.append((code,values))
        (prefix.parent/(prefix.name+'_scored.mol2')).write_text(SCORED)
        if code=='BAD':raise ScientificToolError('dock6 terminou com SIGSEGV (-11)','TOOL_EXIT_FAILED','Revise o composto.')
        return True

    def minimized(self):
        folder=Path(self.engine._Dock6__base_output_path)/'energy_min';folder.mkdir()
        for code in self.codes:(folder/(code+'.lig.min_scored.mol2')).write_text(SCORED)
        return folder

    def test_native_docking_crash_skips_compound_cleans_partial_pose_and_keeps_others(self):
        self.minimized()
        with patch.object(self.engine,'perform_subprocess',side_effect=self.external):self.engine.perform_dock6_evaluation()
        output=Path(self.engine.outputpath)
        self.assertEqual({p.name for p in output.glob('*_scored.mol2')},{c+'_scored.mol2' for c in self.codes if c!='BAD'})
        self.assertEqual(self.report.records[0]['compound_identifier'],'BAD')
        self.assertIn('-11',self.report.records[0]['reason'])
        self.assertEqual(self.report.records[0]['preparation_step'],'Docking DOCK6')
        self.assertTrue((self.root/'output/preparation_failures.csv').exists())
        self.assertTrue(list(output.glob('*.in')))
        for code,values in self.calls:
            self.assertTrue(values['grid_score_grid_prefix'].endswith('/binding_site'))
            self.assertTrue(values['receptor_site_file'].endswith('/binding_site.sph'))
        self.assertTrue(any(v['ligand_outfile_prefix'].endswith('/A.B') for _,v in self.calls))

    def test_minimization_failure_is_not_sent_to_docking(self):
        with patch.object(self.engine,'perform_subprocess',side_effect=self.external):
            self.engine.prepare_minimization();self.calls=[];self.engine.perform_dock6_evaluation()
        self.assertEqual({c for c,_ in self.calls},set(self.codes)-{'BAD'})
        self.assertEqual(len(self.report.records),1)

    def test_success_exit_with_no_usable_pose_is_also_excluded(self):
        self.minimized()
        with patch.object(self.engine,'perform_subprocess',return_value=True):self.engine.perform_dock6_evaluation()
        self.assertEqual(len(self.report.records),len(self.codes))
        self.assertFalse(self.engine._ligand_files(self.engine.ligandpath))

    def test_footprint_failure_is_isolated_and_removes_docking_pose(self):
        self.minimized()
        folder=Path(self.engine._Dock6__base_output_path)/'flex';folder.mkdir()
        for code in self.codes:(folder/(code+'_scored.mol2')).write_text(SCORED)
        def footprint(command,cwd=None,**kwargs):
            config=Path(command[command.index('-i')+1]);values=dict(line.split(None,1) for line in config.read_text().splitlines() if line.strip())
            prefix=Path(values['ligand_outfile_prefix'])
            if prefix.name=='BAD':raise ScientificToolError('SIGSEGV','TOOL_EXIT_FAILED','Revise o composto.')
            (prefix.parent/(prefix.name+'_footprint_scored.txt')).write_text('ALA 1 -1.0 -2.0 x -1.0 -2.0 x\n')
            return True
        with patch.object(self.engine,'perform_subprocess',side_effect=footprint):self.engine.prepare_footprint(docked=True)
        self.assertFalse((folder/'BAD_scored.mol2').exists())
        self.assertTrue((folder/'A.B_scored.mol2').exists())
        self.assertEqual(self.report.records[0]['preparation_step'],'Footprint DOCK6')

    def test_missing_engine_is_global_error_and_does_not_become_bad_compounds(self):
        with self.assertRaises(ScientificToolError):self.engine.validate_setup()
        self.assertEqual(self.report.records,[])

    def test_shared_site_generates_one_selector_box_and_grid_for_all_compounds(self):
        calls=[]
        def surface(source,destination,**options):
            Path(destination).write_text('surface');return {}
        def external(command,cwd=None,**options):
            command=command.split() if isinstance(command,str) else command
            calls.append(command[0]);folder=Path(cwd)
            if command[0]=='sphgen':(folder/'1ABC_A.sph').write_text('spheres')
            if command[0]=='sphere_selector':(folder/'selected_spheres.sph').write_text('cluster 1 number 1\n1 1 2 3 1.5 1 0 0\n')
            if command[0]=='showbox':(folder/'binding_site.box.pdb').write_text('box')
            if command[0]=='grid':
                for suffix in ('.nrg','.bmp'):(folder/('binding_site'+suffix)).write_text('valid grid')
            return True
        with patch('caad.docking.generate_surface',side_effect=surface),patch.object(self.engine,'perform_subprocess',side_effect=external):
            self.engine.prepare_surface();self.engine.prepare_showbox();self.engine.prepare_gridbox()
        self.assertEqual(calls,['sphgen','sphere_selector','showbox','grid'])
        self.assertEqual(len(list(Path(self.engine.outputpath).glob('*.nrg'))),1)

    def test_missing_executable_from_jobs_is_not_suppressed(self):
        self.minimized()
        with patch.object(self.engine,'perform_subprocess',side_effect=ScientificToolError('absent','TOOL_NOT_FOUND','Install DOCK6')):
            with self.assertRaises(ScientificToolError):self.engine.perform_dock6_evaluation()
        self.assertEqual(self.report.records,[])

    def test_mol2_validation_detects_disconnected_invalid_numeric_and_empty_outputs(self):
        path=self.root/'check.mol2';path.write_text(MOL2)
        validate_mol2(path)
        for text in ('',MOL2.replace('0.0','nan'),MOL2.replace('1 0 1 0 0','2 0 1 0 0').replace('@<TRIPOS>BOND','2 C 5 0 0 C.3 1 LIG 0.0\n@<TRIPOS>BOND')):
            path.write_text(text)
            with self.assertRaises(ValueError):validate_mol2(path)
        self.assertEqual(compound_code('A.B.lig.min_scored.mol2'),'A.B')


if __name__=='__main__':unittest.main()
