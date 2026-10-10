"""Per-compound preparation failures never reach docking inputs."""
import csv
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock,patch

from biomolexplorer.docking_data import prepare_ligands,read_compounds,write_csv
from biomolexplorer.docking_preparation import prepare_candidates
from biomolexplorer.preparation_failures import REPORT_NAME
from biomolexplorer.operations import execute_operation
from wrappers.docking import perform_docking_vina,perform_docking_dock6

class PreparationFailureTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.input=self.root/'compounds.csv'
        write_csv(self.input,[dict(molecule_chembl_id='FIRST',canonical_smiles='CCO'),
            dict(molecule_chembl_id='BAD',canonical_smiles='CCN'),dict(molecule_chembl_id='LAST',canonical_smiles='CCC')])

    def test_failed_conversion_continues_and_removes_partial_or_stale_pose(self):
        output=self.root/'out';output.mkdir();(output/'BAD.lig.pdbqt').write_text('stale pose')
        def convert(source,destination,**options):
            destination=Path(destination);destination.write_text('partial pose')
            if destination.stem.startswith('BAD'):raise ValueError('cannot convert this compound')
        with patch('biomolexplorer.docking_data.convert_structure',side_effect=convert):
            outputs=prepare_ligands(self.input,output,'pdbqt')
        self.assertEqual([p.name for p in outputs],['FIRST.lig.pdbqt','LAST.lig.pdbqt'])
        self.assertFalse((output/'BAD.lig.pdbqt').exists())
        failures=read_compounds(output/REPORT_NAME)
        self.assertEqual([r['compound_identifier'] for r in failures],['BAD'])
        self.assertIn('cannot convert',failures[0]['reason'])
        self.assertEqual(len(read_compounds(self.input)),3)

    def test_chimera_error_with_zero_exit_is_excluded_and_both_exports_are_cleaned(self):
        count=0
        def convert(source,destination,**options):Path(destination).write_text('converted')
        def chimera(command,**options):
            nonlocal count
            count+=1;work=Path(options['cwd']);(work/'prepared.mol2').write_text('prepared')
            return SimpleNamespace(returncode=0,stdout='Error while sourcing: bad molecule' if count==2 else '',stderr='')
        with patch('biomolexplorer.docking_preparation.convert_structure',side_effect=convert),patch('biomolexplorer.docking_preparation.subprocess.run',side_effect=chimera):
            result=prepare_candidates(self.input,self.root/'both',engines='both')
        self.assertEqual([r['molecule_chembl_id'] for r in read_compounds(result)],['FIRST','LAST'])
        self.assertFalse((result.parent/'BAD.lig.mol2').exists())
        self.assertEqual(read_compounds(result.parent/REPORT_NAME)[0]['compound_identifier'],'BAD')

    def test_failed_second_export_removes_first_format_and_keeps_next_compound(self):
        count=0
        def convert(source,destination,**options):
            nonlocal count
            Path(destination).write_text('converted')
            if Path(destination).name=='BAD.lig.pdbqt':raise ValueError('export failed')
        def chimera(command,**options):
            (Path(options['cwd'])/'prepared.mol2').write_text('prepared');return SimpleNamespace(returncode=0)
        with patch('biomolexplorer.docking_preparation.convert_structure',side_effect=convert),patch('biomolexplorer.docking_preparation.subprocess.run',side_effect=chimera):
            result=prepare_candidates(self.input,self.root/'exports',engines='both')
        self.assertEqual(len(read_compounds(result)),2)
        self.assertFalse((result.parent/'BAD.lig.mol2').exists());self.assertFalse((result.parent/'BAD.lig.pdbqt').exists())

    def test_all_failed_preparation_writes_empty_table_with_headers(self):
        with patch('biomolexplorer.docking_preparation.convert_structure',side_effect=ValueError('bad compound')):
            result=prepare_candidates(self.input,self.root/'empty',engines='vina')
        self.assertEqual(read_compounds(result),[])
        self.assertIn('canonical_smiles',result.read_text())
        self.assertEqual(len(read_compounds(result.parent/REPORT_NAME)),3)

    def test_missing_engine_and_disk_errors_are_not_silently_skipped(self):
        for error in (FileNotFoundError(2,'engine missing','obabel'),OSError(28,'disk full')):
            with self.subTest(error=error),patch('biomolexplorer.docking_data.convert_structure',side_effect=error):
                with self.assertRaises(type(error)):prepare_ligands(self.input,self.root/'infrastructure','pdbqt')

    def test_operation_registers_separate_report_and_quality_exclusions(self):
        def execute(**options):
            output=Path(options['base_output_path']);prepare_ligands(self.input,output/'Molecules','pdbqt')
        def convert(source,destination,**options):
            if Path(destination).name.startswith('BAD'):raise ValueError('bad conversion')
            Path(destination).write_text('prepared')
        with patch('wrappers.docking.perform_docking_vina',side_effect=execute),patch('biomolexplorer.docking_data.convert_structure',side_effect=convert):
            result=execute_operation('docking_vina',dict(base_input_path=str(self.root),target='Target',base_selected_mols=str(self.root),mol_filename='compounds'),self.root/'worker')
        self.assertEqual(result.details['preparation_failures'],1)
        self.assertEqual(result.details['excluded_records'],1)
        self.assertIn(str(self.root/'worker'/REPORT_NAME),result.artifacts)

    def test_empty_ligand_lists_do_not_start_vina_or_dock6(self):
        for engine in ('vina','dock6'):
            with self.subTest(engine=engine),patch('wrappers.docking.prepare_docking_receptor',return_value=str(self.root)),patch('biomolexplorer.docking_data.prepare_ligands',return_value=[]),patch('biomolexplorer.docking_data.binding_center',return_value=[0,0,0]),patch('wrappers.docking.DockVina') as vina,patch('wrappers.docking.Dock6') as dock6:
                output=self.root/engine
                if engine=='vina':result=perform_docking_vina(str(self.root),'Target',str(output),str(self.root),'compounds',pdb_code=['1ABC','LIG',1,'A'])
                else:result=perform_docking_dock6(str(self.root),'Target',str(output),str(self.root),'/dock6','gas','compounds',pdb_code=['1ABC','LIG',1,'A'])
                self.assertIn('skipped_reason',result)
                vina.return_value.docking.assert_not_called();dock6.assert_not_called()
                self.assertEqual(read_compounds(output/'docking_results.csv'),[])

    def test_mixed_batch_docks_only_successful_compounds_with_each_engine(self):
        for engine in ('vina','dock6'):
            output=self.root/('mixed-'+engine)
            def convert(source,destination,**options):
                Path(destination).write_text('@<TRIPOS>ATOM\n1 C 1.0 2.0 3.0 C.3 1 LIG 0.0\n@<TRIPOS>BOND\n')
                if Path(destination).name.startswith('BAD'):raise ValueError('bad preparation')
            def run_vina(*args):
                poses=output/'Target/Vina';poses.mkdir(parents=True)
                for ligand in (output/'Molecules').glob('*.lig.pdbqt'):
                    (poses/ligand.name).write_text('REMARK VINA RESULT: -7 0 0\n')
            def run_dock6(*args):
                poses=output/'Target/1ABC_A/Dock6/flex';poses.mkdir(parents=True)
                for ligand in (output/'Target/1ABC_A/Molecules').glob('*.lig.mol2'):
                    (poses/(ligand.name.removesuffix('.lig.mol2')+'_scored.mol2')).write_text('Grid_Score: -5\n')
            with self.subTest(engine=engine),patch('wrappers.docking.prepare_docking_receptor',return_value=str(self.root)),patch('biomolexplorer.docking_data.convert_structure',side_effect=convert),patch('biomolexplorer.docking_data.binding_center',return_value=[0,0,0]),patch('wrappers.docking.DockVina') as vina,patch('wrappers.docking.Dock6') as dock6:
                vina.return_value.docking.side_effect=run_vina
                dock6.return_value.perform_dock6_evaluation.side_effect=run_dock6
                if engine=='vina':perform_docking_vina(str(self.root),'Target',str(output),str(self.root),'compounds',pdb_code=['1ABC','LIG',1,'A'])
                else:perform_docking_dock6(str(self.root),'Target',str(output),str(self.root),'/dock6','gas','compounds',pdb_code=['1ABC','LIG',1,'A'])
            self.assertEqual({r['molecule_chembl_id'] for r in read_compounds(output/'docking_results.csv')},{'FIRST','LAST'})
            reports=list(output.rglob(REPORT_NAME));self.assertTrue(reports)
            self.assertEqual(read_compounds(reports[0])[0]['compound_identifier'],'BAD')
