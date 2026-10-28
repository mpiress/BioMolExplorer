"""Regressions from real Alzheimer ligand and receptor preparation failures."""
import ast
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from biomolexplorer.chimera_preparation import (
    ligand_hydrogenation,ligand_charging,validate_ligand_preparation,
    recover_hydrogen_placement,record_recovery,receptor_cleanup)
from biomolexplorer.ligand_mol2 import prepare_chimera_ligand

FIXTURES=Path(__file__).parent/'fixtures/redocking'


class ChimeraPreparationTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)

    def script(self,prefix='1XLU_ISP_1001A',method='am1'):
        path=self.root/('prepare_ligand_'+prefix+'.com');output=self.root/(prefix+'.lig.mol2')
        path.write_text('open input.mol2\naddh\naddcharge all spec #0 chargeModel 14sb method '+method+
            '\nminimize spec #0\nwrite format mol2 #0 '+str(output)+'\nclose all\n')
        return path,output

    def test_charge_comes_from_actual_hydrogenated_graph_not_chimera_heuristic(self):
        for prefix,charge in [('1XLU_ISP_1001A',1),('3VF3_0GS_401A',0)]:
            for method in ('am1','gas'):
                with self.subTest(prefix=prefix,method=method):
                    script,output=self.script(prefix,method);phase=ligand_hydrogenation(script)
                    self.assertNotIn('addcharge',phase['script'].read_text())
                    shutil.copy2(FIXTURES/(prefix+'.hydrogenated.mol2'),phase['intermediate'])
                    ligand_charging(script,phase)
                    helper=script.with_suffix('.charge.py').read_text()
                    self.assertIn('model.residues, '+str(charge),helper)
                    self.assertIn('method='+repr(method),helper)
                    self.assertIn('minimize spec #0',script.read_text())
                    self.assertIn('minimize spec #0 prep false',script.read_text())
                    self.assertEqual(json.loads(output.with_suffix('.preparation.json').read_text())['net_charge'],charge)

    def test_heavy_atom_rmsd_excludes_exceptional_phosphoryl_hydrogens(self):
        from biomolexplorer.dock6_validation import validate_mol2
        from biomolexplorer.redocking_dock6 import heavy_atom_rmsd
        source=FIXTURES/'1XLU_ISP_1001A.hydrogenated.mol2';mol=validate_mol2(source)
        lines=[]
        for i,atom in enumerate(mol.GetAtoms(),1):
            if atom.GetAtomicNum()==1:continue
            point=mol.GetConformer().GetAtomPosition(atom.GetIdx());name=atom.GetProp('_TriposAtomName')
            lines.append(f'ATOM  {i:5d} {name:>4} UNL     1    {point.x:8.3f}{point.y:8.3f}{point.z:8.3f}  0.00  0.00     0.000 {atom.GetSymbol()}')
        pose=self.root/'pose.pdbqt';pose.write_text('\n'.join(lines)+'\n')
        self.assertLess(heavy_atom_rmsd(source,pose),.001)

    def test_final_charge_validation_rejects_zero_charges_for_charged_ligand(self):
        script,output=self.script();phase=ligand_hydrogenation(script)
        source=FIXTURES/'1XLU_ISP_1001A.hydrogenated.mol2'
        shutil.copy2(source,phase['intermediate']);shutil.copy2(source,output)
        ligand_charging(script,phase)
        with self.assertRaisesRegex(ValueError,'cargas parciais'):validate_ligand_preparation(phase)

    def test_final_validation_records_preserved_graph_and_neutral_charge(self):
        script,output=self.script('3VF3_0GS_401A');phase=ligand_hydrogenation(script)
        source=FIXTURES/'3VF3_0GS_401A.hydrogenated.mol2'
        shutil.copy2(source,phase['intermediate']);shutil.copy2(source,output)
        ligand_charging(script,phase);validate_ligand_preparation(phase)
        report=json.loads(output.with_suffix('.preparation.json').read_text())
        self.assertTrue(report['validated_graph_preserved']);self.assertEqual(report['partial_charge_sum'],0)

    def test_known_native_hbond_bug_retries_only_hydrogen_placement_once(self):
        script,output=self.script();before=script.read_text()
        self.assertFalse(recover_hydrogen_placement(script,'unrelated AM1 failure'))
        self.assertEqual(script.read_text(),before)
        error='accPhiPsi() takes exactly 8 arguments (9 given)'
        self.assertTrue(recover_hydrogen_placement(script,error))
        self.assertFalse(recover_hydrogen_placement(script,error))
        self.assertEqual(script.read_text(),before.replace('\naddh\n','\naddh hbond false\n'))
        record_recovery(script)
        self.assertEqual(json.loads(output.with_suffix('.preparation.json').read_text())['charge_and_minimization_settings'],'unchanged')

    def test_ccd_oxygen_and_substituted_pyrrole_types_before_addh(self):
        from biomolexplorer.ligand_mol2 import component_definition
        cache=self.root/'cache/calculations/ccd';cache.mkdir(parents=True)
        for code in ('BUA','S6Q'):shutil.copy2(FIXTURES/(code+'.cif'),cache/(code+'.cif'))
        with patch.dict('os.environ',{'BIOMOL_CACHE_DIR':str(self.root/'cache')}):
            for code in ('BUA','S6Q'):
                data,_=component_definition(code);pdb=self.root/(code+'.complex.pdb');lines=[]
                for index,(name,element,x,y,z) in enumerate(zip(data['_chem_comp_atom.atom_id'],data['_chem_comp_atom.type_symbol'],
                        data['_chem_comp_atom.pdbx_model_Cartn_x_ideal'],data['_chem_comp_atom.pdbx_model_Cartn_y_ideal'],data['_chem_comp_atom.pdbx_model_Cartn_z_ideal']),1):
                    if element=='H':continue
                    lines.append(f'HETATM{index:5d} {name:>4} {code:>3} A   1    {float(x):8.3f}{float(y):8.3f}{float(z):8.3f}  1.00 20.00          {element:>2}')
                pdb.write_text('\n'.join(lines)+'\n');script,output=self.script('TEST_'+code+'_1A')
                script.write_text(f'open {pdb}\nselect :1.A\nselect invert\ndelete selected\ndelete element.H\naddh\nwrite format mol2 #0 {output}\n')
                self.assertTrue(prepare_chimera_ligand(script))
                typing=output.with_suffix('.input.types.py').read_text()
                types=ast.literal_eval(typing.split('atom_types = ')[1].split('\n')[0])
                if code=='BUA':self.assertIn('O3',types.values());self.assertIn('O2',types.values())
                else:
                    self.assertIn('Npl',types.values())

    def test_receptor_hetero_amino_acid_cleanup_preserves_explicit_cofactors(self):
        pdb=self.root/'complex.pdb'
        pdb.write_text('HETATM    1  N   GLY A1531      10.000  10.000  10.000  1.00 20.00           N\n')
        script=self.root/'receptor.com'
        for selection,expected in [('select protein',True),('select protein | :GLY',False)]:
            script.write_text(f'open {pdb}\n{selection}\naddh\nwrite format mol2 #0 {self.root}/receptor.mol2\n')
            receptor_cleanup(script);helper=script.with_suffix('.cleanup.py').read_text()
            candidates=ast.literal_eval(helper.split('isolated_candidates = ')[1].split('\n')[0])
            self.assertEqual(('GLY',1531,'A') in candidates,expected)
            self.assertIn('atom.neighbors',helper)  # Peptide-bound residues survive.


if __name__=='__main__':unittest.main()
