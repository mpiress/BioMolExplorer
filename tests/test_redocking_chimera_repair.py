"""Structural DockPrep must precede charges without implicit recharging."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace

from biomolexplorer.redocking_config import configure_template
from biomolexplorer.paths import SOURCE_ROOT


class ReceptorRepairTests(unittest.TestCase):
    def test_cleanup_precedes_hydrogens_charges_and_minimization(self):
        source=(SOURCE_ROOT/'biomolexplorer/resources/chimera/prepare_receptor.template').read_text()
        script=configure_template(source,'prepare_receptor.template',{'cofactors':['FAD'],'receptor':{'charge_type':'gas'}})
        self.assertLess(script.index('runscript biomol_receptor_cleanup.py FAD'),script.index('addh'))
        self.assertLess(script.index('addh'),script.index('addcharge'))
        self.assertIn('method gas',script)
        self.assertIn('minimize spec #0 prep false',script)

    def test_generated_helper_repairs_structure_without_charging_or_removing_cofactors(self):
        from caad.docking import Docking
        with tempfile.TemporaryDirectory() as folder:
            dock=Docking(complex_input_path=folder,output_path=folder)
            dock.generate_docking_script('src/scripts/chimera/prepare_receptor.template',str(Path(folder)/'prepare.com'),
                input_complex='complex',receptor='rec',_preparation={'cofactors':['FAD']})
            helper=(Path(folder)/'biomol_receptor_cleanup.py').read_text()
            compile(helper,'cleanup','exec')
            self.assertIn('addHFunc=None, addCharges=False',helper)
            self.assertIn('delLigands=False',helper)
            self.assertIn('delAltLocs=True',helper)
            self.assertIn('residue.type in arguments',helper)

    def test_helper_removes_free_amino_acid_additive_but_keeps_polymer_and_explicit_cofactor(self):
        import sys
        from caad.docking import Docking
        class Residue:
            def __init__(self,name,het=True):self.type=name;self.isHet=het;self.atoms=[]
            def oslIdent(self):return self.type
        gly=Residue('GLY');polymer=Residue('GLY',False);cofactor=Residue('ALA');linked=Residue('MSE')
        linked.atoms=[SimpleNamespace(neighbors=[SimpleNamespace(residue=polymer)])]
        model=SimpleNamespace(residues=[gly,polymer,cofactor,linked]);removed=[];model.deleteResidue=removed.append
        fake=SimpleNamespace(Molecule=object,openModels=SimpleNamespace(list=lambda **kw:[model]),replyobj=SimpleNamespace(info=lambda message:None))
        with tempfile.TemporaryDirectory() as folder:
            Docking(complex_input_path=folder,output_path=folder).generate_docking_script(
                'src/scripts/chimera/prepare_receptor.template',str(Path(folder)/'prepare.com'),
                input_complex='complex',receptor='rec',_preparation={})
            source=(Path(folder)/'biomol_receptor_cleanup.py').read_text()
        prep=SimpleNamespace(prep=lambda *a,**kw:None)
        with patch.dict(sys.modules,{'chimera':fake,'chimera.resCode':SimpleNamespace(protein3to1={'ALA':'A','GLY':'G','MSE':'M'}),'DockPrep':prep}):
            exec(compile(source,'cleanup','exec'),{'arguments':['ALA']})
        self.assertEqual(removed,[gly])


class Dock6RedockingTests(unittest.TestCase):
    def test_dock6_pose_updates_shared_metadata_checkpoint_and_viewer_pose(self):
        from biomolexplorer.dock6_redocking import redock
        import pandas as pd
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);prepared=root/'Prepared';prepared.mkdir();output=root/'out';output.mkdir()
            prefix='1ABC_LIG_1A';(prepared/(prefix+'.lig.mol2')).write_text('@<TRIPOS>MOLECULE\nC\n1 0 0 0 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 C 0 0 0 C.3 1 LIG 0.0\n@<TRIPOS>BOND\n')
            (prepared/(prefix+'.lig.pdbqt')).write_text('reference')
            def export():
                pose=output/'dock6'/prefix/'results/flex'/(prefix+'_scored.mol2')
                pose.parent.mkdir(parents=True);pose.write_text('pose')
            def convert(source,target):Path(target).write_text('viewer pose')
            checkpoints=[]
            with patch('caad.docking.Dock6') as dock,patch('biomolexplorer.docking_data.binding_center',return_value=[1,2,3]),patch('biomolexplorer.docking_data.convert_structure',side_effect=convert),patch('biomolexplorer.dock6_redocking.heavy_atom_rmsd',return_value=0.7):
                dock.return_value.export_results.side_effect=export
                redock([['1ABC','LIG',1,'A',2.0]],prepared,output,'/dock6','flex',lambda r,v:checkpoints.append((r,v)))
                dock.return_value.prepare_site.assert_called_once()
                self.assertEqual(dock.call_args.kwargs['binding_site_center'],[1,2,3])
            self.assertTrue((output/(prefix+'.lig.pdbqt')).is_file())
            self.assertEqual(checkpoints[0][1],0.7)
            self.assertEqual(pd.read_csv(root/'pdb_codes.csv')['RMSD'].tolist(),[0.7])

    def test_runtime_manifest_does_not_resolve_external_ambertools(self):
        from biomolexplorer.stage_cache import runtime_manifest
        with patch('biomolexplorer.preparation_tools.amber_tool',side_effect=AssertionError('Unexpected AmberTools')):
            tools=runtime_manifest('redocking',{'preparation_pairs':{'1ABC|LIG|1|A':{'preparation_backend':'auto'}}})
        self.assertIn('chimera',tools);self.assertIn('vina',tools)

    def test_engine_is_validated_and_exposed_in_operation_schema(self):
        from biomolexplorer.catalog import operation_fields
        from biomolexplorer.operations import validate_operation
        fields={f['name']:f for f in operation_fields('redocking')}
        self.assertIn('redocking_engine',fields)
        self.assertEqual(fields['redocking_engine']['default'],'vina')
        with self.assertRaisesRegex(ValueError,'Motor'):
            validate_operation('redocking',{'base_input_path':'/tmp','target':'Target','pdb_codes':[['1ABC','LIG',1,'A']], 'redocking_engine':'bad'})


class ChimeraMol2Tests(unittest.TestCase):
    def test_carboxylate_resonance_keeps_atoms_and_charges_and_is_idempotent(self):
        from biomolexplorer.chimera_mol2 import normalize_carboxylates
        from biomolexplorer.dock6_validation import validate_mol2
        from rdkit import Chem
        text='@<TRIPOS>MOLECULE\nacetate\n7 6 1 0 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 LIG -0.1\n2 C2 1.5 0 0 C.2 1 LIG 0.7\n3 O1 2 1 0 O.co2 1 LIG -0.8\n4 O2 2 -1 0 O.co2 1 LIG -0.8\n5 H1 -1 0 0 H 1 LIG 0.0\n6 H2 0 1 0 H 1 LIG 0.0\n7 H3 0 0 1 H 1 LIG 0.0\n@<TRIPOS>BOND\n1 1 2 1\n2 2 3 2\n3 2 4 2\n4 1 5 1\n5 1 6 1\n6 1 7 1\n'
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'acetate.mol2';path.write_text(text)
            self.assertEqual(normalize_carboxylates(path),1)
            self.assertEqual(normalize_carboxylates(path),0)
            self.assertEqual(text.split('@<TRIPOS>BOND')[0],path.read_text().split('@<TRIPOS>BOND')[0])
            molecule=validate_mol2(path)
            self.assertEqual(Chem.GetFormalCharge(molecule),-1)
            self.assertTrue(all(a.GetNumRadicalElectrons()==0 for a in molecule.GetAtoms()))


class ChimeraChargeRepairTests(unittest.TestCase):
    def test_ccd_confirmed_nitrone_changes_types_without_changing_charge_method(self):
        from biomolexplorer.chimera_charge_repair import charge_commands
        from rdkit import Chem
        molecule=Chem.MolFromSmiles('C=[N+]([O-])C')
        for atom in molecule.GetAtoms():atom.SetProp('_TriposAtomName',atom.GetSymbol()+str(atom.GetIdx()))
        with patch('biomolexplorer.structure_selection.component_path',return_value=Path('/chemistry/LIG.cif')),patch('biomolexplorer.structure_selection.chemical_template',return_value=(molecule,{})):
            commands=charge_commands('/tmp',['1ABC','LIG',5,'A'],{})
        self.assertEqual(len(commands),4)
        self.assertIn('idatmType N2+',commands[2]);self.assertIn('idatmType O3-',commands[3])
        self.assertTrue(all(':5.A@' in c for c in commands))

    def test_chimera_charge_failure_retries_once_with_explicit_types(self):
        from caad.docking import Docking
        from biomolexplorer.processes import ScientificToolError
        with tempfile.TemporaryDirectory() as folder:
            dock=Docking(complex_input_path=folder,output_path=folder)
            name='prepare_ligand_1ABC_LIG_5A.com';path=Path(folder)/name
            path.write_text('open complex.pdb\naddh\naddcharge all spec #0 chargeModel 14sb method am1\nclose all\n')
            dock._chimera_ligand_settings={name:(['1ABC','LIG',5,'A'],{})}
            error=ScientificToolError('Failure running ANTECHAMBER','TOOL_REPORTED_ERROR','Inspect log')
            observed=[]
            def run(*args):
                observed.append(path.read_text())
                if len(observed)==1:raise error
                return True
            with patch.object(dock,'perform_subprocess',side_effect=run),patch('biomolexplorer.chimera_charge_repair.charge_commands',return_value=['setattr a idatmType N2+ :5.A@N1']):
                self.assertTrue(dock.prepare_on_chimera(name))
            self.assertEqual(len(observed),2)
            self.assertIn('idatmType N2+',observed[1])
            self.assertIn('method am1',observed[1])
            self.assertNotIn('method gas',observed[1])

    def test_available_chemical_definition_changes_cache_identity_without_network(self):
        from biomolexplorer.chimera_charge_repair import local_definition
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);(root/'chemistry').mkdir();path=root/'chemistry/LIG.cif'
            self.assertEqual(local_definition(root,'LIG',{}),{})
            path.write_text('first');first=local_definition(root,'LIG',{})
            path.write_text('changed');self.assertNotEqual(first,local_definition(root,'LIG',{}))


class CanonicalExportTests(unittest.TestCase):
    def test_ccd_repairs_nitro_bond_orders_and_preserves_prepared_atom_data(self):
        from rdkit import Chem
        from biomolexplorer.chimera_mol2 import canonicalize_component
        from biomolexplorer.dock6_validation import validate_mol2
        import json
        template=Chem.MolFromSmiles('C[N+](=O)[O-]')
        for atom,name in zip(template.GetAtoms(),('C1','N1','O1','O2')):atom.SetProp('_TriposAtomName',name)
        text='@<TRIPOS>MOLECULE\nnitromethane\n7 6 1 0 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 LIG 0.1\n2 N1 1.5 0 0 N.pl3 1 LIG 0.6\n3 O1 2 1 0 O.2 1 LIG -0.35\n4 O2 2 -1 0 O.2 1 LIG -0.35\n5 H1 -1 0 0 H 1 LIG 0.0\n6 H2 0 1 0 H 1 LIG 0.0\n7 H3 0 0 1 H 1 LIG 0.0\n@<TRIPOS>BOND\n1 1 2 1\n2 2 3 2\n3 2 4 2\n4 1 5 1\n5 1 6 1\n6 1 7 1\n'
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'input.mol2';path.write_text(text);ccd=Path(folder)/'LIG.cif';ccd.write_text('trusted')
            with patch('biomolexplorer.structure_selection.component_path',return_value=ccd),patch('biomolexplorer.structure_selection.chemical_template',return_value=(template,{})):
                canonicalize_component(path,folder,['1ABC','LIG',1,'A'],{})
            self.assertEqual(path.read_text().split('@<TRIPOS>BOND')[0],text.split('@<TRIPOS>BOND')[0])
            molecule=validate_mol2(path)
            self.assertEqual(Chem.GetFormalCharge(molecule),0)
            self.assertTrue(all(a.GetNumRadicalElectrons()==0 for a in molecule.GetAtoms()))
            self.assertEqual(json.loads(path.with_suffix('.export.json').read_text())['scope'],'complete_component')

    def test_oxime_hydroxyl_type_is_recovered_from_ccd(self):
        from rdkit import Chem
        from biomolexplorer.chimera_charge_repair import charge_commands
        molecule=Chem.MolFromSmiles('C=NO')
        for atom in molecule.GetAtoms():atom.SetProp('_TriposAtomName',atom.GetSymbol()+str(atom.GetIdx()))
        with patch('biomolexplorer.structure_selection.component_path',return_value=Path('/chemistry/LIG.cif')),patch('biomolexplorer.structure_selection.chemical_template',return_value=(molecule,{})):
            commands=charge_commands('/tmp',['1ABC','LIG',5,'A'],{})
        self.assertIn('setattr a idatmType O3 :5.A@O2',commands)

    def test_dock6_rmsd_uses_fixed_receptor_frame(self):
        from biomolexplorer.dock6_redocking import heavy_atom_rmsd
        with tempfile.TemporaryDirectory() as folder:
            ref=Path(folder)/'reference.mol2';pose=Path(folder)/'pose.mol2'
            text='@<TRIPOS>MOLECULE\nC\n1 0 0 0 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 C 0 0 0 C.3 1 LIG 0.0\n@<TRIPOS>BOND\n'
            ref.write_text(text);pose.write_text(text.replace('1 C 0 0 0','1 C 1 0 0'))
            self.assertAlmostEqual(heavy_atom_rmsd(ref,pose),1.0)

    def test_ligand_alternates_are_cleaned_before_protonation_without_charge_assignment(self):
        from caad.docking import Docking
        with tempfile.TemporaryDirectory() as folder:
            dock=Docking(complex_input_path=folder,output_path=folder)
            script=Path(folder)/'ligand.com'
            dock.generate_docking_script('src/scripts/chimera/prepare_ligand.template',str(script),input_complex='input',resnum=1,chain='A',charge_type='am1',output_ligand='ligand',_preparation={})
            self.assertLess(script.read_text().index('runscript biomol_ligand_cleanup.py'),script.read_text().index('addh'))
            helper=(Path(folder)/'biomol_ligand_cleanup.py').read_text()
            self.assertIn('delAltLocs=True',helper);self.assertIn('addCharges=False',helper)


class PhosphorylValidationTests(unittest.TestCase):
    def test_neutral_pentavalent_phosphorus_does_not_acquire_formal_charge(self):
        from rdkit import Chem
        from biomolexplorer.dock6_validation import validate_mol2
        molecule=Chem.AddHs(Chem.MolFromSmiles('CP(=O)(C)C'))
        lines=['@<TRIPOS>MOLECULE','phosphoryl',f'{molecule.GetNumAtoms()} {molecule.GetNumBonds()} 1 0 0','SMALL','USER_CHARGES','@<TRIPOS>ATOM']
        for atom in molecule.GetAtoms():
            symbol=atom.GetSymbol();kind={'C':'C.3','O':'O.2','P':'P.3','H':'H'}[symbol]
            lines.append(f'{atom.GetIdx()+1} {symbol}{atom.GetIdx()+1} {atom.GetIdx()} 0 0 {kind} 1 LIG 0.0')
        lines.append('@<TRIPOS>BOND')
        for bond in molecule.GetBonds():lines.append(f'{bond.GetIdx()+1} {bond.GetBeginAtomIdx()+1} {bond.GetEndAtomIdx()+1} {int(bond.GetBondTypeAsDouble())}')
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'phosphoryl.mol2';path.write_text('\n'.join(lines)+'\n')
            self.assertEqual(Chem.GetFormalCharge(validate_mol2(path)),0)


class PreparedPhosphorusStateTests(unittest.TestCase):
    def test_export_keeps_prepared_p_oh_bonds_for_an_observed_fragment(self):
        import json
        from rdkit import Chem
        from biomolexplorer.chimera_mol2 import canonicalize_component
        from biomolexplorer.dock6_validation import validate_mol2
        template=Chem.MolFromSmiles('COP(=O)(O)O')
        for atom,name in zip(template.GetAtoms(),('C1','O1P','P','O4P','O3P','O2P')):atom.SetProp('_TriposAtomName',name)
        text='@<TRIPOS>MOLECULE\nobserved\n10 9 1 0 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 LIG 0.0\n2 O1P 1 0 0 O.3 1 LIG 0.0\n3 P 2 0 0 P.3 1 LIG 0.0\n4 O4P 3 1 0 O.3 1 LIG 0.0\n5 O3P 3 -1 0 O.3 1 LIG 0.0\n6 H1 -1 0 0 H 1 LIG 0.0\n7 H2 0 1 0 H 1 LIG 0.0\n8 H3 0 0 1 H 1 LIG 0.0\n9 H4 4 1 0 H 1 LIG 0.0\n10 H5 4 -1 0 H 1 LIG 0.0\n@<TRIPOS>BOND\n1 1 2 1\n2 2 3 1\n3 3 4 1\n4 3 5 1\n5 1 6 1\n6 1 7 1\n7 1 8 1\n8 4 9 1\n9 5 10 1\n'
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=root/'fragment.mol2';path.write_text(text);ccd=root/'CCD.cif';ccd.write_text('definition')
            with patch('biomolexplorer.structure_selection.component_path',return_value=ccd),patch('biomolexplorer.structure_selection.chemical_template',return_value=(template,{})):
                canonicalize_component(path,folder,['1ABC','LIG',1,'A'],{})
            mol=validate_mol2(path)
            self.assertEqual(Chem.GetFormalCharge(mol),0)
            self.assertFalse(any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()))
            metadata=json.loads(path.with_suffix('.export.json').read_text())
            self.assertEqual(metadata['missing_heavy_atoms'],['O2P'])
            self.assertTrue(metadata['prepared_phosphorus_bonds_preserved'])
            self.assertEqual(mol.GetNumAtoms(),10)


class PrimeAtomNameTests(unittest.TestCase):
    def test_ccd_prime_names_are_preserved_in_chimera_atom_selectors(self):
        from rdkit import Chem
        from biomolexplorer.chimera_charge_repair import charge_commands
        molecule=Chem.MolFromSmiles('CC')
        for atom,name in zip(molecule.GetAtoms(),("C1'","C2'")):atom.SetProp('_TriposAtomName',name)
        with patch('biomolexplorer.structure_selection.component_path',return_value=Path('/chemistry/LIG.cif')),patch('biomolexplorer.structure_selection.chemical_template',return_value=(molecule,{})):
            commands=charge_commands('/tmp',['1ABC','LIG',5,'A'],{})
        self.assertIn("setattr a idatmType C3 :5.A@C1'",commands)
