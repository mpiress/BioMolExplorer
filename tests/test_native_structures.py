"""Native structural identity, authoritative chemistry and Amber protocol checks."""
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from rdkit import Chem

from biomolexplorer.structure_selection import read_selected,subset,write_structure,crystal_molecule,gemmi_module
from biomolexplorer.molecular_core import MolecularPreparationCore
from biomolexplorer.molecular_adapters import Dock6Adapter,AutoDockAdapter
from biomolexplorer.amber_preparation import protonation,MINIMIZATION


class StructuralTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.env=patch.dict(os.environ,BIOMOL_PREPARATION_CACHE=str(self.root/'cache/preparations'))
        self.env.start();self.addCleanup(self.env.stop)
        self.gemmi=gemmi_module()

    def test_pdb_and_mmcif_selection_preserve_ids_water_metals_and_chosen_altloc(self):
        text=('ATOM      1  CA AALA A   2A      1.000   2.000   3.000  0.40 20.00           C\n'
              'ATOM      2  CA BALA A   2A      4.000   5.000   6.000  0.60 20.00           C\n'
              'HETATM    3  O   HOH A 300       7.000   8.000   9.000  1.00 20.00           O\n'
              'HETATM    4 ZN    ZN B 301      10.000  11.000  12.000  1.00 20.00          ZN\nEND\n')
        path=self.root/'source.pdb';path.write_text(text)
        pdb,mapping=read_selected(path)
        self.assertEqual(len(mapping),3);self.assertEqual(mapping[0]['altloc'],'B')
        self.assertEqual(mapping[0]['insertion_code'],'A')
        np.testing.assert_array_equal(mapping[0]['coordinates'],[4,5,6])
        write_structure(pdb,self.root/'selected.cif')
        from biomolexplorer.input_validation import validate_file
        from biomolexplorer.ui.file_selection import compatible_file
        validate_file(self.root/'selected.cif','structures')
        self.assertTrue(compatible_file(self.root/'selected.cif',{'structures'}))
        cif,again=read_selected(self.root/'selected.cif')
        self.assertEqual([(a['chain'],a['residue_number'],a['insertion_code'],a['atom_name'],a['coordinates']) for a in mapping],
                         [(a['chain'],a['residue_number'],a['insertion_code'],a['atom_name'],a['coordinates']) for a in again])
        selected=subset(cif,lambda chain,residue:residue.name in ('HOH','ZN'))
        self.assertEqual(sum(len(r) for c in selected[0] for r in c),2)
        with self.assertRaises(ValueError):read_selected(path,{'altloc':'C'})

    def crystal(self):
        prepared=MolecularPreparationCore(self.root).prepare(Chem.MolFromSmiles('CCO'))
        for atom,name in zip(prepared.molecule.GetAtoms(),prepared.atom_names):
            info=Chem.AtomPDBResidueInfo();info.SetName(name.ljust(4));info.SetResidueName('LIG');info.SetResidueNumber(21);info.SetIsHeteroAtom(True)
            atom.SetMonomerInfo(info)
        template=self.root/'LIG.mol2';Dock6Adapter().export(prepared,template)
        pdb=Chem.MolToPDBBlock(prepared.molecule)
        structure=self.gemmi.read_pdb_string(pdb)
        chain=structure[0][0];chain.name='A';residue=chain[0];residue.name='LIG';residue.seqid=self.gemmi.SeqId(21,'B')
        return prepared,structure,{'chemistry_templates':{'LIG':str(template)},'ligand_insertion_code':'B'}

    def test_crystal_topology_mapping_has_zero_displacement(self):
        prepared,structure,config=self.crystal()
        molecule,mapping,path=crystal_molecule(structure,['TEST','LIG',21,'A'],config,self.root)
        np.testing.assert_allclose(molecule.GetConformer().GetPositions(),prepared.coordinates[:3],atol=.00051)
        self.assertEqual(molecule.GetNumAtoms(),3);self.assertEqual(len(mapping),3)
        self.assertEqual(mapping[0]['original']['insertion_code'],'B')
        self.assertEqual(Chem.GetFormalCharge(molecule),0)

    def test_incomplete_chemical_template_fails_instead_of_guessing(self):
        _,structure,config=self.crystal()
        structure[0][0][0][0].name='UNKNOWN'
        with self.assertRaisesRegex(ValueError,'names/count'):crystal_molecule(structure,['TEST','LIG',21,'A'],config,self.root)

    def test_inter_residue_covalent_bonds_require_dedicated_protocol(self):
        _,structure,config=self.crystal()
        structure[0][0][0][0].serial=42
        residue=self.gemmi.Residue();residue.name='SER';residue.seqid=self.gemmi.SeqId(22,' ')
        atom=self.gemmi.Atom();atom.name='OG';atom.element=self.gemmi.Element('O')
        atom.serial=99;atom.pos=self.gemmi.Position(1.4,0,0);residue.add_atom(atom)
        structure[0][0].add_residue(residue)
        structure.add_conect(42,99,1)
        with self.assertRaisesRegex(ValueError,'Covalent crystal complex'):crystal_molecule(structure,['TEST','LIG',21,'A'],config,self.root)

    def test_protonation_preserves_named_states_and_records_automatic_choices(self):
        structure=self.gemmi.Structure();model=self.gemmi.Model('1');chain=self.gemmi.Chain('A')
        for i,name in enumerate(('HIS','HID','HIP','ASH','GLH','CYM','LYN')):
            residue=self.gemmi.Residue();residue.name=name;residue.seqid=self.gemmi.SeqId(i+1,' ')
            atom=self.gemmi.Atom();atom.name='CA';atom.element=self.gemmi.Element('C');atom.pos=self.gemmi.Position(i*4,0,0)
            residue.add_atom(atom);chain.add_residue(residue)
        model.add_chain(chain);structure.add_model(model)
        normalized,states,bonds=protonation(structure,7.4)
        self.assertEqual([r.name for r in normalized[0][0]],['HIE','HID','HIP','ASH','GLH','CYM','LYN'])
        self.assertTrue(all(s['ph']==7.4 and s['method'] for s in states))

    def test_disulfides_keep_pdb_identity_and_use_dense_amber_addresses(self):
        structure=self.gemmi.Structure();model=self.gemmi.Model('1');chain=self.gemmi.Chain('A')
        for number,x in ((10,0.),(30,2.05)):
            residue=self.gemmi.Residue();residue.name='CYS';residue.seqid=self.gemmi.SeqId(number,' ')
            atom=self.gemmi.Atom();atom.name='SG';atom.element=self.gemmi.Element('S');atom.pos=self.gemmi.Position(x,0,0)
            residue.add_atom(atom);chain.add_residue(residue)
        model.add_chain(chain);structure.add_model(model)
        normalized,states,bonds=protonation(structure,7.4)
        self.assertEqual([residue.name for residue in normalized[0][0]],['CYX','CYX'])
        self.assertEqual(bonds,[{'residues':[{'chain':'A','number':10,'insertion_code':''},
                                              {'chain':'A','number':30,'insertion_code':''}],
                                 'amber_residues':[1,2]}])

    def test_minimization_is_explicitly_minimization_not_a_dynamics_step(self):
        self.assertEqual(MINIMIZATION['imin'],1)
        self.assertEqual((MINIMIZATION['ncyc'],MINIMIZATION['maxcyc']),(100,110))
        self.assertEqual((MINIMIZATION['igb'],MINIMIZATION['ntb']),(0,0))
        from biomolexplorer.molecular_core import charge_strategy
        self.assertEqual(charge_strategy('am1','vina','fast',minimize=True)['effective'],'am1bcc')
        self.assertEqual(charge_strategy('gas','vina','fast',minimize=True)['engine'],'antechamber')

    def test_nested_cache_namespaces_do_not_share_lock_inodes(self):
        from biomolexplorer.calculation_cache import CalculationCache,identity
        first=CalculationCache('structure',self.root);second=CalculationCache('charges',self.root)
        # Deliberately collide under the historical cross-namespace stripe rule.
        key='parent';stripe=int(identity(['structure',key])[:2],16)%64
        child=next(str(i) for i in range(1000) if int(identity(['charges',str(i)])[:2],16)%64==stripe)
        with first.locked(key),second.locked(child):
            self.assertEqual(len(list((first.root/'locks').glob('*/*.lock'))),2)

    def test_pose_recovery_requires_embedded_graph_and_preserves_heavy_coordinates(self):
        prepared=MolecularPreparationCore(self.root).prepare(Chem.MolFromSmiles('CCO'))
        path=self.root/'pose.pdbqt';AutoDockAdapter().export(prepared,path)
        recovered,_=MolecularPreparationCore(self.root).load(path)
        expected=Chem.RemoveHs(prepared.molecule).GetConformer().GetPositions()
        actual=Chem.RemoveHs(recovered).GetConformer().GetPositions()
        from scipy.spatial import cKDTree
        self.assertLess(cKDTree(expected).query(actual)[0].max(),.002)
        path.write_text('\n'.join(r for r in path.read_text().splitlines() if not r.startswith('REMARK SMILES')))
        with self.assertRaisesRegex(ValueError,'trusted chemical graph'):MolecularPreparationCore(self.root).load(path)


@unittest.skipUnless(Path('/tmp/biomol-md-validation/bin/sander').is_file(),'AmberTools validation environment unavailable')
class AmberNativeTests(unittest.TestCase):
    def test_minimization_reuses_charge_topology_and_sander_results(self):
        with tempfile.TemporaryDirectory() as folder,patch.dict(os.environ,BIOMOL_PREPARATION_CACHE=str(Path(folder)/'shared/preparations')):
            core=MolecularPreparationCore(folder,'/tmp/biomol-md-validation')
            initial=core.prepare(Chem.MolFromSmiles('CCO'),minimize=True)
            self.assertLessEqual(initial.metadata['minimization']['final_energy'],initial.metadata['minimization']['initial_energy'])
            with patch('biomolexplorer.molecular_adapters.run_command',side_effect=AssertionError('Repeated parameterization')), \
                 patch('biomolexplorer.amber_preparation.run_command',side_effect=AssertionError('Repeated minimization')):
                again=core.prepare(Chem.MolFromSmiles('CCO'),minimize=True)
            np.testing.assert_array_equal(initial.coordinates,again.coordinates)
            from biomolexplorer import amber_preparation
            with patch('biomolexplorer.molecular_adapters.run_command',side_effect=AssertionError('Repeated topology')), \
                 patch.dict(amber_preparation.MINIMIZATION,maxcyc=111), \
                 patch('biomolexplorer.amber_preparation.run_command',wraps=amber_preparation.run_command) as sander:
                core.prepare(Chem.MolFromSmiles('CCO'),minimize=True)
                self.assertEqual(sander.call_count,1)

    def test_receptor_preserves_altloc_and_reuses_ff14sb_parameters(self):
        from biomolexplorer.amber_preparation import receptor_parameters,export_receptor
        source=Path('/home/michel/progs/dock6/install/test/amber_score_1lgu/1lgu.pdb')
        if not source.exists():self.skipTest('DOCK6 protein fixture unavailable')
        with tempfile.TemporaryDirectory() as folder,patch.dict(os.environ,BIOMOL_PREPARATION_CACHE=str(Path(folder)/'shared/preparations')):
            path=Path(folder);structure,_=read_selected(source)
            structure[0][0][0][0].altloc='B'
            parameters,_=receptor_parameters(structure,path,prefix='/tmp/biomol-md-validation')
            report=json.loads(parameters['receptor.parameters.json'].read_text())
            self.assertEqual(report['forcefield'],'ff14SB')
            with patch('biomolexplorer.amber_preparation.run_command',side_effect=AssertionError('Repeated tleap')):
                receptor_parameters(structure,path,prefix='/tmp/biomol-md-validation')
                result=export_receptor(parameters,structure,path/'rec.mol2',path/'rec.pdb')
            self.assertEqual(result[0][0][0][0].altloc,'B')
            text=(path/'rec.mol2').read_text();self.assertIn('\n\n\n@<TRIPOS>ATOM',text)
            self.assertIn('@<TRIPOS>SUBSTRUCTURE',text)
            self.assertEqual(len(report['states']),162)


if __name__=='__main__':unittest.main()
