"""Native ligand graph/charge/export contracts and optional Amber integration."""
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import numpy.testing
from rdkit import Chem

from biomolexplorer.molecular_core import MolecularPreparationCore,PreparedMolecule,charge_method,charge_strategy,validate_charges
from biomolexplorer.molecular_adapters import Dock6Adapter,AutoDockAdapter,AmberAdapter


class MolecularCoreTests(unittest.TestCase):
    def setUp(self):
        self.folder=tempfile.TemporaryDirectory();self.addCleanup(self.folder.cleanup)
        self.core=MolecularPreparationCore(self.folder.name)

    def test_aliases_and_representative_chemistry_without_external_processes(self):
        self.assertEqual(charge_method('gas'),'gasteiger');self.assertEqual(charge_method('am1'),'am1bcc')
        for smiles in ('CCO','C[NH3+]','CC(=O)[O-]','c1ccccc1','c1ccncc1','CCCl','C[C@H](O)F'):
            with self.subTest(smiles=smiles),patch('subprocess.run',side_effect=AssertionError('Unexpected external process')):
                p=self.core.prepare(Chem.MolFromSmiles(smiles))
                self.assertAlmostEqual(float(p.partial_charges.sum()),Chem.GetFormalCharge(p.molecule),places=6)
                self.assertEqual(p.molecule.GetNumHeavyAtoms(),Chem.MolFromSmiles(smiles).GetNumHeavyAtoms())
                self.assertEqual(Chem.MolToSmiles(Chem.RemoveHs(p.molecule),isomericSmiles=True),Chem.MolToSmiles(Chem.MolFromSmiles(smiles),isomericSmiles=True))
                text=Dock6Adapter().export_string(p)
                self.assertIn('USER_CHARGES',text)
                from openbabel import openbabel as ob
                converter=ob.OBConversion();converter.SetInFormat('mol2');read=ob.OBMol()
                self.assertTrue(converter.ReadString(read,text))
                self.assertEqual(read.NumAtoms(),p.molecule.GetNumAtoms())
                np.testing.assert_allclose([a.GetPartialCharge() for a in ob.OBMolAtomIter(read)],p.partial_charges,atol=1e-7)

    def test_charge_cache_is_shared_across_compound_names_and_respects_disable(self):
        first=self.core.prepare(Chem.MolFromSmiles('CCO'))
        with patch('biomolexplorer.molecular_core.GasteigerEngine.calculate',side_effect=AssertionError('Cache miss')):
            second=self.core.prepare(Chem.MolFromSmiles('CCO'))
        self.assertTrue(second.metadata['charge_cache_hit'])
        np.testing.assert_array_equal(first.partial_charges,second.partial_charges)
        with patch.dict(os.environ,BIOMOL_REUSE_PREPARATIONS='0'):
            p=MolecularPreparationCore(self.folder.name).prepare(Chem.MolFromSmiles('CCO'))
        self.assertFalse(p.metadata['charge_cache_hit'])

    def test_mol2_library_uses_the_selected_graph_not_the_first_molecule(self):
        first=self.core.prepare(Chem.MolFromSmiles('CCO'));second=self.core.prepare(Chem.MolFromSmiles('C[NH3+]'))
        path=Path(self.folder.name)/'library.mol2'
        path.write_text(Dock6Adapter().export_string(first)+Dock6Adapter().export_string(second))
        selected,_=self.core.load(path,index=2)
        self.assertEqual(Chem.GetFormalCharge(selected),1)
        self.assertEqual(selected.GetNumHeavyAtoms(),2)
        with self.assertRaises(ValueError):self.core.load(path,index=3)

    def test_metals_fragments_nonfinite_and_unsupported_minimization_fail_explicitly(self):
        for smiles in ('[Zn+2]','CCO.[Na+]','[CH3]'):
            with self.assertRaises(ValueError):self.core.prepare(Chem.MolFromSmiles(smiles))
        molecule=Chem.MolFromSmiles('CCO')
        for values in ([0,0],[float('nan')]*3,[float('inf')]*3,[1,0,0]):
            with self.assertRaises(ValueError):validate_charges(molecule,values)
        with patch('biomolexplorer.chimera_minimization.minimize_ligand',side_effect=lambda prepared,*args,**kwargs:prepared) as minimize:
            self.core.prepare(molecule,minimize=True)
        minimize.assert_called_once()

    def test_backend_charge_policy_never_substitutes_rdkit_for_dock6(self):
        self.assertEqual(charge_strategy('gas','dock6')['engine'],'antechamber')
        self.assertEqual(charge_strategy('gas','both')['engine'],'antechamber')
        self.assertEqual(charge_strategy('am1','dock6')['effective'],'am1bcc')
        self.assertEqual(charge_strategy('am1','vina')['effective'],'gasteiger')
        self.assertEqual(charge_strategy('am1','vina','preserve')['effective'],'am1bcc')
        path=Path(self.folder.name)/'vina.template';path.write_text('scoring = ad4\n')
        with patch('biomolexplorer.paths.resolve_path',return_value=path):
            strategy=charge_strategy('am1','vina')
        self.assertEqual(strategy['effective'],'am1bcc');self.assertEqual(strategy['engine'],'antechamber')

    def test_carboxylate_bond_orders_and_halogen_names_are_preserved(self):
        for smiles in ('CC(=O)[O-]','CCCl'):
            p=self.core.prepare(Chem.MolFromSmiles(smiles));text=Dock6Adapter().export_string(p)
            self.assertEqual(len(p.atom_names),len(set(p.atom_names)))
            if smiles=='CCCl':self.assertTrue(any(name.startswith('Cl') for name in p.atom_names))
            else:
                bonds=text.split('@<TRIPOS>BOND',1)[1].split('@<TRIPOS>',1)[0]
                self.assertFalse(any(line.split()[-1]=='ar' for line in bonds.splitlines() if line.strip()))
                self.assertTrue(any(line.split()[-1]=='2' for line in bonds.splitlines() if line.strip()))

    def test_meeko_reads_supplied_charges_and_preserves_heavy_geometry(self):
        try:import meeko
        except ImportError:self.skipTest('Meeko is not installed')
        p=self.core.prepare(Chem.MolFromSmiles('CCO'))
        text=AutoDockAdapter().export_string(p)
        atoms=[r for r in text.splitlines() if r.startswith('ATOM')]
        self.assertGreater(len(atoms),p.molecule.GetNumHeavyAtoms())
        self.assertEqual([r[12:16].strip() for r in atoms[:3]],p.atom_names[:3])
        self.assertAlmostEqual(sum(float(r.split()[-2]) for r in atoms),0,places=2)
        np.testing.assert_allclose([[float(r[a:b]) for a,b in ((30,38),(38,46),(46,54))] for r in atoms[:3]],p.coordinates[:3],atol=.00051)


@unittest.skipUnless(__import__('shutil').which('chimera'),'Chimera not installed')
class NativeAmberCoreTests(unittest.TestCase):
    def test_am1bcc_and_gaff_reuse_charges_in_native_tleap(self):
        with tempfile.TemporaryDirectory(prefix='bme-native-core-') as folder:
            root=Path(folder);core=MolecularPreparationCore(root)
            prepared=core.prepare(Chem.MolFromSmiles('CCO'),method='am1')
            self.assertAlmostEqual(float(prepared.partial_charges.sum()),0,places=3)
            with patch('biomolexplorer.molecular_adapters.AM1BCCEngine.calculate',side_effect=AssertionError('Repeated AM1-BCC')):
                reused=core.prepare(Chem.MolFromSmiles('CCO'),method='am1bcc')
            self.assertTrue(reused.metadata['charge_cache_hit'])
            files=AmberAdapter().parameterize(reused,root/'amber',forcefield='gaff')
            self.assertTrue(all(p.is_file() for p in files.values()))
            with patch('biomolexplorer.molecular_adapters.run_command',side_effect=AssertionError('Repeated parameterization')):
                again=AmberAdapter().parameterize(reused,root/'amber',forcefield='gaff')
            self.assertEqual(set(files),set(again))

    @unittest.skipUnless(Path('/home/michel/progs/dock6/install/test/amber_score_1lgu/lig.mol2').is_file(),'DOCK6 example not installed')
    def test_dock6_amber_reuses_verified_native_am1bcc(self):
        from biomolexplorer.md_engine import tools,run_pair
        from biomolexplorer.md_config import DEFAULTS
        from biomolexplorer.stage_cache import file_digest
        with tempfile.TemporaryDirectory(prefix='bme-charge-md-') as folder:
            root=Path(folder)
            with patch.dict(os.environ,BIOMOL_PREPARATION_CACHE=str(root/'cache/preparations')):
                core=MolecularPreparationCore(root)
                settings=dict(DEFAULTS,amber_score_before_md_minimization_cycles=1,amber_score_md_steps=1,amber_score_after_md_minimization_cycles=1)
                dock,found=tools('/home/michel/progs/dock6','',settings);example=dock/'install/test/amber_score_1lgu'
                p=core.prepare(example/'lig.mol2',method='am1bcc');Dock6Adapter().export(p,root/'native.mol2')
                metadata=dict(p.metadata,exports={'mol2':file_digest(root/'native.mol2')})
                result=run_pair(dict(receptor={'path':str(example/'1lgu.pdb')},ligand={'path':str(root/'native.mol2'),
                    'compound':{'preparation_metadata':json.dumps(metadata)}}),root/'md',dock,found,settings)
                protocol=json.loads((root/'md/protocol.json').read_text())
                self.assertEqual(protocol['preparation']['charge_method'],'verified prepared AM1-BCC')
                self.assertTrue(np.isfinite(result['amber_score']))


class NativeCandidateIntegrationTests(unittest.TestCase):
    @unittest.skipUnless(__import__('shutil').which('chimera'),'Chimera not installed')
    def test_both_exports_incremental_reuse_and_isolated_invalid_chemistry(self):
        try:import meeko
        except ImportError:self.skipTest('Meeko is not installed')
        from biomolexplorer.docking_data import write_csv,read_compounds
        from biomolexplorer.docking_preparation import prepare_candidates
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);dataset=root/'input.csv'
            write_csv(dataset,[dict(molecule_chembl_id=k,canonical_smiles=s) for k,s in [('A','CCO'),('B','C[NH3+]'),('BAD','not-a-molecule')]])
            options=dict(preparation_backend='python',ambertools_path='/tmp/biomol-md-validation',ligand={'minimize':False,'charge_type':'gas'})
            with patch.dict(os.environ,BIOMOL_PREPARATION_CACHE=str(root/'cache/preparations')):
                first=read_compounds(prepare_candidates(dataset,root/'first',options,engines='both'))
                self.assertEqual([row['molecule_chembl_id'] for row in first],['A','B'])
                self.assertTrue((root/'first/preparation_failures.csv').is_file())
                with patch('biomolexplorer.molecular_adapters.AntechamberGasteigerEngine.calculate',side_effect=AssertionError('Repeated charges')):
                    second=read_compounds(prepare_candidates(dataset,root/'second',options,engines='dock6'))
                self.assertEqual(len(second),2)
                for row in second:self.assertEqual(json.loads(row['preparation_metadata'])['charge_method'],'gasteiger')
                write_csv(dataset,first[:1])
                changed=dict(options,ligand={'minimize':False,'charge_type':'am1'})
                third=read_compounds(prepare_candidates(dataset,root/'changed-model',changed,engines='dock6'))
                self.assertEqual(json.loads(third[0]['preparation_metadata'])['charge_method'],'am1bcc')



@unittest.skipUnless(__import__('shutil').which('chimera'),'Chimera not installed')
class ChimeraMinimizationIntegrationTests(unittest.TestCase):
    def test_minimization_preserves_graph_charges_and_reuses_completed_result(self):
        with tempfile.TemporaryDirectory(prefix='biomol-chimera-min-') as folder:
            root=Path(folder);core=MolecularPreparationCore(root)
            initial=core.prepare(Chem.MolFromSmiles('CCO'),method='am1bcc')
            minimized=core.prepare(Chem.MolFromSmiles('CCO'),method='am1bcc',minimize=True)
            np.testing.assert_allclose(minimized.partial_charges,initial.partial_charges,atol=1e-8)
            self.assertEqual(minimized.atom_names,initial.atom_names)
            self.assertTrue(np.isfinite(minimized.coordinates).all())
            self.assertEqual(minimized.metadata['minimization']['engine'],'Chimera MMTK')
            with patch('biomolexplorer.chimera_minimization.run_command',side_effect=AssertionError('Minimization repeated')):
                reused=core.prepare(Chem.MolFromSmiles('CCO'),method='am1bcc',minimize=True)
            np.testing.assert_allclose(reused.coordinates,minimized.coordinates,atol=1e-8)


if __name__=='__main__':unittest.main()
