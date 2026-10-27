"""Charge-only migration must never introduce another preparation protocol."""
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from rdkit import Chem

from biomolexplorer.molecular_core import MolecularPreparationCore
from biomolexplorer.molecular_adapters import Dock6Adapter, charge_existing_mol2


class ChargeOnlyTests(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root=Path(self.temporary.name)
        self.addCleanup(patch.stopall)
        patch.dict(os.environ,BIOMOL_PREPARATION_CACHE=str(self.root/'cache'),BIOMOL_REUSE_PREPARATIONS='1').start()
        self.prepared=MolecularPreparationCore(self.root).prepare(Chem.MolFromSmiles('CCO'))
        self.path=self.root/'pose.mol2'
        Dock6Adapter().export(self.prepared,self.path)
        self.original=self.path.read_text()
        patch('biomolexplorer.molecular_adapters.amber_tool',return_value=self.root/'bin/antechamber').start()
        patch('biomolexplorer.molecular_adapters.AM1BCCEngine.identity',return_value={'engine':'test','method':'am1bcc'}).start()

    @staticmethod
    def without_charges(text):
        atoms=False;result=[]
        for line in text.splitlines():
            if line.startswith('@<TRIPOS>'):atoms=line=='@<TRIPOS>ATOM'
            result.append(line.split()[:8]+line.split()[9:] if atoms and not line.startswith('@') and line.strip() else line)
        return result

    def test_preserves_every_record_except_charges_and_reuses_completed_result(self):
        values=np.zeros(self.prepared.molecule.GetNumAtoms())
        values[0]=.12;values[1]=-.12
        with patch('biomolexplorer.molecular_adapters.AM1BCCEngine.calculate',return_value=values) as calculate:
            charge_existing_mol2(self.path)
            self.assertEqual(self.without_charges(self.original),self.without_charges(self.path.read_text()))
            charge_existing_mol2(self.path)
            self.assertEqual(calculate.call_count,1)
            self.path.write_text(self.original)
            charge_existing_mol2(self.path)
            self.assertEqual(calculate.call_count,1)
            self.path.write_text(self.original)
            with patch.dict(os.environ,BIOMOL_REUSE_PREPARATIONS='0'):
                charge_existing_mol2(self.path)
            self.assertEqual(calculate.call_count,2)

    def test_failed_calculation_preserves_original_and_can_resume(self):
        with patch('biomolexplorer.molecular_adapters.AM1BCCEngine.calculate',side_effect=RuntimeError('interrupted')):
            with self.assertRaises(RuntimeError):charge_existing_mol2(self.path)
        self.assertEqual(self.path.read_text(),self.original)
        with patch('biomolexplorer.molecular_adapters.AM1BCCEngine.calculate',return_value=self.prepared.partial_charges) as calculate:
            charge_existing_mol2(self.path)
            self.assertEqual(calculate.call_count,1)

    def test_invalid_charge_vector_is_not_published(self):
        with patch('biomolexplorer.molecular_adapters.AM1BCCEngine.calculate',return_value=[float('nan')]*self.prepared.molecule.GetNumAtoms()):
            with self.assertRaises(ValueError):charge_existing_mol2(self.path)
        self.assertEqual(self.path.read_text(),self.original)

    def test_charge_method_change_cannot_reuse_another_method(self):
        with patch('biomolexplorer.molecular_adapters.AM1BCCEngine.calculate',return_value=self.prepared.partial_charges) as am1:
            charge_existing_mol2(self.path,'am1bcc')
            self.assertEqual(am1.call_count,1)
        with patch('biomolexplorer.molecular_adapters.AntechamberGasteigerEngine.calculate',return_value=self.prepared.partial_charges) as gas:
            charge_existing_mol2(self.path,'gas')
            self.assertEqual(gas.call_count,1)

    def test_duplicate_imported_names_are_mapped_without_renaming_output(self):
        import re
        prefix,atoms=self.original.split('@<TRIPOS>ATOM',1)
        atoms,tail=atoms.split('@<TRIPOS>BOND',1)
        atoms=re.sub(r'(?m)^(\s*[12]\s+)\S+',r'\1DUP',atoms)
        text=prefix+'@<TRIPOS>ATOM'+atoms+'@<TRIPOS>BOND'+tail
        self.assertEqual(text.count('DUP'),2)
        self.path.write_text(text)
        def calculate(molecule,work):
            names=[atom.GetProp('_TriposAtomName') for atom in molecule.GetAtoms()]
            self.assertEqual(len(names),len(set(names)))
            return self.prepared.partial_charges
        with patch('biomolexplorer.molecular_adapters.AM1BCCEngine.calculate',side_effect=calculate):
            charge_existing_mol2(self.path)
        self.assertEqual(self.without_charges(text),self.without_charges(self.path.read_text()))

    def test_changed_geometry_invalidates_only_its_charge_result(self):
        with patch('biomolexplorer.molecular_adapters.AM1BCCEngine.calculate',return_value=self.prepared.partial_charges) as calculate:
            charge_existing_mol2(self.path)
            molecule=self.prepared.molecule
            xyz=molecule.GetConformer().GetAtomPosition(0);xyz.x+=.1
            molecule.GetConformer().SetAtomPosition(0,xyz)
            Dock6Adapter().export(self.prepared,self.path)
            charge_existing_mol2(self.path)
            self.assertEqual(calculate.call_count,2)

    def test_legacy_conversion_routes_am1_to_direct_antechamber(self):
        from biomolexplorer.docking_data import prepare_ligands
        dataset=self.root/'compounds.csv'
        dataset.write_text('molecule_chembl_id,canonical_smiles,conformer_file\nTEST,CCO,'+str(self.path)+'\n')
        def convert(source,destination,**kwargs):
            destination.write_text(self.original)
            return destination
        with patch('biomolexplorer.docking_data.convert_structure',side_effect=convert), \
             patch('biomolexplorer.ligand_cache.preparation_identity',return_value={}), \
             patch('biomolexplorer.molecular_adapters.AM1BCCEngine.calculate',return_value=self.prepared.partial_charges) as calculate:
            outputs=prepare_ligands(dataset,self.root/'output','mol2',charge_type='am1')
            self.assertEqual(len(outputs),1)
            self.assertEqual(calculate.call_count,1)
            prepare_ligands(dataset,self.root/'output','mol2',charge_type='am1')
            self.assertEqual(calculate.call_count,1)
            self.assertFalse(list((self.root/'output').glob('*.com')))


if __name__=='__main__':unittest.main()
