import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from biomolexplorer.complex_classification import classify_complex,COVALENT,NON_COVALENT,CovalentComplexError
from biomolexplorer.docking_capabilities import dock6_capabilities,CovalentDocking
from biomolexplorer.structure_selection import read_selected,crystal_molecule


class ComplexClassificationTests(unittest.TestCase):
    ROOT=Path('/home/michel/Downloads/Alzheimer')

    def structure(self,code):
        try:path=next(self.ROOT.rglob(code+'.pdb'))
        except StopIteration:self.skipTest('Real Alzheimer structure fixture unavailable')
        return path,read_selected(path)[0]

    def test_1p0i_attachment_is_explicit_and_backend_gated(self):
        _,structure=self.structure('1P0I');report=classify_complex(structure,['1P0I','BUA',606,'A'])
        self.assertEqual(report['interaction_type'],COVALENT)
        bond=report['attachments'][0]
        self.assertEqual((bond['receptor_residue'],bond['receptor_number'],bond['receptor_atom']),('SER',198,'OG'))
        self.assertEqual((bond['ligand_residue'],bond['ligand_atom']),('BUA','C4'))
        self.assertIn('PDB LINK',bond['evidence']);self.assertIn('PDB CONECT',bond['evidence'])

    def test_7myi_altloc_conect_is_not_a_covalent_bond(self):
        path,structure=self.structure('7MYI');report=classify_complex(structure,['7MYI','ZQS',401,'A'])
        self.assertEqual(report['interaction_type'],NON_COVALENT)
        chemistry=self.ROOT/'.biomolexplorer/cache/calculations/ccd/ZQS.cif'
        if not chemistry.is_file():self.skipTest('ZQS CCD fixture unavailable')
        molecule,mapping,_=crystal_molecule(structure,['7MYI','ZQS',401,'A'],
            {'chemistry_templates':{'ZQS':str(chemistry)}},path.parent)
        self.assertEqual(molecule.GetNumHeavyAtoms(),21);self.assertEqual(len(mapping),21)

    def test_dock_611_never_enables_covalent_template_keyword(self):
        with patch('biomolexplorer.docking_capabilities.dock6_version',return_value='6.11'):
            capabilities=dock6_capabilities('/home/michel/progs/dock6')
        self.assertTrue(capabilities['standard_docking']);self.assertFalse(capabilities['covalent_attach_and_grow'])
        with self.assertRaisesRegex(ValueError,'DOCK >= 6.12'):CovalentDocking(capabilities).conformer_search_type()

    def test_future_dock_version_enables_adapter_without_affecting_611(self):
        with patch('biomolexplorer.docking_capabilities.dock6_version',return_value='6.12'):
            capabilities=dock6_capabilities('/home/michel/progs/dock6')
        self.assertTrue(capabilities['covalent_attach_and_grow'])
        self.assertEqual(CovalentDocking(capabilities).conformer_search_type(),'covalent')


if __name__=='__main__':unittest.main()
