"""Scientific footprint inputs retain atom identity and reject empty energy tables."""
import tempfile
import unittest
from pathlib import Path
from biomolexplorer.footprints import footprint_receptor, footprint_rows

class FootprintTests(unittest.TestCase):
    def test_atom_groups_preserve_coordinates_charges_and_remap_bonds(self):
        text='''@<TRIPOS>MOLECULE
REC
4 2 2
PROTEIN
USER_CHARGES
@<TRIPOS>ATOM
1 C1 1 2 3 C.3 1 ALA1 0.15
2 C2 4 5 6 C.3 2 GLY2 -0.25
3 H1 7 8 9 H 1 ALA1 0.1
4 H2 10 11 12 H 2 GLY2 0.0
@<TRIPOS>BOND
1 1 3 1
2 2 4 1
@<TRIPOS>SUBSTRUCTURE
1 ALA1 1 RESIDUE 0 A ROOT
2 GLY2 2 RESIDUE 0 A ROOT
'''
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'source.mol2';dest=Path(temp)/'footprint.mol2';source.write_text(text)
            footprint_receptor(source,dest);result=dest.read_text()
            self.assertEqual(source.read_text(),text)
            atomrows=result.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>BOND')[0].splitlines()
            self.assertEqual([r.split()[1] for r in atomrows],['C1','H1','C2','H2'])
            original=text.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>BOND')[0].splitlines()
            expected=[original[i].split()[2:] for i in (0,2,1,3)]
            self.assertEqual([row.split()[2:] for row in atomrows],expected)
            self.assertIn('1 1 2 1\n2 3 4 1',result)
            self.assertIn('2 GLY2 3 RESIDUE',result)

    def test_empty_or_nonfinite_energies_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'energies.txt'
            for text in ('','resname resid vdw_ref es_ref hb_ref vdw_pose es_pose hb_pose\n','ALA1 1 0 0 0 nan 0 0\n'):
                path.write_text(text)
                with self.assertRaises(ValueError):footprint_rows(path)
            path.write_text('resname resid vdw_ref es_ref hb_ref vdw_pose es_pose hb_pose\nALA1 1 -1 2 0 -3 4 0\n')
            self.assertEqual(len(footprint_rows(path)),1)
