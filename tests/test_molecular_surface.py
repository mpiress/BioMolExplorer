"""Compare the native SES against retained outputs of the official UCSF DMS.

The oracle is independent C code, not a second implementation of the port.
"""
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from biomolexplorer.molecular_surface import generate_surface

FIXTURES = Path(__file__).parent/'fixtures/dms'


def canonical(text):
    # Server completion order is unspecified. Negative zero has no geometry.
    return Counter(line.replace('-0.000', ' 0.000') for line in text.splitlines())


class MolecularSurfaceTests(unittest.TestCase):
    def test_official_dms_golden_surfaces(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder)/'native.dms'
            for case in json.loads((FIXTURES/'manifest.json').read_text()):
                with self.subTest(case=case['name']):
                    expected = gzip.decompress((FIXTURES/f"{case['name']}.dms.gz").read_bytes())
                    self.assertEqual(hashlib.sha256(expected).hexdigest(), case['sha256'])
                    summary = generate_surface(FIXTURES/f"{case['name']}.pdb", output,
                                               density=case['density'], probe_radius=case['probe_radius'])
                    self.assertEqual(canonical(expected.decode()), canonical(output.read_text()))
                    self.assertEqual(summary.atoms+summary.contact_points+summary.saddle_points+summary.concave_points,
                                     len(expected.decode().splitlines()))
                    self.assertGreater(summary.area, 0)

    def test_collinear_and_coincident_atoms_produce_finite_surfaces(self):
        line = (FIXTURES/'single.pdb').read_text().splitlines()[0]
        with tempfile.TemporaryDirectory() as folder:
            pdb, output = Path(folder)/'input.pdb', Path(folder)/'surface.dms'
            for coords in [(0., 2., 4.), (0., 0., 0.)]:
                pdb.write_text(''.join(line[:30]+f'{x:8.3f}'+line[38:]+'\n' for x in coords)+'END\n')
                generate_surface(pdb, output)
                for row in output.read_text().splitlines():
                    values = [float(row[13:21]), float(row[22:30]), float(row[31:39])]
                    self.assertTrue(np.isfinite(values).all())
                self.assertNotIn('nan', output.read_text().lower())

    def test_validation_and_failed_write_preserve_existing_output(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder)/'surface.dms'
            output.write_text('previous successful surface')
            for kwargs in ({'density': 0}, {'density': float('nan')}, {'density': 11},
                           {'probe_radius': .9}, {'probe_radius': float('inf')}):
                with self.assertRaises(ValueError):
                    generate_surface(FIXTURES/'single.pdb', output, **kwargs)
                self.assertEqual(output.read_text(), 'previous successful surface')
            with patch('biomolexplorer.molecular_surface.os.replace', side_effect=OSError('disk')):
                with self.assertRaises(OSError):
                    generate_surface(FIXTURES/'single.pdb', output)
            self.assertEqual(output.read_text(), 'previous successful surface')
            self.assertEqual(list(Path(folder).glob('*.tmp')), [])

    def test_atom_selection_custom_radii_and_normals(self):
        with tempfile.TemporaryDirectory() as folder:
            pdb, output, radii = (Path(folder)/n for n in ('input.pdb', 'surface.dms', 'radii'))
            line = (FIXTURES/'single.pdb').read_text().splitlines()[0]
            pdb.write_text(line.replace('ALA', 'LIG')+'\n'+line.replace('ATOM  ', 'HETATM').replace('ALA', 'HOH')+'\nEND\n')
            with self.assertRaisesRegex(ValueError, 'no eligible'):
                generate_surface(pdb, output)
            radii.write_text('C 2.0\ndefault 1.9\n')
            summary = generate_surface(pdb, output, include_hetero=True, normals=False, radii_path=radii)
            self.assertEqual(summary.atoms, 2)
            rows = output.read_text().splitlines()
            self.assertTrue(any('1A*' in row for row in rows))
            self.assertTrue(any('LIG' in row for row in rows))
            surface = next(row for row in rows if 'SC0' in row)
            self.assertEqual(len(surface), 50)
            self.assertAlmostEqual(float(surface[31:39]), 2.0, places=3)


if __name__ == '__main__':
    unittest.main()
