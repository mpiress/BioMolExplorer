"""Isolated hydrogens retain identity without being mistaken for failures."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors

from biomolexplorer.molecular_parsing import mol_from_smiles
from biomolexplorer.molecule_quality import QualityReport, sanitize_csv
from biomolexplorer.operations import execute_operation


class MolecularParsingTests(unittest.TestCase):
    def test_parser_preserves_default_identity_stereo_isotopes_and_fingerprints(self):
        for smiles in ('[H]','CC.[HH]','CC.[H-].[Na+]','[2H]O','[3H]',
                       '[H][H]','[H]C([H])([H])[H]','[H][C@](F)(Cl)Br','c1cc[nH]c1','C[Hg+]'):
            with self.subTest(smiles=smiles):
                with rdBase.BlockLogs():expected=Chem.MolFromSmiles(smiles)
                actual=mol_from_smiles(smiles)
                self.assertEqual(Chem.MolToSmiles(actual),Chem.MolToSmiles(expected))
                self.assertEqual(actual.GetNumAtoms(),expected.GetNumAtoms())
                self.assertEqual(rdMolDescriptors.CalcMolFormula(actual),rdMolDescriptors.CalcMolFormula(expected))
                self.assertEqual(Chem.RDKFingerprint(actual).ToBitString(),Chem.RDKFingerprint(expected).ToBitString())

    def test_isolated_hydrogens_are_not_repeated_as_warnings_and_invalid_smiles_still_log(self):
        script="""
from biomolexplorer.molecular_parsing import mol_from_smiles
from biomolexplorer.visualizations import molecule_image
from biomolexplorer.network_scaffolds import individual
for smiles in ('CC.[HH]', 'CC.[H-].[Na+]', '[3H]'):
    assert mol_from_smiles(smiles) is not None
    assert molecule_image(smiles)
    assert individual(smiles)['image']
assert mol_from_smiles('invalid') is None
"""
        result=subprocess.run([sys.executable,'-c',script],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertNotIn('not removing hydrogen atom without neighbors',result.stderr)
        self.assertIn('SMILES Parse Error',result.stderr)

    def test_quality_report_explains_retained_hydrogens_without_excluding_compounds(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);source=root/'source.csv'
            source.write_text('molecule_chembl_id,canonical_smiles\nA,CC.[H-].[Na+]\nB,CC.[HH]\n')
            original=source.read_bytes();report=QualityReport('fingerprints')
            cleaned=sanitize_csv(source,root/'clean.csv',report=report)
            self.assertEqual(cleaned.read_text(),source.read_text())
            self.assertEqual(source.read_bytes(),original)
            path=report.write(root)
            saved=json.loads(path.read_text())
            self.assertEqual(saved['excluded_records'],0)
            self.assertEqual(len(saved['notices']),2)
            self.assertEqual(saved['notices'][0]['retained_hydrogen_atoms'],1)
            self.assertIn('não excluído',saved['notices'][0]['reason'])

    def test_real_execution_failure_is_preserved_alongside_hydrogen_notice(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);inputs=root/'inputs';inputs.mkdir()
            (inputs/'compounds.csv').write_text('molecule_chembl_id,canonical_smiles\nA,CC.[H-].[Na+]\n')
            with patch('biomolexplorer.operations._execute_operation',side_effect=RuntimeError('engine failed')):
                with self.assertRaisesRegex(RuntimeError,'engine failed'):
                    execute_operation('fingerprints',{'base_input_path':str(inputs)},root/'out')
            saved=json.loads((root/'out'/'molecule_exclusions.json').read_text())
            self.assertEqual(saved['excluded_records'],0)
            self.assertTrue(saved['notices'])
