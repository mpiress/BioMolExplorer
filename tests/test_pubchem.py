"""Offline checks of cross-database exclusions and resumable HTTP retrieval."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import pandas as pd
import requests
from crawlers.pubchem import PubChemSimilarMols


class PubChemTests(unittest.TestCase):
    def test_exclusions_shared_hits_and_cache(self):
        with tempfile.TemporaryDirectory() as directory:
            crawler = PubChemSimilarMols(directory)
            crawler.session = Mock()
            calls = []

            def post(url, data, params, timeout):
                calls.append((url, data))
                if '/fastsimilarity_2d/' in url:
                    result = {'IdentifierList': {'CID': [1, 2, 3, 4]}}
                elif '/smiles/cids/' in url:
                    result = {'IdentifierList': {'CID': [1 if data['smiles'] == 'CCO' else 2]}}
                else:
                    self.assertEqual(data['cid'], '3,4')  # ChEMBL properties never downloaded
                    result = {'PropertyTable': {'Properties': [
                        {'CID': 3, 'SMILES': 'CCC', 'MolecularWeight': '44', 'InChIKey': 'new'},
                        # Different CID but an existing ChEMBL structure
                        {'CID': 4, 'SMILES': 'OCC', 'MolecularWeight': '46', 'InChIKey': 'duplicate'},
                    ]}}
                response = Mock(status_code=200)
                response.json.return_value = result
                return response

            crawler.session.post.side_effect = post
            chembl = pd.DataFrame([
                {'molecule_chembl_id': 'CHEMBL1', 'canonical_smiles': 'CCO'},
                {'molecule_chembl_id': 'CHEMBL2', 'canonical_smiles': 'CCN'},
            ])
            with patch('crawlers.pubchem.time.sleep'):
                compounds = crawler.search(chembl)
                count = len(calls)
                repeated = crawler.search(chembl)
            self.assertEqual(compounds['PubChem_CID'].tolist(), [3])
            self.assertEqual(len(pd.read_csv(Path(directory) / 'matches.csv')), 2)
            self.assertEqual(len(calls), count)
            pd.testing.assert_frame_equal(compounds, repeated)

    def test_no_hits_and_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            for threshold in [0, 101, 75.5]:
                with self.assertRaises(ValueError):
                    PubChemSimilarMols(directory, threshold=threshold)
            crawler = PubChemSimilarMols(directory)
            crawler._request = Mock(return_value={})
            result = crawler.search(pd.DataFrame([{'molecule_chembl_id': 'CHEMBL1', 'canonical_smiles': 'CCO'}]))
            self.assertTrue(result.empty)
            self.assertEqual(list(pd.read_csv(Path(directory) / 'compounds.csv').columns), crawler.COLUMNS)

    def test_http_failure_not_cached(self):
        with tempfile.TemporaryDirectory() as directory:
            crawler = PubChemSimilarMols(directory)
            crawler.session = Mock()
            response = Mock(status_code=503)
            response.raise_for_status.side_effect = requests.HTTPError('service unavailable')
            crawler.session.post.return_value = response
            with self.assertRaisesRegex(RuntimeError, 'consulta PubChem'):
                crawler._request('smiles/cids/JSON', {'smiles': 'CCO'})
            self.assertEqual(list(crawler.cache_path.glob('*.json')), [])

    def test_read_nested_chembl_csv(self):
        with tempfile.TemporaryDirectory() as directory:
            pd.DataFrame([{'molecule_chembl_id': 'CHEMBL1', 'molecule_structures':
                           {'canonical_smiles': 'CCO', 'standard_inchi_key': 'key'}}]).to_csv(
                               Path(directory) / 'CHEMBL1.csv', index=False)
            result = PubChemSimilarMols.read_chembl([directory])
            self.assertEqual(result.iloc[0]['canonical_smiles'], 'CCO')
            self.assertEqual(result.iloc[0]['InChIKey'], 'key')


if __name__ == '__main__':
    unittest.main()
