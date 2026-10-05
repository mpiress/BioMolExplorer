"""Regression checks using the real wrapper module, without ChEMBL requests."""
import runpy
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import pandas as pd
from wrappers import crawlers
from crawlers.pubchem import PubChemSimilarMols


class CompoundRetrievalTests(unittest.TestCase):
    def test_workflow_calls_new_entry_point(self):
        with patch.object(crawlers, 'retrieve_compounds') as retrieve:
            runpy.run_path(str(ROOT / 'workflow/1-InformationRetrieval/retrieve_compounds.py'), run_name='__main__')
        retrieve.assert_called_once_with(search_term='CHEMBL220', base_output_path='/datasets',
            include_pubchem=True, pubchem_threshold=75, pubchem_max_records=1000)

    def test_retrieval_dispatches_expansion(self):
        names = ['Targets', 'Bioactivity', 'Molecule', 'SimilarMols', 'fileHandling', 'read_filters']
        with patch.multiple(crawlers, **{name: Mock(return_value={}) if name == 'read_filters' else Mock() for name in names}), \
                patch.object(crawlers, 'expand_similar_compounds') as expand:
            crawlers.retrieve_compounds('CHEMBL220', '/datasets', include_pubchem=True)
            expand.assert_called_once_with('CHEMBL220', '/datasets', 75, 1000)

    def test_expansion_saves_neutral_dataset(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            chembl_path = root / 'downloads/ChEMBL/molecules/CHEMBL220'
            chembl_path.mkdir(parents=True)
            pd.DataFrame([{'molecule_chembl_id': 'CHEMBL1', 'molecule_properties': {'full_mwt': 46},
                           'molecule_structures': {'canonical_smiles': 'CCO'}}]).to_csv(
                               chembl_path / 'CHEMBL1.csv', index=False)
            extra = pd.DataFrame([{'molecule_chembl_id': 'PUBCHEM3', 'canonical_smiles': 'CCC',
                                  'molecule_properties': {'full_mwt': 44}, 'source': 'PubChem'}])
            with patch.object(PubChemSimilarMols, 'search', return_value=extra):
                result = crawlers.expand_similar_compounds('CHEMBL220', str(root / 'output'),
                                                          base_input_path=str(root / 'downloads'))
            combined = pd.read_csv(root / 'output/compounds/CHEMBL220/compounds.csv')
            self.assertEqual(combined['molecule_chembl_id'].tolist(), ['CHEMBL1', 'PUBCHEM3'])
            self.assertEqual(combined['source'].tolist(), ['ChEMBL', 'PubChem'])
            self.assertEqual(result['source'].tolist(), ['ChEMBL', 'PubChem'])

    def test_expansion_uses_only_the_curated_input_table(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'selected.csv').write_text('molecule_chembl_id,canonical_smiles\nCHEMBL2,CCC\n')
            (root/'unselected.csv').write_text('molecule_chembl_id,canonical_smiles\nCHEMBL1,CCO\n')
            extra=pd.DataFrame([{'molecule_chembl_id':'PUBCHEM3','canonical_smiles':'CCCC','source':'PubChem'}])
            with patch.object(PubChemSimilarMols,'search',return_value=extra) as search:
                result=crawlers.expand_similar_compounds('CHEMBL220',str(root/'out'),
                    base_input_path=str(root),input_file='selected.csv')
            self.assertEqual(search.call_args.args[0]['molecule_chembl_id'].tolist(),['CHEMBL2'])
            self.assertEqual(result['molecule_chembl_id'].tolist(),['CHEMBL2','PUBCHEM3'])


if __name__ == '__main__':
    unittest.main()
