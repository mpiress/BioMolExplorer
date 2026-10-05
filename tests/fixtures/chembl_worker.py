"""Fake only the external provider; run the real Python worker and business code."""
import os
import time
from unittest.mock import Mock, patch
import requests


def response(key, rows):
    return Mock(ok=True, status_code=200, json=lambda: {key: rows, 'page_meta': {'next': None}})


def get(url, params=None, **kwargs):
    time.sleep(float(os.environ.get('BIOMOL_TEST_CHEMBL_DELAY', '0')))
    if '/target.json' in url:
        return response('targets', [{'target_chembl_id': 'CHEMBL220', 'pref_name': 'Test target',
                                     'organism': 'Homo sapiens', 'target_type': 'SINGLE PROTEIN', 'target_components': []}])
    if '/activity.json' in url:
        assert params.get('standard_type') == 'IC50' or params.get('standard_type__in') == 'Ki,IC50', params
        return response('activities', [{'activity_id': 1, 'molecule_chembl_id': 'CHEMBL4087364',
            'canonical_smiles': 'CCO', 'standard_value': '10', 'standard_units': 'nM', 'standard_type': 'IC50'}])
    if '/molecule.json' in url:
        if os.environ.get('BIOMOL_TEST_CHEMBL_FAILURE') == 'timeout':
            raise requests.ReadTimeout('Simulated provider read timeout')
        if os.environ.get('BIOMOL_TEST_CHEMBL_FAILURE') == '500':
            return Mock(ok=False, status_code=500)
        return response('molecules', [{'molecule_chembl_id': 'CHEMBL4087364', 'molecule_type': 'Small molecule',
            'natural_product': 1, 'molecule_properties': {'full_mwt': '46'},
            'molecule_structures': {'canonical_smiles': 'CCO'}}])
    if '/similarity/' in url:
        return response('molecules', [])
    raise AssertionError('Unexpected provider request: ' + url)


def post(url, data=None, params=None, **kwargs):
    """Serve PubChem through the same real HTTP boundary as the application."""
    time.sleep(float(os.environ.get('BIOMOL_TEST_CHEMBL_DELAY', '0')))
    if '/fastsimilarity_2d/' in url:
        assert params == {'Threshold': 75, 'MaxRecords': 1000}, params
        body = {'IdentifierList': {'CID': [1, 3]}}
    elif '/smiles/cids/' in url:
        body = {'IdentifierList': {'CID': [1]}}
    elif '/cid/property/' in url:
        assert data['cid'] == '3', data
        body = {'PropertyTable': {'Properties': [{'CID': 3, 'SMILES': 'CCC',
            'InChIKey': 'ATUOYWHBWRKTHZ-UHFFFAOYSA-N', 'MolecularFormula': 'C3H8', 'MolecularWeight': '44.10'}]}}
    else:
        raise AssertionError('Unexpected PubChem request: ' + url)
    return Mock(ok=True, status_code=200, json=lambda: body)


if __name__ == '__main__':
    from biomolexplorer.worker import main
    session = Mock()
    session.headers = {}
    session.get.side_effect = get
    session.post.side_effect = post
    with patch('crawlers.chembl_client.requests.Session', return_value=session):
        raise SystemExit(main())
