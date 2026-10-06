"""ChEMBL REST and CHEMBL220/IC50 regressions without relying on live availability."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock,patch
import pandas as pd
from crawlers.chembl_client import ChEMBLClient

class ChEMBLRestTests(unittest.TestCase):
    def response(self,resource,records,next_page=None):
        key={'target':'targets','activity':'activities','molecule':'molecules','similarity':'molecules'}[resource]
        return Mock(ok=True,status_code=200,json=Mock(return_value={key:records,'page_meta':{'next':next_page}}))
    def test_pagination_only_and_no_spore(self):
        session=Mock();session.get.side_effect=[self.response('activity',[{'activity_id':1}],'/chembl/api/data/activity.json?offset=1'),self.response('activity',[{'activity_id':2}])]
        with patch('crawlers.chembl_client.requests.Session',return_value=session):
            query=ChEMBLClient().activity.filter(target_chembl_id='CHEMBL220',standard_type__in=['IC50']).only(['activity_id'])
            self.assertEqual(pd.DataFrame.from_records(query)['activity_id'].tolist(),[1,2])
            self.assertEqual(len(query),2)
        self.assertEqual(session.get.call_count,2)
        self.assertEqual(session.get.call_args_list[0].kwargs['params']['standard_type'],'IC50')
        self.assertNotIn('only', session.get.call_args_list[0].kwargs['params'])
        self.assertFalse(any('spore' in c.args[0] for c in session.get.call_args_list))
    def test_activity_names_with_literal_commas_use_exact_filters_and_shared_limit(self):
        session=Mock();session.get.side_effect=[
            self.response('activity',[{'activity_id':1,'standard_type':'Ki'}]),
            self.response('activity',[{'activity_id':2,'standard_type':'K(p,uu,brain)'}])]
        with patch('crawlers.chembl_client.requests.Session',return_value=session):
            query=ChEMBLClient().activity.filter(target_chembl_id='CHEMBL220',standard_units='nM',
                standard_type__in=['Ki','K(p,uu,brain)','K(p,uu,CSF)']).take(2).only(['activity_id'])
            self.assertEqual(list(query),[{'activity_id':1},{'activity_id':2}])
        self.assertEqual(session.get.call_count,2)
        first,second=session.get.call_args_list
        self.assertEqual(first.kwargs['params']['standard_type'],'Ki')
        self.assertEqual(second.kwargs['params']['standard_type'],'K(p,uu,brain)')
        self.assertEqual(second.kwargs['params']['limit'],1)
        self.assertEqual(second.kwargs['params']['standard_units'],'nM')
        self.assertNotIn('standard_type__in',second.kwargs['params'])

    def test_server_error_is_clear(self):
        session=Mock();session.get.return_value=Mock(ok=False,status_code=500)
        with patch('crawlers.chembl_client.requests.Session',return_value=session),self.assertRaisesRegex(RuntimeError,'HTTP 500'):
            list(ChEMBLClient().target.filter(target_chembl_id='CHEMBL220'))
        session.close.assert_called_once()

    def test_target_identity_query_preserves_organism_and_type_without_server_joins(self):
        session = Mock()
        session.get.return_value = self.response('target', [{'target_chembl_id': 'CHEMBL220',
            'pref_name': 'Acetylcholinesterase', 'organism': 'Homo sapiens', 'target_type': 'SINGLE PROTEIN'}])
        with patch('crawlers.chembl_client.requests.Session', return_value=session):
            client = ChEMBLClient()
            query = client.target.filter(target_chembl_id='CHEMBL220', organism='Homo sapiens',
                target_type__in=['SINGLE PROTEIN']).only(['target_chembl_id', 'pref_name'])
            self.assertEqual(list(query), [{'target_chembl_id': 'CHEMBL220', 'pref_name': 'Acetylcholinesterase'}])
            self.assertEqual(session.get.call_args.kwargs['params'], {'target_chembl_id': 'CHEMBL220'})
            self.assertEqual(list(client.target.filter(target_chembl_id='CHEMBL220', organism='Mus musculus')), [])
            self.assertEqual(list(client.target.filter(target_chembl_id='CHEMBL220', target_type__in=['PROTEIN FAMILY'])), [])

    def test_empty_identifier_selection_never_requests_entire_database(self):
        session = Mock()
        with patch('crawlers.chembl_client.requests.Session', return_value=session):
            self.assertEqual(list(ChEMBLClient().molecule.filter(molecule_chembl_id__in=[])), [])
        session.get.assert_not_called()

    def test_intersecting_identifier_filters_are_not_overwritten(self):
        session = Mock()
        session.get.return_value = self.response('molecule', [])
        with patch('crawlers.chembl_client.requests.Session', return_value=session):
            list(ChEMBLClient().molecule.filter(molecule_chembl_id='CHEMBL1', molecule_chembl_id__in=['CHEMBL2']))
        self.assertEqual(session.get.call_args.kwargs['params']['molecule_chembl_id'], 'CHEMBL1')
        self.assertEqual(session.get.call_args.kwargs['params']['molecule_chembl_id__in'], 'CHEMBL2')

    def test_malformed_provider_response_has_clear_error(self):
        session = Mock()
        session.get.return_value = Mock(ok=True, status_code=200, json=lambda: [])
        with patch('crawlers.chembl_client.requests.Session', return_value=session):
            with self.assertRaisesRegex(RuntimeError, 'Resposta inválida'):
                list(ChEMBLClient().activity.filter(target_chembl_id='CHEMBL220'))

    def test_batched_molecules_preserve_all_filters_and_per_compound_files(self):
        from crawlers.molecules import Molecule
        rows = [{'molecule_chembl_id': f'CHEMBL{index}', 'molecule_type': 'Small molecule',
            'natural_product': index % 2, 'molecule_properties': {'full_mwt': str(index + 45)},
            'molecule_structures': {'canonical_smiles': 'CCO'}} for index in range(1, 28)]
        session = Mock()
        def get(url, params=None, **kwargs):
            ids = str(params.get('molecule_chembl_id__in') or params['molecule_chembl_id']).split(',')
            self.assertLessEqual(len(ids), 25)
            return self.response('molecule', [row for row in rows if row['molecule_chembl_id'] in ids])
        session.get.side_effect = get
        with tempfile.TemporaryDirectory() as temp, patch('crawlers.chembl_client.requests.Session', return_value=session):
            root = Path(temp); (root / 'activity').mkdir()
            pd.DataFrame({'molecule_chembl_id': [row['molecule_chembl_id'] for row in rows] + ['CHEMBL1']}).to_csv(root / 'activity/CHEMBL220.csv', index=False)
            crawler = Molecule(path=str(root / 'molecules'), bioactivity_path=str(root / 'activity'))
            crawler.search({'natural_product': 0, 'molecule_type': 'small molecule', 'molecule_weight': 60})
            saved = {path.stem for path in (root / 'molecules').glob('*.csv')}
            self.assertEqual(saved, {f'CHEMBL{index}' for index in range(2, 16, 2)})
            self.assertEqual(session.get.call_count, 2)
    def test_sample_limit_and_cross_host_pagination(self):
        session=Mock();session.get.return_value=self.response('activity',[{'activity_id':1},{'activity_id':2}],'/chembl/api/data/activity.json?offset=2')
        with patch('crawlers.chembl_client.requests.Session',return_value=session):
            self.assertEqual(len(ChEMBLClient().activity.filter().take(1)),1)
        session.get.assert_called_once()
        session.get.return_value=self.response('activity',[], 'https://example.org/data')
        with patch('crawlers.chembl_client.requests.Session',return_value=session),self.assertRaisesRegex(RuntimeError,'paginação'):
            list(ChEMBLClient().activity.filter())
    def test_natural_product_filter_applies_to_both_choices(self):
        from crawlers.molecules import Molecule
        row={'molecule_chembl_id':'CHEMBL1','molecule_type':'Small molecule','natural_product':0,'molecule_properties':{'full_mwt':'46'},'molecule_structures':{'canonical_smiles':'CCO'}}
        session=Mock();session.get.return_value=self.response('molecule',[row])
        for natural_product in (0,1):
            with self.subTest(natural_product=natural_product),tempfile.TemporaryDirectory() as folder,patch('crawlers.chembl_client.requests.Session',return_value=session):
                crawler=Molecule(path=folder)
                crawler._Molecule__search_mol('CHEMBL1',{},natural_product,None,None)
                self.assertEqual((Path(folder)/'CHEMBL1.csv').exists(),natural_product==0)
    def test_chembl220_ic50_full_retrieval(self):
        from wrappers.crawlers import retrieve_compounds
        activities=[{'activity_id':1,'molecule_chembl_id':'CHEMBL1','canonical_smiles':'CCO','value':'0.01','units':'uM','standard_value':'10','standard_units':'nM','standard_type':'IC50'},
                    {'activity_id':2,'molecule_chembl_id':'CHEMBL2','canonical_smiles':'CCC','value':'8','units':'uM','standard_value':'8000','standard_units':'nM','standard_type':'IC50'}]
        molecule={'molecule_chembl_id':'CHEMBL1','molecule_type':'Small molecule','natural_product':1,'molecule_properties':{'full_mwt':'46.0'},'molecule_structures':{'canonical_smiles':'CCO'}}
        def get(url,params=None,**kwargs):
            if '/target.json' in url: return self.response('target',[{'target_chembl_id':'CHEMBL220','pref_name':'Test target','organism':'Homo sapiens','target_type':'SINGLE PROTEIN','target_components':[]}])
            if '/activity.json' in url:
                self.assertEqual(params['standard_type'],'IC50');return self.response('activity',activities)
            if '/molecule.json' in url:
                self.assertEqual(params['molecule_chembl_id'],'CHEMBL1');return self.response('molecule',[molecule])
            if '/similarity/' in url:return self.response('similarity',[])
            self.fail('Unexpected API request: '+url)
        session=Mock();session.get.side_effect=get
        with tempfile.TemporaryDirectory() as folder,patch('crawlers.chembl_client.requests.Session',return_value=session):
            result=retrieve_compounds('CHEMBL220',folder,include_pubchem=False,chembl_filters={'bioactivity':{'standard_type__in':['IC50'],'standard_units':'nM','max_value_ref':5000}})
            self.assertEqual(result['molecule_chembl_id'].tolist(),['CHEMBL1'])
            self.assertTrue((Path(folder)/'compounds/CHEMBL220/compounds.csv').exists())
            activity=pd.read_csv(Path(folder)/'ChEMBL/bioactivity/CHEMBL220/CHEMBL220.csv')
            self.assertEqual(activity['standard_type'].tolist(),['IC50'])
            self.assertEqual(activity['standard_value'].tolist(),[10])
