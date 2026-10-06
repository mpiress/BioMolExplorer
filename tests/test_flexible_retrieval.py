"""Search contracts, downloads and output handoffs without live providers."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch
from types import SimpleNamespace

import pandas as pd
from biomolexplorer.retrieval import identifiers, collection_name, target_criteria
from biomolexplorer.operations import validate_operation
from biomolexplorer.catalog import new_stage
from crawlers.chembl_client import ChEMBLClient
from crawlers.complex import PDBComplex
from wrappers.crawlers import retrieve_compounds, load_pdb

ATOM='ATOM      1  C   ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n'
LIGAND='HETATM    2  C   LIG A   2       3.000   0.000   0.000  1.00  0.00           C\n'
MOLECULE={'molecule_chembl_id':'CHEMBL25','pref_name':'ASPIRIN','natural_product':0,
          'molecule_type':'Small molecule','molecule_properties':{'full_mwt':'180.16'},
          'molecule_structures':{'canonical_smiles':'CC(=O)Oc1ccccc1C(=O)O'}}


class FlexibleRetrievalTests(unittest.TestCase):
    def test_identifier_normalization_and_uniprot_detection(self):
        self.assertEqual(identifiers('1crn,1CRN; 4hhb','pdb'),['1CRN','4HHB'])
        self.assertEqual(target_criteria('P00533'),{'target_components__accession__in':['P00533']})
        self.assertEqual(target_criteria('kinase2'),{'pref_name__icontains':'kinase2'})
        self.assertEqual(target_criteria('CHEMBL220, CHEMBL240','target_id'),{'target_chembl_id__in':['CHEMBL220','CHEMBL240']})
        for value,kind in [('evil/path','pdb'),('notuniprot','uniprot'),('CHEMBLX','chembl')]:
            with self.assertRaises(ValueError):identifiers(value,kind)
        self.assertNotIn('/',collection_name('F/C=C/F'))
        self.assertEqual(collection_name('Meu Alvo'),'MeuAlvo')

    def test_pdb_validation_requires_a_search_criterion_not_an_enzyme(self):
        validate_operation('retrieve_structures',{'pdb_ids':'1CRN'})
        validate_operation('retrieve_structures',{'pdb_query':'kinase'})
        validate_operation('retrieve_structures',{'organism':['Mus musculus']})
        with self.assertRaises(ValueError):validate_operation('retrieve_structures',{})
        with self.assertRaises(ValueError):validate_operation('retrieve_structures',{'pdb_ids':'../bad'})
        with self.assertRaises(ValueError):validate_operation('retrieve_structures',{'max_records':0,'pdb_query':'kinase'})
        validate_operation('retrieve_compounds',{'search_term':'F/C=C/F','search_mode':'substructure'})
        with self.assertRaises(ValueError):validate_operation('retrieve_compounds',{'search_term':'../bad'})
        with self.assertRaises(ValueError):validate_operation('retrieve_compounds',{'search_term':'CHEMBLX','search_mode':'molecule_id'})

    def test_pdb_query_combines_text_ids_and_optional_filters_and_limits_iterator(self):
        crawler=object.__new__(PDBComplex);crawler.logger=Mock()
        with patch('crawlers.complex.TextQuery') as text, patch('crawlers.complex.AttributeQuery') as attr:
            query=Mock();query.__iand__=Mock(return_value=query)
            crawler._search=Mock(return_value=['1CRN','4HHB'])
            text.return_value=query;attr.return_value=query
            result=crawler.get_pdb_ids_with_filters({'pdb_query':'kinase','pdb_ids':['1CRN'], 'max_records':2})
            self.assertEqual(result,['1CRN','4HHB'])
            attr.assert_called_once_with('rcsb_id','in',['1CRN'])
            crawler._search.assert_called_once_with(query,2)

    def test_rcsb_pagination_has_timeouts_and_honors_limit(self):
        crawler=object.__new__(PDBComplex)
        query=Mock();query.to_json.return_value='{"type":"terminal"}'
        session=Mock();session.post.side_effect=[
            Mock(status_code=200,json=lambda:{'total_count':3,'result_set':[{'identifier':'1CRN'}]}),
            Mock(status_code=200,json=lambda:{'total_count':3,'result_set':[{'identifier':'4HHB'}]})]
        with patch('requests.Session',return_value=session):
            self.assertEqual(crawler._search(query,2),['1CRN','4HHB'])
        self.assertEqual(session.post.call_args.kwargs['json']['request_options']['paginate'],{'start':1,'rows':1})
        self.assertEqual(session.post.call_args.kwargs['timeout'],(5,30))
        self.assertTrue(crawler._last_search['limit_reached'])
        session.close.assert_called_once()

    def test_pdb_download_reports_partial_failure_and_retains_valid_metadata(self):
        with tempfile.TemporaryDirectory() as folder:
            crawler=PDBComplex(folder)
            crawler.get_pdb_ids_with_filters=Mock(return_value=['1CRN','4HHB'])
            session=Mock()
            def response(url,**kwargs):
                if '4HHB' in url:raise RuntimeError('unavailable')
                return SimpleNamespace(text=ATOM+LIGAND+'END\n',raise_for_status=lambda:None)
            session.get.side_effect=response
            with patch('requests.Session',return_value=session):
                report=crawler.get_pdb_files({'pdb_ids':['1CRN','4HHB']})
            self.assertEqual([r['status'] for r in report['downloads']],['downloaded','failed'])
            self.assertTrue((Path(folder)/'1CRN.pdb').is_file())
            rows=pd.read_csv(Path(folder)/'pdb_codes.csv')
            self.assertEqual(rows['LIGAND'].tolist(),['LIG'])
            self.assertEqual(json.loads((Path(folder)/'retrieval_report.json').read_text())['matched'],['1CRN','4HHB'])

    def test_empty_pdb_search_is_a_clear_failure(self):
        with tempfile.TemporaryDirectory() as folder:
            crawler=PDBComplex(folder);crawler.get_pdb_ids_with_filters=Mock(return_value=[])
            with self.assertRaisesRegex(ValueError,'Nenhuma estrutura'):crawler.get_pdb_files({'pdb_query':'unknown'})
            self.assertTrue((Path(folder)/'retrieval_report.json').is_file())

    def test_direct_molecule_search_bypasses_bioactivities_and_preserves_compound_contract(self):
        with tempfile.TemporaryDirectory() as folder:
            resource=Mock();resource.filter.return_value.take.return_value=[MOLECULE]
            client=SimpleNamespace(molecule=resource)
            with patch('crawlers.chembl_client.ChEMBLClient',return_value=client),patch('wrappers.crawlers.Bioactivity') as activity:
                result=retrieve_compounds('CHEMBL25',folder,search_mode='molecule_id',max_records=5)
            activity.assert_not_called()
            resource.filter.assert_called_once_with(molecule_chembl_id__in=['CHEMBL25'])
            resource.filter.return_value.take.assert_called_once_with(5)
            self.assertEqual(result['molecule_chembl_id'].tolist(),['CHEMBL25'])
            self.assertTrue((Path(folder)/'compounds/CHEMBL25/compounds.csv').is_file())
            self.assertFalse(json.loads((Path(folder)/'retrieval_report.json').read_text())['activity_evidence'])

    def test_structural_routes_encode_smiles_and_apply_record_limit(self):
        session=Mock();session.get.return_value=Mock(ok=True,json=lambda:{'molecules':[MOLECULE],'page_meta':{'next':None}})
        with patch('crawlers.chembl_client.requests.Session',return_value=session):
            self.assertEqual(list(ChEMBLClient().substructure.filter(smiles='F/C=C/F').take(1)),[MOLECULE])
        self.assertIn('substructure/F%2FC%3DC%2FF.json',session.get.call_args.args[0])
        self.assertEqual(session.get.call_args.kwargs['params']['limit'],1)
        session.close.assert_called_once()

    def test_activity_cutoff_requires_units_and_unrestricted_search_retains_missing_values(self):
        from crawlers.bioactivities import Bioactivity
        with tempfile.TemporaryDirectory() as folder:
            crawler=Bioactivity(path=folder)
            query=Mock();query.filter.return_value.only.return_value.take.return_value=[{
                'canonical_smiles':'CCO','molecule_chembl_id':'CHEMBL1','standard_value':None}]
            crawler._Bioactivity__bioactivity=query
            with self.assertRaisesRegex(ValueError,'unidade'):
                crawler._Bioactivity__search_bioactivity('CHEMBL220',{'max_value_ref':5000})
            query.filter.assert_not_called()
            crawler._Bioactivity__search_bioactivity('CHEMBL220',{})
            result=pd.read_csv(Path(folder)/'CHEMBL220.csv')
            self.assertEqual(len(result),1)
            self.assertTrue(result['standard_value'].isna().all())

    def test_zinc_table_parser_accepts_whitespace_and_rejects_bad_columns(self):
        from crawlers.molecules import ZincMols
        crawler=ZincMols();session=Mock()
        response=Mock(status_code=200,text='smiles zinc_id\nCCO    ZINC1\nCCC\tZINC2\n')
        session.get.return_value=response
        with patch('crawlers.molecules.requests.Session',return_value=session):
            frame=crawler._ZincMols__search_in_zinc(1,'https://example.org/data',False)
            self.assertEqual(frame['zinc_id'].tolist(),['ZINC1','ZINC2'])
            response.text='smiles zinc_id\nCCO ZINC1 unexpected\n'
            with self.assertRaisesRegex(ValueError,'Tabela ZINC inválida'):
                crawler._ZincMols__search_in_zinc(1,'https://example.org/data',False)

    def test_pdb_wrapper_can_search_without_collection_name_or_ec(self):
        with patch('wrappers.crawlers.PDBComplex') as crawler:
            load_pdb(base_output_path='/tmp/test',pdb_query='kinase')
            params=crawler.return_value.get_pdb_files.call_args.kwargs['filters']
            self.assertEqual(params['pdb_query'],'kinase')
            self.assertNotIn('ec_target',params)


if __name__=='__main__':unittest.main()
