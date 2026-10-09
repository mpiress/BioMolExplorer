"""ZINC URI/script lists, compressed structures and molecular handoffs."""
import gzip
import bz2
import shutil
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import Mock,patch

import test_pipeline_execution as execution
from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.zinc_retrieval import read_download_list,retrieve_tranches
from biomolexplorer.docking_data import read_compounds,input_records,prepare_ligands,write_csv
from biomolexplorer.ui.file_selection import compatible_file

MOL2='''@<TRIPOS>MOLECULE
ZINC2
3 2 1 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 C.3 1 LIG -0.1
2 C2 1.5 0.0 0.0 C.3 1 LIG 0.1
3 O1 2.8 0.0 0.0 O.3 1 LIG -0.3
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
@<TRIPOS>SUBSTRUCTURE
1 LIG 1
'''


def response(data=b'',status=200,location=None):
    result=Mock(status_code=status,headers={'Location':location} if location else {})
    result.__enter__=Mock(return_value=result);result.__exit__=Mock(return_value=False)
    result.iter_content.return_value=[data]
    return result


class ZincTrancheTests(unittest.TestCase):
    def test_parallel_downloads_preserve_list_order_and_bound_prefetch(self):
        from biomolexplorer.zinc_retrieval import _downloads
        with tempfile.TemporaryDirectory() as folder:
            urls=['https://files.docking.org/'+str(i)+'.smi' for i in range(6)]
            second_finished=threading.Event()
            lock=threading.Lock();started=[];finished=[]
            def download(index,url,temporary):
                with lock:started.append(index)
                if index==0:self.assertTrue(second_finished.wait(5))
                path=Path(temporary)/(str(index)+'.smi');path.write_text('CCO ZINC1\n')
                with lock:finished.append(index)
                if index==1:second_finished.set()
                return path,url
            with patch('biomolexplorer.zinc_retrieval._download_job',side_effect=download):
                with _downloads(urls,folder,2) as downloads:
                    first=next(downloads)
                    self.assertEqual(first[0],urls[0])
                    self.assertEqual(sorted(started),[0,1])
                    self.assertEqual(finished,[1,0])
                    ordered=[first[0]]+[item[0] for item in downloads]
                    self.assertEqual(ordered,urls)
            self.assertEqual(list(Path(folder).iterdir()),[])

    def test_parallel_and_sequential_outputs_are_identical(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=root/'list.uri'
            path.write_text('https://files.docking.org/a.smi\nhttps://files.docking.org/b.smi\n')
            def download(index,url,temporary):
                target=Path(temporary)/(str(index)+'.smi')
                target.write_text('CCO ZINC1\n' if index==0 else 'CCO ZINC1\nCCN ZINC2\n')
                return target,url
            with patch('biomolexplorer.zinc_retrieval._download_job',side_effect=download):
                sequential=retrieve_tranches(path,root/'one',download_workers=1)
                parallel=retrieve_tranches(path,root/'two',download_workers=2)
            self.assertEqual((root/'one/compounds.csv').read_bytes(),(root/'two/compounds.csv').read_bytes())
            self.assertEqual(sequential['downloads'],parallel['downloads'])
            self.assertEqual(parallel['download_workers'],2)
            self.assertEqual(parallel['duplicates'],1)

    def test_failed_parallel_download_does_not_publish_table(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=root/'list.uri'
            path.write_text('https://files.docking.org/a.smi\nhttps://files.docking.org/b.smi\n')
            def download(index,url,temporary):
                if index==1:raise ConnectionError('download failed')
                target=Path(temporary)/'a.smi';target.write_text('CCO ZINC1\n')
                return target,url
            with patch('biomolexplorer.zinc_retrieval._download_job',side_effect=download):
                with self.assertRaises(ConnectionError):retrieve_tranches(path,root/'out')
            self.assertFalse((root/'out/compounds.csv').exists())
            self.assertFalse((root/'out/retrieval_report.json').exists())

    def test_worker_limits_are_validated_before_download(self):
        from biomolexplorer.operations import validate_operation
        for value in (0,17,True,2.5,'4',None):
            with self.subTest(value=value),patch('requests.Session') as session:
                with self.assertRaisesRegex(ValueError,'1 e 16'):
                    retrieve_tranches('unused.uri','unused',download_workers=value)
                with self.assertRaisesRegex(ValueError,'1 e 16'):
                    validate_operation('retrieve_zinc',{'base_input_path':'input','download_workers':value})
                session.assert_not_called()

    def test_uri_and_script_formats_deduplicate_and_upgrade_http(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'zinc-download.sh'
            path.write_text('#!/bin/bash\nmkdir -p AA\ncurl -o AA/a.smi http://files.docking.org/2D/AA/a.smi\n'
                'wget "https://files2.docking.org/3D/AA/a.mol2.gz"\nhttps://files.docking.org/2D/AA/a.smi\n')
            self.assertEqual(read_download_list(path),['https://files.docking.org/2D/AA/a.smi',
                'https://files2.docking.org/3D/AA/a.mol2.gz'])
            self.assertTrue(compatible_file(path,{'zinc_urls'}))

    def test_wrong_hosts_and_unsupported_formats_fail_before_requests(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'list.uri'
            for value in ('https://localhost/a.smi','https://files.docking.org.evil.org/a.smi',
                          'https://files.docking.org/a.db2.gz','https://user@files.docking.org/a.smi','empty'):
                path.write_text(value)
                with self.subTest(value=value),patch('requests.Session') as session:
                    with self.assertRaises(ValueError):retrieve_tranches(path,Path(folder)/'out')
                    session.assert_not_called()

    def test_smi_without_header_and_mol2_gzip_produce_standard_bundle(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=root/'list.uri'
            path.write_text('https://files.docking.org/a.smi\nhttps://files.docking.org/b.mol2.gz\n')
            with patch('requests.Session') as factory:
                session=factory.return_value.__enter__.return_value
                session.get.side_effect=lambda url,**kwargs:response(b'CCN ZINC1\n' if url.endswith('.smi') else gzip.compress(MOL2.encode()))
                report=retrieve_tranches(path,root/'out')
            rows=read_compounds(root/'out/compounds.csv')
            self.assertEqual({r['molecule_chembl_id'] for r in rows},{'ZINC1','ZINC2'})
            molecule=next(r for r in rows if r['molecule_chembl_id']=='ZINC2')
            self.assertEqual(Path(molecule['conformer_file']).read_text(),MOL2)
            self.assertEqual(report['conformations'],1)
            self.assertTrue((root/'out/retrieval_report.json').exists())
            files=[p for p in (root/'out').rglob('*') if p.is_file()]
            self.assertEqual(len(input_records([root/'out/compounds.csv'],files)),2)

    def test_compressed_smi_header_order_and_duplicates(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=root/'list.uri';path.write_text('https://files.docking.org/a.smi.bz2\n')
            data=bz2.compress(b'zinc_id smiles\n123 CCO\n123 CCO\n')
            with patch('requests.Session') as factory:
                factory.return_value.__enter__.return_value.get.return_value=response(data)
                report=retrieve_tranches(path,root/'out')
            self.assertEqual(report['compounds'],1);self.assertEqual(report['duplicates'],1)
            self.assertEqual(read_compounds(root/'out/compounds.csv')[0]['molecule_chembl_id'],'ZINC000000000123')

    def test_redirects_are_validated_before_following(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=root/'list.uri';path.write_text('https://files.docking.org/a.smi\n')
            with patch('requests.Session') as factory:
                session=factory.return_value.__enter__.return_value
                session.get.return_value=response(status=302,location='http://localhost/secret.smi')
                with self.assertRaisesRegex(ValueError,'não autorizado'):retrieve_tranches(path,root/'out')
                session.get.assert_called_once()

    def test_invalid_smiles_and_identifier_conflicts_do_not_publish_table(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=root/'list.uri';path.write_text('https://files.docking.org/a.smi\n')
            for index,data in enumerate((b'invalid ZINC1\n',b'CCO ZINC1\nCCC ZINC1\n')):
                output=root/str(index)
                with patch('requests.Session') as factory:
                    factory.return_value.__enter__.return_value.get.return_value=response(data)
                    with self.assertRaises(ValueError):retrieve_tranches(path,output)
                self.assertFalse((output/'compounds.csv').exists())

    def test_library_conformations_are_centered_for_dock6_without_changing_source(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);source=root/'zinc.mol2';source.write_text(MOL2)
            dataset=root/'compounds.csv'
            write_csv(dataset,[dict(molecule_chembl_id='ZINC2',canonical_smiles='CCO',
                conformer_file=str(source),conformer_origin='library')])
            def convert(source,destination,**kwargs):shutil.copy2(source,destination)
            with patch('biomolexplorer.docking_data.convert_structure',side_effect=convert):
                outputs=prepare_ligands(dataset,root/'direct','mol2',center=[10.,20.,30.])
            def center(path):
                atoms=Path(path).read_text().split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0].strip().splitlines()
                return [sum(float(line.split()[k]) for line in atoms)/len(atoms) for k in (2,3,4)]
            for actual,expected in zip(center(outputs[0]),[10.,20.,30.]):self.assertAlmostEqual(actual,expected,places=3)
            self.assertEqual(source.read_text(),MOL2)
            write_csv(dataset,[dict(molecule_chembl_id='ZINC2',canonical_smiles='CCO',
                prepared_mol2=str(source),prepared_origin='library',docking_engines='dock6')])
            with patch('biomolexplorer.docking_data.convert_structure') as convert:
                outputs=prepare_ligands(dataset,root/'reused','mol2',center=[5.,6.,7.])
                convert.assert_not_called()
            for actual,expected in zip(center(outputs[0]),[5.,6.,7.]):self.assertAlmostEqual(actual,expected,places=3)


class ZincPipelineTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    upload=execution.PipelineExecutionTests.upload

    def test_uploaded_uri_file_resolves_directly_and_multiple_lists_merge(self):
        first=self.upload('download.uri','http://files.docking.org/2D/a.smi\n','zinc_urls')
        second=self.upload('download.sh','wget https://files.docking.org/3D/b.mol2.gz\n','zinc_urls')
        stage=new_stage('retrieve_zinc');stage['bindings']={'base_input_path':{'sources':[{'asset':first},{'asset':second}]}}
        service=PipelineService(self.store)
        try:
            params=service._resolve(self.project_id,self.store.user(self.token)['id'],stage,{})
            urls=read_download_list(Path(params['base_input_path'])/params['filename'])
            self.assertEqual(len(urls),2);self.assertTrue(all(u.startswith('https://') for u in urls))
        finally:service.close()

    def test_invalid_list_upload_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'Padrão esperado'):
            self.upload('download.uri','https://example.org/a.smi\n','zinc_urls')

    def test_mixed_import_routes_only_download_lists_to_zinc(self):
        csv=self.upload()
        uri=self.upload('download.uri','https://files.docking.org/a.smi\n','zinc_urls')
        origin=new_stage('import_results');origin['parameters'].update(kind='other',asset_ids=[csv,uri],
            asset_types={csv:'compounds',uri:'zinc_urls'})
        stage=new_stage('retrieve_zinc');stage['bindings']={'base_input_path':{'stage':origin['id']}}
        service=PipelineService(self.store)
        try:
            files=service._import(self.project_id,origin['parameters'],self.store.project_dir(self.project_id)/'imported')
            params=service._resolve(self.project_id,self.store.user(self.token)['id'],stage,{origin['id']:files})
            self.assertEqual(read_download_list(Path(params['base_input_path'])/params['filename']),['https://files.docking.org/a.smi'])
        finally:service.close()
