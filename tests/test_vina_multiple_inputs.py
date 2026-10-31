"""Multi-file Vina preparation remains bounded, incremental and cancellable."""
import copy
import csv
import tempfile
import threading
import time
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock,patch

import flet as ft
from biomolexplorer.catalog import new_stage
from biomolexplorer.docking_data import input_records,read_compounds
from biomolexplorer.docking_inputs import compound_identifiers
from biomolexplorer.ui.compound_choice import CompoundChoice
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.workspace import WorkspaceStore
import test_pipeline_execution as execution


def write(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['molecule_chembl_id','canonical_smiles']);writer.writeheader();writer.writerows(rows)


class MultipleInputTests(unittest.TestCase):
    def test_csv_selection_never_reads_unrelated_pose_tables(self):
        with tempfile.TemporaryDirectory() as temp:
            selected=Path(temp)/'BBB+.csv';write(selected,[dict(molecule_chembl_id='A',canonical_smiles='CCO')])
            other=Path(temp)/'other.csv';other.write_text('molecule_chembl_id,canonical_smiles,conformer_file\nB,CCN,/unauthorized.pdbqt\n')
            with patch('biomolexplorer.docking_data.read_compounds',wraps=read_compounds) as reader:
                rows=input_records([selected],[selected,other])
            self.assertEqual([r['molecule_chembl_id'] for r in rows],['A']);reader.assert_called_once_with(selected)

    def test_selected_table_still_rejects_unauthorized_pose_aliases(self):
        with tempfile.TemporaryDirectory() as temp:
            for field in ('conformer_file','pose_file','vina_pose','dock6_pose','prepared_pdbqt','prepared_mol2'):
                path=Path(temp)/'selected.csv';path.write_text('molecule_chembl_id,canonical_smiles,'+field+'\nA,CCO,/unauthorized.pose\n')
                with self.subTest(field=field),self.assertRaisesRegex(ValueError,'autorizados'):
                    input_records([path],[path])

    def test_identifiers_are_lightweight_and_chemistry_is_still_validated_at_execution(self):
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'BBB+.csv';write(source,[dict(molecule_chembl_id='A',canonical_smiles='invalid SMILES')])
            with patch('biomolexplorer.molecule_quality._validated_row',side_effect=AssertionError('chemistry in UI')):
                self.assertEqual(compound_identifiers([source]),{'A':''})
            from biomolexplorer.molecule_quality import sanitize_csv
            with patch('rdkit.Chem.MolFromSmiles',return_value=None):sanitize_csv(source,Path(temp)/'clean.csv')
            self.assertEqual(read_compounds(Path(temp)/'clean.csv'),[])

    def test_compound_choice_keeps_all_ids_reachable_with_bounded_controls(self):
        choice=ft.Dropdown(value='M49999');picker=CompoundChoice(choice,SimpleNamespace(update=Mock()))
        picker.set_records({f'M{i:05d}':'' for i in range(50000)})
        self.assertLessEqual(len(choice.options),102);self.assertIn('M49999',[o.key for o in choice.options])
        picker.query.value='M49876';picker.search()
        self.assertIn('M49876',[o.key for o in choice.options]);self.assertEqual(choice.value,'M49999')
        picker.query.value='';picker.search();picker.move(1)
        self.assertIn('M00100',[o.key for o in choice.options]);self.assertEqual(choice.options[0].key,'')

    def test_normalized_tables_reuse_across_receptors_and_reject_changed_or_corrupt_data(self):
        with tempfile.TemporaryDirectory() as temp:
            store=WorkspaceStore(temp);token=store.register('Tester','multi@local.test','temporary-test-password');project=store.create_project(token,'Test');pid=project['id'];root=store.project_dir(pid)
            a=root/'BBB+.csv';b=root/'second'/'BBB+.csv'
            write(a,[dict(molecule_chembl_id='A',canonical_smiles='CCO')]);write(b,[dict(molecule_chembl_id='B',canonical_smiles='CCN')])
            source=new_stage('admet');stage=new_stage('docking_vina');refs=[dict(stage=source['id'],selector='BBB+.csv')]
            results={source['id']:[str(a)]};service=PipelineService(store);self.addCleanup(service.close)
            from biomolexplorer.molecule_quality import sanitize_csv
            with patch('biomolexplorer.molecule_quality.sanitize_csv',wraps=sanitize_csv) as normalizer:
                first,_,_=service._materialize_docking_compounds(pid,stage,'base_selected_mols',refs,results,None)
                stage['bindings']['base_input_path']={'stage':'another-receptor'}
                second,_,_=service._materialize_docking_compounds(pid,stage,'base_selected_mols',refs,results,root/'another-job')
                self.assertEqual(normalizer.call_count,1)
                self.assertEqual((first/'compounds.csv').read_bytes(),(second/'compounds.csv').read_bytes())
                (first/'compounds.csv').write_text('tampered')
                service._materialize_docking_compounds(pid,stage,'base_selected_mols',refs,results,None)
                self.assertEqual(normalizer.call_count,1)
                self.assertEqual(read_compounds(first/'compounds.csv')[0]['molecule_chembl_id'],'A')
                write(a,[dict(molecule_chembl_id='A',canonical_smiles='CCC')])
                changed,_,_=service._materialize_docking_compounds(pid,stage,'base_selected_mols',refs,results,None)
                self.assertNotEqual(first,changed);self.assertEqual(normalizer.call_count,2)
                from biomolexplorer.calculation_cache import CalculationCache
                cache=CalculationCache('docking-compound-tables',root/'.input-cache')
                with cache.connect() as db:sha=db.execute('select sha256 from artifacts where namespace=? order by experiment_key',('docking-compound-tables',)).fetchall()
                for row in sha:cache.object_path(row[0]).write_text('corrupt')
                service._materialize_docking_compounds(pid,stage,'base_selected_mols',refs,results,None)
                self.assertEqual(normalizer.call_count,3)
            self.assertEqual(read_compounds(a)[0]['canonical_smiles'],'CCC')

    def test_source_edit_during_normalization_cannot_publish_a_stale_checkpoint(self):
        with tempfile.TemporaryDirectory() as temp:
            store=WorkspaceStore(temp);token=store.register('Tester','change@local.test','temporary-test-password');project=store.create_project(token,'Change');pid=project['id'];root=store.project_dir(pid)
            source=root/'BBB+.csv';write(source,[dict(molecule_chembl_id='A',canonical_smiles='CCO')])
            producer=new_stage('admet');stage=new_stage('docking_vina');refs=[dict(stage=producer['id'],selector='BBB+.csv')];results={producer['id']:[str(source)]}
            service=PipelineService(store);self.addCleanup(service.close);original=service._write_docking_compounds
            def changed(selections,destination):
                path=original(selections,destination)
                write(source,[dict(molecule_chembl_id='B',canonical_smiles='CCN')]);return path
            with patch.object(service,'_write_docking_compounds',side_effect=changed):
                with self.assertRaisesRegex(ValueError,'resultados foram alterados'):
                    service._materialize_docking_compounds(pid,stage,'base_selected_mols',refs,results,None)
            folder,_,_=service._materialize_docking_compounds(pid,stage,'base_selected_mols',refs,results,None)
            self.assertEqual([r['molecule_chembl_id'] for r in read_compounds(folder/'compounds.csv')],['B'])


class ResumeLockTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def test_cancel_remains_available_during_slow_resume_validation(self):
        source=new_stage('retrieve_compounds');stage=new_stage('admet');stage['bindings']['base_input_path']={'stage':source['id']};self.save([source,stage])
        manager,_=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store);self.addCleanup(service.close)
            run=self.wait(service.submit(self.token,self.project_id));configuration=copy.deepcopy(stage)
            configuration['bindings']['base_input_path']['selector']=Path(run['stages'][0]['artifacts'][0]).name
            entered=threading.Event();release=threading.Event()
            def slow(*args,**kwargs):entered.set();self.assertTrue(release.wait(5));return {}
            with patch.object(service,'_resolve',side_effect=slow),ThreadPoolExecutor(max_workers=1) as pool:
                future=pool.submit(service.resume,self.token,run['id'],configuration)
                try:
                    self.assertTrue(entered.wait(5));started=time.monotonic();service.cancel(self.token,run['id']);cancelled=self.store.get_run(self.token,run['id'])
                    self.assertLess(time.monotonic()-started,1);self.assertEqual(cancelled['status'],'cancelled')
                finally:release.set()
                with self.assertRaises(ValueError):future.result(timeout=5)
            self.assertEqual(self.store.get_run(self.token,run['id'])['status'],'cancelled')
