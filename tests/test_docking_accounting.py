import csv
import json
from pathlib import Path
import tempfile
import unittest

from biomolexplorer.docking_accounting import reconcile
from biomolexplorer.docking_data import write_csv
from biomolexplorer.ligand_cache import validate_prepared_ligand
from biomolexplorer.result_files import ResultFiles
from biomolexplorer.artifacts import index_run
from biomolexplorer.artifact_registry import register_run
import test_clear_all_results as cleanup
from test_preparation_failures import PDBQT


class AccountingTests(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root=Path(temporary.name)
        self.dataset=self.root/'compounds.csv'
        write_csv(self.dataset,[dict(molecule_chembl_id=code,canonical_smiles='CCO') for code in ['A','A','B','C']])

    def test_union_counts_and_failure_reasons_reconcile_every_pair(self):
        rows=[dict(molecule_chembl_id='A',receptor_id='4M0E_A')]
        failures=[dict(compound_identifier='B',reason='Invalid charge')]
        execution=[dict(compound_identifier='C',configuration='4M0E_LIG_1A_C.vina',error='Engine rejected input')]
        report=reconcile('vina',self.dataset,[['4M0E','LIG',1,'A']],rows,self.root,failures,execution)
        self.assertEqual(report['input_rows'],4)
        self.assertEqual(report['unique_compounds'],3)
        self.assertEqual(report['expected_pairs'],3)
        self.assertEqual(report['successful_pairs']+report['unsuccessful_pairs'],3)
        self.assertEqual(report['unaccounted_pairs'],0)
        with (self.root/'docking_unprocessed.csv').open() as stream:
            missing=list(csv.DictReader(stream))
        self.assertEqual([row['status'] for row in missing],['preparation_failed','docking_failed'])

    def test_duplicates_and_foreign_results_are_rejected(self):
        row=dict(molecule_chembl_id='A',receptor_id='4M0E_A')
        for rows in ([row,row],[dict(row,receptor_id='OTHER_A')]):
            with self.assertRaises(ValueError):
                reconcile('dock6',self.dataset,[['4M0E','LIG',1,'A']],rows,self.root)

    def test_distinct_receptor_poses_are_not_reported_as_duplicate_inputs(self):
        rows=[dict(molecule_chembl_id='A',canonical_smiles='CCO',receptor_id=receptor,conformer_file='pose.pdbqt')
              for receptor in ('4M0E_A','2WSL_A')]
        write_csv(self.dataset,rows)
        report=reconcile('dock6',self.dataset,[['4M0E','LIG',1,'A'],['2WSL','LIG',1,'A']],rows,self.root)
        self.assertEqual(report['duplicate_input_rows'],0)
        self.assertEqual(report['successful_pairs'],2)

    def test_multifragment_pdbqt_is_rejected_before_native_engine(self):
        path=self.root/'salt.pdbqt'
        path.write_text(PDBQT+PDBQT)
        with self.assertRaisesRegex(ValueError,'fragmentos ROOT'):
            validate_prepared_ligand(path,'pdbqt')
        path.write_text('MODEL 1\n'+PDBQT+'ENDMDL\nMODEL 2\n'+PDBQT+'ENDMDL\n')
        validate_prepared_ligand(path,'pdbqt')

    def test_optional_footprint_failure_preserves_completed_docking_and_minimization(self):
        from caad.docking import Dock6
        from biomolexplorer.preparation_failures import PreparationFailures
        engine=object.__new__(Dock6)
        engine.failure_report=PreparationFailures()
        engine.compound_dataset=self.dataset
        engine.compound_rows={'A':(2,dict(molecule_chembl_id='A',canonical_smiles='CCO'))}
        engine.export_path=self.root
        engine._Dock6__base_output_path=str(self.root)
        files=['A.lig.min_scored.mol2','A_scored.mol2','A_footprint_scored.txt','A.pdf']
        for name in files:(self.root/name).write_text('scientific output')
        with engine._compound_guard('A','Footprint DOCK6'):
            raise ValueError('Optional footprint failed')
        self.assertTrue((self.root/files[0]).is_file())
        self.assertTrue((self.root/files[1]).is_file())
        self.assertFalse((self.root/files[2]).exists())
        self.assertFalse((self.root/files[3]).exists())
        (self.root/files[2]).write_text('completed footprint')
        with engine._compound_guard('A','Gráfico do footprint DOCK6'):
            raise ValueError('Optional plot failed')
        self.assertTrue((self.root/files[2]).is_file())

    def test_recorded_files_metric_excludes_missing_files_and_duplicate_references(self):
        from biomolexplorer.ui.run_presentation import outcome_overview
        path=self.root/'pose.pdbqt'
        path.write_text(PDBQT)
        stage=dict(status='succeeded',artifacts=[str(path),str(path),str(self.root/'deleted.pdbqt')])
        overview=outcome_overview(dict(stages=[stage],status='succeeded',created=0,updated=1),[])
        metric=overview.content.controls[1].controls[1]
        self.assertEqual(metric.content.controls[0].value,'1')


class DeletionIndexTests(unittest.TestCase):
    setUp=cleanup.ClearAllResultsTests.setUp

    def test_incomplete_docking_is_not_reused_as_a_complete_stage_or_batch(self):
        from biomolexplorer.pipeline import PipelineService
        stage=self.stages[2]
        stage.update(operation='docking_vina',docking_failures=1)
        with self.store.connect() as db:
            db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(self.stages),self.rid))
        service=object.__new__(PipelineService)
        service.store=self.store
        self.assertIsNone(service._cached_stage(self.pid,self.store.user(self.token)['id'],
            dict(stage['configuration'],operation='docking_vina'),'key','software'))
        self.assertIsNone(service._cached_batch(self.pid,stage['id'],'key'))

    def test_stage_deletion_updates_index_and_manifest_atomically(self):
        index_run(self.stages,self.rid)
        with self.store.connect() as db:
            register_run(db,self.pid,self.rid,self.stages)
            db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(self.stages),self.rid))
        stage=self.stages[2]
        ResultFiles(self.store).remove_all(self.token,self.pid,self.rid,stage['id'])
        with self.store.connect() as db:
            self.assertEqual(db.execute('SELECT count(*) FROM artifact_instances WHERE project_id=? AND stage_id=?',
                (self.pid,stage['id'])).fetchone()[0],0)
            manifest=json.loads(db.execute('SELECT manifest FROM execution_manifests WHERE run_id=? AND stage_id=?',
                (self.rid,stage['id'])).fetchone()[0])
        self.assertEqual(manifest['outputs'],{})
        self.assertTrue(self.original.is_file())


if __name__=='__main__':
    unittest.main()
