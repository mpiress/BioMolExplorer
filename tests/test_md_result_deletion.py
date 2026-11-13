"""AMBER result deletion uses one write transaction and permits recomputation."""
import json
from pathlib import Path
import sqlite3
import tempfile
import time
import unittest
from unittest.mock import patch
from uuid import uuid4

from biomolexplorer.catalog import new_stage,operation_fields
from biomolexplorer.md_results import MDResults
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.result_files import ResultFiles
from biomolexplorer.stage_cache import artifact_manifest
from biomolexplorer.workspace import AccessDenied,WorkspaceStore


class MDResultDeletionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.store=WorkspaceStore(self.temp.name)
        self.token=self.store.register('Owner','owner@example.org','owner-password')
        self.pid=self.store.create_project(self.token,'AMBER cleanup')['id'];self.rid=uuid4().hex
        self.stage=new_stage('dynamics_dock6');self.key='a'*64
        self.store.save_pipeline(self.token,self.pid,[self.stage],0)
        root=self.store.project_dir(self.pid)/'runs'/self.rid/self.stage['id']/'artifacts';root.mkdir(parents=True)
        self.table=root/'dynamics_results.json';self.pose=root/'pairs'/self.key/'final.pdb';self.pose.parent.mkdir(parents=True)
        self.pose.write_text('scientific result')
        self.table.write_text(json.dumps({'pairs':[dict(id=self.key,receptor='2WSL',ligand='CHEMBL92663',status='succeeded',amber_score=2)]}))
        self.files=[str(self.table),str(self.pose)]
        item=dict(self.stage,configuration=self.stage,status='succeeded',artifacts=self.files,artifact_manifest=artifact_manifest(self.files))
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',(self.rid,self.pid,self.store.user(self.token)['id'],'succeeded',json.dumps([item]),time.time(),time.time(),None))
        self.service=ResultFiles(self.store)

    def args(self):return self.token,self.pid,self.rid,self.stage['id']

    def test_bulk_delete_avoids_nested_projection_writer_and_revisions_cache(self):
        with patch('biomolexplorer.result_catalog.read_run',side_effect=sqlite3.OperationalError('database is locked')):
            self.assertEqual(self.service.remove_all(*self.args()),2)
        self.assertFalse(self.pose.exists());self.assertFalse(self.table.exists())
        stage=self.store.project(self.token,self.pid)['pipeline'][0]
        self.assertTrue(stage['result_cache_epoch'])
        runtime=PipelineService._runtime_parameters(type('Service',(),{'dock6_path':Path('/dock6')})(),stage)
        self.assertEqual(runtime['md_cache_epoch'],stage['result_cache_epoch'])
        self.assertNotIn('md_cache_epoch',{field['name'] for field in operation_fields('dynamics_dock6')})
        result=self.store.get_result_run(self.token,self.rid)['stages'][0]
        self.assertTrue(result['outputs_removed']);self.assertEqual(result['artifacts'],[])

    def test_single_file_delete_and_rollback_restore_metadata_and_bytes(self):
        with patch('biomolexplorer.result_catalog.read_run',side_effect=AssertionError('Nested result reader')):
            self.service.remove(*self.args(),str(self.pose))
        epoch=self.store.project(self.token,self.pid)['pipeline'][0]['result_cache_epoch']
        self.assertTrue(epoch);self.assertFalse(self.pose.exists());self.assertTrue(self.table.exists())
        history=self.store.history(self.token,self.pid)[0]
        self.store.rollback(self.token,self.pid,history['id'])
        self.assertEqual(self.pose.read_text(),'scientific result')
        self.assertNotIn('result_cache_epoch',self.store.project(self.token,self.pid)['pipeline'][0])

    def test_failed_transaction_restores_files_and_cache_generation(self):
        with patch('biomolexplorer.project_state.record',side_effect=OSError('snapshot failure')):
            with self.assertRaises(OSError):self.service.remove_all(*self.args())
        self.assertTrue(self.pose.exists());self.assertTrue(self.table.exists())
        self.assertNotIn('result_cache_epoch',self.store.project(self.token,self.pid)['pipeline'][0])

    def test_pair_exclusion_does_not_enter_result_projection_inside_transaction(self):
        with patch('biomolexplorer.result_catalog.read_run',side_effect=sqlite3.OperationalError('database is locked')):
            project=MDResults(self.store).remove(*self.args(),self.key)
        self.assertIn(self.key,project['pipeline'][0]['parameters']['excluded_pairs'])
        self.assertTrue(self.pose.exists())

    def test_cross_project_result_cannot_be_deleted(self):
        other=self.store.create_project(self.token,'Other')['id']
        with self.assertRaises(AccessDenied):self.service.remove_all(self.token,other,self.rid,self.stage['id'])
        self.assertTrue(self.pose.exists())


if __name__=='__main__':unittest.main()
