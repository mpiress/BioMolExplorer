"""Restarted pending runs must cancel atomically even under SQLite spill locks."""
import json
import sqlite3
import unittest
from contextlib import contextmanager
from unittest.mock import patch
from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.result_files import ResultFiles
import test_pipeline_execution as execution

class ExclusiveWriteConnection(sqlite3.Connection):
    def execute(self,sql,*args,**kwargs):
        # SQLite also takes an exclusive lock when a large stages blob spills
        # dirty pages. Reproduce that state without giant scientific fixtures.
        if sql=='BEGIN IMMEDIATE':sql='BEGIN EXCLUSIVE'
        return super().execute(sql,*args,**kwargs)

class RestartedCancellationTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    @contextmanager
    def exclusive_writes(self):
        db=sqlite3.connect(self.store.database,timeout=.1,factory=ExclusiveWriteConnection)
        db.row_factory=sqlite3.Row
        try:
            with db:yield db
        finally:db.close()

    def paused_pipeline(self):
        source=new_stage('retrieve_compounds');pending=new_stage('admet');later=new_stage('fingerprints')
        pending['bindings']={'base_input_path':{'stage':source['id']}}
        later['bindings']={'base_input_path':{'stage':pending['id']}}
        self.save([source,pending,later])
        service=PipelineService(self.store)
        try:
            run=self.wait(service.submit(self.token,self.project_id))
            self.assertEqual(run['status'],'awaiting_input',run.get('error'))
        finally:service.close()
        return run,source

    def test_restarted_run_cancels_preserves_completed_work_and_unlocks_delete_and_rerun(self):
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            run,source=self.paused_pipeline()
            service=PipelineService(self.store);self.addCleanup(service.close)
            self.assertFalse(service._active)
            original=self.store.get_run(self.token,run['id'])['stages'][0]['artifacts']
            with patch.object(self.store,'connect',self.exclusive_writes):
                service.cancel(self.token,run['id'])
            cancelled=self.store.get_run(self.token,run['id'])
            self.assertEqual(cancelled['status'],'cancelled')
            self.assertEqual([s['status'] for s in cancelled['stages']],['succeeded','cancelled','skipped'])
            self.assertEqual(cancelled['stages'][0]['artifacts'],original)
            events=self.store.history(self.token,self.project_id)
            self.assertEqual(events[0]['action'],'execution')
            self.assertIn('cancelled',events[0]['summary'])
            # Repeated clicks must not create a second terminal history record.
            service.cancel(self.token,run['id'])
            self.assertEqual(len(self.store.history(self.token,self.project_id)),len(events))
            removed=ResultFiles(self.store).remove_all(self.token,self.project_id,run['id'],source['id'])
            self.assertGreater(removed,0)
            again=self.wait(service.submit(self.token,self.project_id))
            self.assertNotEqual(again['id'],run['id'])
            self.assertEqual(again['status'],'awaiting_input',again.get('error'))
            self.assertEqual(len(calls),2)
