"""Persisted, confirmed file choices can complete a restarted pipeline."""
import copy
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import test_pipeline_execution as execution
from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.workspace import WorkspaceStore


class RestartReuseTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    finish_selected=execution.PipelineExecutionTests.finish_selected
    upload=execution.PipelineExecutionTests.upload
    imported=execution.PipelineExecutionTests.imported
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def connected(self):
        source,consumer=new_stage('retrieve_compounds'),new_stage('admet')
        consumer['bindings']['base_input_path']={'stage':source['id']}
        self.save([source,consumer])
        return source,consumer

    def first_run(self):
        service=PipelineService(self.store,worker_python=sys.executable)
        try:
            run=self.finish_selected(service,service.submit(self.token,self.project_id))
            self.assertEqual(run['status'],'succeeded',run['error'])
            return run
        finally:service.close()

    def reopen(self):
        self.store=WorkspaceStore(self.store.root)
        return PipelineService(self.store,worker_python=sys.executable)

    def test_restart_reuses_completed_chain_without_repeating_file_selection(self):
        self.connected()
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager),patch('biomolexplorer.pipeline.implementation_digest',return_value='release-a'):
            first=self.first_run()
            service=self.reopen()
            try:second=self.wait(service.submit(self.token,self.project_id,reuse_results=True))
            finally:service.close()
        self.assertEqual(second['status'],'succeeded',second['error'])
        self.assertEqual(len(calls),2)
        for original,reused in zip(first['stages'],second['stages']):
            self.assertTrue(reused['reused'])
            self.assertEqual(reused['artifacts'],original['artifacts'])
        self.assertEqual(second['stages'][1]['configuration']['bindings']['base_input_path']['selector'],'compounds.csv')

    def test_restart_after_app_upgrade_verifies_and_reuses_saved_provenance(self):
        self.connected()
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            with patch('biomolexplorer.pipeline.implementation_digest',return_value='release-a'):
                first=self.first_run()
            with patch('biomolexplorer.pipeline.implementation_digest',return_value='release-b'):
                service=self.reopen()
                try:second=self.wait(service.submit(self.token,self.project_id,reuse_results=True))
                finally:service.close()
        self.assertEqual(second['status'],'succeeded',second['error'])
        self.assertEqual(len(calls),2)
        self.assertTrue(all(stage.get('reused') for stage in second['stages']))
        self.assertTrue(all(stage['cache_key']!=old['cache_key'] for stage,old in zip(second['stages'],first['stages'])))
        self.assertTrue(all(stage['cache_software'].startswith('release-b:') for stage in second['stages']))

    def test_declining_reuse_runs_collection_and_waits_for_new_file_selection(self):
        self.connected()
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager),patch('biomolexplorer.pipeline.implementation_digest',return_value='release-a'):
            self.first_run()
            service=self.reopen()
            try:
                second=self.wait(service.submit(self.token,self.project_id,reuse_results=False))
                self.assertEqual(second['status'],'awaiting_input')
                self.assertEqual(len(calls),3)
                self.assertFalse(any(stage.get('reused') for stage in second['stages']))
                finished=self.finish_selected(service,second)
            finally:service.close()
        self.assertEqual(finished['status'],'succeeded',finished['error'])
        self.assertEqual(len(calls),4)

    def test_changed_processing_mode_after_restart_requires_new_selection(self):
        source,consumer=self.connected()
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager),patch('biomolexplorer.pipeline.implementation_digest',return_value='release-a'):
            self.first_run()
            consumer['input_processing']='individual'
            self.save([source,consumer])
            service=self.reopen()
            try:second=self.wait(service.submit(self.token,self.project_id,reuse_results=True))
            finally:service.close()
        self.assertEqual(second['status'],'awaiting_input')
        self.assertTrue(second['stages'][0]['reused'])
        self.assertEqual(second['stages'][1]['status'],'awaiting_input')
        self.assertEqual(len(calls),2)

    def test_modified_consumer_output_cannot_be_reused_after_upgrade(self):
        self.connected()
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            with patch('biomolexplorer.pipeline.implementation_digest',return_value='release-a'):
                first=self.first_run()
            Path(first['stages'][1]['artifacts'][0]).write_text('molecule_chembl_id,canonical_smiles\nCHANGED,O\n')
            with patch('biomolexplorer.pipeline.implementation_digest',return_value='release-b'):
                service=self.reopen()
                try:second=self.wait(service.submit(self.token,self.project_id,reuse_results=True))
                finally:service.close()
        self.assertEqual(second['status'],'awaiting_input')
        self.assertTrue(second['stages'][0]['reused'])
        self.assertFalse(second['stages'][1].get('reused',False))
        self.assertEqual(len(calls),2)

    def test_changed_uploaded_input_invalidates_consumer_after_upgrade(self):
        asset=self.upload()
        source=self.imported(asset)
        consumer=new_stage('admet')
        consumer['bindings']['base_input_path']={'stage':source['id']}
        self.save([source,consumer])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            with patch('biomolexplorer.pipeline.implementation_digest',return_value='release-a'):
                self.first_run()
            self.store.asset_path(self.token,self.project_id,asset).write_text('name,smiles\nWATER,O\n')
            with patch('biomolexplorer.pipeline.implementation_digest',return_value='release-b'):
                service=self.reopen()
                try:second=self.wait(service.submit(self.token,self.project_id,reuse_results=True))
                finally:service.close()
        self.assertEqual(second['status'],'awaiting_input')
        self.assertFalse(any(stage.get('reused') for stage in second['stages']))
        self.assertEqual(len(calls),1)

    def test_tampered_original_cache_key_is_not_adopted_after_upgrade(self):
        stage=new_stage('retrieve_compounds')
        self.save([stage])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            with patch('biomolexplorer.pipeline.implementation_digest',return_value='release-a'):
                first=self.first_run()
            first['stages'][0]['cache_key']='0'*64
            with self.store.connect() as db:
                db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(first['stages']),first['id']))
            with patch('biomolexplorer.pipeline.implementation_digest',return_value='release-b'):
                service=self.reopen()
                try:second=self.wait(service.submit(self.token,self.project_id,reuse_results=True))
                finally:service.close()
        self.assertEqual(second['status'],'succeeded',second['error'])
        self.assertFalse(second['stages'][0].get('reused',False))
        self.assertEqual(len(calls),2)

    def test_restart_restores_only_the_previously_selected_file(self):
        source,consumer=self.connected()
        base,calls=self.immediate_manager()
        class MultipleFiles(base):
            def submit(manager,operation,parameters):
                job=super().submit(operation,parameters)
                if operation=='retrieve_compounds':
                    extra=manager.config.workspace/'extra.csv'
                    extra.write_text('molecule_chembl_id,canonical_smiles\nEXTRA,CCN\n')
                    job['result']['artifacts'].append(str(extra))
                return job
        with patch('biomolexplorer.pipeline.JobManager',MultipleFiles),patch('biomolexplorer.pipeline.implementation_digest',return_value='release-a'):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:
                pending=self.wait(service.submit(self.token,self.project_id))
                self.assertEqual(pending['status'],'awaiting_input')
                selected=copy.deepcopy(consumer)
                selected['bindings']['base_input_path']['selector']='extra.csv'
                first=self.wait(service.resume(self.token,pending['id'],selected))
                self.assertEqual(first['status'],'succeeded',first['error'])
            finally:service.close()
            service=self.reopen()
            try:second=self.wait(service.submit(self.token,self.project_id,reuse_results=True))
            finally:service.close()
        self.assertEqual(second['status'],'succeeded',second['error'])
        self.assertEqual(len(calls),2)
        self.assertEqual(second['stages'][1]['configuration']['bindings']['base_input_path']['selector'],'extra.csv')
        self.assertEqual(second['stages'][1]['input_files']['base_input_path'],first['stages'][1]['input_files']['base_input_path'])


if __name__=='__main__':unittest.main()
