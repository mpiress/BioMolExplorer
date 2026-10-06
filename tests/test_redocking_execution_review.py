"""Real project/pipeline handoffs with external engines replaced at the job boundary."""
import copy
import unittest
from pathlib import Path
from unittest.mock import patch
import test_pipeline_execution as execution
from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService

PDB=('ATOM      1  CA  ALA A   2       0.000   0.000   0.000  1.00  0.00           C\n'
     'HETATM    2  C   LIG A   1       1.000   0.000   0.000  1.00  0.00           C\n')


class RedockingExecutionReviewTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    upload=execution.PipelineExecutionTests.upload
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def configuration(self,selected=True):
        asset=self.upload('1ABC.pdb',PDB,'structures')
        source=new_stage('import_results')
        source['parameters'].update(kind='structures',asset_ids=[asset])
        stage=new_stage('redocking');stage['input_processing']='individual'
        stage['bindings']['base_input_path']={'stage':source['id'],'selector':'auto'}
        if selected:stage['parameters']['pdb_codes']=[['1ABC','LIG',1,'A']]
        return source,stage

    def test_preconfigured_pairs_execute_without_another_selection_popup(self):
        source,stage=self.configuration()
        self.save([source,stage]);manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager), patch('biomolexplorer.redocking_config.validate_redocking_tools'):
            service=PipelineService(self.store)
            try:run=self.wait(service.submit(self.token,self.project_id))
            finally:service.close()
        self.assertEqual(run['status'],'succeeded',run['error'])
        self.assertEqual([op for op,_ in calls],['redocking'])
        self.assertEqual(calls[0][1]['pdb_codes'],[['1ABC','LIG',1,'A']])
        self.assertFalse(run['stages'][1]['requires_curation'])
        self.assertEqual(len(run['stages'][1]['batches']),1)

    def test_unconfigured_pairs_still_require_explicit_selection(self):
        source,stage=self.configuration(False);self.save([source,stage])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store)
            try:
                run=self.wait(service.submit(self.token,self.project_id))
                self.assertEqual(run['status'],'awaiting_input')
                with self.assertRaisesRegex(ValueError,'Selecione pelo menos um par'):
                    service.resume(self.token,run['id'],copy.deepcopy(run['stages'][1]['configuration']))
                self.assertEqual(calls,[])
            finally:service.close()

    def test_invalid_residue_fails_before_submitting_a_scientific_job(self):
        source,stage=self.configuration();stage['parameters']['pdb_codes'][0][2]=99
        self.save([source,stage]);manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store)
            try:run=self.wait(service.submit(self.token,self.project_id))
            finally:service.close()
        self.assertEqual(run['status'],'failed')
        self.assertIn('Ligante ou resíduo',run['error'])
        self.assertEqual(calls,[])

    def test_missing_chimera_fails_before_submitting_a_scientific_job(self):
        source,stage=self.configuration();self.save([source,stage]);manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager),patch('shutil.which',return_value=None):
            service=PipelineService(self.store)
            try:run=self.wait(service.submit(self.token,self.project_id))
            finally:service.close()
        self.assertEqual(run['status'],'failed')
        self.assertIn('chimera',run['error'])
        self.assertEqual(calls,[])

