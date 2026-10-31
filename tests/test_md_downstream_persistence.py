"""Downstream edits preserve upstream history; only verified work is reused."""
import copy
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.stage_cache import file_digest,operation_digest
import test_pipeline_execution as execution


class DownstreamPersistenceTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    finish_selected=execution.PipelineExecutionTests.finish_selected
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def test_restart_and_add_md_reuses_upstream_without_touching_files(self):
        stage=new_stage('retrieve_compounds');self.save([stage])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:first=self.finish_selected(service,service.submit(self.token,self.project_id))
            finally:service.close()
            before={p:(file_digest(p),Path(p).stat().st_mtime_ns) for p in first['stages'][0]['artifacts']}
            md=new_stage('dynamics_dock6');md['depends_on']=[stage['id']]
            receptor=self.store.project_dir(self.project_id)/'receptor.pdb';receptor.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n')
            ligand=self.store.project_dir(self.project_id)/'ligand.mol2';ligand.write_text('@<TRIPOS>MOLECULE\nLIG\n1 0 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 LIG 0\n@<TRIPOS>BOND\n')
            md['parameters'].update(base_input_path=str(receptor),base_selected_mols=str(ligand))
            self.assertEqual(self.store.result_removal_plan(self.token,self.project_id,[stage,md])['affected'],0)
            self.save([stage,md])
            # Reopen the service exactly as after a desktop restart. The MD
            # block awaits configuration; no upstream job may be submitted.
            service=PipelineService(self.store,worker_python=sys.executable)
            try:second=self.wait(service.submit(self.token,self.project_id))
            finally:service.close()
        self.assertEqual(second['status'],'awaiting_input',second.get('error'))
        self.assertTrue(second['stages'][0]['reused'])
        self.assertEqual(len(calls),1)
        self.assertEqual(before,{p:(file_digest(p),Path(p).stat().st_mtime_ns) for p in before})

    def test_downstream_edit_does_not_reconcile_an_unrelated_historical_request(self):
        stage=new_stage('retrieve_compounds');self.save([stage])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:first=self.finish_selected(service,service.submit(self.token,self.project_id))
            finally:service.close()
        # A historic producer request predates the current saved baseline.
        # Editing only MD is not authorization to clean that upstream history.
        with self.store.connect() as db:
            items=copy.deepcopy(first['stages']);items[0]['configuration']['parameters']['search_term']='OTHER'
            db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(items),first['id']))
        md=new_stage('dynamics_dock6');md['depends_on']=[stage['id']]
        proposal=self.store.result_removal_plan(self.token,self.project_id,[stage,md])
        self.assertEqual(proposal['affected'],0)
        self.save([stage,md]);md['parameters']['amber_score_temperature']=310
        self.assertEqual(self.store.result_removal_plan(self.token,self.project_id,[stage,md])['affected'],0)
        self.assertTrue(all(Path(p).is_file() for p in first['stages'][0]['artifacts']))

    def test_md_source_code_does_not_change_retrieval_contract(self):
        import biomolexplorer.stage_cache as cache
        original=cache.SOURCE_ROOT
        real_read=Path.read_bytes
        def read(path):
            value=real_read(path)
            return value+b'\n# MD-only revision\n' if path==original/'biomolexplorer/md_inputs.py' else value
        operation_digest.cache_clear();before=operation_digest('retrieve_compounds')
        with patch.object(Path,'read_bytes',read):
            operation_digest.cache_clear();self.assertEqual(operation_digest('retrieve_compounds'),before)
        operation_digest.cache_clear()
