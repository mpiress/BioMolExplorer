"""Incremental execution through the real project store and pipeline scheduler."""
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import test_pipeline_execution as execution
from biomolexplorer.catalog import new_stage,template_names
from biomolexplorer.compound_tables import CompoundTables
from biomolexplorer.pipeline import PipelineService,select_input
from biomolexplorer.templates import RESOURCE_ROOT


class StageCacheTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    finish_selected=execution.PipelineExecutionTests.finish_selected
    upload=execution.PipelineExecutionTests.upload
    imported=execution.PipelineExecutionTests.imported
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def run_twice(self,stages,change,force=False):
        self.save(stages)
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:
                first=self.finish_selected(service,service.submit(self.token,self.project_id))
                self.assertEqual(first['status'],'succeeded',first['error'])
                change(stages,first)
                self.save(stages)
                second=self.finish_selected(service,service.submit(self.token,self.project_id,reuse_results=not force))
                self.assertEqual(second['status'],'succeeded',second['error'])
                return first,second,calls
            finally:service.close()

    def test_adding_connected_block_runs_only_new_stage(self):
        source=new_stage('retrieve_compounds')
        def add(stages,run):
            downstream=new_stage('admet')
            downstream['bindings']['base_input_path']={'stage':source['id']}
            stages.append(downstream)
        first,second,calls=self.run_twice([source],add)
        self.assertEqual([c[0] for c in calls],['retrieve_compounds','admet'])
        self.assertTrue(second['stages'][0]['reused'])
        self.assertEqual(second['stages'][0]['reused_from_run'],first['id'])
        self.assertEqual(second['stages'][0]['artifacts'],first['stages'][0]['artifacts'])

    def test_adding_unconnected_block_and_moving_canvas_preserves_existing_results(self):
        source=new_stage('retrieve_compounds')
        def add(stages,run):
            source.update(name='Renamed block',position={'x':750,'y':300})
            stages.append(new_stage('retrieve_compounds'))
        _,second,calls=self.run_twice([source],add)
        self.assertEqual(len(calls),2)
        self.assertTrue(second['stages'][0]['reused'])
        self.assertNotIn('reused',second['stages'][1])

    def test_configuration_change_reruns_changed_block(self):
        def change(stages,run):stages[0]['parameters']['search_term']='CHEMBL221'
        _,second,calls=self.run_twice([new_stage('retrieve_compounds')],change)
        self.assertEqual(len(calls),2)
        self.assertNotIn('reused',second['stages'][0])

    def test_template_change_invalidates_previous_output(self):
        source=new_stage('retrieve_compounds')
        def change(stages,run):
            name=template_names(source['operation'])[0]
            source['templates'][name]=(RESOURCE_ROOT/name).read_text()+'\n'
        _,second,calls=self.run_twice([source],change)
        self.assertEqual(len(calls),2)
        self.assertNotIn('reused',second['stages'][0])

    def test_explicit_compounds_filename_selects_integrated_dataset(self):
        root=self.store.project_dir(self.project_id)
        pubchem=root/'artifacts'/'PubChem'/'similars'/'CHEMBL220'/'compounds.csv'
        integrated=root/'artifacts'/'compounds'/'CHEMBL220'/'compounds.csv'
        for path in (pubchem,integrated):
            path.parent.mkdir(parents=True);path.write_text('molecule_chembl_id,canonical_smiles\nCHEMBL1,CCO\n')
        self.assertEqual(select_input([str(pubchem),str(integrated)],'base_input_path','compounds.csv'),integrated.parent)
        self.assertEqual(select_input([str(pubchem),str(integrated)],'base_input_path','PubChem/similars/CHEMBL220/compounds.csv'),pubchem.parent)

    def test_missing_or_modified_output_is_never_reused(self):
        for delete in (True,False):
            with self.subTest(delete=delete):
                def change(stages,run):
                    path=Path(run['stages'][0]['artifacts'][0])
                    if delete:path.unlink()
                    else:path.write_text('molecule_chembl_id,canonical_smiles\nCHANGED,O\n')
                _,second,calls=self.run_twice([new_stage('retrieve_compounds')],change)
                self.assertEqual(len(calls),2)
                self.assertNotIn('reused',second['stages'][0])

    def test_forced_execution_bypasses_valid_results(self):
        _,second,calls=self.run_twice([new_stage('retrieve_compounds')],lambda stages,run:None,force=True)
        self.assertEqual(len(calls),2)
        self.assertNotIn('reused',second['stages'][0])

    def test_existing_successful_legacy_results_are_adopted(self):
        def legacy(stages,run):
            for stage in run['stages']:
                stage.pop('cache_key',None);stage.pop('artifact_manifest',None)
            with self.store.connect() as db:
                db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(run['stages']),run['id']))
        first,second,calls=self.run_twice([new_stage('retrieve_compounds')],legacy)
        self.assertEqual(len(calls),1)
        self.assertTrue(second['stages'][0]['reused'])
        self.assertEqual(second['stages'][0]['reused_from_run'],first['id'])

    def test_changed_uploaded_input_reruns_import_and_consumers(self):
        asset=self.upload()
        source=self.imported(asset)
        downstream=new_stage('admet');downstream['bindings']['base_input_path']={'stage':source['id']}
        def change(stages,run):
            self.store.asset_path(self.token,self.project_id,asset).write_text('name,smiles\nWATER,O\n')
        first,second,calls=self.run_twice([source,downstream],change)
        self.assertEqual(len(calls),2)
        self.assertTrue(all(not s.get('reused') for s in second['stages']))
        self.assertNotEqual(second['stages'][0]['cache_key'],first['stages'][0]['cache_key'])

    def test_identical_connected_pipeline_reuses_every_stage(self):
        source=new_stage('retrieve_compounds');downstream=new_stage('admet')
        downstream['bindings']['base_input_path']={'stage':source['id']}
        _,second,calls=self.run_twice([source,downstream],lambda stages,run:None)
        self.assertEqual(len(calls),2)
        self.assertTrue(all(s.get('reused') for s in second['stages']))

    def test_curation_reuses_collection_but_recalculates_consumers(self):
        source=new_stage('retrieve_compounds');consumer=new_stage('admet')
        consumer['bindings']['base_input_path']={'stage':source['id']}
        self.save([source,consumer])
        base,calls=self.immediate_manager()
        class Manager(base):
            def submit(manager,operation,parameters):
                job=super().submit(operation,parameters)
                if operation=='retrieve_compounds':
                    original=Path(job['result']['artifacts'][0])
                    output=manager.config.workspace/'artifacts'/'compounds'/'CHEMBL220'/'compounds.csv'
                    output.parent.mkdir(parents=True);original.replace(output)
                    job['result']['artifacts']=[str(output)]
                return job
        with patch('biomolexplorer.pipeline.JobManager',Manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:
                first=self.finish_selected(service,service.submit(self.token,self.project_id))
                self.assertEqual(first['status'],'succeeded',first['error'])
                path=first['stages'][0]['artifacts'][0]
                tables=CompoundTables(self.store)
                page=tables.page(self.token,self.project_id,first['id'],source['id'],path)
                tables.remove(self.token,self.project_id,first['id'],source['id'],path,0,page['version'])
                second=self.finish_selected(service,service.submit(self.token,self.project_id))
            finally:service.close()
        self.assertEqual(second['status'],'succeeded',second['error'])
        self.assertTrue(second['stages'][0]['reused'])
        self.assertNotIn('reused',second['stages'][1])
        self.assertEqual([c[0] for c in calls],['retrieve_compounds','admet','admet'])
        self.assertEqual(tables.page(self.token,self.project_id,second['id'],source['id'],path)['total'],0)

    def test_legacy_consumer_with_recently_edited_input_is_not_adopted(self):
        asset=self.upload();source=self.imported(asset);downstream=new_stage('admet')
        downstream['bindings']['base_input_path']={'stage':source['id']}
        def change(stages,run):
            for stage in run['stages']:stage.pop('cache_key',None);stage.pop('artifact_manifest',None)
            with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(run['stages']),run['id']))
            self.store.asset_path(self.token,self.project_id,asset).write_text('name,smiles\nWATER,O\n')
        _,second,calls=self.run_twice([source,downstream],change)
        self.assertEqual(len(calls),2)
        self.assertNotIn('reused',second['stages'][1])
