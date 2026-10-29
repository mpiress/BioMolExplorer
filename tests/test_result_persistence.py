"""Durable incremental execution and confirmed, dependency-aware cleanup."""
import copy
import json
import os
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import test_pipeline_execution as execution
from test_clear_all_results import ClearAllResultsTests
from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.workspace import WorkspaceStore
from biomolexplorer.result_files import ResultFiles
from biomolexplorer import preparation_cache


class CrashRecoveryTests(unittest.TestCase):
    setUp=ClearAllResultsTests.setUp
    def exercise_crash(self,committed):
        script='''
import json,os,sys
from pathlib import Path
from biomolexplorer.workspace import WorkspaceStore
from biomolexplorer.result_transactions import ResultTransaction
root,token,pid,rid,filename,committed=sys.argv[1:]
store=WorkspaceStore(root);source=Path(filename)
with store.change(token,pid,'editor','remove_all_results','crash recovery test') as db:
    journal=ResultTransaction(store,db,pid)
    backup=store.project_dir(pid)/'.history'/('crash-'+journal.id)
    journal.prepare([(source,backup)]);journal.move(source,backup);journal.flush()
    stages=json.loads(db.execute('SELECT stages FROM runs WHERE id=?',(rid,)).fetchone()[0])
    for stage in stages:stage['artifacts']=[]
    db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),rid))
    db.execute("UPDATE projects SET pipeline='[]' WHERE id=?",(pid,))
    from biomolexplorer.project_state import mirror
    mirror(store,db,pid)
    if committed=='False':os._exit(17)
os._exit(17)
'''
        result=subprocess.run([sys.executable,'-c',script,str(self.store.root),self.token,self.pid,self.rid,str(self.generated),str(committed)],capture_output=True,text=True)
        self.assertEqual(result.returncode,17,result.stderr)
        reopened=WorkspaceStore(self.store.root)
        stages=reopened.get_run(self.token,self.rid)['stages']
        self.assertEqual(self.generated.exists(),not committed)
        self.assertEqual(stages,[dict(s,artifacts=[]) for s in self.stages] if committed else self.stages)
        self.assertTrue(self.original.is_file())
        expected=[] if committed else self.pipeline
        self.assertEqual(reopened.project(self.token,self.pid)['pipeline'],expected)
        self.assertEqual(json.loads((self.root/'project.json').read_text())['pipeline'],expected)
        self.assertEqual(list((self.root/'.history/.result-transactions').glob('*.json')),[])
        self.assertEqual(list((self.root/'.history').glob('crash-*')),[])
        with reopened.connect() as db:self.assertEqual(db.execute('SELECT count(*) FROM result_file_transactions').fetchone()[0],0)
        WorkspaceStore(self.store.root) # Recovery is idempotent.

    def test_death_before_commit_restores_files_and_metadata(self):self.exercise_crash(False)

    def test_death_after_commit_finishes_file_cleanup(self):self.exercise_crash(True)


class PersistentBatchTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    finish_selected=execution.PipelineExecutionTests.finish_selected
    upload=execution.PipelineExecutionTests.upload
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def test_removing_one_imported_source_prunes_only_its_copies_and_dependent_batch(self):
        a=self.upload('a.csv');b=self.upload('b.csv','name,smiles\nB,CCC\n')
        source=new_stage('import_results');source['parameters']['asset_ids']=[a,b]
        stage=new_stage('fingerprints');stage['input_processing']='individual'
        stage['bindings']['base_input_path']={'stage':source['id']};self.save([source,stage])
        manager,calls=self.immediate_manager()
        from test_input_processing import InputProcessingTests
        def finish(service,run):
            run=self.wait(run)
            while run['status']=='awaiting_input':run=InputProcessingTests.confirm_all(self,service,run,'individual')
            return run
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:
                first=finish(service,service.submit(self.token,self.project_id))
                self.assertEqual(first['status'],'succeeded',first['error']);self.assertEqual(len(calls),2)
                source['parameters']['asset_ids']=[b]
                proposal=self.store.result_removal_plan(self.token,self.project_id,[source,stage])
                self.assertGreater(proposal['affected'],0)
                self.assertTrue(all(Path(p).exists() for p in first['stages'][0]['artifacts']))
                project=self.store.project(self.token,self.project_id)
                self.store.save_pipeline(self.token,self.project_id,[source,stage],project['revision'],prune_results=True)
                retained=self.store.get_run(self.token,first['id'])['stages']
                self.assertEqual([Path(p).name for p in retained[0]['artifacts']],['b.csv'])
                self.assertEqual(len(retained[1]['batches']),1)
                for asset in (a,b):self.assertTrue(self.store.asset_path(self.token,self.project_id,asset).is_file())
                second=finish(service,service.submit(self.token,self.project_id))
                self.assertEqual(second['status'],'succeeded',second['error']);self.assertEqual(len(calls),2)
            finally:service.close()

    def test_verified_v1_cache_survives_migration_and_growth_without_recomputation(self):
        self._exercise_v1_migration(reuse_first=True)

    def test_v1_batches_are_adopted_when_new_sources_are_added_immediately(self):
        self._exercise_v1_migration(reuse_first=False)

    def _exercise_v1_migration(self,reuse_first):
        from biomolexplorer.templates import RESOURCE_ROOT
        from biomolexplorer.pipeline import execution_stage,execution_context
        source_build=json.loads((RESOURCE_ROOT/'legacy_cache_compatibility.json').read_text())['source_builds'][0]
        a=self.upload('a.csv');b=self.upload('b.csv','name,smiles\nB,CCC\n')
        stage=new_stage('fingerprints');stage['input_processing']='individual';stage['bindings']['base_input_path']={'asset':a};self.save([stage])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:
                first=self.finish_selected(service,service.submit(self.token,self.project_id))
                old=copy.deepcopy(first['stages']);user=self.store.user(self.token)['id'];provenance=source_build+':'+sys.executable
                config=execution_stage(old[0]['configuration'],old[0].get('input_mode'))
                old[0]['cache_key']=service._legacy_cache_key(self.project_id,user,config,{},provenance,pipeline=execution_context(old))
                old[0]['cache_software']=provenance
                for field in ('processing_version','processing_runtime','producer_software'):old[0].pop(field,None)
                for batch in old[0]['batches']:
                    for field in ('cache_key','status','artifact_manifest','processing_version','processing_runtime','producer_software'):batch.pop(field,None)
                with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(old),first['id']))
            finally:service.close()
            service=PipelineService(WorkspaceStore(self.store.root),worker_python=sys.executable)
            try:
                if reuse_first:
                    second=self.finish_selected(service,service.submit(self.token,self.project_id))
                    self.assertEqual(second['status'],'succeeded',second['error']);self.assertEqual(len(calls),1)
                    self.assertTrue(second['stages'][0]['reused']);self.assertIsNone(second['stages'][0]['processing_version'])
                stage['bindings']['base_input_path']={'sources':[{'asset':a},{'asset':b}]};self.save([stage])
                third=self.finish_selected(service,service.submit(self.token,self.project_id))
                self.assertEqual(third['status'],'succeeded',third['error']);self.assertEqual(len(calls),2)
                self.assertTrue(third['stages'][0]['batches'][0]['reused'])
            finally:service.close()

    def test_restart_growth_removal_and_selective_configuration_invalidation(self):
        a=self.upload('a.csv');b=self.upload('b.csv','name,smiles\nB,CCC\n')
        stage=new_stage('fingerprints');stage['input_processing']='individual'
        stage['bindings']['base_input_path']={'asset':a};self.save([stage])
        manager,calls=self.immediate_manager()
        def run():
            service=PipelineService(WorkspaceStore(self.store.root),worker_python=sys.executable)
            try:return self.finish_selected(service,service.submit(self.token,self.project_id))
            finally:service.close()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            first=run();self.assertEqual(first['status'],'succeeded',first['error'])
            paths=first['stages'][0]['artifacts'];self.assertEqual(len(calls),1)
            second=run();self.assertTrue(second['stages'][0]['reused']);self.assertEqual(len(calls),1)
            stage['bindings']['base_input_path']={'sources':[{'asset':a},{'asset':b}]};self.save([stage])
            third=run();self.assertEqual(third['status'],'succeeded',third['error'])
            self.assertEqual(len(calls),2);self.assertTrue(third['stages'][0]['batches'][0]['reused'])
            self.assertTrue(all(Path(p).is_file() for p in paths))
            stage['bindings']['base_input_path']={'asset':b}
            proposal=self.store.result_removal_plan(self.token,self.project_id,[stage]);self.assertGreater(proposal['affected'],0)
            # A destructive edit cannot become effective before confirmation.
            with self.assertRaisesRegex(ValueError,'Confirme'):self.save([stage])
            self.assertTrue(all(Path(p).is_file() for p in paths))
            project=self.store.project(self.token,self.project_id)
            self.store.save_pipeline(self.token,self.project_id,[stage],project['revision'],prune_results=True)
            fourth=run();self.assertEqual(fourth['status'],'succeeded',fourth['error'])
            self.assertEqual(len(calls),2);self.assertEqual(len(fourth['stages'][0]['artifacts']),1)
            stage['parameters']['radius']=3
            project=self.store.project(self.token,self.project_id)
            self.store.save_pipeline(self.token,self.project_id,[stage],project['revision'],prune_results=True)
            fifth=run()
            self.assertEqual(fifth['status'],'succeeded',fifth['error']);self.assertEqual(len(calls),3)

    def test_completed_batches_survive_a_later_failed_batch_and_service_restart(self):
        a=self.upload('a.csv');b=self.upload('b.csv','name,smiles\nB,CCC\n')
        stage=new_stage('fingerprints');stage['input_processing']='individual'
        stage['bindings']['base_input_path']={'sources':[{'asset':a},{'asset':b}]};self.save([stage])
        base,calls=self.immediate_manager();fail={'enabled':True}
        class Manager(base):
            def submit(manager,operation,parameters):
                job=super().submit(operation,parameters)
                if fail['enabled'] and parameters.get('files')==['b.csv']:
                    job.update(status='failed',error='fixture failure')
                return job
        with patch('biomolexplorer.pipeline.JobManager',Manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:first=self.finish_selected(service,service.submit(self.token,self.project_id))
            finally:service.close()
            self.assertEqual(first['status'],'failed');self.assertEqual(len(first['stages'][0]['batches']),1)
            fail['enabled']=False
            service=PipelineService(WorkspaceStore(self.store.root),worker_python=sys.executable)
            try:second=self.finish_selected(service,service.submit(self.token,self.project_id))
            finally:service.close()
            self.assertEqual(second['status'],'succeeded',second['error']);self.assertEqual(len(calls),3)
            self.assertTrue(second['stages'][0]['batches'][0]['reused'])

    def test_real_changed_dataset_reuses_other_fingerprints_similarity_and_graphs(self):
        from test_input_processing import InputProcessingTests
        from biomolexplorer.jobs import JobManager
        assets=[self.upload('alpha.csv','name,smiles\nA,CCO\nB,CCCO\n'),
                self.upload('beta.csv','name,smiles\nC,c1ccccc1\nD,c1ccncc1\n')]
        source=new_stage('import_results');source['parameters']['asset_ids']=assets
        fp,sim,graph=[new_stage(op) for op in ('fingerprints','similarity','graphs')]
        fp['bindings']['base_input_path']={'stage':source['id']}
        sim['bindings']['base_input_path']={'stage':fp['id']}
        graph['bindings']['similarity_path']={'stage':sim['id']}
        self.save([source,fp,sim,graph]);calls=[]
        class Counting(JobManager):
            def submit(manager,operation,parameters):
                calls.append(operation);return super().submit(operation,parameters)
        def finish(service,run):
            run=self.wait(run)
            while run['status']=='awaiting_input':run=InputProcessingTests.confirm_all(self,service,run,'individual')
            self.assertEqual(run['status'],'succeeded',run['error']);return run
        with patch('biomolexplorer.pipeline.JobManager',Counting):
            service=PipelineService(self.store,worker_python=sys.executable,cpu_workers=1)
            try:first=finish(service,service.submit(self.token,self.project_id))
            finally:service.close()
            self.assertEqual(len(calls),6)
            self.store.asset_path(self.token,self.project_id,assets[0]).write_text('name,smiles\nA,CCN\nB,CCCO\n')
            project=self.store.project(self.token,self.project_id)
            self.store.save_pipeline(self.token,self.project_id,project['pipeline'],project['revision'],prune_results=True)
            service=PipelineService(WorkspaceStore(self.store.root),worker_python=sys.executable,cpu_workers=1)
            try:second=finish(service,service.submit(self.token,self.project_id))
            finally:service.close()
        self.assertEqual(len(calls),9,calls)
        for item in second['stages'][1:]:self.assertEqual(sum(b['reused'] for b in item['batches']),1)
        self.assertTrue(all(Path(p).is_file() for item in first['stages'][1:] for p in item['batches'][1]['artifacts']))

    def test_removing_one_batch_output_preserves_reuse_of_other_batches(self):
        a=self.upload('a.csv');b=self.upload('b.csv','name,smiles\nB,CCC\n')
        stage=new_stage('fingerprints');stage['input_processing']='individual'
        stage['bindings']['base_input_path']={'sources':[{'asset':a},{'asset':b}]};self.save([stage])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,worker_python=sys.executable)
            try:
                first=self.finish_selected(service,service.submit(self.token,self.project_id))
                removed=first['stages'][0]['batches'][0]['artifacts'][0]
                other=first['stages'][0]['batches'][1]['artifacts'][0]
                ResultFiles(self.store).remove(self.token,self.project_id,first['id'],stage['id'],removed)
                self.assertTrue(Path(other).is_file())
                second=self.finish_selected(service,service.submit(self.token,self.project_id))
                self.assertEqual(second['status'],'succeeded',second['error']);self.assertEqual(len(calls),3)
                self.assertTrue(second['stages'][0]['batches'][1]['reused'])
            finally:service.close()


class CleanupContractTests(unittest.TestCase):
    setUp=ClearAllResultsTests.setUp

    def test_remove_all_preserves_every_configuration_and_source(self):
        before=self.store.project(self.token,self.pid)
        cached=self.root/'.biomolexplorer/cache/preparation/receptors/key/manifest.json';cached.parent.mkdir(parents=True);cached.write_text('{}')
        molecular=self.root/'.biomolexplorer/cache/calculations/molecular/fingerprints.sqlite3';molecular.parent.mkdir(parents=True);molecular.write_text('computed cache')
        ResultFiles(self.store).remove_project_results(self.token,self.pid)
        self.assertEqual(self.store.project(self.token,self.pid)['pipeline'],before['pipeline'])
        self.assertFalse(self.generated.exists());self.assertTrue(self.original.is_file())
        self.assertFalse(cached.exists());self.assertFalse(molecular.exists());self.assertTrue(self.ccd.exists())
        reopened=WorkspaceStore(self.store.root)
        self.assertEqual(reopened.project(self.token,self.pid)['pipeline'],before['pipeline'])

    def test_clear_all_resets_scientific_defaults_and_import_selections_preserving_graph(self):
        before=self.store.project(self.token,self.pid)
        ResultFiles(self.store).clear_project_results(self.token,self.pid)
        after=self.store.project(self.token,self.pid)
        self.assertEqual([s['id'] for s in after['pipeline']],[s['id'] for s in before['pipeline']])
        for stage in after['pipeline']:
            self.assertEqual(stage['parameters'],new_stage(stage['operation'])['parameters'])
            self.assertEqual(stage['templates'],{})
        self.assertEqual(after['pipeline'][2]['bindings']['base_input_path']['stage'],before['pipeline'][1]['id'])
        self.assertTrue(self.original.exists());self.assertTrue(self.log.exists())

    def test_confirmed_block_removal_preserves_shared_artifact_used_elsewhere(self):
        pipeline=copy.deepcopy(self.pipeline);removed=pipeline.pop(2)
        # The fingerprint stage has an independent, still valid reference to the bytes.
        self.stages[3]['input_files']={'base_input_path':[str(self.generated)]}
        with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(self.stages),self.rid))
        project=self.store.project(self.token,self.pid)
        proposal=self.store.result_removal_plan(self.token,self.pid,pipeline)
        self.assertGreater(proposal['shared_preserved'],0)
        self.store.save_pipeline(self.token,self.pid,pipeline,project['revision'],prune_results=True)
        self.assertTrue(self.generated.exists());self.assertTrue(self.original.exists())
        old=next(s for s in self.store.get_run(self.token,self.rid)['stages'] if s['id']==removed['id'])
        self.assertEqual(old['artifacts'],[])

    def test_failed_confirmed_removal_rolls_back_files_and_pipeline(self):
        before=self.store.project(self.token,self.pid);pipeline=copy.deepcopy(self.pipeline);pipeline.pop(2)
        with patch('biomolexplorer.project_state.record',side_effect=RuntimeError('audit failed')):
            with self.assertRaises(RuntimeError):self.store.save_pipeline(self.token,self.pid,pipeline,before['revision'],prune_results=True)
        self.assertTrue(self.generated.exists())
        self.assertEqual(self.store.project(self.token,self.pid)['pipeline'],before['pipeline'])

    def test_source_removal_cascades_only_through_its_actual_artifact_lineage(self):
        import time
        from uuid import uuid4
        a=self.asset
        ticket=self.store.prepare_upload(self.token,self.pid,'other.csv','compounds');(self.store.staging/ticket).write_text('name,smiles\nB,CCC\n')
        b=self.store.finish_upload(self.token,ticket)
        fp,sim,graph=[new_stage(op) for op in ('fingerprints','similarity','graphs')]
        fp['bindings']['base_input_path']={'sources':[{'asset':a},{'asset':b}]}
        paths=[];items=[];previous=None
        for stage,field in [(fp,'base_input_path'),(sim,'base_input_path'),(graph,'similarity_path')]:
            directory=self.root/'lineage'/stage['id'];directory.mkdir(parents=True)
            outputs=[directory/'a.csv',directory/'b.csv']
            for path in outputs:path.write_text('persisted '+str(path))
            if previous:stage['bindings'][field]={'sources':[{'stage':previous[0]['id'],'selector':p.name} for p in previous[1]]}
            batches=[]
            for index,path in enumerate(outputs):
                config=copy.deepcopy(stage)
                config['bindings'][field]={'asset':(a,b)[index]} if previous is None else {'stage':previous[0]['id'],'selector':previous[1][index].name}
                inputs=[str(self.store.asset_path(self.token,self.pid,(a,b)[index]))] if previous is None else [str(previous[1][index])]
                batches.append(dict(configuration=config,status='succeeded',artifacts=[str(path)],input_files={field:inputs}))
            items.append(dict(stage,configuration=copy.deepcopy(stage),status='succeeded',artifacts=list(map(str,outputs)),batches=batches))
            paths.append(outputs);previous=(stage,outputs)
        project=self.store.project(self.token,self.pid);pipeline=[fp,sim,graph]
        self.store.save_pipeline(self.token,self.pid,pipeline,project['revision'],prune_results=True)
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES(?,?,?,?,?,?,?,?)',(uuid4().hex,self.pid,self.store.user(self.token)['id'],'succeeded',json.dumps(items),time.time(),time.time(),None))
        pipeline=copy.deepcopy(pipeline);pipeline[0]['bindings']['base_input_path']={'asset':b}
        project=self.store.project(self.token,self.pid)
        self.store.save_pipeline(self.token,self.pid,pipeline,project['revision'],prune_results=True)
        for first,second in paths:self.assertFalse(first.exists());self.assertTrue(second.exists())
        saved=self.store.project(self.token,self.pid)['pipeline']
        self.assertEqual(saved[1]['bindings']['base_input_path']['selector'],'b.csv')
        self.assertEqual(saved[2]['bindings']['similarity_path']['selector'],'b.csv')


class SharedPreparationTests(unittest.TestCase):
    def test_shared_receptor_is_verified_and_independent_of_ligand_and_engine(self):
        with tempfile.TemporaryDirectory() as folder,patch.dict(os.environ,{'BIOMOL_CACHE_DIR':folder}):
            root=Path(folder);a=root/'first';b=root/'second';a.mkdir();b.mkdir()
            for suffix in ('.complex.pdb','.dockprep.mol2','.noH.pdb'):(a/('1ABC_A'+suffix)).write_text('validated native bytes '+suffix)
            preparation_cache.publish('key',a,'1ABC_A')
            self.assertTrue(preparation_cache.restore('key',b,'1ABC_A'))
            self.assertEqual((a/'1ABC_A.dockprep.mol2').read_bytes(),(b/'1ABC_A.dockprep.mol2').read_bytes())
            entry=preparation_cache.root()/'key'/'receptor.dockprep.mol2';entry.write_text('tampered')
            self.assertFalse(preparation_cache.restore('key',b,'1ABC_A'))
            preparation_cache.publish('key',a,'1ABC_A')
            self.assertTrue(preparation_cache.restore('key',b,'1ABC_A'))
            with patch.dict(os.environ,{'BIOMOL_REUSE_RESULTS':'0'}):self.assertIsNone(preparation_cache.root())

    def test_receptor_key_changes_only_for_receptor_requirements(self):
        with tempfile.TemporaryDirectory() as folder:
            source=Path(folder)/'source.pdb';source.write_text('crystal')
            a=['1ABC','LIG',1,'A'];b=['1ABC','OTHER',2,'A']
            key=preparation_cache.receptor_key(source,a,{'ligand':{'charge_type':'am1'}},7.4)
            self.assertEqual(key,preparation_cache.receptor_key(source,b,{'ligand':{'charge_type':'gas'}},7.4))
            self.assertNotEqual(key,preparation_cache.receptor_key(source,a,{'receptor':{'minimize':False}},7.4))
            self.assertNotEqual(key,preparation_cache.receptor_key(source,a,{'cofactors':['FAD']},7.4))
            source.write_text('modified crystal');self.assertNotEqual(key,preparation_cache.receptor_key(source,a,{},7.4))


class MolecularCacheTests(unittest.TestCase):
    def test_new_molecules_and_algorithm_changes_are_selective_and_force_is_honored(self):
        import pandas as pd
        from kernel.descriptors import Descriptors
        from biomolexplorer.molecular_cache import fingerprints
        with tempfile.TemporaryDirectory() as folder,patch.dict(os.environ,{'BIOMOL_CACHE_DIR':folder,'BIOMOL_REUSE_RESULTS':'1'}):
            calls=[]
            def compute(frame,radius=2):
                calls.extend(frame['smiles'].tolist())
                return pd.DataFrame([Descriptors.morgan_worker(row.smiles,row.molecule_chembl_id,radius,128) for row in frame.itertuples()])
            first=pd.DataFrame(dict(smiles=['CCO','OCC'],molecule_chembl_id=['A','B']))
            expected=compute(first);calls.clear()
            values=fingerprints(first,compute,'morgan',2,128)
            self.assertEqual(values['fingerprint'].tolist(),expected['fingerprint'].tolist());self.assertEqual(len(calls),1)
            grown=pd.concat([first,pd.DataFrame(dict(smiles=['CCC'],molecule_chembl_id=['C']))],ignore_index=True)
            fingerprints(grown,compute,'morgan',2,128);self.assertEqual(len(calls),2)
            fingerprints(grown,compute,'morgan',2,128);self.assertEqual(len(calls),2)
            fingerprints(grown,lambda rows:compute(rows,3),'morgan',3,128);self.assertEqual(len(calls),4)
            with patch.dict(os.environ,{'BIOMOL_REUSE_RESULTS':'0'}):fingerprints(grown,compute,'morgan',2,128)
            self.assertEqual(len(calls),7)


class ScientificVersionTests(unittest.TestCase):
    def test_rmsd_changes_do_not_invalidate_fingerprints_but_fingerprint_changes_do(self):
        import ast
        import shutil
        from biomolexplorer import stage_cache
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            files=['wrappers/molecular_analyzer.py','kernel/descriptors.py']+[
                'biomolexplorer/'+name for name in ('molecular_parsing.py','molecule_quality.py','input_validation.py','storage.py','molecular_cache.py')]
            for name in files:
                destination=root/name;destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(stage_cache.SOURCE_ROOT/name,destination)
            try:
                with patch.object(stage_cache,'SOURCE_ROOT',root):
                    stage_cache.operation_digest.cache_clear();before=stage_cache.operation_digest('fingerprints')
                    path=root/'kernel/descriptors.py';tree=ast.parse(path.read_text())
                    cls=next(n for n in tree.body if getattr(n,'name',None)=='Descriptors')
                    method=next(n for n in cls.body if getattr(n,'name',None)=='calcRMSD')
                    method.body=[ast.Raise(exc=ast.Call(func=ast.Name(id='RuntimeError',ctx=ast.Load()),args=[ast.Constant('changed RMSD')],keywords=[]))]
                    path.write_text(ast.unparse(ast.fix_missing_locations(tree)))
                    stage_cache.operation_digest.cache_clear();self.assertEqual(stage_cache.operation_digest('fingerprints'),before)
                    method=next(n for n in cls.body if getattr(n,'name',None)=='morgan_worker');method.body=[ast.Return(value=ast.Dict(keys=[],values=[]))]
                    path.write_text(ast.unparse(ast.fix_missing_locations(tree)))
                    stage_cache.operation_digest.cache_clear();self.assertNotEqual(stage_cache.operation_digest('fingerprints'),before)
            finally:stage_cache.operation_digest.cache_clear()

    def test_unused_engine_options_do_not_invalidate_redocking(self):
        from biomolexplorer.stage_cache import stage_key
        stage=new_stage('redocking');params=copy.deepcopy(stage['parameters'])
        first=stage_key(stage,params,{},'code:python')
        params['density']=.8;self.assertEqual(stage_key(stage,params,{},'code:python'),first)
        params['exhaustiveness']=99;self.assertNotEqual(stage_key(stage,params,{},'code:python'),first)
        params['redocking_engine']='dock6';second=stage_key(stage,params,{},'code:python')
        params['sizeof_box']=[50,50,50];self.assertEqual(stage_key(stage,params,{},'code:python'),second)
        params['radius']=2.;self.assertNotEqual(stage_key(stage,params,{},'code:python'),second)

    def test_unrelated_operation_change_does_not_invalidate_other_science(self):
        from biomolexplorer import stage_cache,operations
        from biomolexplorer.contracts import OperationSpec
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);(root/'wrappers').mkdir()
            (root/'wrappers/a.py').write_text('from .shared import constant\n')
            (root/'wrappers/b.py').write_text('constant=1\n')
            (root/'wrappers/shared.py').write_text('constant=1\n')
            specs={'a':OperationSpec('wrappers.a','run',()),'b':OperationSpec('wrappers.b','run',())}
            try:
                with patch.object(stage_cache,'SOURCE_ROOT',root),patch.object(operations,'OPERATIONS',specs):
                    stage_cache.operation_digest.cache_clear();a=stage_cache.operation_digest('a');b=stage_cache.operation_digest('b')
                    (root/'wrappers/b.py').write_text('constant=2\n')
                    stage_cache.operation_digest.cache_clear()
                    self.assertEqual(stage_cache.operation_digest('a'),a);self.assertNotEqual(stage_cache.operation_digest('b'),b)
                    (root/'wrappers/shared.py').write_text('constant=2\n')
                    stage_cache.operation_digest.cache_clear();self.assertNotEqual(stage_cache.operation_digest('a'),a)
            finally:stage_cache.operation_digest.cache_clear()
