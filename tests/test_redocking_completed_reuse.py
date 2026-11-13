"""Reuse complete redocking reports without retrying unchanged excluded pairs."""
import copy
import csv
import json
import time
import unittest
from pathlib import Path
from unittest.mock import patch
from uuid import uuid4

from biomolexplorer.catalog import new_stage
from biomolexplorer.docking_data import write_csv
from biomolexplorer.redocking_reuse import accounted_preparation_exclusions
from biomolexplorer.stage_cache import artifact_manifest
import test_pipeline_execution as execution


class RedockingCompletedReuseTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait

    def fixture(self):
        root=self.store.project_dir(self.project_id)
        source=root/'inputs';(source/'Target').mkdir(parents=True)
        pdb=source/'Target/1ABC.pdb';pdb.write_text('original receptor and ligand structure')
        stage=new_stage('redocking')
        stage['parameters'].update(base_input_path=str(source),target='Target',
            pdb_codes=[['1ABC','LIG',1,'A'],['1ABC','BAD',2,'A']])
        self.save([stage])
        output=root/'completed';output.mkdir()
        rows=[dict(PDB_CODE='1ABC',LIGAND='LIG',RESNUM=1,CHAIN='A',ENGINE='vina',RMSD=.25,
                   STATUS='succeeded',FAILED_STEP='',ERROR=''),
              dict(PDB_CODE='1ABC',LIGAND='BAD',RESNUM=2,CHAIN='A',ENGINE='vina',RMSD='',
                   STATUS='preparation_failed',FAILED_STEP='Preparação do complexo',ERROR='Invalid prepared charges')]
        write_csv(output/'pdb_codes.csv',rows)
        write_csv(output/'redocking_failures.csv',[{key:value for key,value in rows[1].items() if key not in ('ENGINE','RMSD')}])
        artifacts=[str(p) for p in sorted(output.iterdir())]
        entry=dict(id=stage['id'],operation='redocking',name=stage['name'],configuration=copy.deepcopy(stage),
                   status='succeeded',redocking_failures=1,artifacts=artifacts,artifact_manifest=artifact_manifest(artifacts),
                   input_files={'base_input_path':[str(pdb)]},
                   input_manifest={'selected:base_input_path':list(artifact_manifest([pdb]).items())},cache_key='previous-key')
        rid=uuid4().hex
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',(rid,self.project_id,self.store.user(self.token)['id'],
                'succeeded',json.dumps([entry]),time.time(),time.time(),None))
        return stage,entry,pdb

    def test_complete_unchanged_stage_reuses_after_restart_without_any_worker(self):
        from biomolexplorer.pipeline import PipelineService
        stage,entry,pdb=self.fixture()
        original={p:Path(p).read_bytes() for p in entry['artifacts']}
        with patch('biomolexplorer.pipeline.JobManager',side_effect=AssertionError('Redocking worker restarted')):
            for attempt in range(2):
                service=PipelineService(self.store)
                try:
                    with patch.object(service,'_resolve',side_effect=AssertionError('Complement inputs regenerated')):
                        run=self.wait(service.submit(self.token,self.project_id))
                finally:service.close()
                self.assertEqual(run['status'],'succeeded',run.get('error'))
                self.assertTrue(run['stages'][0]['reused'])
                self.assertEqual(run['stages'][0]['redocking_failures'],1)
        self.assertEqual(original,{p:Path(p).read_bytes() for p in original})

    def test_pair_settings_input_changes_and_missing_reports_invalidate_reuse(self):
        from biomolexplorer.published_results import verify
        stage,entry,pdb=self.fixture()
        self.assertTrue(verify(self.store,self.project_id,stage,entry,[entry],{}))
        changed=copy.deepcopy(stage)
        changed['parameters']['preparation_pairs']={'1ABC|LIG|1|A':{'ligand':{'minimize':False}}}
        self.assertFalse(verify(self.store,self.project_id,changed,entry,[entry],{}))
        changed=copy.deepcopy(stage);changed['parameters']['pdb_codes'].append(['1ABC','NEW',3,'A'])
        self.assertFalse(verify(self.store,self.project_id,changed,entry,[entry],{}))
        pdb.write_text('changed receptor')
        self.assertFalse(verify(self.store,self.project_id,stage,entry,[entry],{}))
        Path(entry['artifacts'][0]).unlink()
        self.assertFalse(accounted_preparation_exclusions(entry,self.store.project_dir(self.project_id)))

    def test_engine_failure_cannot_be_adopted_as_preparation_exclusion(self):
        stage,entry,pdb=self.fixture()
        for name in ('pdb_codes.csv','redocking_failures.csv'):
            path=Path(entry['artifacts'][0]).parent/name
            with path.open() as stream:rows=list(csv.DictReader(stream))
            for row in rows:
                if row['STATUS']!='succeeded':row['FAILED_STEP']='Redocking Vina / RMSD'
            write_csv(path,rows)
        entry['artifact_manifest']=artifact_manifest(entry['artifacts'])
        self.assertFalse(accounted_preparation_exclusions(entry,self.store.project_dir(self.project_id)))

    def test_exact_stage_and_batch_checkpoints_use_the_same_validation(self):
        from biomolexplorer.pipeline import PipelineService
        stage,entry,pdb=self.fixture()
        service=object.__new__(PipelineService);service.store=self.store
        self.assertIsNotNone(service._cached_stage(self.project_id,'user',stage,'previous-key','software'))
        self.assertIsNotNone(service._cached_batch(self.project_id,stage['id'],'previous-key',unit=stage))


if __name__=='__main__':
    unittest.main()
