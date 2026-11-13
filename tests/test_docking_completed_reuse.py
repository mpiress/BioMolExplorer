"""Completed docking exclusions must not restart preparation or native engines."""
import copy
import csv
import json
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch
from uuid import uuid4

from biomolexplorer.catalog import new_stage
from biomolexplorer.docking_accounting import reconcile
from biomolexplorer.docking_data import write_csv
from biomolexplorer.docking_reuse import accounted_preparation_exclusions
from biomolexplorer.stage_cache import artifact_manifest
import test_pipeline_execution as execution


def completed(root, engine, codes=('A','B'), successful=('A',), failed='B'):
    root.mkdir(parents=True,exist_ok=True)
    dataset=root/'input.csv'
    write_csv(dataset,[dict(molecule_chembl_id=code,canonical_smiles='CCO') for code in codes])
    rows=[dict(molecule_chembl_id=code,canonical_smiles='CCO',engine=engine,
               receptor_id='1ABC_A',score=-7 if engine=='vina' else -20) for code in successful]
    write_csv(root/'docking_results.csv',rows)
    accounting=reconcile(engine,dataset,[['1ABC','LIG',1,'A']],rows,root,
              [dict(compound_identifier=failed,reason='Invalid MOL2 topology')] if failed else [])
    files=[str(root/name) for name in ('docking_results.csv','docking_accounting.json','docking_unprocessed.csv')]
    return dict(operation='docking_'+engine,status='succeeded',docking_failures=int(bool(failed)),
                excluded_records=int(bool(failed)),docking_accounting=accounting,
                artifacts=files,artifact_manifest=artifact_manifest(files)),dataset


class CompletedDockingTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)

    def test_both_engines_reuse_accounted_preparation_exclusions(self):
        for engine in ('vina','dock6'):
            entry,_=completed(self.root/engine,engine)
            self.assertTrue(accounted_preparation_exclusions(entry,'docking_'+engine,self.root))

    def test_engine_failures_unaccounted_pairs_and_tampering_are_rejected(self):
        for scenario in ('engine','unaccounted','tampered','missing'):
            entry,_=completed(self.root/scenario,'vina')
            folder=self.root/scenario
            if scenario in ('engine','unaccounted'):
                path=folder/'docking_unprocessed.csv'
                with path.open() as stream:rows=list(csv.DictReader(stream))
                rows[0]['status']='docking_failed' if scenario=='engine' else 'unaccounted'
                write_csv(path,rows)
                entry['artifact_manifest']=artifact_manifest(entry['artifacts'])
            elif scenario=='tampered':
                (folder/'docking_results.csv').write_text('corrupted')
            else:
                (folder/'docking_accounting.json').unlink()
            self.assertFalse(accounted_preparation_exclusions(entry,'docking_vina',self.root))

    def test_all_batch_reports_are_required_and_failed_stages_never_reuse(self):
        a,_=completed(self.root/'a','dock6')
        b,_=completed(self.root/'b','dock6',codes=('C',),successful=('C',),failed=None)
        entry=dict(a,batches=[a,b])
        self.assertTrue(accounted_preparation_exclusions(entry,'docking_dock6',self.root))
        path=self.root/'b/docking_accounting.json'
        report=json.loads(path.read_text());report['unaccounted_pairs']=1;path.write_text(json.dumps(report))
        b['artifact_manifest']=artifact_manifest(b['artifacts'])
        self.assertFalse(accounted_preparation_exclusions(entry,'docking_dock6',self.root))
        a['status']='failed'
        self.assertFalse(accounted_preparation_exclusions(a,'docking_dock6',self.root))

    def test_quality_exclusions_alone_do_not_invalidate_intact_completed_results(self):
        from biomolexplorer.published_results import intact
        entry,_=completed(self.root/'quality','vina',successful=('A','B'),failed=None)
        entry['excluded_records']=1
        self.assertTrue(intact(None,'project',entry,self.root))


class CompletedPipelineReuseTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait

    def fixture(self):
        root=self.store.project_dir(self.project_id)
        receptors=root/'receptors';receptors.mkdir()
        receptor=receptors/'receptor.pdb';receptor.write_text('receptor bytes')
        vina=new_stage('docking_vina');dock6=new_stage('docking_dock6')
        first,dataset=completed(root/'vina','vina',codes=('A','B','C'),successful=('A','C'),failed='B')
        second,_=completed(root/'dock6','dock6',codes=('A','C'),successful=('A',),failed='C')
        ligands=root/'ligands';ligands.mkdir()
        destination=ligands/'input.csv';dataset.rename(destination);dataset=destination
        for stage in (vina,dock6):
            stage['parameters'].update(base_input_path=str(receptors),target='Target',
                                       base_selected_mols=str(dataset.parent),mol_filename='input')
        dock6['bindings']['base_selected_mols']={'stage':vina['id'],'selector':'auto'}
        installation=root/'tools';installation.mkdir()
        dock6['parameters']['dock6_app_path']=str(installation)
        self.save([vina,dock6])
        records=[]
        for stage,entry,input_file in ((vina,first,dataset),(dock6,second,root/'vina/docking_results.csv')):
            files={'base_input_path':[str(receptor)],'base_selected_mols':[str(input_file)]}
            entry.update(id=stage['id'],name=stage['name'],configuration=copy.deepcopy(stage),
                         input_files=files,input_manifest={'selected:'+field:list(artifact_manifest(paths).items()) for field,paths in files.items()},
                         finished_at=time.time())
            records.append(entry)
        rid=uuid4().hex
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',
                (rid,self.project_id,self.store.user(self.token)['id'],'succeeded',json.dumps(records),time.time(),time.time(),None))
        return [vina,dock6],records,installation

    def test_completed_vina_to_dock6_pipeline_starts_no_worker_after_restart(self):
        from biomolexplorer.pipeline import PipelineService
        stages,records,installation=self.fixture()
        expected={p:Path(p).read_bytes() for item in records for p in item['artifacts']}
        with patch('biomolexplorer.pipeline.JobManager',side_effect=AssertionError('Docking restarted')):
            for attempt in range(2):
                service=PipelineService(self.store,dock6_path=installation)
                try:
                    with patch.object(service,'_resolve',side_effect=AssertionError('Inputs or preparation regenerated')):
                        run=self.wait(service.submit(self.token,self.project_id))
                finally:service.close()
        self.assertEqual(run['status'],'succeeded',run.get('error'))
        self.assertTrue(all(item.get('reused') for item in run['stages']))
        self.assertEqual([item['docking_failures'] for item in run['stages']],[1,1])
        self.assertEqual([item['docking_accounting'] for item in run['stages']],
                         [item['docking_accounting'] for item in records])
        self.assertEqual(expected,{path:Path(path).read_bytes() for path in expected})

    def test_new_inputs_or_changed_configuration_invalidate_completed_request(self):
        from biomolexplorer.published_results import verify
        stages,records,installation=self.fixture()
        results={entry['id']:entry['artifacts'] for entry in records}
        self.assertTrue(verify(self.store,self.project_id,stages[0],records[0],records,results))
        changed=copy.deepcopy(stages[0]);changed['parameters']['exhaustiveness']+=1
        self.assertFalse(verify(self.store,self.project_id,changed,records[0],records,results))
        path=Path(records[0]['input_files']['base_selected_mols'][0])
        original=path.read_bytes()
        with path.open('a') as stream:stream.write('D,CCC\n')
        self.assertFalse(verify(self.store,self.project_id,stages[0],records[0],records,results))
        path.write_bytes(original)
        self.assertTrue(verify(self.store,self.project_id,stages[0],records[0],records,results))
        receptor=Path(records[0]['input_files']['base_input_path'][0])
        receptor.write_text('changed receptor coordinates')
        self.assertFalse(verify(self.store,self.project_id,stages[0],records[0],records,results))

    def test_exact_stage_and_batch_cache_accept_the_same_complete_experiment(self):
        from biomolexplorer.pipeline import PipelineService
        stages,records,installation=self.fixture()
        for entry in records:entry['cache_key']='key-'+entry['id']
        with self.store.connect() as db:
            db.execute('UPDATE runs SET stages=? WHERE project_id=?',(json.dumps(records),self.project_id))
        service=object.__new__(PipelineService);service.store=self.store
        for stage in stages:
            key='key-'+stage['id']
            cached=service._cached_stage(self.project_id,self.store.user(self.token)['id'],stage,key,'software')
            self.assertIsNotNone(cached)
            self.assertEqual(cached['docking_failures'],1)
            batch=service._cached_batch(self.project_id,stage['id'],key,unit=stage)
            self.assertIsNotNone(batch)
            self.assertEqual(batch['excluded_records'],1)
