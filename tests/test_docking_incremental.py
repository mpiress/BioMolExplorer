"""Scientific identity, partial recovery and indexed artifact isolation."""
import json
import os
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch

from biomolexplorer.calculation_cache import CalculationCache
from biomolexplorer.storage import write_json
from biomolexplorer.vina_batch import execute, run_resilient

ATOM='ATOM      1  C   UNL A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\n'
POSE='REMARK VINA RESULT: -7.0 0.0 0.0\nROOT\n'+ATOM+'ENDROOT\nTORSDOF 0\n'


class IndexedCacheTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.env=patch.dict(os.environ,{'BIOMOL_PREPARATION_CACHE':str(self.root/'shared/ligands'),
                                      'BIOMOL_REUSE_PREPARATIONS':'1'})
        self.env.start();self.addCleanup(self.env.stop)

    def test_content_addressing_and_independent_materializations(self):
        cache=CalculationCache('test',self.root)
        source=self.root/'source';source.write_text('same content')
        cache.remember('one',{'pose':source},{'compound':'A'})
        cache.remember('two',{'pose':source},{'compound':'B'})
        with sqlite3.connect(cache.database) as db:
            self.assertEqual(db.execute('SELECT count(*) FROM experiments').fetchone()[0],2)
            self.assertEqual(db.execute('SELECT count(*) FROM objects').fetchone()[0],1)
        destination=self.root/'restored'
        self.assertIsNotNone(cache.restore('one',{'pose':destination}))
        destination.write_text('user edit')
        self.assertIsNotNone(cache.restore('two',{'pose':destination}))
        self.assertEqual(destination.read_text(),'same content')

    def test_corrupted_object_and_missing_artifact_are_not_hits(self):
        cache=CalculationCache('test',self.root)
        source=self.root/'source';source.write_text('complete')
        cache.remember('one',{'pose':source})
        with cache.connect() as db:
            digest=db.execute('SELECT sha256 FROM objects').fetchone()[0]
        cache.object_path(digest).write_text('partial')
        self.assertIsNone(cache.restore('one',{'pose':self.root/'out'}))
        self.assertFalse((self.root/'out').exists())
        self.assertIsNone(cache.restore('missing',{'pose':self.root/'out'}))

    def test_automatic_docking_binding_does_not_restore_an_outdated_file_selection(self):
        from biomolexplorer.catalog import new_stage
        from biomolexplorer.pipeline import PipelineService
        import copy
        for operation in ('docking_vina','docking_dock6','redocking'):
            stage=new_stage(operation)
            stage['bindings']={'base_input_path':{'stage':'source','selector':'auto'}}
            previous=copy.deepcopy(stage)
            previous['bindings']['base_input_path']['selector']='old_compounds.csv'
            restored=PipelineService._reuse_configuration(stage,previous)
            self.assertEqual(restored['bindings'],stage['bindings'])


class VinaIncrementalTests(IndexedCacheTests):
    def setUp(self):
        super().setUp()
        (self.root/'rec.pdbqt').write_text(ATOM)
        self.calls=[];self.fail=set();self.interrupt=set()
        def dock(session,configuration,directory):
            options=dict(line.split(' = ',1) for line in configuration.read_text().splitlines())
            compound=Path(options['ligand']).stem;self.calls.append(compound)
            if compound in self.interrupt:raise KeyboardInterrupt()
            if compound in self.fail:raise RuntimeError('controlled compound failure')
            (directory/options['out']).write_text('REMARK '+compound+'\n'+POSE)
            return True
        self.runner=patch('biomolexplorer.vina_session.VinaSessions.run',dock)
        self.runner.start();self.addCleanup(self.runner.stop)

    def config(self,code,exhaustiveness=8):
        ligand=self.root/(code+'.pdbqt')
        if not ligand.exists():ligand.write_text('REMARK '+code+'\nROOT\n'+ATOM+'ENDROOT\nTORSDOF 0\n')
        config=self.root/(code+'.vina')
        config.write_text(f'receptor = rec.pdbqt\nligand = {code}.pdbqt\nout = {code}.out.pdbqt\n'
                          f'exhaustiveness = {exhaustiveness}\nseed = 73\n')
        return config

    def batch(self,codes,exhaustiveness=8):
        configurations=[self.config(code,exhaustiveness).name for code in codes]
        request=self.root/'request.json';state=self.root/'state.json'
        write_json(dict(directory=str(self.root),configurations=configurations),request)
        return execute(request,state)

    def test_unchanged_batch_and_new_compound_only_run_missing_work(self):
        self.batch(['A','B']);self.assertEqual(self.calls,['A','B'])
        result=self.batch(['A','B']);self.assertEqual(self.calls,['A','B'])
        self.assertTrue(all(row['reused'] for row in result['completed']))
        result=self.batch(['A','B','C']);self.assertEqual(self.calls,['A','B','C'])
        self.assertEqual(len(result['completed']),3)
        with patch('biomolexplorer.vina_batch.subprocess.run') as child:
            result=run_resilient([self.root/(code+'.vina') for code in ('A','B','C')],self.root)
        child.assert_not_called();self.assertEqual(len(result['completed']),3)

    def test_changed_ligand_affects_only_that_compound_and_parameters_invalidate_docking(self):
        self.batch(['A','B'])
        with (self.root/'A.pdbqt').open('a') as stream:stream.write('REMARK new conformer\n')
        self.batch(['A','B']);self.assertEqual(self.calls,['A','B','A'])
        self.batch(['A','B'],exhaustiveness=20);self.assertEqual(self.calls,['A','B','A','A','B'])
        with (self.root/'rec.pdbqt').open('a') as stream:stream.write('REMARK changed receptor\n')
        self.batch(['A','B'],exhaustiveness=20);self.assertEqual(self.calls[-2:],['A','B'])

    def test_partial_failure_and_interruption_preserve_completed_checkpoint(self):
        self.fail={'B'};result=self.batch(['A','B'])
        self.assertEqual(len(result['completed']),1);self.assertEqual(len(result['failed']),1)
        self.fail.clear();self.batch(['A','B']);self.assertEqual(self.calls,['A','B','B'])
        # A new coordinator and output directory still find committed objects.
        for code in ('A','B'):(self.root/(code+'.out.pdbqt')).unlink()
        self.batch(['A','B']);self.assertEqual(self.calls,['A','B','B'])

    def test_forced_rerun_recalculates_and_missing_ligand_is_not_reused(self):
        self.batch(['A','B'])
        with patch.dict(os.environ,{'BIOMOL_REUSE_PREPARATIONS':'0'}):self.batch(['A','B'])
        self.assertEqual(self.calls,['A','B','A','B'])
        (self.root/'A.pdbqt').unlink()
        write_json(dict(directory=str(self.root),configurations=['A.vina']),self.root/'request.json')
        result=execute(self.root/'request.json',self.root/'state.json')
        self.assertEqual(len(result['failed']),1)
        self.assertIn('ausente',result['failed'][0]['error'])

    def test_interruption_commits_completed_items_and_concurrent_requests_share_checkpoint(self):
        from concurrent.futures import ThreadPoolExecutor
        self.interrupt={'B'}
        with self.assertRaises(KeyboardInterrupt):self.batch(['A','B'])
        self.interrupt.clear();self.batch(['A','B'])
        self.assertEqual(self.calls,['A','B','B'])
        config=self.config('C')
        requests=[]
        for number in range(2):
            request=self.root/('concurrent-'+str(number)+'.json')
            state=self.root/('state-'+str(number)+'.json')
            write_json(dict(directory=str(self.root),configurations=[config.name]),request)
            requests.append((request,state))
        with ThreadPoolExecutor(max_workers=2) as pool:
            results=list(pool.map(lambda args:execute(*args),requests))
        self.assertEqual(self.calls.count('C'),1)
        self.assertTrue(all(len(result['completed'])==1 for result in results))


if __name__=='__main__':unittest.main()
