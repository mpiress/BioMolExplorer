"""Incremental runs reuse real prepared files and poses, never unfinished data."""
import csv
import copy
import json
from pathlib import Path
import shutil
import sqlite3
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
from biomolexplorer.operations import execute_operation
from biomolexplorer.redocking_cache import NAME, stage_sources
from biomolexplorer.stage_cache import artifact_manifest
from caad.docking import Docking
from test_redocking_partial_failures import ATOM, LIG, BAD


class Store:
    def __init__(self, root):
        self.root=root
        self.db=sqlite3.connect(':memory:');self.db.row_factory=sqlite3.Row
        self.db.execute('CREATE TABLE runs (project_id TEXT, status TEXT, created INT, stages TEXT)')
    def connect(self):return self.db
    def scoped_path(self,pid,path):
        result=Path(path).resolve();result.relative_to(self.root.resolve());return result


class IncrementalRedockingTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.store=Store(self.root);self.addCleanup(self.store.db.close)
        self.source=self.root/'input'/'Target';self.source.mkdir(parents=True)
        for code in ('1ABC','2ABC'):(self.source/(code+'.pdb')).write_text(ATOM+LIG+BAD)
        self.scripts=[];self.docked=[];self.index=0
        def chimera(dock,filename):
            self.scripts.append(filename);directory=Path(dock.outputpath)
            if filename.startswith('prepare_complex_'):
                (directory/(filename[16:-4]+'.complex.pdb')).write_text(ATOM+LIG+BAD)
            elif filename.startswith('prepare_receptor_'):
                (directory/(filename[17:-4]+'.dockprep.mol2')).write_text('mol2')
            else:
                prefix=filename[15:-4]
                (directory/(prefix+'.lig.mol2')).write_text('mol2')
                (directory/(prefix+'.lig.pdb')).write_text(LIG)
        def obabel(dock,inputfile,outputfile,*args,**kwargs):(Path(dock.outputpath)/outputfile).write_text(ATOM)
        def vina(factory):
            args=factory.call_args.kwargs;frame=args['pdb_codes'].copy();frame['RMSD']=.25
            for row in frame.itertuples():
                self.docked.append((row.PDB_CODE,row.LIGAND))
                output=Path(args['output_path']);output.mkdir(parents=True,exist_ok=True)
                (output/f'{row.PDB_CODE}_{row.LIGAND}_{row.RESNUM}{row.CHAIN}.lig.pdbqt').write_text(ATOM)
            frame.to_csv(Path(args['complex_input_path'])/'pdb_codes.csv',index=False)
        self.patches=[patch.object(Docking,'prepare_on_chimera',chimera),patch.object(Docking,'prepare_on_obabel',obabel),
            patch.object(Docking,'calculate_ligand_centerofmass',return_value=[0,0,0]),
            patch('biomolexplorer.redocking_config.validate_redocking_tools'),patch('wrappers.redocking.DockVina')]
        started=[p.start() for p in self.patches]
        for p in self.patches:self.addCleanup(p.stop)
        started[-1].return_value.redocking.side_effect=lambda **kw:vina(started[-1])

    def run_pairs(self, records, status='succeeded', reuse=True, only_selected_inputs=False, **params):
        self.index+=1
        input_root=self.root/f'inputs-{self.index}';shutil.copytree(self.source,input_root/'Target')
        if only_selected_inputs:
            for path in (input_root/'Target').glob('*.pdb'):
                if path.stem not in {r[0] for r in records}:path.unlink()
        if reuse:stage_sources(self.store,'project','redock',input_root/'Target')
        configuration=dict(operation='redocking',parameters=dict(base_input_path=str(input_root),target='Target',pdb_codes=records,**params))
        output=self.root/f'run-{self.index}'
        result=execute_operation('redocking',configuration['parameters'],output)
        stage=dict(id='redock',status='succeeded',artifacts=result.artifacts,
                   artifact_manifest=artifact_manifest(result.artifacts),configuration=configuration)
        self.store.db.execute('INSERT INTO runs VALUES (?,?,?,?)',('project',status,self.index,json.dumps([stage])))
        with (output/'structures/Target/pdb_codes.csv').open() as f:rows=list(csv.DictReader(f))
        self.last_scripts=self.scripts[:];self.scripts.clear();docked=self.docked[:];self.docked.clear()
        return output,stage,rows,docked

    def test_new_ligand_reuses_receptor_and_accumulates_without_repeating_old_pair(self):
        first,_,_,_=self.run_pairs([['1ABC','LIG',1,'A']])
        self.assertEqual(len(json.loads((first/NAME).read_text())['entries']),1)
        with patch('biomolexplorer.redocking_preparation.prepare_pairs', wraps=__import__('biomolexplorer.redocking_preparation',fromlist=['prepare_pairs']).prepare_pairs) as prepare:
            output,_,rows,docked=self.run_pairs([['1ABC','BAD',3,'A']])
        self.assertEqual(prepare.call_args.kwargs['ready'],{('1ABC','A')})
        self.assertEqual(docked,[('1ABC','BAD')])
        self.assertEqual(self.last_scripts,['prepare_ligand_1ABC_BAD_3A.com'])
        self.assertEqual({r['LIGAND'] for r in rows},{'LIG','BAD'})
        self.assertEqual({r['RMSD'] for r in rows},{'0.25'})
        self.assertTrue((output/'Target/1ABC_LIG_1A.lig.pdbqt').exists())
        _,_,rows,docked=self.run_pairs([['1ABC','LIG',1,'A'],['1ABC','BAD',3,'A']])
        self.assertEqual(docked,[]);self.assertEqual(len(rows),2)

    def test_native_preparation_provenance_survives_pair_reuse(self):
        prepare=Docking.prepare_on_chimera
        metadata={'net_charge':1,'charge_method':'am1','validated_graph_preserved':True}
        def annotated(dock,filename):
            prepare(dock,filename)
            if filename.startswith('prepare_ligand_'):
                name=filename[15:-4]+'.lig.preparation.json'
                (Path(dock.outputpath)/name).write_text(json.dumps(metadata))
        with patch.object(Docking,'prepare_on_chimera',annotated):
            first,_,_,_=self.run_pairs([['1ABC','LIG',1,'A']])
        entry=json.loads((first/NAME).read_text())['entries'][0]
        self.assertTrue(any(name.endswith('.lig.preparation.json') for name in entry['files']))
        output,_,_,docked=self.run_pairs([['1ABC','LIG',1,'A']])
        self.assertEqual(docked,[])
        self.assertEqual(json.loads((output/'structures/Target/Prepared/1ABC_LIG_1A.lig.preparation.json').read_text()),metadata)

    def test_new_receptor_prepares_and_docks_only_new_pair(self):
        self.run_pairs([['1ABC','LIG',1,'A']])
        _,_,rows,docked=self.run_pairs([['1ABC','LIG',1,'A'],['2ABC','LIG',1,'A']])
        self.assertEqual(docked,[('2ABC','LIG')]);self.assertEqual(len(rows),2)
        self.assertEqual(self.last_scripts,['prepare_complex_2ABC_A.com','prepare_receptor_2ABC_A.com','prepare_ligand_2ABC_LIG_1A.com'])

    def test_new_receptor_selection_keeps_previous_receptors_in_final_table(self):
        self.run_pairs([['1ABC','LIG',1,'A']])
        output,_,rows,docked=self.run_pairs([['2ABC','LIG',1,'A']],only_selected_inputs=True)
        self.assertEqual(docked,[('2ABC','LIG')])
        self.assertEqual({r['PDB_CODE'] for r in rows},{'1ABC','2ABC'})
        self.assertTrue((output/'structures/Target/1ABC.pdb').exists())
        self.assertTrue((output/'Target/1ABC_LIG_1A.lig.pdbqt').exists())

    def test_removed_original_is_not_resurrected_by_accumulation(self):
        self.run_pairs([['1ABC','LIG',1,'A']])
        (self.root/'inputs-1/Target/1ABC.pdb').unlink()
        _,_,rows,docked=self.run_pairs([['2ABC','LIG',1,'A']],only_selected_inputs=True)
        self.assertEqual(docked,[('2ABC','LIG')])
        self.assertEqual([r['PDB_CODE'] for r in rows],['2ABC'])

    def test_individual_batches_partition_previous_results_without_duplicates(self):
        self.run_pairs([['1ABC','LIG',1,'A'],['2ABC','LIG',1,'A']])
        first=self.root/'first'/'Target';second=self.root/'second'/'Target'
        self.assertEqual(stage_sources(self.store,'project','redock',first,exclude_codes={'2ABC'}),1)
        self.assertEqual(stage_sources(self.store,'project','redock',second,only_codes={'2ABC'}),1)
        first_codes={e['record'][0] for p in first.rglob(NAME) for e in json.loads(p.read_text())['entries']}
        second_codes={e['record'][0] for p in second.rglob(NAME) for e in json.loads(p.read_text())['entries']}
        self.assertEqual(first_codes,{'1ABC'});self.assertEqual(second_codes,{'2ABC'})

    def test_changing_input_source_does_not_mix_previous_collection(self):
        _,stage,_,_=self.run_pairs([['1ABC','LIG',1,'A']])
        configuration=copy.deepcopy(stage['configuration'])
        configuration['bindings']={'base_input_path':{'stage':'different_collection'}}
        self.assertEqual(stage_sources(self.store,'project','redock',self.root/'other'/'Target',
                                      configuration=configuration),0)

    def test_cancelled_and_failed_runs_are_never_reused(self):
        for status in ('cancelled','failed','running','awaiting_input'):
            with self.subTest(status=status):
                self.store.db.execute('DELETE FROM runs')
                self.run_pairs([['1ABC','LIG',1,'A']],status=status)
                _,_,_,docked=self.run_pairs([['1ABC','LIG',1,'A']])
                self.assertEqual(docked,[('1ABC','LIG')])

    def test_truncated_pose_invalidates_reuse(self):
        output,_,_,_=self.run_pairs([['1ABC','LIG',1,'A']])
        (output/'Target/1ABC_LIG_1A.lig.pdbqt').write_text('partial')
        _,_,_,docked=self.run_pairs([['1ABC','LIG',1,'A']])
        self.assertEqual(docked,[('1ABC','LIG')])

    def test_changed_docking_options_reuse_receptor_but_recalculate_pose(self):
        self.run_pairs([['1ABC','LIG',1,'A']])
        with patch('biomolexplorer.redocking_preparation.prepare_pairs', wraps=__import__('biomolexplorer.redocking_preparation',fromlist=['prepare_pairs']).prepare_pairs) as prepare:
            _,_,_,docked=self.run_pairs([['1ABC','LIG',1,'A']],exhaustiveness=40)
        self.assertEqual(prepare.call_args.kwargs['ready'],{('1ABC','A')})
        self.assertEqual(docked,[('1ABC','LIG')])

    def test_changed_structure_or_receptor_settings_require_preparation(self):
        self.run_pairs([['1ABC','LIG',1,'A']])
        with (self.source/'1ABC.pdb').open('a') as f:f.write('REMARK changed source\n')
        with patch('biomolexplorer.redocking_preparation.prepare_pairs', wraps=__import__('biomolexplorer.redocking_preparation',fromlist=['prepare_pairs']).prepare_pairs) as prepare:
            _,_,_,docked=self.run_pairs([['1ABC','LIG',1,'A']])
        self.assertEqual(prepare.call_args.kwargs['ready'],set())
        self.assertEqual(docked,[('1ABC','LIG')])
        with patch('biomolexplorer.redocking_preparation.prepare_pairs', wraps=__import__('biomolexplorer.redocking_preparation',fromlist=['prepare_pairs']).prepare_pairs) as prepare:
            _,_,_,docked=self.run_pairs([['1ABC','LIG',1,'A']],
                preparation_pairs={'1ABC|LIG|1|A':{'receptor':{'minimize':False}}})
        self.assertEqual(prepare.call_args.kwargs['ready'],set())
        self.assertEqual(docked,[('1ABC','LIG')])

    def test_forced_rerun_ignores_cache(self):
        self.run_pairs([['1ABC','LIG',1,'A']])
        _,_,_,docked=self.run_pairs([['1ABC','LIG',1,'A']],reuse=False)
        self.assertEqual(docked,[('1ABC','LIG')])

    def test_completed_legacy_native_results_are_adopted(self):
        output,stage,_,_=self.run_pairs([['1ABC','LIG',1,'A']])
        stage['artifacts'].remove(str(output/NAME));stage['artifact_manifest'].pop(str(output/NAME))
        (output/NAME).unlink()
        self.store.db.execute('UPDATE runs SET stages=?',(json.dumps([stage]),))
        _,_,rows,docked=self.run_pairs([['1ABC','BAD',3,'A']])
        self.assertEqual(docked,[('1ABC','BAD')]);self.assertEqual(len(rows),2)


class SchedulerIncrementalTests(unittest.TestCase):
    """Exercise the real store and scheduler, with only scientific engines mocked."""
    def test_scheduler_passes_verified_seeds_and_results_table_accumulates(self):
        import test_pipeline_execution as execution
        from biomolexplorer.catalog import new_stage
        from biomolexplorer.pipeline import PipelineService
        from biomolexplorer.redocking_results import RedockingResults
        fixture=execution.PipelineExecutionTests()
        fixture.setUp();self.addCleanup(fixture.tearDown)
        science=IncrementalRedockingTests();science.setUp();self.addCleanup(science.doCleanups)
        source=fixture.store.project_dir(fixture.project_id)/'input'/'Target'
        shutil.copytree(science.source,source)
        stage=new_stage('redocking')
        stage['parameters'].update(base_input_path=str(source.parent),target='Target',pdb_codes=[['1ABC','LIG',1,'A']])
        fixture.save([stage]);submitted=[]
        class Manager:
            def __init__(self,config):self.config=config
            def submit(self,operation,parameters):
                submitted.append(copy.deepcopy(parameters))
                result=execute_operation(operation,parameters,self.config.workspace/'artifacts')
                self.job=dict(id='worker',status='succeeded',log_path=str(self.config.workspace/'worker.log'),
                    error=None,progress=None,result=result.to_dict())
                return self.job
            def get(self,job):return self.job
            def close(self):pass
            def cancel(self,job):self.job['status']='cancelled'
        with patch('biomolexplorer.pipeline.JobManager',Manager):
            service=PipelineService(fixture.store);self.addCleanup(service.close)
            first=fixture.finish_selected(service,service.submit(fixture.token,fixture.project_id))
            self.assertEqual(first['status'],'succeeded',first.get('error'))
            stage['parameters']['pdb_codes']=[['1ABC','BAD',3,'A']];fixture.save([stage])
            second=fixture.finish_selected(service,service.submit(fixture.token,fixture.project_id))
        self.assertEqual(second['status'],'succeeded',second.get('error'))
        self.assertIn('redocking-inputs',submitted[-1]['base_input_path'])
        self.assertFalse(Path(submitted[-1]['base_input_path']).exists())
        self.assertEqual(science.docked,[('1ABC','LIG'),('1ABC','BAD')])
        rows=RedockingResults(fixture.store).simulations(fixture.token,fixture.project_id,second['id'],stage['id'])
        self.assertEqual({r['ligand'] for r in rows},{'LIG','BAD'})
        self.assertEqual([r['rmsd'] for r in rows],[.25,.25])
        # A whole-stage match must also reject a run that did not commit.
        with fixture.store.connect() as db:
            db.execute("UPDATE runs SET status='interrupted' WHERE id=?",(second['id'],))
        with patch('biomolexplorer.pipeline.JobManager',Manager):
            third=fixture.finish_selected(service,service.submit(fixture.token,fixture.project_id))
        self.assertEqual(third['status'],'succeeded',third.get('error'))
        self.assertEqual(science.docked,[('1ABC','LIG'),('1ABC','BAD'),('1ABC','BAD')])
        self.assertFalse(third['stages'][0].get('reused'))
        # Returning to the first selection must keep later accumulated results,
        # rather than returning the first run's smaller, whole-stage table.
        stage['parameters']['pdb_codes']=[['1ABC','LIG',1,'A']];fixture.save([stage])
        with patch('biomolexplorer.pipeline.JobManager',Manager):
            fourth=fixture.finish_selected(service,service.submit(fixture.token,fixture.project_id))
        self.assertEqual(fourth['status'],'succeeded',fourth.get('error'))
        self.assertEqual(science.docked,[('1ABC','LIG'),('1ABC','BAD'),('1ABC','BAD')])
        rows=RedockingResults(fixture.store).simulations(fixture.token,fixture.project_id,fourth['id'],stage['id'])
        self.assertEqual({r['ligand'] for r in rows},{'LIG','BAD'})


if __name__=='__main__':unittest.main()
