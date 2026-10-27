"""MD contracts, exact pair reuse, failure isolation, and real DOCK integration."""
from contextlib import ExitStack
import copy
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from biomolexplorer.catalog import new_stage,operation_fields
from biomolexplorer.flow import compatible
from biomolexplorer.md_config import DEFAULTS,validate_settings
from biomolexplorer.md_inputs import MDInputs,pair_id,validate_structure
from biomolexplorer.md_results import MDResults
from biomolexplorer.stage_cache import artifact_manifest,stage_key
from biomolexplorer.result_visibility import project_run
from wrappers.dynamics import perform_dynamics_dock6
import test_pipeline_execution as execution

PDB='ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n'
def mol2(name='LIG'):
    return '@<TRIPOS>MOLECULE\n'+name+'\n1 0 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n1 C1 0.0 1.0 2.0 C.3 1 LIG 0.0\n@<TRIPOS>BOND\n'


class MDContractTests(unittest.TestCase):
    def test_node_has_typed_ports_and_every_documented_setting(self):
        stage=new_stage('dynamics_dock6')
        self.assertTrue(set(DEFAULTS)<=stage['parameters'].keys())
        self.assertTrue(set(DEFAULTS)<={f['name'] for f in operation_fields('dynamics_dock6')})
        self.assertTrue(compatible(new_stage('docking_dock6'),stage,'base_selected_mols'))
        self.assertTrue(compatible(new_stage('fingerprints'),stage,'base_selected_mols'))
        self.assertTrue(compatible(new_stage('docking_vina'),stage,'base_selected_mols'))
        self.assertTrue(compatible(new_stage('prepare_structures'),stage,'base_input_path'))
    def test_invalid_settings_rejected_including_nonfinite_and_empty(self):
        for key,value in [('amber_score_md_steps',0),('amber_score_temperature',math.nan),
                ('amber_score_minimization_rmsgrad',''),('amber_score_gb_model',3),('md_time_step',.01),
                ('amber_score_md_steps',True),('md_movable_mask',"':1\n"),('selected_pairs',['invalid'])]:
            with self.subTest(key=key),self.assertRaises(ValueError):validate_settings({key:value})
        with self.assertRaises(ValueError):validate_settings(dict(amber_score_type='dock6_sander',amber_score_md_steps=1,trajectory_stride=10))
    def test_structure_types_presence_and_finite_coordinates(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);receptor=root/'rec.pdb';ligand=root/'lig.mol2'
            receptor.write_text(PDB);ligand.write_text(mol2())
            validate_structure(receptor,'receptor');validate_structure(ligand,'ligand')
            ligand.write_text(mol2().replace('0.0 1.0 2.0','nan 1.0 2.0'))
            with self.assertRaises(ValueError):validate_structure(ligand,'ligand')
            with self.assertRaises(ValueError):validate_structure(receptor,'ligand')
            with self.assertRaises(ValueError):validate_structure(root/'missing.pdb','receptor')


class MDPairCacheTests(unittest.TestCase):
    def setUp(self):
        self.stack=ExitStack();self.addCleanup(self.stack.close)
        self.root=Path(self.stack.enter_context(tempfile.TemporaryDirectory()))
        receptor=self.root/'rec.pdb';receptor.write_text(PDB)
        self.pairs=[]
        for code in ('A','B','C'):
            file=self.root/(code+'.mol2');file.write_text(mol2(code))
            r=dict(path=str(receptor),name=receptor.name,identity='rec')
            l=dict(path=str(file),name=file.name,identity=code)
            self.pairs.append(dict(id=pair_id(r,l),receptor=r,ligand=l))
        self.stack.enter_context(patch('biomolexplorer.md_engine.tools',return_value=(self.root,{})))
        self.stack.enter_context(patch('biomolexplorer.md_engine.tool_identity',return_value={'version':'test'}))
        self.runner=self.stack.enter_context(patch('biomolexplorer.md_engine.run_pair',side_effect=self.compute))
    def compute(self,pair,directory,*args):
        directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
        (directory/'energies.json').write_text('{}');(directory/'final.pdb').write_text(PDB)
        (directory/'topology.prmtop').write_text('topology')
        return dict(amber_score=-2.,complex_energy=-10.,receptor_contribution=5.,ligand_contribution=3.)
    def run_batch(self,name,pairs,seeds=None,**kwargs):
        output=self.root/name
        details=perform_dynamics_dock6(str(self.root),str(self.root),str(output),md_inputs=pairs,md_seed_path=seeds,**kwargs)
        return details,json.loads((output/'dynamics_results.json').read_text())
    def seed(self,name):
        root=self.root/('seed-'+name);shutil.copytree(self.root/name,root/'0');return str(root)
    def test_addition_removal_and_changed_science_only_run_needed_pairs(self):
        self.run_batch('first',self.pairs[:2]);self.assertEqual(self.runner.call_count,2)
        details,_=self.run_batch('add',self.pairs,self.seed('first'))
        self.assertEqual(details['md_reused'],2);self.assertEqual(self.runner.call_count,3)
        details,table=self.run_batch('remove',self.pairs[1:],self.seed('add'))
        self.assertEqual(details['md_reused'],2);self.assertEqual(len(table['pairs']),2)
        self.assertNotIn(self.pairs[0]['id'],[p['id'] for p in table['pairs']])
        self.run_batch('changed',self.pairs[1:],self.seed('remove'),amber_score_temperature=310.)
        self.assertEqual(self.runner.call_count,5)
    def test_failed_pair_is_red_state_and_only_failure_retried(self):
        def fail(pair,*args):
            if pair['ligand']['name']=='B.mol2':raise ValueError('topology problem')
            return self.compute(pair,*args)
        self.runner.side_effect=fail
        details,table=self.run_batch('failed',self.pairs[:2])
        self.assertEqual(details['md_failures'],1);self.assertEqual(table['pairs'][1]['status'],'failed')
        self.runner.reset_mock();self.runner.side_effect=self.compute
        details,_=self.run_batch('retry',self.pairs[:2],self.seed('failed'))
        self.assertEqual(details['md_reused'],1);self.runner.assert_called_once()
    def test_project_index_reuses_without_seed_and_respects_opt_out(self):
        with patch.dict(os.environ,BIOMOL_PREPARATION_CACHE=str(self.root/'cache'/'preparations')):
            self.run_batch('shared-a',self.pairs[:1])
            details,_=self.run_batch('shared-b',self.pairs[:2])
            self.assertEqual(details['md_reused'],1);self.assertEqual(self.runner.call_count,2)
            with patch.dict(os.environ,BIOMOL_REUSE_PREPARATIONS='0'):
                details,_=self.run_batch('shared-off',self.pairs[:2])
                self.assertEqual(details['md_reused'],0);self.assertEqual(self.runner.call_count,4)
    def test_tampered_cache_and_changed_pose_are_not_reused(self):
        self.run_batch('a',self.pairs[:1]);seed=self.seed('a')
        next(Path(seed).rglob('final.pdb')).write_text('tampered')
        details,_=self.run_batch('b',self.pairs[:1],seed);self.assertEqual(details['md_reused'],0)
        Path(self.pairs[0]['ligand']['path']).write_text(mol2('CHANGED'))
        details,_=self.run_batch('c',self.pairs[:1],self.seed('b'));self.assertEqual(details['md_reused'],0)


class MDInputAndResultsTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    upload=execution.PipelineExecutionTests.upload
    save=execution.PipelineExecutionTests.save
    def test_csv_preserves_same_compound_poses_for_distinct_receptors(self):
        from biomolexplorer.md_inputs import molecular_identity
        root=self.store.project_dir(self.project_id)/'md-input-test';root.mkdir()
        pose=root/'lig.mol2';pose.write_text(mol2())
        table=root/'docking_results.csv'
        table.write_text('molecule_chembl_id,receptor_id,conformer_file\nA,REC1,lig.mol2\nA,REC2,lig.mol2\n')
        stage=new_stage('dynamics_dock6');stage['parameters']['base_selected_mols']=str(root)
        inputs=MDInputs(self.store,self.project_id,{})
        rows=inputs.files(stage,'base_selected_mols')
        self.assertEqual(len(rows),2);self.assertNotEqual(rows[0]['identity'],rows[1]['identity'])
        original=molecular_identity(rows[0])
        table.write_text(table.read_text()+'B,REC1,lig.mol2\n')
        self.assertEqual(original,molecular_identity(inputs.files(stage,'base_selected_mols')[0]))
        pose.write_text(mol2('changed'))
        self.assertNotEqual(original,molecular_identity(rows[0]))
    def test_asset_pairs_exclusions_and_empty_selection(self):
        receptor=self.upload('rec.pdb',PDB,'receptor_pdb')
        ligands=[self.upload(code+'.mol2',mol2(code),'molecular_poses') for code in ('A','B')]
        stage=new_stage('dynamics_dock6');stage['bindings']={'base_input_path':{'asset':receptor},
            'base_selected_mols':{'sources':[{'asset':key} for key in ligands]}}
        inputs=MDInputs(self.store,self.project_id,{})
        with self.assertRaises(ValueError):inputs.resolve(stage)
        stage['parameters']['selected_pairs']=[pair_id(r,l) for r in inputs.files(stage,'base_input_path') for l in inputs.files(stage,'base_selected_mols')]
        pairs=inputs.resolve(stage)['md_inputs'];self.assertEqual(len(pairs),2)
        stage['parameters']['excluded_pairs']=[pairs[0]['id']]
        self.assertEqual([p['id'] for p in inputs.resolve(stage)['md_inputs']],[pairs[1]['id']])
        stage['parameters']['selected_pairs']=[]
        with self.assertRaises(ValueError):inputs.resolve(stage)
    def test_run_paths_do_not_change_scientific_cache_key(self):
        stage=new_stage('dynamics_dock6')
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);files=[]
            for name in ('one','two'):
                path=root/name;path.write_text('identical');files.append(path)
            def params(file):return dict(stage['parameters'],md_inputs=[dict(id='a'*64,
                receptor=dict(path=str(file),identity='rec'),ligand=dict(path=str(file),identity='lig'))])
            self.assertEqual(stage_key(stage,params(files[0]),{},'software'),stage_key(stage,params(files[1]),{},'software'))


class SelectionVisibilityTests(unittest.TestCase):
    def test_removed_batch_and_descendants_disappear_but_history_is_unchanged(self):
        source=new_stage('retrieve_compounds');stage=new_stage('fingerprints');child=new_stage('similarity')
        stage['bindings']={'base_input_path':{'sources':[{'stage':source['id'],'selector':'A.csv'},{'stage':source['id'],'selector':'B.csv'}]}}
        batches=[]
        for name in ('A','B'):
            config=copy.deepcopy(stage);config['bindings']['base_input_path']={'stage':source['id'],'selector':name+'.csv'}
            batches.append(dict(configuration=config,artifacts=['/output/'+name+'.csv'],input_files={'base_input_path':['/input/'+name+'.csv']}))
        run=dict(status='succeeded',stages=[dict(id=stage['id'],operation='fingerprints',configuration=stage,status='succeeded',
            artifacts=['/output/A.csv','/output/B.csv'],batches=batches),
            dict(id=child['id'],operation='similarity',configuration=child,status='succeeded',artifacts=['/child/B.csv'],input_files={'base_input_path':['/output/B.csv']})])
        active=copy.deepcopy(stage);active['bindings']['base_input_path']={'stage':source['id'],'selector':'A.csv'}
        projected=project_run(run,[source,active,child])
        self.assertEqual(projected['stages'][0]['artifacts'],['/output/A.csv'])
        self.assertEqual(projected['stages'][1]['artifacts'],[])
        self.assertEqual(len(run['stages'][0]['artifacts']),2)


@unittest.skipUnless(Path('/home/michel/progs/dock6/install/test/amber_score_1lgu/lig.mol2').is_file(),'DOCK6 example not installed')
class NativeMDTests(unittest.TestCase):
    def test_native_am1bcc_parameterization(self):
        from biomolexplorer.md_engine import tools,run_pair
        from biomolexplorer.md_analysis import topology
        settings=dict(DEFAULTS,amber_score_before_md_minimization_cycles=1,amber_score_md_steps=1,amber_score_after_md_minimization_cycles=1)
        root,found=tools('/home/michel/progs/dock6','',settings);example=root/'install/test/amber_score_1lgu'
        with tempfile.TemporaryDirectory(prefix='bme-md-am1bcc-') as folder:
            result=run_pair(dict(receptor={'path':str(example/'1lgu.pdb')},ligand={'path':str(example/'lig.mol2')}),Path(folder),root,found,settings)
            charges=topology(Path(folder)/'lig.1.prmtop')['charges']
            self.assertLess(abs(float(charges.sum())-round(float(charges.sum()))),.01)
            self.assertTrue(math.isfinite(result['amber_score']))

    @unittest.skipUnless(Path('/tmp/biomol-md-validation/bin/sander').is_file(),'isolated AmberTools not installed')
    def test_native_sander_and_separate_interface_analysis(self):
        from biomolexplorer.md_engine import tools,run_pair
        settings=dict(DEFAULTS,amber_score_type='dock6_sander',amber_score_before_md_minimization_cycles=1,
            amber_score_md_steps=20,trajectory_stride=5,amber_score_after_md_minimization_cycles=1,use_existing_ligand_charges=True)
        root,found=tools('/home/michel/progs/dock6','/tmp/biomol-md-validation',settings)
        example=root/'install/test/amber_score_1lgu'
        with tempfile.TemporaryDirectory(prefix='bme-md-sander-') as folder:
            result=run_pair(dict(receptor={'path':str(example/'1lgu.pdb')},ligand={'path':str(example/'lig.mol2')}),Path(folder),root,found,settings)
            data=json.loads((Path(folder)/'analysis.json').read_text())
            self.assertEqual(data['summary']['frames'],4)
            self.assertAlmostEqual(data['summary']['simulation_duration_ps'],.02)
            for key in ('protein_backbone_rmsd','ligand_rmsd','hydrogen_bonds','sasa'):
                self.assertTrue(all(math.isfinite(r[key]) for r in data['series']))
            self.assertTrue(all(0<=r['occupancy']<=1 for r in data['hbonds']))
            self.assertTrue(math.isfinite(result['sander_mean_potential_energy']))

    def test_native_amber_with_preserved_charges_generates_finite_energies(self):
        from biomolexplorer.md_engine import tools,run_pair
        settings=dict(DEFAULTS,amber_score_before_md_minimization_cycles=1,amber_score_md_steps=1,amber_score_after_md_minimization_cycles=1,use_existing_ligand_charges=True)
        root,found=tools('/home/michel/progs/dock6','',settings);example=root/'install/test/amber_score_1lgu'
        with tempfile.TemporaryDirectory(prefix='bme-md-native-') as folder:
            result=run_pair(dict(receptor={'path':str(example/'1lgu.pdb')},ligand={'path':str(example/'lig.mol2')}),Path(folder),root,found,settings)
            self.assertTrue(all(math.isfinite(v) for v in result.values()))
            self.assertTrue((Path(folder)/'energies.png').is_file())
            self.assertTrue((Path(folder)/'rec.lig.1.final_pose.amber.restart').is_file())
    def test_native_parameter_checkpoints_survive_temperature_change(self):
        from biomolexplorer.md_engine import tools,run_pair
        settings=dict(DEFAULTS,amber_score_before_md_minimization_cycles=1,amber_score_md_steps=1,
            amber_score_after_md_minimization_cycles=1,use_existing_ligand_charges=True)
        root,found=tools('/home/michel/progs/dock6','',settings);example=root/'install/test/amber_score_1lgu'
        pair=dict(receptor={'path':str(example/'1lgu.pdb')},ligand={'path':str(example/'lig.mol2')})
        with tempfile.TemporaryDirectory(prefix='bme-md-cache-') as folder:
            base=Path(folder)
            with patch.dict(os.environ,BIOMOL_PREPARATION_CACHE=str(base/'cache/preparations')):
                run_pair(pair,base/'first',root,found,settings)
                run_pair(pair,base/'second',root,found,dict(settings,amber_score_temperature=301.))
            for name in ('amberize_receptor.out','amberize_ligand.lig.1.out','amberize_complex.lig.1.out'):
                self.assertGreater((base/'first'/name).stat().st_size,0)
                self.assertEqual((base/'second'/name).stat().st_size,0)


if __name__=='__main__':unittest.main()
