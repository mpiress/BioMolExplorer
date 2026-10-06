"""Selected datasets stay separate across workers unless merge is requested."""
import csv
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import test_pipeline_execution as execution
from biomolexplorer.artifact_choices import choices
from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService,validate_pipeline
from biomolexplorer.ui.file_selection import FileSelection


class InputProcessingTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    upload=execution.PipelineExecutionTests.upload
    imported=execution.PipelineExecutionTests.imported
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def confirm_all(self,service,run,mode):
        pending=next(s for s in run['stages'] if s['status']=='awaiting_input')
        form=FileSelection(run,pending)
        self.assertEqual(form.processing.value,pending['configuration'].get('input_processing','individual'))
        form.processing.value=mode
        for rows in form.rows.values():
            for check,_ in rows:check.value=True
        return self.wait(service.resume(self.token,run['id'],form.read()))

    def test_real_fingerprints_similarity_and_graphs_keep_each_dataset_separate(self):
        assets=[self.upload('alpha.csv','name,smiles\nA,CCO\nB,CCCO\n'),
                self.upload('beta.csv','name,smiles\nC,c1ccccc1\nD,c1ccncc1\n')]
        source=self.imported(assets[0]);source['parameters']['asset_ids']=assets
        fingerprint,similarity,graph=new_stage('fingerprints'),new_stage('similarity'),new_stage('graphs')
        fingerprint['bindings']['base_input_path']={'stage':source['id']}
        similarity['bindings']['base_input_path']={'stage':fingerprint['id']}
        graph['bindings']['similarity_path']={'stage':similarity['id']}
        self.save([source,fingerprint,similarity,graph])
        service=PipelineService(self.store,worker_python=sys.executable,cpu_workers=1)
        try:
            run=self.wait(service.submit(self.token,self.project_id))
            for _ in range(3):
                self.assertEqual(run['status'],'awaiting_input',run['error'])
                run=self.confirm_all(service,run,'individual')
            self.assertEqual(run['status'],'succeeded',run['error'])
            for item in run['stages'][1:]:self.assertEqual(len(item['batches']),2)
            fingerprints=[Path(p) for p in run['stages'][1]['artifacts'] if p.endswith('.csv')]
            self.assertEqual({p.name for p in fingerprints},{'morgan_alpha.csv','morgan_beta.csv'})
            datasets=[]
            for path in fingerprints:
                with path.open() as stream:datasets.append({r['molecule_chembl_id'] for r in csv.DictReader(stream)})
            self.assertEqual(datasets,[{'A','B'},{'C','D'}])
            graph_inputs=[b['input_files']['base_input_path'] for b in run['stages'][3]['batches']]
            self.assertEqual([len(paths) for paths in graph_inputs],[1,1])
            self.assertEqual({Path(paths[0]).name for paths in graph_inputs},{'morgan_alpha.csv','morgan_beta.csv'})
        finally:service.close()

    def test_real_fingerprint_merge_creates_one_combined_dataset(self):
        assets=[self.upload('alpha.csv','name,smiles\nA,CCO\n'),self.upload('beta.csv','name,smiles\nB,CCC\n')]
        source=self.imported(assets[0]);source['parameters']['asset_ids']=assets
        stage=new_stage('fingerprints');stage['bindings']['base_input_path']={'stage':source['id']}
        self.save([source,stage])
        service=PipelineService(self.store,worker_python=sys.executable,cpu_workers=1)
        try:
            run=self.confirm_all(service,self.wait(service.submit(self.token,self.project_id)),'merge')
            self.assertEqual(run['status'],'succeeded',run['error'])
            paths=[Path(p) for p in run['stages'][1]['artifacts'] if p.endswith('.csv')]
            self.assertEqual(len(paths),1)
            with paths[0].open() as stream:self.assertEqual({r['molecule_chembl_id'] for r in csv.DictReader(stream)},{'A','B'})
            self.assertEqual(run['stages'][1]['configuration']['input_processing'],'merge')
        finally:service.close()

    def test_duplicate_names_remain_selectable_cache_reuses_individual_batches_and_mode_change_recalculates(self):
        source,stage=new_stage('retrieve_compounds'),new_stage('fingerprints')
        stage['bindings']['base_input_path']={'stage':source['id']}
        self.save([source,stage])
        base,calls=self.immediate_manager()
        class Multiple(base):
            def submit(manager,operation,parameters):
                job=super().submit(operation,parameters)
                if operation=='retrieve_compounds':
                    paths=[]
                    for index in (1,2):
                        path=manager.config.workspace/str(index)/'compounds.csv';path.parent.mkdir()
                        path.write_text(f'molecule_chembl_id,canonical_smiles\nM{index},CCO\n');paths.append(str(path))
                    job['result']['artifacts']=paths
                return job
        with patch('biomolexplorer.pipeline.JobManager',Multiple):
            service=PipelineService(self.store)
            try:
                run=self.confirm_all(service,self.wait(service.submit(self.token,self.project_id)),'individual')
                self.assertEqual(run['status'],'succeeded',run['error'])
                self.assertEqual(len(calls),3)
                self.assertEqual(len(choices(run['stages'][1])),2)
                for _,params in calls[1:]:
                    with (Path(params['base_input_path'])/params['files'][0]).open() as stream:
                        self.assertEqual(len(list(csv.DictReader(stream))),1)
                again=self.wait(service.submit(self.token,self.project_id))
                self.assertEqual(again['status'],'succeeded',again['error'])
                self.assertEqual(len(calls),3)
                self.assertTrue(again['stages'][1]['reused'])
                self.assertEqual(len(again['stages'][1]['batches']),2)
                stage['input_processing']='merge';self.save([source,stage])
                merged=self.confirm_all(service,self.wait(service.submit(self.token,self.project_id)),'merge')
                self.assertEqual(merged['status'],'succeeded',merged['error'])
                self.assertEqual(len(calls),4)
                self.assertFalse(merged['stages'][1].get('reused',False))
            finally:service.close()

    def test_multiple_input_ports_create_file_combinations_and_invalid_modes_are_rejected(self):
        stage=new_stage('docking_vina')
        stage['bindings']={field:{'sources':[{'asset':self.upload(name+'.csv')} for name in names]}
            for field,names in (('base_input_path',('a','b')),('base_selected_mols',('c','d')))}
        service=PipelineService(self.store)
        try:
            variants=list(service._variants(self.project_id,stage))
            self.assertEqual(len(variants),4)
            for variant,_ in variants:
                self.assertTrue(all('asset' in g for g in variant['bindings'].values()))
        finally:service.close()
        stage['input_processing']='invalid'
        with self.assertRaisesRegex(ValueError,'processar individualmente'):validate_pipeline([stage])

    def test_graph_merge_combines_selected_similarity_files(self):
        assets=[self.upload('first.csv','source,target,value\nA,B,0.8\n','similarity'),
                self.upload('second.csv','source,target,value\nC,D,0.9\n','similarity')]
        stage=new_stage('graphs');stage['input_processing']='merge'
        stage['bindings']['similarity_path']={'sources':[{'asset':asset} for asset in assets]}
        service=PipelineService(self.store)
        try:
            params=service._resolve(self.project_id,self.store.user(self.token)['id'],stage,{},pipeline=[stage])
            self.assertEqual(len(params['graph_inputs']),1)
            with Path(params['graph_inputs'][0]['file']).open() as stream:
                self.assertEqual({(r['source'],r['target']) for r in csv.DictReader(stream)},{('A','B'),('C','D')})
        finally:service.close()

    def test_individual_structure_preparation_limits_configured_records_to_each_file(self):
        source,stage=new_stage('retrieve_structures'),new_stage('prepare_structures')
        root=self.store.project_dir(self.project_id)/'structures';root.mkdir()
        atom='ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n'
        paths=[]
        for code in ('1ABC','2ABC'):
            path=root/(code+'.pdb');path.write_text(atom);paths.append(str(path))
        stage['parameters']['pdb_codes']=[['1ABC','LIG',1,'A'],['2ABC','LIG',1,'A']]
        stage['bindings']['base_input_path']={'sources':[{'stage':source['id'],'selector':Path(p).name} for p in paths]}
        service=PipelineService(self.store)
        try:
            for (variant,_),code in zip(service._variants(self.project_id,stage),('1ABC','2ABC')):
                params=service._resolve(self.project_id,self.store.user(self.token)['id'],variant,{source['id']:paths})
                self.assertEqual(params['pdb_codes'],[[code,'LIG',1,'A']])
        finally:service.close()

    def test_redocking_individual_batches_preserve_curated_pairs_and_settings(self):
        source,stage=new_stage('retrieve_structures'),new_stage('redocking')
        root=self.store.project_dir(self.project_id)/'structures';root.mkdir()
        paths=[]
        for code in ('1ABC','2ABC','3ABC'):
            path=root/(code+'.pdb')
            path.write_text('ATOM      1  C   ALA A   2       0.000   0.000   0.000  1.00  0.00           C\nHETATM    2  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n')
            paths.append(str(path))
        stage['parameters']['pdb_codes']=[['1ABC','LIG',1,'A'],['2ABC','LIG',1,'A']]
        stage['parameters']['preparation_pairs']={f'{code}|LIG|1|A':{'cofactors':[]} for code in ('1ABC','2ABC')}
        stage['bindings']['base_input_path']={'sources':[{'stage':source['id'],'selector':Path(p).name} for p in paths]}
        service=PipelineService(self.store)
        try:
            variants=list(service._variants(self.project_id,stage))
            self.assertEqual(len(variants),2)
            for (variant,_),code in zip(variants,('1ABC','2ABC')):
                params=service._resolve(self.project_id,self.store.user(self.token)['id'],variant,{source['id']:paths})
                self.assertEqual(params['pdb_codes'],[[code,'LIG',1,'A']])
                self.assertEqual(set(params['preparation_pairs']),{f'{code}|LIG|1|A'})
        finally:service.close()

    def test_merge_prepared_outputs_combines_their_metadata_and_centers(self):
        source,stage=new_stage('prepare_structures'),new_stage('docking_vina')
        atom='ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n'
        files=[];refs=[]
        for code in ('1ABC','2ABC'):
            root=self.store.project_dir(self.project_id)/code/'MeuAlvo'
            prepared=root/'Prepared';prepared.mkdir(parents=True)
            (prepared/f'{code}_A.dockprep.pdbqt').write_text(atom)
            (root/'pdb_codes.csv').write_text(f'PDB_CODE,LIGAND,RESNUM,CHAIN\n{code},LIG,1,A\n')
            (prepared/'centers.csv').write_text(f'{code}_LIG_1_A\n0\n0\n0\n')
            files.extend(str(p) for p in root.rglob('*') if p.is_file())
            refs.append({'stage':source['id'],'selector':f'{code}_A.dockprep.pdbqt'})
        service=PipelineService(self.store)
        try:
            path,_,_=service._materialize_inputs(self.project_id,stage,'base_input_path',refs,{source['id']:files})
            with (path/'MeuAlvo'/'pdb_codes.csv').open() as stream:
                self.assertEqual({r['PDB_CODE'] for r in csv.DictReader(stream)},{'1ABC','2ABC'})
            with (path/'MeuAlvo'/'Prepared'/'centers.csv').open() as stream:
                self.assertEqual(next(csv.reader(stream)),['1ABC_LIG_1_A','2ABC_LIG_1_A'])
        finally:service.close()
