"""Execution contracts for independent branches, selected blocks and recovery."""
import copy
import json
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest.mock import patch
from uuid import uuid4

from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService, validate_pipeline, execution_modes
from biomolexplorer.workspace import AccessDenied, WorkspaceStore


class PipelineExecutionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='biomol pipeline execution ')
        self.store = WorkspaceStore(self.temporary.name)
        self.token = self.store.register('Owner','execution@example.org','execution-password')
        self.project = self.store.create_project(self.token,'Execution review')
        self.project_id = self.project['id']

    def tearDown(self):
        self.temporary.cleanup()

    def upload(self, name='compounds.csv', content='name,smiles\nETHANOL,CCO\n', kind='compounds'):
        ticket = self.store.prepare_upload(self.token,self.project_id,name,kind)
        (self.store.staging / ticket).write_text(content)
        return self.store.finish_upload(self.token,ticket)

    def imported(self, asset):
        stage = new_stage('import_results')
        stage['parameters']['asset_ids'] = [asset]
        return stage

    def save(self, stages):
        revision = self.store.project(self.token,self.project_id)['revision']
        self.store.save_pipeline(self.token,self.project_id,stages,revision)

    def wait(self, run):
        deadline = time.monotonic() + 40
        while time.monotonic() < deadline:
            run = self.store.get_run(self.token,run['id'])
            if run['status'] not in ('queued','running'):
                return run
            time.sleep(.025)
        self.fail('The pipeline did not finish in 40 seconds')

    def finish_selected(self, service, run):
        """Explicitly confirm fixture outputs at every stage boundary."""
        from biomolexplorer.bindings import sources
        run=self.wait(run)
        while run['status']=='awaiting_input':
            pending=next(s for s in run['stages'] if s['status']=='awaiting_input')
            configuration=copy.deepcopy(pending['configuration'])
            completed={s['id']:s for s in run['stages'] if s['status']=='succeeded'}
            for group in configuration.get('bindings',{}).values():
                for ref in sources(group):
                    if 'stage' in ref and ref.get('selector','auto')=='auto':
                        ref['selector']=Path(completed[ref['stage']]['artifacts'][0]).name
            run=self.wait(service.resume(self.token,run['id'],configuration))
        return run

    def immediate_manager(self, fail_operations=()):
        """Provider-independent supervisor boundary, with real pipeline/store calls."""
        calls = []
        class ImmediateManager:
            def __init__(self, config):
                self.config = config
                self.closed = False

            def submit(self, operation, parameters):
                calls.append((operation,dict(parameters)))
                artifact = self.config.workspace / 'compounds.csv'
                artifact.write_text('molecule_chembl_id,canonical_smiles\nCHEMBL1,CCO\n')
                failed = operation in fail_operations
                self.job = {'id':uuid4().hex,'status':'failed' if failed else 'succeeded',
                    'log_path':str(self.config.workspace / 'execution.log'), 'progress':None,
                    'error':'HTTP 500 do provedor' if failed else None,
                    'result':{'artifacts':[str(artifact)]}}
                return self.job

            def get(self, job_id):
                return self.job

            def cancel(self, job_id):
                self.job['status'] = 'cancelled'

            def close(self):
                self.closed = True
        return ImmediateManager,calls

    def test_single_default_chembl_block_reaches_supervisor_with_correct_configuration(self):
        stage = new_stage('retrieve_compounds')
        self.save([stage])
        manager,calls = self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service = PipelineService(self.store,worker_python=sys.executable)
            try:
                run = self.finish_selected(service,service.submit(self.token,self.project_id))
            finally:
                service.close()
        self.assertEqual(run['status'],'succeeded',run['error'])
        self.assertEqual(len(calls),1)
        self.assertEqual(calls[0],('retrieve_compounds',stage['parameters']))
        self.assertEqual(calls[0][1]['search_term'],'CHEMBL220')
        self.assertTrue(calls[0][1]['include_pubchem'])
        self.assertEqual(run['stages'][0]['status'],'succeeded')

    def test_independent_default_blocks_all_execute(self):
        stages = [new_stage('retrieve_compounds') for _ in range(3)]
        self.save(stages)
        manager,calls = self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service = PipelineService(self.store)
            try:
                run = self.wait(service.submit(self.token,self.project_id))
            finally:
                service.close()
        self.assertEqual(run['status'],'succeeded',run['error'])
        self.assertEqual(len(calls),3)
        self.assertEqual([s['status'] for s in run['stages']],['succeeded']*3)
        self.assertEqual(len({s['artifacts'][0] for s in run['stages']}),3)

    def test_provider_failure_only_blocks_its_descendants(self):
        source = new_stage('retrieve_compounds')
        dependent = new_stage('admet')
        dependent['bindings']['base_input_path'] = {'stage':source['id']}
        descendant = new_stage('fingerprints')
        descendant['bindings']['base_input_path'] = {'stage':dependent['id']}
        imported = self.imported(self.upload())
        independent = new_stage('admet')
        independent['bindings']['base_input_path'] = {'stage':imported['id']}
        self.save([source,dependent,descendant,imported,independent])
        manager,calls = self.immediate_manager({'retrieve_compounds'})
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service = PipelineService(self.store)
            try:
                run = self.finish_selected(service,service.submit(self.token,self.project_id))
            finally:
                service.close()
        states = {s['id']:s for s in run['stages']}
        self.assertEqual(run['status'],'failed')
        self.assertEqual(states[source['id']]['status'],'failed')
        self.assertEqual(states[dependent['id']]['status'],'skipped')
        self.assertIn(source['name'],states[dependent['id']]['error'])
        self.assertEqual(states[descendant['id']]['status'],'skipped')
        self.assertEqual(states[independent['id']]['status'],'succeeded')
        self.assertEqual([operation for operation,_ in calls],['retrieve_compounds','admet'])
        self.assertIn('HTTP 500',run['error'])

    def test_invalid_uploaded_input_is_rejected_before_any_experiment(self):
        with self.assertRaisesRegex(ValueError,'Padrão esperado'):
            self.upload('invalid.csv','name,value\nINVALID,10\n')
        self.assertEqual(self.store.assets(self.token,self.project_id),[])
        self.assertEqual(self.store.list_runs(self.token,self.project_id),[])

    def test_selected_block_runs_ancestors_and_excludes_unrelated_blocks(self):
        source = self.imported(self.upload())
        selected = new_stage('admet')
        selected['bindings']['base_input_path'] = {'stage':source['id']}
        unrelated = new_stage('retrieve_compounds')
        self.save([selected,unrelated,source])
        manager,calls = self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service = PipelineService(self.store)
            try:
                run = self.finish_selected(service,service.submit(self.token,self.project_id,[selected['id']]))
            finally:
                service.close()
        self.assertEqual(run['status'],'succeeded',run['error'])
        self.assertEqual([s['id'] for s in run['stages']],[source['id'],selected['id']])
        self.assertEqual([operation for operation,_ in calls],['admet'])
        self.assertEqual(calls[0][1]['input_file'],'selected_compounds.csv')

    def test_cancelling_one_stage_stops_the_entire_run_including_independent_blocks(self):
        stages = [new_stage('retrieve_compounds') for _ in range(2)]
        self.save(stages)
        manager,calls = self.immediate_manager()
        entered = threading.Event()
        class RunningManager(manager):
            def submit(self, operation, parameters):
                job = super().submit(operation,parameters)
                job['status'] = 'running'
                return job

            def get(self, job_id):
                entered.set()
                return super().get(job_id)
        with patch('biomolexplorer.pipeline.JobManager',RunningManager):
            service = PipelineService(self.store)
            try:
                submitted = service.submit(self.token,self.project_id)
                self.assertTrue(entered.wait(5))
                service.cancel(self.token,submitted['id'])
                run = self.wait(submitted)
            finally:
                service.close()
        self.assertEqual(run['status'],'cancelled')
        self.assertEqual([s['status'] for s in run['stages']],['cancelled','skipped'])
        self.assertEqual(len(calls),1)
        self.assertIn('cancelada',run['stages'][1]['error'])

    def test_startup_recovery_normalizes_all_unfinished_stage_statuses(self):
        stages = []
        for status in ('succeeded','running','queued'):
            stage = new_stage('retrieve_compounds')
            stages.append({'id':stage['id'],'name':stage['name'],'operation':stage['operation'],
                           'status':status,'configuration':stage,'artifacts':[]})
        run_id = uuid4().hex
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,NULL)',(run_id,self.project_id,
                self.store.user(self.token)['id'],'running',json.dumps(stages),time.time(),time.time()))
        service = PipelineService(self.store)
        try:
            run = self.store.get_run(self.token,run_id)
        finally:
            service.close()
        self.assertEqual(run['status'],'interrupted')
        self.assertEqual([s['status'] for s in run['stages']],['succeeded','interrupted','skipped'])
        self.assertTrue(run['stages'][1]['finished_at'])
        self.assertTrue(run['stages'][2]['error'])

    def test_invalid_inputs_are_reported_before_a_run_is_created(self):
        stage = new_stage('admet')
        self.save([stage])
        service = PipelineService(self.store)
        try:
            with self.assertRaisesRegex(ValueError,'Avaliar ADMET.*Conecte'):
                service.submit(self.token,self.project_id)
            stage['parameters']['base_input_path'] = 'missing-input-directory'
            self.save([stage])
            with self.assertRaisesRegex(ValueError,'não está disponível'):
                service.submit(self.token,self.project_id)
            self.assertEqual(self.store.list_runs(self.token,self.project_id),[])
        finally:
            service.close()

    def test_effective_graph_expansion_and_consensus_inputs_are_validated_early(self):
        service = PipelineService(self.store)
        root = self.store.project_dir(self.project_id) / 'local-input'
        root.mkdir()
        try:
            for operation in ('graphs','expand_similar_compounds','consensus'):
                with self.subTest(operation=operation):
                    stage = new_stage(operation)
                    if operation != 'expand_similar_compounds':
                        stage['parameters']['base_input_path'] = str(root)
                    with self.assertRaisesRegex(ValueError,'Conecte um bloco'):
                        service._validate_parameters(self.project_id,stage,partial=True)
            (root / 'Similarity').mkdir()
            stage = new_stage('graphs')
            stage['parameters']['similarity_path'] = str(root / 'Similarity')
            service._validate_parameters(self.project_id,stage,partial=True)
            for folder in ('Vina','Dock6'):
                (root / folder).mkdir()
            stage = new_stage('consensus')
            stage['parameters']['base_input_path'] = str(root)
            service._validate_parameters(self.project_id,stage,partial=True)
        finally:
            service.close()

    def test_imported_docking_results_connect_through_the_two_visual_consensus_ports(self):
        vina = self.imported(self.upload('MOL1.lig.pdbqt','REMARK VINA RESULT: -6.0\nATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n','vina'))
        vina['parameters']['kind'] = 'vina'
        dock6 = self.imported(self.upload('MOL1_scored.mol2','Grid_Score: -20.0\n@<TRIPOS>MOLECULE\nMOL1\n@<TRIPOS>ATOM\n1 C 0 0 0 C.3\n','dock6'))
        dock6['parameters']['kind'] = 'dock6'
        consensus = new_stage('consensus')
        consensus['bindings'] = {'base_vina_path':{'stage':vina['id']},
                                 'base_dock6_path':{'stage':dock6['id']}}
        self.save([consensus,vina,dock6])
        manager,calls = self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service = PipelineService(self.store)
            try:
                run = self.finish_selected(service,service.submit(self.token,self.project_id))
            finally:
                service.close()
        self.assertEqual(run['status'],'succeeded',run['error'])
        self.assertEqual(len(calls),1)
        operation,params = calls[0]
        self.assertEqual(operation,'consensus')
        self.assertTrue((Path(params['base_vina_path']) / 'MOL1.lig.pdbqt').is_file())
        self.assertTrue((Path(params['base_dock6_path']) / 'MOL1_scored.mol2').is_file())
        self.assertTrue(Path(params['base_input_path']).is_dir())

    def test_empty_selection_and_disabled_pipeline_do_not_create_successful_runs(self):
        stage = new_stage('retrieve_compounds')
        stage['enabled'] = False
        self.save([stage])
        service = PipelineService(self.store)
        try:
            for selection in (None,[],[stage['id']]):
                with self.assertRaises(ValueError):
                    service.submit(self.token,self.project_id,selection)
            self.assertEqual(self.store.list_runs(self.token,self.project_id),[])
        finally:
            service.close()

    def test_project_deletion_after_validation_is_rechecked_before_insertion(self):
        stage = new_stage('retrieve_compounds')
        self.save([stage])
        service = PipelineService(self.store)
        original = service._validate_parameters
        def delete_after_validation(*args,**kwargs):
            original(*args,**kwargs)
            self.store.delete_project(self.token,self.project_id)
        try:
            with patch.object(service,'_validate_parameters',delete_after_validation):
                with self.assertRaises(AccessDenied):
                    service.submit(self.token,self.project_id)
            with self.store.connect() as db:
                self.assertEqual(db.execute('SELECT COUNT(*) FROM runs').fetchone()[0],0)
        finally:
            service.close()

    def test_malformed_connections_raise_validation_errors(self):
        stage = new_stage('admet')
        mutations = [({'id':{}},'identificadores'),({'depends_on':None},'dependências'),
            ({'name':''},'nome válido'),({'bindings':{'base_input_path':{'asset':[]}}},'identificador'),
            ({'bindings':{'base_input_path':{'asset':uuid4().hex,'selector':12}}},'nome de arquivo')]
        for changes,message in mutations:
            with self.subTest(changes=changes):
                broken = copy.deepcopy(stage)
                broken.update(changes)
                with self.assertRaisesRegex(ValueError,message):
                    validate_pipeline([broken])
        unsupported = new_stage('retrieve_compounds')
        unsupported['bindings']['base_input_path'] = {'asset':uuid4().hex}
        with self.assertRaisesRegex(ValueError,'entrada inválida'):
            validate_pipeline([unsupported])
        structures = new_stage('retrieve_structures')
        analysis = new_stage('admet')
        analysis['bindings']['base_input_path'] = {'stage':structures['id']}
        with self.assertRaisesRegex(ValueError,'dados incompatíveis'):
            validate_pipeline([structures,analysis])

    def test_legacy_automatic_pipeline_requires_confirmation_and_merges_only_selected_files(self):
        source=new_stage('retrieve_compounds')
        source['process_all']=True
        analysis=new_stage('admet')
        analysis['bindings']['base_input_path']={'stage':source['id'],'selector':'compounds.csv'}
        self.save([source,analysis])
        manager,calls=self.immediate_manager()
        class MultipleFiles(manager):
            def submit(self,operation,parameters):
                job=super().submit(operation,parameters)
                if operation=='retrieve_compounds':
                    extra=self.config.workspace/'extra.csv'
                    extra.write_text('molecule_chembl_id,canonical_smiles\nCHEMBL2,CCC\n')
                    unrelated=self.config.workspace/'report.csv'
                    unrelated.write_text('status\nok\n')
                    job['result']['artifacts'] += [str(extra),str(unrelated)]
                return job
        with patch('biomolexplorer.pipeline.JobManager',MultipleFiles):
            service=PipelineService(self.store)
            try:
                run=self.wait(service.submit(self.token,self.project_id))
                self.assertEqual(run['status'],'awaiting_input')
                self.assertEqual(len(calls),1)
                analysis['bindings']['base_input_path']={'sources':[
                    {'stage':source['id'],'selector':'compounds.csv'},
                    {'stage':source['id'],'selector':'extra.csv'}]}
                run=self.wait(service.resume(self.token,run['id'],analysis))
            finally:service.close()
        self.assertEqual(run['status'],'succeeded',run['error'])
        params=calls[1][1]
        import csv
        with (Path(params['base_input_path'])/params['input_file']).open() as stream:
            self.assertEqual({r['molecule_chembl_id'] for r in csv.DictReader(stream)},{'CHEMBL1','CHEMBL2'})
        self.assertEqual(len(run['stages'][1]['input_files']['base_input_path']),2)

    def test_curation_pauses_at_each_descendant_and_resumes_same_run(self):
        source=new_stage('retrieve_compounds')
        first,second=new_stage('admet'),new_stage('admet')
        first['bindings']['base_input_path']={'stage':source['id']}
        second['bindings']['base_input_path']={'stage':first['id']}
        self.save([second,first,source])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store)
            try:
                run=self.wait(service.submit(self.token,self.project_id))
                run_id=run['id']
                self.assertEqual(run['status'],'awaiting_input')
                self.assertEqual([s['status'] for s in run['stages']],['succeeded','awaiting_input','queued'])
                self.assertIn(first['name'],run['error'])
                self.assertEqual(len(calls),1)
                with self.assertRaisesRegex(ValueError,'explicitamente'):
                    service.resume(self.token,run_id,first)
                with self.assertRaisesRegex(ValueError,'execução ativa'):
                    service.submit(self.token,self.project_id)
                with self.assertRaises(ValueError):self.store.delete_project(self.token,self.project_id)
                first['bindings']['base_input_path']['selector']='compounds.csv'
                run=self.wait(service.resume(self.token,run_id,first))
                self.assertEqual(run['id'],run_id)
                self.assertEqual(run['status'],'awaiting_input')
                self.assertEqual([s['status'] for s in run['stages']],['succeeded','succeeded','awaiting_input'])
                second['bindings']['base_input_path']['selector']='compounds.csv'
                run=self.wait(service.resume(self.token,run_id,second))
                self.assertEqual(run['status'],'succeeded',run['error'])
                self.assertEqual([op for op,_ in calls],['retrieve_compounds','admet','admet'])
                with self.assertRaisesRegex(ValueError,'não está aguardando'):
                    service.resume(self.token,run_id,second)
            finally:service.close()

    def test_selected_inputs_clean_bad_records_before_merging_and_publish_report(self):
        source=new_stage('retrieve_compounds');analysis=new_stage('admet')
        analysis['bindings']['base_input_path']={'stage':source['id'],'selector':'auto'}
        self.save([source,analysis])
        manager,calls=self.immediate_manager()
        class BadCachedData(manager):
            def submit(self,operation,parameters):
                job=super().submit(operation,parameters)
                if operation=='retrieve_compounds':
                    path=Path(job['result']['artifacts'][0])
                    path.write_text(path.read_text()+'BAD,invalid\nBLANK,\n')
                return job
        with patch('biomolexplorer.pipeline.JobManager',BadCachedData):
            service=PipelineService(self.store)
            try:run=self.finish_selected(service,service.submit(self.token,self.project_id))
            finally:service.close()
        self.assertEqual(run['status'],'succeeded',run['error'])
        item=run['stages'][1]
        self.assertEqual(item['excluded_records'],2)
        report=next(Path(p) for p in item['artifacts'] if Path(p).name=='molecule_exclusions.json')
        self.assertEqual(json.loads(report.read_text())['excluded_records'],2)
        selected=Path(calls[1][1]['base_input_path'])/calls[1][1]['input_file']
        self.assertNotIn('BAD',selected.read_text())

    def test_pending_run_survives_service_restart_and_can_be_cancelled(self):
        source=new_stage('retrieve_compounds');source['process_all']=False
        stage=new_stage('admet');stage['bindings']['base_input_path']={'stage':source['id']}
        self.save([source,stage])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store)
            try:run=self.wait(service.submit(self.token,self.project_id))
            finally:service.close()
            service=PipelineService(self.store)
            try:
                self.assertEqual(self.store.get_run(self.token,run['id'])['status'],'awaiting_input')
                stage['bindings']['base_input_path']['selector']='compounds.csv'
                resumed=self.wait(service.resume(self.token,run['id'],stage))
                self.assertEqual(resumed['status'],'succeeded',resumed['error'])
                run=self.wait(service.submit(self.token,self.project_id,reuse_results=False))
                service.cancel(self.token,run['id'])
                cancelled=self.store.get_run(self.token,run['id'])
                self.assertEqual(cancelled['status'],'cancelled')
                self.assertEqual(cancelled['stages'][1]['status'],'cancelled')
            finally:service.close()

    def test_curation_selects_only_requested_file_and_invalid_choice_stays_pending(self):
        source=new_stage('retrieve_compounds');source['process_all']=False
        stage=new_stage('admet');stage['bindings']['base_input_path']={'stage':source['id']}
        self.save([source,stage])
        manager,calls=self.immediate_manager()
        class MultipleFiles(manager):
            def submit(self,operation,parameters):
                job=super().submit(operation,parameters)
                if operation=='retrieve_compounds':
                    extra=self.config.workspace/'extra.csv'
                    extra.write_text('molecule_chembl_id,canonical_smiles\nCHEMBL2,CCC\n')
                    job['result']['artifacts'].append(str(extra))
                return job
        with patch('biomolexplorer.pipeline.JobManager',MultipleFiles):
            service=PipelineService(self.store)
            try:
                run=self.wait(service.submit(self.token,self.project_id))
                stage['bindings']['base_input_path']['selector']='missing.csv'
                with self.assertRaisesRegex(ValueError,'não está nos resultados'):
                    service.resume(self.token,run['id'],stage)
                self.assertEqual(self.store.get_run(self.token,run['id'])['status'],'awaiting_input')
                stage['bindings']['base_input_path']['selector']='extra.csv'
                run=self.wait(service.resume(self.token,run['id'],stage))
                self.assertEqual(run['status'],'succeeded',run['error'])
                self.assertEqual([Path(p).name for p in run['stages'][1]['input_files']['base_input_path']],['extra.csv'])
            finally:service.close()

    def test_curated_run_does_not_hold_an_executor_slot(self):
        source=new_stage('retrieve_compounds');source['process_all']=False
        stage=new_stage('admet');stage['bindings']['base_input_path']={'stage':source['id']}
        self.save([source,stage])
        other=self.store.create_project(self.token,'Independent')
        self.store.save_pipeline(self.token,other['id'],[new_stage('retrieve_compounds')],0)
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store,max_runs=1)
            try:
                waiting=self.wait(service.submit(self.token,self.project_id))
                independent=self.wait(service.submit(self.token,other['id']))
                self.assertEqual(waiting['status'],'awaiting_input')
                self.assertEqual(independent['status'],'succeeded',independent['error'])
            finally:service.close()

    def test_all_modes_require_curation_including_legacy_automatic_stages(self):
        automatic,curated,unrelated=[new_stage('retrieve_compounds') for _ in range(3)]
        curated['process_all']=False
        join=new_stage('admet')
        join['bindings']['base_input_path']={'sources':[{'stage':automatic['id']},{'stage':curated['id']}]}
        downstream=new_stage('fingerprints');downstream['bindings']['base_input_path']={'stage':join['id']}
        modes=execution_modes([downstream,join,unrelated,curated,automatic])
        self.assertEqual(modes[join['id']],'curated')
        self.assertEqual(modes[downstream['id']],'curated')
        self.assertEqual(modes[unrelated['id']],'curated')
        automatic.pop('process_all')
        self.assertEqual(execution_modes([automatic])[automatic['id']],'curated')
        automatic['process_all']='yes'
        with self.assertRaisesRegex(ValueError,'marcada ou desmarcada'):validate_pipeline([automatic])

    def test_structure_selection_preserves_metadata_for_only_the_selected_complex(self):
        source=new_stage('retrieve_structures');stage=new_stage('prepare_structures')
        root=self.store.project_dir(self.project_id)/'structures'/'MeuAlvo';root.mkdir(parents=True)
        atom='ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n'
        for code in ('1ABC','2ABC'):(root/(code+'.pdb')).write_text(atom)
        (root/'pdb_codes.csv').write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RESOLUTION\n1ABC,LIG,1,A,1.5\n2ABC,LIG,1,A,1.5\n')
        results={source['id']:[str(p) for p in root.iterdir()]}
        stage['bindings']['base_input_path']={'stage':source['id'],'selector':'1ABC.pdb'}
        service=PipelineService(self.store)
        try:
            params=service._resolve(self.project_id,self.store.user(self.token)['id'],stage,results)
            selected=Path(params['base_input_path'])/'MeuAlvo'
            self.assertTrue((selected/'1ABC.pdb').exists())
            self.assertFalse((selected/'2ABC.pdb').exists())
            self.assertNotIn('2ABC',(selected/'pdb_codes.csv').read_text())
        finally:service.close()

    def test_pending_block_can_be_configured_and_resumed_through_the_ui(self):
        try:
            import flet as ft
            from biomolexplorer.ui.app import WorkspaceUI
        except ImportError:self.skipTest('Install the ui extra')
        import asyncio
        from types import SimpleNamespace
        from unittest.mock import AsyncMock
        source=new_stage('retrieve_compounds');source['process_all']=False
        stage=new_stage('admet');stage['bindings']['base_input_path']={'stage':source['id']}
        self.save([source,stage])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store)
            try:
                run=self.wait(service.submit(self.token,self.project_id))
                dialogs=[]
                page=SimpleNamespace(width=1440,height=1000,update=lambda:None,show_dialog=dialogs.append,
                    run_task=lambda *args:SimpleNamespace(done=lambda:False))
                ui=WorkspaceUI(page,self.store,service)
                ui.token=self.token;ui.current=self.store.project(self.token,self.project_id)
                ui.base_pipeline=copy.deepcopy(ui.current['pipeline'])
                ui.draw_project=AsyncMock()
                async def local_call(function,*args,**kwargs):return function(*args,**kwargs)
                ui.call=local_call
                asyncio.run(ui.open_stage_dialog(stage['id'],waiting_run_id=run['id']))
                dialog=dialogs[-1]
                def walk(control):
                    yield control
                    for child in getattr(control,'controls',[]) or []:yield from walk(child)
                    child=getattr(control,'content',None)
                    if isinstance(child,ft.Control):yield from walk(child)
                selector=next(c for c in walk(dialog.content) if isinstance(c,ft.Dropdown) and c.label=='Resultado usado nesta entrada')
                self.assertIn('compounds.csv',[o.key for o in selector.options])
                selector.value='compounds.csv'
                self.assertEqual(dialog.actions[-1].content,'Aplicar e continuar pipeline')
                asyncio.run(dialog.actions[-1].on_click(None))
                finished=self.wait(run)
                self.assertEqual(finished['status'],'succeeded',finished['error'])
                self.assertFalse(dialog.open)
                saved=self.store.project(self.token,self.project_id)
                configured=next(s for s in saved['pipeline'] if s['id']==stage['id'])
                self.assertEqual(configured['bindings']['base_input_path']['selector'],'compounds.csv')
            finally:service.close()


if __name__ == '__main__':
    unittest.main()
