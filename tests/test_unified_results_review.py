"""Read-only consolidation, stable shared forms and alert information events."""
import asyncio
import copy
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock,patch
import flet as ft
from biomolexplorer.catalog import new_stage
from biomolexplorer.result_visibility import consolidated_run
from biomolexplorer.ui.guided import GuidedForm
from biomolexplorer.ui.localization import LocalizedPage,Translator
from biomolexplorer.ui.settings_layout import fit_dialog
from test_frontend_navigation import Page,walk


class ConsolidatedViewTests(unittest.TestCase):
    def test_cancelled_latest_uses_original_producer_once_without_editing_history(self):
        config=new_stage('retrieve_compounds');old=dict(id=config['id'],name=config['name'],operation=config['operation'],status='succeeded',artifacts=['molecules.csv'],publication_run_id='original',configuration=config)
        latest=dict(id='cancelled',status='cancelled',stages=[dict(old,status='cancelled',artifacts=[])])
        before=copy.deepcopy((latest,old));view=consolidated_run([latest],[config],[old])
        self.assertEqual(len(view['stages']),1);self.assertEqual(view['stages'][0]['_result_run_id'],'original')
        self.assertEqual(view['stages'][0]['artifacts'],['molecules.csv']);self.assertEqual((latest,old),before)

    def test_new_completed_result_supersedes_old_and_removed_block_is_hidden(self):
        config=new_stage('retrieve_compounds');stage=dict(id=config['id'],name=config['name'],operation=config['operation'],status='succeeded',artifacts=['new.csv'])
        removed=dict(stage,id='removed',artifacts=['historical.csv'])
        view=consolidated_run([dict(id='new',stages=[stage,removed])],[config],[dict(stage,publication_run_id='old',artifacts=['old.csv'])])
        self.assertEqual([s['id'] for s in view['stages']],[config['id']]);self.assertEqual(view['stages'][0]['_result_run_id'],'new')


class SharedFormReviewTests(unittest.TestCase):
    def ui(self,stages,language='en'):
        return SimpleNamespace(current={'pipeline':stages},service=SimpleNamespace(dock6_path=None),page=LocalizedPage(Page(),language),tr=Translator(language))

    def test_repeated_layout_does_not_accumulate_padding_or_change_configuration(self):
        for operation in ('retrieve_compounds','admet','fingerprints','redocking'):
            stage=new_stage(operation);ui=self.ui([stage]);form=GuidedForm(ui,stage,[],True)
            before=form.read();first=form.layout();count=len(list(walk(first)))
            for _ in range(8):
                layout=form.layout()
                self.assertEqual(len(list(walk(layout))),count,operation)
                self.assertEqual(form.read(),before,operation)

    def test_every_admet_alert_click_opens_its_own_description_in_both_languages(self):
        from biomolexplorer.admet_config import ALERT_CATALOGS
        for language in ('pt','en'):
            stage=new_stage('admet');ui=self.ui([stage],language);form=GuidedForm(ui,stage,[],True);layout=form.layout();before=form.read()
            buttons=[c for c in walk(layout) if isinstance(c,ft.IconButton) and c.icon==ft.Icons.INFO_OUTLINE and c.tooltip in ALERT_CATALOGS.values()]
            self.assertEqual(len(buttons),len(ALERT_CATALOGS))
            for name,button in zip(ALERT_CATALOGS,buttons):
                button.on_click(None);dialog=ui.page._page.dialogs[-1]
                self.assertEqual(dialog.title.value,name)
                expected=ui.tr(ALERT_CATALOGS[name])
                self.assertIn(expected,[c.value for c in walk(dialog.content) if isinstance(c,ft.Text)])
                self.assertEqual(form.read(),before)
                dialog.actions[0].on_click(None)
                self.assertFalse(dialog.open)


class ChemistryBoundaryReviewTests(unittest.TestCase):
    def test_minimize_has_one_explicit_prep_false_even_after_repeat_configuration(self):
        from biomolexplorer.chimera_preparation import ligand_charging
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'input.mol2';source.write_text('unused fixture')
            molecule=SimpleNamespace(GetAtoms=lambda:[])
            phase=dict(intermediate=source,output=root/'output.mol2',model='14sb',method='am1',tail='\nminimize spec #0 nsteps 100 prep false\n')
            with patch('biomolexplorer.ligand_mol2.ensure_ligand_mol2'),patch('biomolexplorer.dock6_validation.validate_mol2',return_value=molecule):
                script=root/'ligand.com';ligand_charging(script,phase)
                self.assertEqual(script.read_text().count('prep false'),1)
                self.assertIn('nsteps 100 prep false',script.read_text())

    def test_syntax_only_intermediate_defers_sanitization_but_final_export_remains_strict(self):
        from biomolexplorer.dock6_validation import validate_mol2
        from rdkit import Chem
        # Pyrrole N needs its explicit H before MOL2 aromatic sanitization.
        text='@<TRIPOS>MOLECULE\nPYRROLE\n5 5 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n'+''.join(f'{i+1} A{i+1} {i}.0 0.0 0.0 '+('N.ar' if i==0 else 'C.ar')+' 1 LIG 0.0\n' for i in range(5))+'@<TRIPOS>BOND\n'+''.join(f'{i+1} {i+1} {(i+1)%5+1} ar\n' for i in range(5))+'@<TRIPOS>SUBSTRUCTURE\n1 LIG 1\n'
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'intermediate.mol2';path.write_text(text)
            with patch.object(Chem,'SanitizeMol',wraps=Chem.SanitizeMol) as sanitize:
                validate_mol2(path,sanitize=False);sanitize.assert_not_called()
            with self.assertRaises(ValueError):validate_mol2(path)

class ReviewHistoryTests(unittest.TestCase):
    def test_order_reversal_uses_existing_history_without_copying_artifact_blobs(self):
        import test_pipeline_execution as fixture
        f=fixture.PipelineExecutionTests();f.setUp();self.addCleanup(f.tearDown)
        asset=f.upload();a,b=new_stage('retrieve_compounds'),new_stage('retrieve_structures');f.save([a,b])
        f.store.save_execution_order(f.token,f.project_id,[a['id'],b['id']])
        blobs=f.store.project_dir(f.project_id)/'.history/blobs';before={p.name for p in blobs.iterdir()}
        f.store.save_execution_order(f.token,f.project_id,[b['id'],a['id']])
        event=f.store.history(f.token,f.project_id)[0]
        self.assertEqual({p.name for p in blobs.iterdir()},before)
        self.assertEqual(event['action'],'execution_order')
        f.store.rollback(f.token,f.project_id,event['id'])
        self.assertEqual(f.store.execution_order(f.token,f.project_id)['order'],[a['id'],b['id']])
        self.assertTrue(f.store.asset_path(f.token,f.project_id,asset).is_file())

class CheckpointConsolidationTests(unittest.TestCase):
    def test_cancel_restart_reuses_batch_and_keeps_checkpoint_visible_under_latest_run(self):
        import json,threading,time
        from uuid import uuid4
        import test_pipeline_execution as fixture
        from biomolexplorer.pipeline import PipelineService
        f=fixture.PipelineExecutionTests();f.setUp();self.addCleanup(f.tearDown)
        a,b=f.upload('a.csv'),f.upload('b.csv',content='name,smiles\nPROPANE,CCC\n')
        stage=new_stage('admet');stage['input_processing']='individual';stage['bindings']['base_input_path']={'sources':[{'asset':a},{'asset':b}]};f.save([stage])
        original,calls=f.immediate_manager();entered,release=threading.Event(),threading.Event()
        class BlockFirst(original):
            def get(manager,job_id):
                if len(calls)==1:entered.set();release.wait(8)
                return super().get(job_id)
        with patch('biomolexplorer.pipeline.JobManager',BlockFirst):
            service=PipelineService(f.store)
            try:
                run=service.submit(f.token,f.project_id)
                self.assertTrue(entered.wait(8));service.pause(f.token,run['id']);release.set();paused=f.wait(run)
                self.assertEqual(paused['status'],'awaiting_input',paused.get('error'))
                self.assertEqual(len(paused['stages'][0]['batches']),1)
                service.cancel(f.token,run['id'])
                checkpoint=f.store.available_results(f.token,f.project_id,include_checkpoints=True)
                self.assertEqual(len(checkpoint),1);self.assertTrue(checkpoint[0]['coverage_pending'])
                self.assertEqual(checkpoint[0]['reuse_provenance']['policy'],'verified_completed_units_v1')
                self.assertEqual(len(checkpoint[0]['reuse_provenance']['units']),1)
                key=checkpoint[0]['publication_run_id'];self.assertEqual(key,run['id'])
                files=checkpoint[0]['artifacts'][:]
                with f.store.connect() as db:
                    empty=dict(paused['stages'][0],status='cancelled',artifacts=[],batches=[])
                    db.execute('INSERT INTO runs VALUES(?,?,?,?,?,?,?,?)',(uuid4().hex,f.project_id,f.store.user(f.token)['id'],'cancelled',json.dumps([empty]),time.time()+1,time.time()+1,'cancelled'))
                available=f.store.available_results(f.token,f.project_id,include_checkpoints=True)
                self.assertEqual(available[0]['artifacts'],files)
                view=consolidated_run(f.store.list_runs(f.token,f.project_id),[stage],available)
                self.assertEqual(view['stages'][0]['_result_run_id'],key)
            finally:release.set();service.close()
            service=PipelineService(f.store)
            try:
                finished=f.wait(service.submit(f.token,f.project_id))
                self.assertEqual(finished['status'],'succeeded',finished.get('error'))
                self.assertEqual(len(calls),2)
                self.assertEqual(len(finished['stages'][0]['batches']),2)
                self.assertTrue(set(files)<=set(finished['stages'][0]['artifacts']))
            finally:service.close()

class PreparationCacheReviewTests(unittest.TestCase):
    def test_ligand_protocol_change_does_not_change_receptor_identity(self):
        from biomolexplorer.structure_cache import StructureCache
        from biomolexplorer.stage_cache import file_digest
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'1ABC.pdb').write_text('END\n')
            row=['1ABC','LIG',1,'A'];settings={'1ABC|LIG|1|A':{'preparation_backend':'legacy_chimera'}}
            with patch('biomolexplorer.structure_cache.preparation_identity',return_value={'installed_chimera':'unchanged'}):
                before=StructureCache(root,root/'prepared',settings,'gas',7.4).keys(row)
                def changed(path):return 'changed-ligand-protocol' if Path(path).name=='chimera_preparation.py' else file_digest(path)
                with patch('biomolexplorer.structure_cache.file_digest',side_effect=changed):
                    after=StructureCache(root,root/'prepared',settings,'gas',7.4).keys(row)
            self.assertEqual(before[0],after[0]);self.assertNotEqual(before[1],after[1])

class ComplementaryCheckpointTests(unittest.TestCase):
    def test_disjoint_historical_checkpoints_form_one_block_and_are_reused_without_jobs(self):
        import json,time
        from uuid import uuid4
        import test_pipeline_execution as fixture
        from biomolexplorer.pipeline import PipelineService
        from biomolexplorer.ui.app import WorkspaceUI
        f=fixture.PipelineExecutionTests();f.setUp();self.addCleanup(f.tearDown)
        a,b=f.upload('first.csv'),f.upload('second.csv',content='name,smiles\nPROPANE,CCC\n')
        stage=new_stage('admet');stage['input_processing']='individual';stage['bindings']['base_input_path']={'sources':[{'asset':a},{'asset':b}]};f.save([stage])
        manager,calls=f.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(f.store)
            try:completed=f.wait(service.submit(f.token,f.project_id))
            finally:service.close()
            self.assertEqual(completed['status'],'succeeded',completed.get('error'))
            original=completed['stages'][0];expected=set(original['artifacts'])
            # Model independently persisted checkpoints from interrupted runs.
            with f.store.connect() as db:
                for index,batch in enumerate(original['batches']):
                    part=copy.deepcopy(original);part.update(status='cancelled',batches=[batch],artifacts=batch['artifacts'],artifact_manifest=batch['artifact_manifest'])
                    if index==0:db.execute('UPDATE runs SET status=?,stages=? WHERE id=?',('cancelled',json.dumps([part]),completed['id']))
                    else:db.execute('INSERT INTO runs VALUES(?,?,?,?,?,?,?,?)',(uuid4().hex,f.project_id,f.store.user(f.token)['id'],'cancelled',json.dumps([part]),time.time()+index,time.time()+index,'cancelled'))
            available=f.store.available_results(f.token,f.project_id,include_checkpoints=True)
            self.assertEqual(len(available),1);self.assertEqual(len(available[0]['_result_sources']),2)
            self.assertEqual(set(available[0]['artifacts']),expected)
            # A new paused run's completed unit is visible alongside the older
            # checkpoint, before cancellation or completion of the whole stage.
            with f.store.connect() as db:
                latest=db.execute('SELECT id,stages FROM runs WHERE project_id=? ORDER BY created DESC LIMIT 1',(f.project_id,)).fetchone()
                active=json.loads(latest[1]);active[0]['status']='awaiting_input'
                db.execute('UPDATE runs SET status=?,stages=? WHERE id=?',('awaiting_input',json.dumps(active),latest[0]))
            live=f.store.available_results(f.token,f.project_id,include_checkpoints=True)
            self.assertEqual(set(live[0]['artifacts']),expected)
            self.assertEqual(len(live[0]['_result_sources']),2)
            with f.store.connect() as db:db.execute('UPDATE runs SET status=?,stages=? WHERE id=?',('cancelled',latest[1],latest[0]))
            async def display():
                ui=WorkspaceUI(Page(),f.store,SimpleNamespace(dock6_path=None));ui.token=f.token;ui.current=f.store.project(f.token,f.project_id)
                self.addCleanup(ui.structure_viewer.close)
                async def call(fn,*args,**kwargs):return fn(*args,**kwargs)
                ui.call=call
                controls=await ui.runs_view(True)
                self.assertEqual(len(ui.stage_results),1)
                result=next(iter(ui.stage_results.values()));await result.load()
                self.assertEqual(len(result.checkpoint_children),2)
                paths=[file['path'] for child in result.checkpoint_children for file in child.child.files]
                self.assertEqual(set(paths),expected);self.assertEqual(len(paths),len(set(paths)))
                self.assertFalse(any(isinstance(c,ft.Text) and c.value=='Resultados disponíveis do pipeline' for c in walk(controls)))
                result.close();self.assertTrue(all(not child.valid() for child in result.checkpoint_children))
            asyncio.run(display())
            service=PipelineService(f.store)
            try:
                resumed=f.wait(service.submit(f.token,f.project_id));self.assertEqual(resumed['status'],'succeeded',resumed.get('error'))
                self.assertEqual(len(calls),2)
                self.assertEqual(set(resumed['stages'][0]['artifacts']),expected)
            finally:service.close()

class ResultIntegrityReviewTests(unittest.TestCase):
    def test_tampered_latest_output_is_hidden_without_deleting_history(self):
        import test_pipeline_execution as fixture
        from biomolexplorer.pipeline import PipelineService
        from biomolexplorer.ui.app import WorkspaceUI
        f=fixture.PipelineExecutionTests();f.setUp();self.addCleanup(f.tearDown)
        stage=new_stage('retrieve_compounds');f.save([stage]);manager,_=f.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(f.store)
            try:run=f.wait(service.submit(f.token,f.project_id))
            finally:service.close()
        path=Path(run['stages'][0]['artifacts'][0]);path.write_text('altered outside BioMolExplorer')
        async def display():
            ui=WorkspaceUI(Page(),f.store,SimpleNamespace(dock6_path=None));ui.token=f.token;ui.current=f.store.project(f.token,f.project_id)
            self.addCleanup(ui.structure_viewer.close)
            async def call(fn,*args,**kwargs):return fn(*args,**kwargs)
            ui.call=call;await ui.runs_view(True)
            result=next(iter(ui.stage_results.values()))
            self.assertTrue(result.stage['results_stale']);self.assertEqual(result.stage['artifacts'],[])
            await result.load();self.assertIsNone(result.child)
        asyncio.run(display())
        self.assertTrue(path.is_file())
        original=f.store.get_run(f.token,run['id'])['stages'][0]
        self.assertEqual(original['status'],'succeeded');self.assertIn(str(path),original['artifacts'])
