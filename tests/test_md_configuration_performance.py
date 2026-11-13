"""MD configuration reuses input catalogs and preserves target-specific poses."""
import json
import shutil
from pathlib import Path
import unittest
from unittest.mock import patch

import test_md_pair_selection as pair_fixtures
import test_consensus_target_actions as consensus_fixtures
from test_docking_scene import atom
from biomolexplorer.md_inputs import MDInputs,pair_id,normalize_receptor_variants
from biomolexplorer.catalog import new_stage
from biomolexplorer.ui.localization import Translator


class MDConfigurationPerformanceTests(unittest.TestCase):
    setUp=pair_fixtures.AmberPairSelectionTests.setUp
    form=pair_fixtures.AmberPairSelectionTests.form

    def test_refresh_and_pair_selection_reuse_unchanged_input_catalog(self):
        form=self.form();settings=form.md_settings
        with patch.object(MDInputs,'files',side_effect=AssertionError('Unchanged catalog reread')):
            settings.refresh()
            settings.search.value='lig';settings.filter(None)
            settings.all(None)
            form.read()

    def test_changed_input_invalidates_dialog_catalog(self):
        form=self.form();settings=form.md_settings
        path=Path(form.stage['parameters']['base_selected_mols'])
        path.write_text(path.read_text().replace('LIG','NEW'))
        original=MDInputs.files
        with patch.object(MDInputs,'files',autospec=True,side_effect=original) as files:
            settings.refresh()
        self.assertEqual([call.args[2] for call in files.call_args_list],['base_selected_mols'])

    def test_hidden_md_selector_never_builds_pose_catalog(self):
        form=self.form();editor=form.input_editors['base_input_path']
        with patch('biomolexplorer.artifacts.records',side_effect=AssertionError('Hidden dropdown indexed artifacts')):
            self.assertEqual(len(editor.selectors('stage:'+self.producer['id'])),1)

    def test_parameter_text_and_values_follow_language_without_changing_api_keys(self):
        form=self.form('en');self.ui.page.add(form.layout())
        region=form.md_settings.controls['amber_score_movable_region']
        self.assertEqual(region.options[0].key,'ligand')
        self.assertEqual(region.options[0].text,'Ligand only')
        self.assertEqual(Translator('en')('DINÂMICA MOLECULAR'),'MOLECULAR DYNAMICS')
        self.ui.page.set_language('pt');self.ui.page.update()
        self.assertEqual(region.options[0].text,'Somente o ligante')
        self.assertEqual(form.read()['parameters']['amber_score_movable_region'],'ligand')

    def test_execution_progress_translates_without_changing_molecular_labels(self):
        tr=Translator('en')
        self.assertEqual(tr('Dinâmica: 1/2 · 4M0E / CHEMBL92663'),'Amber Score: 1/2 · 4M0E / CHEMBL92663')
        self.assertEqual(tr('Dinâmica concluída: 1 sucessos; 0 falhas.'),'Amber Score completed: 1 succeeded; 0 failed.')
        self.assertEqual(tr('Etapa da dinâmica: input_validation'),'Amber Score phase: input_validation')


class MDConsensusPoseTests(unittest.TestCase):
    setUp=consensus_fixtures.ConsensusTargetActionsTests.setUp
    fixture=consensus_fixtures.ConsensusTargetActionsTests.fixture

    def test_resolve_selects_consensus_pose_for_actual_receptor_and_keeps_pair_id(self):
        service,selection,tables=self.fixture()
        original=Path(selection['table']);folder=self.root/'runs'/self.rid/self.sid
        folder.mkdir(parents=True);table=folder/original.name;shutil.copy2(original,table)
        with self.store.connect() as db:
            stages=json.loads(db.execute('SELECT stages FROM runs WHERE id=?',(self.rid,)).fetchone()[0])
            db.execute('INSERT INTO runs(id,project_id,user_id,status,stages,created,updated) VALUES(?,?,?,?,?,?,?)',
                ('origin-run',self.pid,self.store.user(self.token)['id'],'succeeded',json.dumps(stages[:2]),1,1))
            stages[-1]['artifacts']=[str(table)]
            db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),self.rid))
        receptor=self.root/'2BBB_A.noH.pdb';receptor.write_text(atom(1,'CA','ALA','A',2,10))
        stage=new_stage('dynamics_dock6');stage['parameters'].update(base_input_path=str(receptor),base_selected_mols=str(table))
        inputs=MDInputs(self.store,self.pid,{})
        rec=inputs.files(stage,'base_input_path')[0];lig=inputs.files(stage,'base_selected_mols')[0]
        identity=pair_id(rec,lig);stage['parameters']['selected_pairs']=[identity]
        result=inputs.resolve(stage)['md_inputs'][0]
        self.assertEqual(result['id'],identity)
        self.assertEqual(result['ligand']['compound']['receptor_id'],'2BBB_A')
        self.assertEqual(Path(result['ligand']['compound']['conformer_file']).name,'2BBB_A.pdbqt')
        self.assertEqual(Path(result['ligand']['compound']['conformer_file']).parent.name,'dock6')

    def test_legacy_pdb_filename_is_accepted_only_for_a_unique_consensus_target(self):
        from biomolexplorer.md_consensus import ConsensusPoses
        resolver=ConsensusPoses(self.store,self.pid)
        ligand={'compound':{'receptor_id':'1AAA_A','target_scores':json.dumps([{'receptor_id':'1AAA_A'}])}}
        self.assertIs(resolver.ligand({'path':self.root/'1AAA.pdb'},ligand),ligand)
        ligand['compound']['target_scores']=json.dumps([{'receptor_id':'1AAA_A'},{'receptor_id':'1AAA_B'}])
        with self.assertRaisesRegex(ValueError,'alvo correspondente'):resolver.ligand({'path':self.root/'1AAA.pdb'},ligand)

    def test_explicit_disulfide_and_backbone_gaps_keep_coordinates(self):
        source=self.root/'source.pdb';destination=self.root/'prepared.pdb'
        source.write_text(atom(1,'SG','CYS','A',1,0,element='S')+atom(2,'C','CYS','A',1,1)+
            atom(3,'SG','CYS','A',2,2,element='S')+atom(4,'N','CYS','A',2,10)+'CONECT    1    3\nCONECT    3    1\n')
        changes=normalize_receptor_variants(source,destination)
        before=[line for line in source.read_text().splitlines() if line.startswith('ATOM')]
        after=[line for line in destination.read_text().splitlines() if line.startswith('ATOM')]
        self.assertEqual([line[30:] for line in before],[line[30:] for line in after])
        self.assertTrue(all(line[17:20]=='CYX' for line in after))
        self.assertIn('TER\n',destination.read_text())
        self.assertEqual(len(changes),3)

    def test_conflicting_thiol_hydrogen_is_rejected(self):
        source=self.root/'source.pdb';destination=self.root/'prepared.pdb'
        source.write_text(atom(1,'SG','CYS','A',1,0,element='S')+atom(2,'SG','CYS','A',2,2,element='S')+
            atom(3,'HG','CYS','A',1,1,element='H')+'CONECT    1    2\n')
        with self.assertRaisesRegex(ValueError,'dissulfeto'):normalize_receptor_variants(source,destination)
        self.assertFalse(destination.exists())


class MDEnergyPrecisionTests(unittest.TestCase):
    def test_native_float32_total_accepts_only_its_actual_rounding_bound(self):
        import tempfile
        from biomolexplorer.md_engine import energies
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'amber_scored.mol2'
            text='Amber_Score: 2391.533691\nAmber_complex_energy: -17393.140872\nAmber_receptor_energy: 19769.696901\nAmber_ligand_energy: 14.977599\n'
            path.write_text(text)
            self.assertEqual(energies(path)['amber_score'],2391.533691)
            path.write_text(text.replace('2391.533691','2391.543691'))
            with self.assertRaisesRegex(ValueError,'Termos AMBER'):energies(path)


class MDDialogTests(unittest.TestCase):
    from test_pipeline_editor_regressions import PipelineEditorRegressionTests
    setUp=PipelineEditorRegressionTests.setUp

    def test_unchanged_md_apply_closes_without_reading_result_history_or_saving(self):
        import asyncio
        import copy
        from unittest.mock import AsyncMock
        async def scenario():
            stage=new_stage('dynamics_dock6')
            self.ui.current['pipeline']=[stage];self.ui.dirty=False
            with patch('biomolexplorer.ui.guided.GuidedForm.read',return_value=copy.deepcopy(stage)):
                with patch.object(self.store,'list_runs',side_effect=AssertionError('Full execution history loaded')):
                    await self.ui.open_stage_dialog(stage['id'])
                dialog=self.page.dialogs[-1]
                self.ui.save_draft=AsyncMock(side_effect=AssertionError('Unchanged configuration saved'))
                await dialog.actions[-1].on_click(None)
            self.assertFalse(dialog.open)
            self.ui.save_draft.assert_not_awaited()
        asyncio.run(scenario())


class MDResolutionIsolationTests(unittest.TestCase):
    from test_molecular_dynamics import MDPairCacheTests
    setUp=MDPairCacheTests.setUp
    compute=MDPairCacheTests.compute
    run_batch=MDPairCacheTests.run_batch

    def test_missing_consensus_target_does_not_execute_wrong_pose_or_block_other_pairs(self):
        broken=dict(self.pairs[0],input_error='Target-specific pose unavailable')
        details,records=self.run_batch('resolved-inputs',[broken,self.pairs[1]])
        self.assertEqual(details['md_failures'],1)
        self.assertEqual(self.runner.call_count,1)
        self.assertEqual(self.runner.call_args.args[0]['id'],self.pairs[1]['id'])


class MDAuditTests(unittest.TestCase):
    from test_consensus_history import ConsensusHistoryTests
    setUp=ConsensusHistoryTests.setUp
    tearDown=ConsensusHistoryTests.tearDown
    completed=ConsensusHistoryTests.completed

    def test_md_only_save_preserves_native_execution_in_compact_rollback_snapshot(self):
        import copy
        from biomolexplorer.project_state import read_snapshot
        rid,item=self.completed();original=item['configuration']
        with self.store.connect() as db:
            db.execute('UPDATE projects SET pipeline=? WHERE id=?',(json.dumps([original]),self.pid))
        project=self.store.project(self.owner,self.pid)
        md=new_stage('dynamics_dock6')
        saved=self.store.save_pipeline(self.owner,self.pid,[original,md],project['revision'])
        event=self.store.history(self.owner,self.pid)[0]['id']
        raw=json.loads((self.folder/'.history/snapshots'/f'{event}-after.json').read_text())
        self.assertTrue(raw['_compact_runs'])
        self.assertIsInstance(raw['runs'][0]['stages'],dict)
        before=read_snapshot(self.store,self.pid,event,'before')
        after=read_snapshot(self.store,self.pid,event,'after')
        self.assertEqual(before['project']['pipeline'],[original])
        self.assertEqual(after['project']['pipeline'],[original,md])
        from biomolexplorer.project_state import translate
        self.assertEqual(after['runs'][0]['stages'],translate([item],str(self.folder),'@project'))
        self.assertEqual(before['runs'],after['runs'])
        blob=next((self.folder/'.history/run-stages').glob('*.json'))
        blobs={p.name:p.stat().st_mtime_ns for p in blob.parent.glob('*.json') if p.name!='index.json'}
        updated=copy.deepcopy(md);updated['parameters']['amber_score_temperature']=301.
        self.store.save_pipeline(self.owner,self.pid,[original,updated],saved['revision'])
        self.assertEqual(blobs,{p.name:p.stat().st_mtime_ns for p in blob.parent.glob('*.json') if p.name!='index.json'})

    def test_other_block_saves_keep_existing_audit_format(self):
        rid,item=self.completed();original=item['configuration']
        with self.store.connect() as db:db.execute('UPDATE projects SET pipeline=? WHERE id=?',(json.dumps([original]),self.pid))
        project=self.store.project(self.owner,self.pid)
        updated=dict(original,name='Renamed consensus')
        self.store.save_pipeline(self.owner,self.pid,[updated],project['revision'])
        event=self.store.history(self.owner,self.pid)[0]['id']
        raw=json.loads((self.folder/'.history/snapshots'/f'{event}-after.json').read_text())
        self.assertNotIn('_compact_runs',raw)
        self.assertIsInstance(raw['runs'][0]['stages'],list)

    def test_pruning_md_results_hydrates_native_records_and_preserves_upstream_provenance(self):
        import copy
        from biomolexplorer.stage_cache import file_digest
        rid,item=self.completed();original=item['configuration'];md=new_stage('dynamics_dock6')
        from biomolexplorer.artifacts import records
        from biomolexplorer.artifact_registry import execution_manifest
        item['artifact_records']=records(item,rid)
        item['artifact_records'][0]['custom_provenance']={'preserve':'native'}
        item['execution_manifest']=execution_manifest(item)
        output=self.folder/'dynamics_results.json';output.write_text('{}')
        native=dict(id=md['id'],name=md['name'],operation='dynamics_dock6',configuration=copy.deepcopy(md),status='succeeded',
            artifacts=[str(output)],artifact_manifest={str(output):file_digest(output)})
        native['artifact_records']=records(native,rid)
        native['execution_manifest']=execution_manifest(native)
        with self.store.connect() as db:
            db.execute('UPDATE projects SET pipeline=? WHERE id=?',(json.dumps([original,md]),self.pid))
            db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps([item,native]),rid))
        project=self.store.project(self.owner,self.pid)
        updated=copy.deepcopy(md);updated['parameters']['amber_score_temperature']=301.
        proposal=self.store.result_removal_plan(self.owner,self.pid,[original,updated])
        self.assertGreater(proposal['affected'],0)
        self.store.save_pipeline(self.owner,self.pid,[original,updated],project['revision'],prune_results=True,confirmation=proposal['confirmation'])
        with self.store.connect() as db:
            stored=json.loads(db.execute('SELECT stages FROM runs WHERE id=?',(rid,)).fetchone()[0])
        self.assertEqual(stored[0],item)
        self.assertTrue(stored[1]['outputs_removed'])
        self.assertFalse(output.exists())

    def test_run_revision_change_falls_back_to_authoritative_native_plan(self):
        from biomolexplorer.result_catalog import read_proofs
        from biomolexplorer.result_sync import plan
        rid,item=self.completed();original=item['configuration']
        with self.store.connect() as db:db.execute('UPDATE projects SET pipeline=? WHERE id=?',(json.dumps([original]),self.pid))
        project=self.store.project(self.owner,self.pid);md=new_stage('dynamics_dock6')
        def changed(store,run_id):
            result=read_proofs(store,run_id)
            with store.connect() as db:
                latest=dict(item,concurrent_diagnostic='preserved')
                db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps([latest]),run_id))
            return result
        with patch('biomolexplorer.result_catalog.read_proofs',side_effect=changed),patch('biomolexplorer.result_sync.plan',wraps=plan) as plans:
            self.store.save_pipeline(self.owner,self.pid,[original,md],project['revision'])
        self.assertTrue(all(call.kwargs.get('history') is None for call in plans.call_args_list))
        with self.store.connect() as db:latest=json.loads(db.execute('SELECT stages FROM runs WHERE id=?',(rid,)).fetchone()[0])
        self.assertEqual(latest[0]['concurrent_diagnostic'],'preserved')
