"""Apply must not reset a curated plan when auto bindings become explicit files."""
import asyncio
import copy
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import AsyncMock,Mock
import unittest
from biomolexplorer.catalog import new_stage
from biomolexplorer.redocking_config import validate_pairs
from biomolexplorer.ui.file_selection import FileSelection
import test_pipeline_editor_regressions as fixtures

class RedockingApplyTests(unittest.TestCase):
    setUp=fixtures.PipelineEditorRegressionTests.setUp

    def prepare(self):
        self.ui.pdb_view_url=Mock(return_value=None)
        source=self.stage;source['parameters']['pdb_ids']='1ABC,2XYZ'
        stage=new_stage('redocking');stage['bindings']['base_input_path']={'stage':source['id']}
        self.ui.current['pipeline']=[source,stage]
        root=self.store.project_dir(self.pid)/'retrieved';root.mkdir()
        paths=[]
        for code in ('1ABC','2XYZ'):
            path=root/(code+'.pdb');path.write_text('END\n');paths.append(str(path))
        metadata=root/'pdb_codes.csv';metadata.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RESOLUTION\n1ABC,ATP,102,A,1.7\n2XYZ,LIG,23,B,2.0\n')
        item=dict(source,configuration=copy.deepcopy(source),status='succeeded',artifacts=paths+[str(metadata)])
        pending=dict(stage,configuration=copy.deepcopy(stage),status='awaiting_input',artifacts=[])
        self.run=dict(id='pending-run',project_id=self.pid,status='awaiting_input',stages=[item,pending])
        self.store.list_runs=Mock(return_value=[self.run]);self.store.get_run=Mock(return_value=self.run)
        form=FileSelection(self.run,pending,page=self.ui.page)
        form.pickers['base_input_path'].add_all_sources()
        selected=form.read();selected['parameters']['exhaustiveness']=17
        selected['parameters']['sizeof_box']=[29,29,29]
        self.ui.draw_project=AsyncMock()
        return stage,selected

    def test_apply_and_continue_preserves_pairs_and_preparation_at_the_service_boundary(self):
        async def scenario():
            stage,selected=self.prepare()
            def resume(token,run_id,configuration):
                validate_pairs(configuration['parameters'].get('pdb_codes'),configuration['parameters'].get('preparation_pairs') or {})
                return dict(id=self.run['id'],status='queued')
            self.ui.service.resume=Mock(side_effect=resume)
            await self.ui.open_stage_dialog(stage['id'],waiting_run_id=self.run['id'],configuration=selected)
            dialog=self.page.dialogs[-1]
            await dialog.actions[-1].on_click(SimpleNamespace(control=dialog.actions[-1]))
            self.ui.service.resume.assert_called_once()
            sent=self.ui.service.resume.call_args.args[-1]
            self.assertEqual(sent['parameters']['pdb_codes'],selected['parameters']['pdb_codes'])
            self.assertEqual(sent['parameters']['exhaustiveness'],17)
            self.assertEqual(sent['parameters']['sizeof_box'],[29,29,29])
            self.assertEqual(len(sent['parameters']['preparation_pairs']),2)
            self.assertFalse(dialog.open)
            saved=self.store.project(self.token,self.pid)['pipeline'][-1]
            self.assertEqual(saved['parameters'],sent['parameters'])
        asyncio.run(scenario())

    def test_apply_configuration_without_a_pending_run_persists_the_curated_plan(self):
        async def scenario():
            stage,selected=self.prepare()
            await self.ui.open_stage_dialog(stage['id'],configuration=selected)
            dialog=self.page.dialogs[-1]
            await dialog.actions[-1].on_click(SimpleNamespace(control=dialog.actions[-1]))
            saved=self.store.project(self.token,self.pid)['pipeline'][-1]
            self.assertEqual(saved['parameters']['pdb_codes'],selected['parameters']['pdb_codes'])
            self.assertEqual(saved['parameters']['exhaustiveness'],17)
            self.assertFalse(dialog.open)
        asyncio.run(scenario())

    def test_finishing_redocking_closes_all_flow_routes_and_invalidates_old_buttons(self):
        import flet as ft
        from unittest.mock import patch
        async def scenario():
            stage,_=self.prepare()
            self.ui.service.resume=Mock(return_value={'id':self.run['id'],'status':'queued'})
            created=[]
            def capture(*args,**kwargs):
                form=FileSelection(*args,**kwargs);created.append(form);return form
            with patch('biomolexplorer.ui.file_selection.FileSelection',side_effect=capture):
                await self.ui.configure_pending(self.run['id'])
                triage=self.ui.selection_dialog;old_configure=triage.actions[1].on_click
                created[0].pickers['base_input_path'].add_all_sources()
                await triage.actions[-1].on_click(SimpleNamespace(control=triage.actions[-1]))
                settings=self.ui.settings_dialog
                # An older information route can still be mounted during its
                # dismiss animation. Completion must target every route.
                triage.open=True
                extra=ft.AlertDialog(actions=[ft.TextButton('Configurar')]);self.page.show_dialog(extra)
                self.ui._track_configuration_dialog(extra,self.run['id'],stage['id'])
                await settings.actions[-1].on_click(SimpleNamespace(control=settings.actions[-1]))
                self.assertTrue(all(not d.open for d in (triage,extra,settings)))
                self.assertIsNone(self.ui.selection_dialog)
                self.assertIsNone(self.ui.settings_dialog)
                self.assertNotIn((self.run['id'],stage['id']),self.ui.pending_file_selections)
                count=len(self.page.dialogs)
                await old_configure(SimpleNamespace(control=triage.actions[1]))
                self.assertEqual(len(self.page.dialogs),count)
                self.ui.service.resume.assert_called_once()
        asyncio.run(scenario())

    def test_finishing_clears_progress_reopen_request_without_closing_unrelated_dialog(self):
        import flet as ft
        stage,_=self.prepare()
        progress=SimpleNamespace(run={'id':self.run['id']},is_open=True,closing=False,
            reopen_requested=True,dialog=ft.AlertDialog())
        unrelated=ft.AlertDialog();self.page.show_dialog(unrelated);self.page.show_dialog(progress.dialog)
        self.ui.run_progress=progress
        self.ui._finish_configuration_dialogs(self.run['id'],stage['id'])
        self.assertFalse(progress.dialog.open)
        self.assertFalse(progress.is_open)
        self.assertFalse(progress.reopen_requested)
        self.assertTrue(unrelated.open)

    def test_failed_confirmation_keeps_settings_for_retry(self):
        async def scenario():
            stage,selected=self.prepare()
            self.ui.service.resume=Mock(side_effect=ValueError('Configuração inválida'))
            await self.ui.open_stage_dialog(stage['id'],waiting_run_id=self.run['id'],configuration=selected)
            settings=self.ui.settings_dialog
            await settings.actions[-1].on_click(SimpleNamespace(control=settings.actions[-1]))
            self.assertTrue(settings.open)
            self.assertIs(self.ui.settings_dialog,settings)
            self.assertIn((self.run['id'],stage['id']),self.ui._configuration_dialogs)
            self.ui.service.resume.side_effect=None
            self.ui.service.resume.return_value={'id':self.run['id'],'status':'queued'}
            await settings.actions[-1].on_click(SimpleNamespace(control=settings.actions[-1]))
            self.assertFalse(settings.open)
            self.assertIsNone(self.ui.settings_dialog)
        asyncio.run(scenario())

    def test_finishing_one_block_keeps_another_blocks_settings_open(self):
        import flet as ft
        stage,_=self.prepare()
        own=ft.AlertDialog();other=ft.AlertDialog()
        self.page.show_dialog(own);self.page.show_dialog(other)
        self.ui._track_configuration_dialog(own,None,stage['id'])
        self.ui._track_configuration_dialog(other,None,'another-block')
        self.ui.settings_dialog=other;self.ui.settings_run_id=None
        self.ui._finish_configuration_dialogs(None,stage['id'])
        self.assertFalse(own.open)
        self.assertTrue(other.open)
        self.assertIs(self.ui.settings_dialog,other)
        self.assertTrue(self.ui.editing_stage)
