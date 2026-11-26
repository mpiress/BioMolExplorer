"""Concurrent autosave, PDB dialogs and bounded canvas update work."""
import asyncio
import copy
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import AsyncMock, Mock, patch

import flet as ft
from biomolexplorer.catalog import new_stage
from biomolexplorer.ui.app import WorkspaceUI
from biomolexplorer.ui.flow_canvas import FlowCanvas
from biomolexplorer.ui.localization import Translator
from biomolexplorer.workspace import WorkspaceStore
from test_frontend_navigation import Page, walk


class PipelineEditorRegressionTests(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.store=WorkspaceStore(temporary.name)
        self.token=self.store.register('Editor','editor@example.org','editor-password')
        project=self.store.create_project(self.token,'Pipeline')
        self.pid=project['id']
        self.page=Page()
        self.ui=WorkspaceUI(self.page,self.store,SimpleNamespace(dock6_path=None))
        self.addCleanup(self.ui.structure_viewer.close)
        self.ui.token=self.token;self.ui.current=project
        self.stage=new_stage('retrieve_structures')
        self.ui.current['pipeline']=[self.stage];self.ui.dirty=True
        async def call(function,*args,**kwargs):return function(*args,**kwargs)
        self.ui.call=call

    def test_debounce_does_not_cancel_a_write_or_create_a_false_collaborator_conflict(self):
        async def scenario():
            entered=asyncio.Event();release=asyncio.Event();tasks=[];writes=0
            async def call(function,*args,**kwargs):
                nonlocal writes
                result=function(*args,**kwargs)
                if function==self.store.save_pipeline:
                    writes+=1
                    if writes==1:
                        entered.set();await release.wait()
                return result
            self.ui.call=call
            async def guard(action):await action()
            self.ui.guard=guard
            def run_task(handler,*args):
                task=asyncio.create_task(handler(*args));tasks.append(task);return task
            self.page.run_task=run_task
            self.ui.schedule_save()
            await asyncio.wait_for(entered.wait(),3)
            self.stage['position']={'x':500,'y':200}
            self.ui.schedule_save()
            await asyncio.sleep(0)
            self.assertFalse(tasks[0].cancelling())
            release.set()
            await tasks[0]
            # The first acknowledgement reschedules the latest draft.
            await asyncio.gather(*tasks[1:],return_exceptions=True)
            self.assertFalse(self.ui.dirty)
            self.assertEqual(self.ui.current['revision'],2)
            self.assertEqual(self.store.project(self.token,self.pid)['pipeline'][0]['position'],{'x':500,'y':200})
        asyncio.run(scenario())

    def test_pdb_dialog_opened_during_autosave_still_applies_after_acknowledgement(self):
        async def scenario():
            entered=asyncio.Event();release=asyncio.Event();first=True
            async def call(function,*args,**kwargs):
                nonlocal first
                result=function(*args,**kwargs)
                if function==self.store.save_pipeline and first:
                    first=False;entered.set();await release.wait()
                return result
            self.ui.call=call
            project=self.ui.current
            save=asyncio.create_task(self.ui.save_draft(silent=True))
            await entered.wait()
            await self.ui.open_stage_dialog(self.stage['id'])
            dialog=self.page.dialogs[-1]
            field=next(c for c in walk(dialog) if isinstance(c,ft.TextField) and c.label=='Nome da coleção PDB (opcional)')
            field.value='Minha coleção'
            release.set();await save
            self.assertIs(self.ui.current,project)
            await dialog.actions[-1].on_click(SimpleNamespace(control=dialog.actions[-1]))
            self.assertFalse(self.ui.editing_stage)
            saved=self.store.project(self.token,self.pid)
            self.assertEqual(saved['pipeline'][0]['parameters']['target'],'Minha coleção')
        asyncio.run(scenario())

    def test_empty_redocking_can_be_saved_while_pipeline_is_being_built(self):
        async def scenario():
            stage=new_stage('redocking');self.ui.current['pipeline']=[stage];self.ui.dirty=False
            await self.ui.open_stage_dialog(stage['id'])
            dialog=self.page.dialogs[-1]
            await dialog.actions[-1].on_click(SimpleNamespace(control=dialog.actions[-1]))
            self.assertFalse(dialog.open);self.assertFalse(self.ui.editing_stage)
            saved=self.store.project(self.token,self.pid)['pipeline'][0]
            self.assertEqual(saved['parameters'].get('pdb_codes'),[])
        asyncio.run(scenario())

    def test_remote_merge_and_edits_made_during_save_are_both_preserved(self):
        async def scenario():
            await self.ui.save_draft(silent=True)
            remote=copy.deepcopy(self.ui.current['pipeline'])
            remote[0]['name']='Shared PDB'
            self.store.save_pipeline(self.token,self.pid,remote,1)
            self.stage['parameters']['target']='Local collection';self.ui.dirty=True
            async def call(function,*args,**kwargs):
                result=function(*args,**kwargs)
                if function==self.store.save_pipeline:
                    self.stage['position']={'x':750,'y':400}
                return result
            self.ui.call=call
            await self.ui.save_draft(silent=True)
            self.ui.call=self._direct_call
            await self.ui.save_draft(silent=True)
            saved=self.store.project(self.token,self.pid)['pipeline'][0]
            self.assertEqual(saved['name'],'Shared PDB')
            self.assertEqual(saved['parameters']['target'],'Local collection')
            self.assertEqual(saved['position'],{'x':750,'y':400})
        asyncio.run(scenario())

    async def _direct_call(self,function,*args,**kwargs):return function(*args,**kwargs)

    def test_poll_does_not_replace_a_newer_save_with_an_older_response(self):
        async def scenario():
            await self.ui.save_draft(silent=True)
            stale=copy.deepcopy(self.ui.current)
            newer=copy.deepcopy(stale['pipeline'])
            newer[0]['position']={'x':800,'y':500}
            saved=self.store.save_pipeline(self.token,self.pid,newer,1)
            async def delayed_response(*args):
                self.ui.current.update(saved)
                self.ui.base_pipeline=copy.deepcopy(saved['pipeline'])
                return stale
            self.ui.call=delayed_response
            self.ui.draw_project=AsyncMock()
            with patch('biomolexplorer.ui.project_tools.asyncio.sleep',
                       side_effect=[None,asyncio.CancelledError()]):
                with self.assertRaises(asyncio.CancelledError):await self.ui.poll_project(self.pid)
            self.ui.draw_project.assert_not_awaited()
            self.assertEqual(self.ui.current['revision'],2)
            self.assertEqual(self.ui.base_pipeline,newer)
        asyncio.run(scenario())

    def test_translation_cache_is_bounded_and_cleared_when_language_changes(self):
        translate=Translator('en')
        self.assertEqual(translate('Pipeline salvo.'),'Pipeline saved.')
        self.assertEqual(translate('Pipeline salvo.'),'Pipeline saved.')
        self.assertGreater(translate._translate.cache_info().hits,0)
        self.assertEqual(translate._translate.cache_info().maxsize,4096)
        translate.set_language('pt')
        self.assertEqual(translate('Pipeline salvo.'),'Pipeline salvo.')

    def test_drag_preserves_controls_and_updates_only_node_and_edges(self):
        self.ui.page=SimpleNamespace(update=Mock())
        editor=FlowCanvas(self.ui,True)
        node=editor.nodes[self.stage['id']]
        options=editor.selection.options
        edges=editor.edges
        header=node.content.controls[0]
        self.assertEqual(header.drag_interval,32)
        header.on_pan_update(SimpleNamespace(local_delta=SimpleNamespace(x=10,y=20)))
        self.assertIs(editor.nodes[self.stage['id']],node)
        self.assertIs(editor.selection.options,options)
        self.assertIs(editor.edges,edges)
        self.ui.page.update.assert_called_once_with(node,edges)
        self.assertEqual((node.left,node.top),(70,84))
        asyncio.run(header.on_tap(None))
        self.assertIs(editor.nodes[self.stage['id']],node)
        self.assertEqual(editor.selection.value,self.stage['id'])


class InvalidationDecisionTests(unittest.TestCase):
    """An invalidation question must never freeze the buttons that ask it."""
    setUp = PipelineEditorRegressionTests.setUp

    def project_with_results(self):
        """A saved pipeline whose recorded results an edit would invalidate."""
        import json, time
        from uuid import uuid4
        from biomolexplorer.stage_cache import artifact_manifest
        retrieval = new_stage('retrieve_compounds'); admet = new_stage('admet')
        admet['bindings'] = {'base_input_path': {'stage': retrieval['id'], 'selector': 'auto'}}
        self.store.save_pipeline(self.token, self.pid, [retrieval, admet], 0)
        run_id = uuid4().hex
        folder = self.store.project_dir(self.pid)/'runs'/run_id/'artifacts'/'compounds'/'CHEMBL220'
        folder.mkdir(parents=True)
        table = folder/'compounds.csv'
        table.write_text('molecule_chembl_id,canonical_smiles\nCHEMBL1,CCO\n')
        stages = [dict(id=retrieval['id'], name=retrieval['name'], operation='retrieve_compounds',
                       configuration=retrieval, status='succeeded', artifacts=[str(table)],
                       artifact_manifest=artifact_manifest([str(table)]), finished_at=time.time())]
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',
                       (run_id, self.pid, self.store.user(self.token)['id'], 'succeeded',
                        json.dumps(stages), time.time(), time.time(), None))
        return retrieval

    async def edited(self):
        retrieval = self.project_with_results()
        await self.ui.open_project(self.pid)
        stage = next(s for s in self.ui.current['pipeline'] if s['id'] == retrieval['id'])
        stage['parameters']['target'] = 'OutroAlvo'
        self.ui.dirty = True
        self.assertTrue(self.store.result_removal_plan(self.token, self.pid, self.ui.current['pipeline'])['affected'])

    def toolbar(self, label):
        return next(control for control in walk(self.page.controls) if isinstance(control, ft.Button)
                    and (control.content.value if isinstance(control.content, ft.Text) else control.content) == label)

    def question(self):
        return next((dialog for dialog in self.page.dialogs if dialog.open and
                     getattr(getattr(dialog, 'title', None), 'value', '') == 'Confirmar alteração e invalidar resultados?'), None)

    def covering(self, dialog):
        return [other for other in self.page.dialogs if other is not dialog and getattr(other, 'open', False)
                and str(getattr(other, 'key', '') or '').startswith('processing')]

    async def answer(self, label=None):
        """Reach the decision the way a person does, or dismiss it."""
        for _ in range(500):
            dialog = self.question()
            if dialog is not None:
                covering = self.covering(dialog)
                if label is None:
                    dialog.on_dismiss(None); dialog.open = False
                else:
                    next(action for action in dialog.actions if
                         (action.content.value if isinstance(action.content, ft.Text) else action.content) == label).on_click(None)
                return covering
            await asyncio.sleep(0)
        raise AssertionError('A decisão nunca foi apresentada.')

    def test_the_decision_is_reachable_and_the_execution_continues(self):
        async def scenario():
            await self.edited()
            pressed = asyncio.ensure_future(self.answer('Confirmar remoção'))
            await asyncio.wait_for(self.toolbar('Executar pipeline').on_click(None), timeout=30)
            # A modal splash keeps its barrier above later dialogs: a question
            # under one could never be answered and the action waited forever.
            self.assertEqual(await pressed, [])
            self.assertFalse(self.ui.submitting)
            self.assertFalse(self.ui.dirty)
        asyncio.run(scenario())

    def test_dismissing_the_decision_keeps_every_other_button_alive(self):
        async def scenario():
            await self.edited()
            pressed = asyncio.ensure_future(self.answer())
            await asyncio.wait_for(self.toolbar('Executar pipeline').on_click(None), timeout=30)
            await pressed
            self.assertFalse(self.ui.submitting)
            # The change was restored, and the other actions still answer.
            await asyncio.wait_for(self.toolbar('Arquivos').on_click(None), timeout=30)
            self.assertEqual(self.ui.tab, 'Arquivos')
            self.ui.tab = 'Pipeline'
            await self.ui.draw_project()
            await asyncio.wait_for(self.toolbar('Executar pipeline').on_click(None), timeout=30)
            self.assertFalse(self.ui.submitting)
        asyncio.run(scenario())

    def test_a_pending_decision_never_blocks_another_action(self):
        async def scenario():
            await self.edited()
            click = asyncio.ensure_future(self.toolbar('Executar pipeline').on_click(None))
            for _ in range(500):
                if self.question() is not None: break
                await asyncio.sleep(0)
            self.assertIsNotNone(self.question())
            # While the question waits, every other action must still answer.
            self.assertIs(await asyncio.wait_for(self.ui.save_draft(silent=True), timeout=10), False)
            # A second execution declines the duplicate instead of waiting.
            await asyncio.wait_for(self.toolbar('Executar pipeline').on_click(None), timeout=30)
            await asyncio.wait_for(self.toolbar('Arquivos').on_click(None), timeout=30)
            self.assertEqual(self.ui.tab, 'Arquivos')
            self.assertIsNotNone(self.question())
            await self.answer('Cancelar')
            await asyncio.wait_for(click, timeout=30)
            self.assertFalse(self.ui.submitting)
        asyncio.run(scenario())

    def test_an_automatic_draft_never_asks_and_keeps_the_change_pending(self):
        async def scenario():
            await self.edited()
            opened = [dialog for dialog in self.page.dialogs if dialog.open]
            self.assertIs(await asyncio.wait_for(self.ui.save_draft(silent=True, confirm=False), timeout=10), False)
            self.assertTrue(self.ui.dirty)
            self.assertIsNone(self.question())
            self.assertEqual([dialog for dialog in self.page.dialogs if dialog.open], opened)
        asyncio.run(scenario())


if __name__=='__main__':unittest.main()
