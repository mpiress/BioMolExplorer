"""Immediate, shared, exception-safe feedback for foreground processing."""
import asyncio
import inspect
import threading
import unittest
from types import SimpleNamespace
from unittest.mock import AsyncMock,Mock

import flet as ft
from biomolexplorer.ui.app import WorkspaceUI
from biomolexplorer.ui.pdb_results import PDBActions
from biomolexplorer.ui.processing import processing,background,save_file
from test_frontend_navigation import Page,walk


class ProcessingFeedbackTests(unittest.IsolatedAsyncioTestCase):
    def ui(self):return SimpleNamespace(page=Page())

    async def test_nested_reads_use_one_splash_until_entire_action_finishes(self):
        ui=self.ui()
        def read_file():
            self.assertTrue(ui._configuration_loading.dialog.open)
            return 'data'
        async with processing(ui,immediate=False):
            self.assertEqual(ui.page.dialogs,[])
            self.assertEqual(await WorkspaceUI.call(ui,read_file),'data')
            dialog=ui._configuration_loading.dialog
            async with processing(ui,'Validando arquivos…'):
                self.assertEqual(await WorkspaceUI.call(ui,read_file),'data')
                self.assertIs(ui._configuration_loading.dialog,dialog)
            self.assertTrue(dialog.open)
            self.assertEqual(len(ui.page.dialogs),1)
        self.assertFalse(dialog.open)
        self.assertEqual(ui._configuration_loading.requests,set())

    async def test_concurrent_actions_and_cancellation_do_not_hide_other_work(self):
        ui=self.ui();entered=asyncio.Event();release=asyncio.Event()
        async def first():
            async with processing(ui,'Importando projeto…'):
                entered.set();await release.wait()
        pending=asyncio.create_task(first());await asyncio.wait_for(entered.wait(),1)
        async with processing(ui,'Salvando pipeline…'):
            self.assertEqual(len(ui._configuration_loading.requests),2)
            pending.cancel()
            with self.assertRaises(asyncio.CancelledError):await pending
            self.assertTrue(ui._configuration_loading.dialog.open)
            self.assertEqual(ui._configuration_loading.title.value,'Salvando pipeline…')
        self.assertFalse(ui._configuration_loading.requests)

    async def test_background_poll_does_not_show_or_change_foreground_feedback(self):
        ui=self.ui()
        @background
        async def poll():return await WorkspaceUI.call(ui,lambda:'updated')
        self.assertEqual(await poll(),'updated')
        self.assertFalse(hasattr(ui,'_configuration_loading'))
        async with processing(ui,'Salvando pipeline…'):
            dialog=ui._configuration_loading.dialog
            await asyncio.create_task(poll())
            self.assertIs(ui._configuration_loading.dialog,dialog)
            self.assertEqual(ui._configuration_loading.title.value,'Salvando pipeline…')
            self.assertEqual(len(ui._configuration_loading.requests),1)

    async def test_native_file_dialog_pauses_splash_and_next_write_resumes_it(self):
        ui=self.ui()
        async def choose(**kwargs):
            self.assertFalse(ui._configuration_loading.requests)
            return '/tmp/result.csv'
        ui.picker=SimpleNamespace(save_file=AsyncMock(side_effect=choose))
        async with processing(ui,'Preparando exportação…'):
            first=ui._configuration_loading.dialog;card=first.content
            self.assertEqual(await save_file(ui,file_name='result.csv'),'/tmp/result.csv')
            self.assertFalse(first.open)
            await WorkspaceUI.call(ui,lambda:None)
            second=ui._configuration_loading.dialog
            self.assertTrue(second.open)
            self.assertIsNot(second.content,card)
        self.assertFalse(second.open)

    async def test_slow_analysis_keeps_cancel_available_on_the_standard_splash(self):
        ui=self.ui();stop=threading.Event();entered=asyncio.Event();loop=asyncio.get_running_loop()
        def common(rows,timeout,limit,cancelled):
            loop.call_soon_threadsafe(entered.set)
            return cancelled.wait(2)
        pending=asyncio.create_task(WorkspaceUI.call(ui,common,[],10,100,stop))
        await asyncio.wait_for(entered.wait(),2)
        dialog=ui._configuration_loading.dialog
        self.assertTrue(dialog.actions)
        dialog.actions[0].on_click(None)
        self.assertTrue(await asyncio.wait_for(pending,2))
        self.assertFalse(dialog.open)


class PDBSavingFeedbackTests(unittest.IsolatedAsyncioTestCase):
    async def asyncSetUp(self):
        self.errors=[];self.page=Page()
        async def guard(action):
            try:await action()
            except Exception as error:self.errors.append(error)
        async def call(function,*args,**kwargs):
            result=function(*args,**kwargs)
            return await result if inspect.isawaitable(result) else result
        self.ui=SimpleNamespace(page=self.page,call=call,guard=guard,notify=Mock(),editing_stage=True)
        self.context=dict(pdb_id='1XYZ',revision='version',resolution=2.,ligands=[],
            available_ligands=[dict(LIGAND='ATP',RESNUM='102',CHAIN='B',RECOMMENDED=False)],recommended_ligands=[])
        self.saved=dict(self.context,reviewed=True,ligands=[dict(LIGAND='ATP',RESNUM='102',CHAIN='B')])
        self.service=SimpleNamespace(pdb_ligands=Mock(return_value=self.context),set_pdb_ligands=AsyncMock(return_value=self.saved))
        result=SimpleNamespace(ui=self.ui,service=self.service,token='token',project_id='project',run_id='run',stage={'id':'stage'},valid=lambda:True)
        self.callback=Mock();self.actions=PDBActions(result,True,on_saved=self.callback)
        self.actions.refresh_pipeline=AsyncMock()
        await self.actions.open_ligands({'name':'1XYZ.pdb','path':'/project/1XYZ.pdb'})
        self.editor=self.page.dialogs[-1]
        controls=self.editor.content.content.controls
        controls[2].value='0';controls[3].on_click(None)
        self.save=self.editor.actions[-1].on_click

    async def test_manual_ligand_save_shows_splash_before_draft_and_keeps_it_until_refresh(self):
        entered=asyncio.Event();release=asyncio.Event();refresh=asyncio.Event();finish=asyncio.Event()
        async def draft(**kwargs):
            splash=self.ui._configuration_loading
            self.assertTrue(splash.dialog.open)
            self.assertEqual(splash.title.value,'Salvando ligantes PDB…')
        self.ui.save_draft=AsyncMock(side_effect=draft)
        async def commit(*args):entered.set();await release.wait();return self.saved
        async def rebuild():refresh.set();await finish.wait()
        self.service.set_pdb_ligands.side_effect=commit
        self.actions.refresh_pipeline.side_effect=rebuild
        pending=asyncio.create_task(self.save(None))
        await asyncio.wait_for(entered.wait(),1)
        dialog=self.ui._configuration_loading.dialog
        self.assertIsNot(dialog,self.editor)
        self.assertTrue(dialog.modal)
        self.assertTrue(any(isinstance(c,ft.ProgressRing) and c.value is None for c in walk(dialog)))
        self.assertTrue(any(isinstance(c,ft.Image) for c in walk(dialog)))
        self.assertTrue(all(c.disabled for c in self.editor.actions))
        await self.save(None)
        self.service.set_pdb_ligands.assert_awaited_once()
        self.assertEqual(self.service.set_pdb_ligands.call_args.args[-2],[dict(LIGAND='ATP',RESNUM='102',CHAIN='B')])
        release.set();await asyncio.wait_for(refresh.wait(),1)
        self.assertTrue(dialog.open)
        self.callback.assert_not_called()
        finish.set();await pending
        self.assertFalse(dialog.open);self.assertFalse(self.editor.open)
        self.callback.assert_called_once()
        await self.save(None)
        self.service.set_pdb_ligands.assert_awaited_once()

    async def test_save_failure_cleans_splash_preserves_fields_and_allows_retry(self):
        self.service.set_pdb_ligands.side_effect=[ValueError('invalid ligand'),self.saved]
        await self.save(None)
        self.assertEqual(len(self.errors),1)
        self.assertTrue(self.editor.open)
        self.assertFalse(self.ui._configuration_loading.requests)
        self.assertTrue(all(not c.disabled for c in self.editor.actions))
        row=self.editor.content.content.controls[1].controls[0]
        self.assertEqual([c.value for c in row.controls[:3]],['ATP','102','B'])
        self.callback.assert_not_called()
        await self.save(None)
        self.assertFalse(self.editor.open)
        self.assertFalse(self.ui._configuration_loading.requests)
        self.callback.assert_called_once()


if __name__=='__main__':unittest.main()
