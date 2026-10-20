"""Slow removals remain visible and accept only one confirmation."""
import asyncio
import unittest
from types import SimpleNamespace
from unittest.mock import AsyncMock, Mock, patch

import flet as ft
from biomolexplorer.ui.app import WorkspaceUI
from biomolexplorer.ui.file_table import FileTable
from biomolexplorer.ui.removal import RemovalProgress


class RemovalProgressTests(unittest.IsolatedAsyncioTestCase):
    def ui(self):
        async def guard(action):
            try:
                await action()
            except Exception as error:
                errors.append(error)
        async def call(function, *args):
            return await function(*args)
        dialogs, errors = [], []
        return SimpleNamespace(token='session', current={'id':'project'},
            page=SimpleNamespace(update=Mock(), show_dialog=dialogs.append),
            guard=guard, call=call, notify=Mock(), dialogs=dialogs, errors=errors)

    async def exercise_slow_removal(self, dialog, remove, entered, release, mutation):
        pending = asyncio.create_task(remove(None))
        await entered.wait()
        content = dialog.content.content if isinstance(dialog.content, ft.Container) else dialog.content
        indicator = content.controls[-1]
        self.assertTrue(indicator.visible)
        self.assertIsInstance(indicator.controls[0], ft.ProgressBar)
        self.assertIsNone(indicator.controls[0].value)
        self.assertTrue(all(control.disabled for control in dialog.actions))
        await remove(None)
        self.assertEqual(mutation.await_count, 1)
        release.set()
        await pending
        await remove(None)  # A queued click after completion must also be ignored.
        self.assertEqual(mutation.await_count, 1)
        self.assertFalse(indicator.visible)

    async def test_bulk_stage_removal_is_visible_and_executes_once(self):
        ui=self.ui();ui.store=object();ui.draw_project=AsyncMock()
        entered, release=asyncio.Event(), asyncio.Event()
        async def slow(*args):
            entered.set();await release.wait()
        mutation=AsyncMock(side_effect=slow)
        with patch('biomolexplorer.result_files.ResultFiles') as service:
            service.return_value.remove_all=mutation
            await WorkspaceUI.remove_stage_results(ui,'project','run','stage')
            dialog=ui.dialogs[-1]
            await self.exercise_slow_removal(dialog, dialog.actions[-1].on_click, entered, release, mutation)
        ui.draw_project.assert_awaited_once()
        self.assertFalse(ui.errors)

    async def test_file_removal_uses_same_single_flight_feedback(self):
        ui=self.ui();entered,release=asyncio.Event(),asyncio.Event()
        async def slow(item):
            entered.set();await release.wait()
        mutation=AsyncMock(side_effect=slow)
        item={'name':'results.csv','path':'results.csv'}
        table=FileTable(ui,'project',[item],writable=True,remove=mutation)
        await table.confirm_remove(item)
        dialog=ui.dialogs[-1]
        await self.exercise_slow_removal(dialog,dialog.actions[-1].on_click,entered,release,mutation)
        self.assertEqual(table.files,[])

    async def test_failure_restores_buttons_and_allows_retry(self):
        ui=self.ui()
        dialog=ft.AlertDialog(content=ft.Text('Confirmar'), actions=[ft.TextButton('Cancelar'),ft.TextButton('Remover')])
        progress=RemovalProgress(ui,dialog)
        mutation=AsyncMock(side_effect=[OSError('disk failure'),None])
        await progress.run(mutation)
        self.assertEqual(len(ui.errors),1)
        self.assertFalse(progress.indicator.visible)
        self.assertTrue(all(not button.disabled for button in dialog.actions))
        await progress.run(mutation)
        await progress.run(mutation)
        self.assertEqual(mutation.await_count,2)

    async def test_failure_preserves_previously_disabled_controls(self):
        ui=self.ui()
        dialog=ft.AlertDialog(content=ft.Column([]),actions=[ft.TextButton('Cancelar'),ft.TextButton('Remover',disabled=True)])
        progress=RemovalProgress(ui,dialog)
        await progress.run(AsyncMock(side_effect=OSError('failure')))
        self.assertFalse(dialog.actions[0].disabled)
        self.assertTrue(dialog.actions[1].disabled)
