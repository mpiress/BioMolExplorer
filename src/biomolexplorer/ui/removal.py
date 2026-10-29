"""Visible, single-flight feedback for destructive confirmation dialogs."""
import asyncio
import flet as ft


class RemovalProgress:
    def __init__(self, ui, dialog):
        self.ui, self.dialog = ui, dialog
        self.running = False
        self.completed = False
        self.indicator = ft.Column([
            ft.ProgressBar(),
            ft.Text('Removendo dados… Aguarde a conclusão.', size=13),
        ], visible=False, tight=True, spacing=8)
        container = dialog.content if isinstance(dialog.content, ft.Container) else dialog
        content = container.content
        if isinstance(content, ft.Column):
            content.controls.append(self.indicator)
        else:
            container.content = ft.Column([content, self.indicator], tight=True, spacing=16)

    async def run(self, action):
        # Acquire before yielding: even already queued click events cannot start twice.
        if self.running or self.completed:
            return
        self.running = True
        disabled = [control.disabled for control in self.dialog.actions]
        for control in self.dialog.actions:
            control.disabled = True
        self.indicator.visible = True
        self.ui.page.update()

        async def execute():
            # Flush the busy state before any deletion or expensive UI rebuild.
            await asyncio.sleep(0)
            result=await action()
            self.completed = result is not False

        try:
            await self.ui.guard(execute)
        finally:
            self.running = False
            self.indicator.visible = False
            if not self.completed:
                for control, previous in zip(self.dialog.actions, disabled):
                    control.disabled = previous
            self.ui.page.update()
