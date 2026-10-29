"""Visible, single-flight feedback for destructive confirmation dialogs."""
from .processing import processing


class RemovalProgress:
    def __init__(self, ui, dialog,title='Removendo dados…'):
        self.ui, self.dialog = ui, dialog
        self.title=title
        self.running = False
        self.completed = False

    async def run(self, action):
        # Acquire before yielding: even already queued click events cannot start twice.
        if self.running or self.completed:
            return
        self.running = True
        disabled = [control.disabled for control in self.dialog.actions]
        for control in self.dialog.actions:
            control.disabled = True
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
            if not self.completed:
                for control, previous in zip(self.dialog.actions, disabled):
                    control.disabled = previous
            self.ui.page.update()
