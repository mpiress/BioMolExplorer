"""Nonblocking update notifications shared by desktop and web sessions."""
import asyncio
import flet as ft
from .feedback import close_dialog
from .processing import background
from .localization import verbatim
from biomolexplorer.workspace import AccessDenied


class UpdateNotice:
    def __init__(self, ui, checker):
        self.ui, self.checker = ui, checker
        self.control = None
        self.task = None
        self.downloading = False

    def user(self):
        try:
            return self.ui.store.user(self.ui.token)['id'] if self.ui.token else 'installation'
        except AccessDenied:
            return 'installation'

    def button(self):
        self.control = ft.TextButton('Atualizações', icon=ft.Icons.SYSTEM_UPDATE_ALT,
            on_click=self.open)
        self.refresh()
        return self.control

    def refresh(self):
        if self.control:
            available=self.checker.visible(self.user())
            label=self.ui.tr('Atualização disponível' if available else 'Atualizações')
            compact=(getattr(self.ui.page,'width',None) or 1440)<650
            self.control.content = '' if compact else label
            self.control.icon = ft.Icons.NEW_RELEASES_OUTLINED if available else ft.Icons.SYSTEM_UPDATE_ALT
            self.control.tooltip = label

    @background
    async def poll(self):
        while True:
            await self.ui.call(self.checker.check)
            self.refresh()
            self.ui.page.update()
            await asyncio.sleep(60)

    async def open(self, event):
        latest = await self.ui.call(self.checker.check, True)
        self.refresh()
        if latest is None:
            self.ui.notify('Não foi possível consultar atualizações. Tente novamente mais tarde.')
            return
        if not latest['available']:
            self.ui.notify('Nenhuma atualização nova na branch master.')
            return
        async def later(e):
            await self.ui.call(self.checker.remind_later, self.user())
            close_dialog(self.ui.page, dialog)
            self.refresh(); self.ui.page.update()
        async def download(e):
            if self.downloading:
                return
            self.downloading = True
            button=e.control if e is not None else None
            if button is not None:
                button.disabled=True
                button.content=self.ui.tr('Instalando atualização...')
                button.update()
            try:
                await self.ui.call(self.checker.install,latest)
            except Exception as exc:
                self.ui.notify('Não foi possível instalar a atualização: '+str(exc))
                if button is not None:
                    button.disabled=False
                    button.content=self.ui.tr('Instalar atualização')
                    button.update()
            else:
                await later(e)
                self.ui.notify('Atualização instalada com segurança. Reinicie o BioMolExplorer para usar a nova versão.')
            finally:
                self.downloading=False
        dialog = ft.AlertDialog(modal=False, title=ft.Text('Atualização disponível'),
            content=ft.Column([verbatim(ft.Text(latest['message'])), verbatim(ft.Text(latest['sha'][:12])),
                ft.Text('A atualização será validada e instalada automaticamente no ambiente atual. Reinicie o aplicativo após a conclusão.')],
                tight=True, width=480),
            actions=[ft.TextButton('Lembrar em 24 horas', on_click=later),
                ft.TextButton('Instalar atualização', on_click=download)])
        self.ui.page.show_dialog(dialog)
