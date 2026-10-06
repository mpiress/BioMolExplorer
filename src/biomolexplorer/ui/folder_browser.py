"""Visible-directory navigation and creation on the application host."""
import os
import inspect
import stat
from pathlib import Path

import flet as ft

from .feedback import close_dialog
from .localization import verbatim


class FolderBrowser:
    def __init__(self,ui,on_select,initial=None,allow_nonempty=False):
        self.ui,self.on_select,self.token=ui,on_select,ui.token
        self.current=Path(initial).expanduser() if initial else Path.home()
        self.path=verbatim(ft.Text('',size=12,selectable=True))
        self.entries=ft.Column(scroll=ft.ScrollMode.AUTO,expand=True,spacing=4)
        self.allow_nonempty=allow_nonempty
        self.selected=False
        self.name=ft.TextField(label='Nome da nova pasta',expand=True)
        self.error=ft.Text('',color='#B91C1C',size=12)
        self.up=ft.TextButton('Pasta anterior',icon=ft.Icons.ARROW_UPWARD,
            on_click=self.ui.event(self.parent))
        roots=[Path(drive) for drive in os.listdrives()] if hasattr(os,'listdrives') else [Path('/')]
        locations=ft.Row([ft.TextButton('Pasta pessoal',icon=ft.Icons.HOME_OUTLINED,
            on_click=self.ui.event(self.navigate,Path.home())),
            *[ft.TextButton(str(root),icon=ft.Icons.STORAGE,on_click=self.ui.event(self.navigate,root)) for root in roots]],wrap=True)
        self.dialog=ft.AlertDialog(modal=True,title=ft.Text('Escolha a pasta do projeto',size=24),
            content=ft.Container(width=560,height=420,content=ft.Column([
                ft.Text('No navegador, as pastas são do computador onde o BioMolExplorer está em execução.',size=12,color='#64748B'),
                locations,ft.Row([self.up,ft.Container(self.path,expand=True)]),self.entries,
                ft.TextButton('Criar pasta',icon=ft.Icons.CREATE_NEW_FOLDER_OUTLINED,
                    on_click=self.ui.event(self.prompt_create)),self.error],spacing=12)),
            actions=[ft.TextButton('Cancelar',on_click=lambda e:close_dialog(self.ui.page,self.dialog)),
                ft.Button('Selecionar esta pasta',icon=ft.Icons.CHECK,on_click=self.ui.event(self.select))])

    async def active(self):
        if self.ui.token!=self.token:return False
        await self.ui.call(self.ui.store.user,self.token)
        return self.ui.token==self.token

    @staticmethod
    def directories(path):
        path=path.resolve(strict=True)
        return path,sorted((p for p in path.iterdir() if not p.name.startswith('.') and p.is_dir()
                           and not (getattr(p.stat(),'st_file_attributes',0)&getattr(stat,'FILE_ATTRIBUTE_HIDDEN',0))
                           and not (getattr(p.stat(),'st_flags',0)&getattr(stat,'UF_HIDDEN',0))),key=lambda p:p.name.casefold())

    async def open(self):
        if not await self.active():return
        await self.navigate(self.current)
        if self.ui.token==self.token:self.ui.page.show_dialog(self.dialog)

    async def navigate(self,path):
        if not await self.active():return
        try:
            current,folders=await self.ui.call(self.directories,Path(path))
        except OSError:
            self.error.value='Não foi possível acessar esta pasta. Verifique as permissões.'
            self.ui.page.update();return
        if self.ui.token!=self.token:return
        self.current=current;self.path.value=str(current);self.error.value=''
        self.up.disabled=current.parent==current
        self.entries.controls=[ft.TextButton(content=verbatim(ft.Text(folder.name,max_lines=1,
            overflow=ft.TextOverflow.ELLIPSIS)),icon=ft.Icons.FOLDER_OUTLINED,
            tooltip=str(folder),on_click=self.ui.event(self.navigate,folder)) for folder in folders]
        if not folders:self.entries.controls=[ft.Text('Esta pasta está vazia.',size=12,color='#64748B')]
        self.ui.page.update()

    async def parent(self):
        await self.navigate(self.current.parent)

    async def prompt_create(self):
        if not await self.active():return
        self.name.value='';self.error.value=''
        parent=self.current
        async def submit(e):
            if self.current!=parent:return
            await self.create()
            if self.selected:close_dialog(self.ui.page,dialog)
            else:
                notice.value=self.error.value;self.ui.page.update()
        notice=ft.Text('',color='#B91C1C',size=12)
        dialog=ft.AlertDialog(modal=True,title=ft.Text('Criar pasta'),
            content=ft.Column([verbatim(ft.Text(str(parent),size=12)),self.name,notice],tight=True),
            actions=[ft.TextButton('Cancelar',on_click=lambda e:close_dialog(self.ui.page,dialog)),
                     ft.Button('Criar e selecionar',on_click=submit)])
        self.ui.page.show_dialog(dialog)

    async def create(self):
        if not await self.active():return
        name=(self.name.value or '').strip()
        if not name or name in ('.','..') or any(c in name for c in '/\\\x00') or Path(name).anchor:
            self.error.value='Informe um nome de pasta válido, sem separadores de caminho.'
            self.ui.page.update();return
        try:
            destination=self.current/name
            await self.ui.call(destination.mkdir,mode=0o700)
        except OSError:
            self.error.value='Não foi possível criar a pasta. Verifique o nome e as permissões.'
            self.ui.page.update();return
        self.name.value=''
        if self.ui.token!=self.token:return
        await self.navigate(destination)
        await self.select()

    async def select(self):
        if self.selected or not await self.active():return
        try:
            empty=await self.ui.call(lambda:self.current.is_dir() and not any(self.current.iterdir()))
        except OSError:
            empty=False
        if self.ui.token!=self.token:return
        if not empty and not self.allow_nonempty:
            self.error.value='Escolha uma pasta nova ou vazia para evitar sobrescrever arquivos existentes.'
            self.ui.page.update();return
        self.selected=True
        try:
            result=self.on_select(str(self.current))
            if inspect.isawaitable(result):await result
        except Exception:
            self.selected=False
            raise
        close_dialog(self.ui.page,self.dialog)
