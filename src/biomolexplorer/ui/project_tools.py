"""Workspace portability, version history and collaboration controls."""
import asyncio
import copy
import shutil
from datetime import datetime
from pathlib import Path
from uuid import uuid4
import flet as ft
from .localization import verbatim

from .feedback import close_dialog
from .folder_browser import FolderBrowser


class ProjectTools:
    def folder_controls(self, value=None, locked=False):
        field=ft.TextField(label='Pasta do projeto',value=value or '',
            hint_text='Selecione uma pasta',read_only=True,multiline=True,min_lines=1,max_lines=3)
        def selected(path):
            field.value=str(path);field.error=None;self.page.update()
        async def choose(e):
            if locked:return
            token=self.token
            async def action():
                if self.page.web:
                    await FolderBrowser(self,selected,field.value).open()
                else:
                    parent=await self.picker.get_directory_path(dialog_title=self.tr('Escolha a pasta do projeto'))
                    if parent and self.token==token:selected(parent)
            await self.guard(action)
        field.suffix_icon=ft.Semantics(container=True,label='Escolher pasta',button=True,
            exclude_semantics=True,disabled=locked,on_tap=None if locked else choose,
            content=ft.IconButton(icon=ft.Icons.FOLDER_OPEN,tooltip='Escolher pasta',
                on_click=choose,disabled=locked))
        return field,ft.Column([field,
            ft.Text('Escolha uma pasta nova ou vazia. Configuração, entradas, resultados e histórico ficarão nela.',size=12,color='#64748B'),
            ft.Text('No navegador, as pastas são do computador onde o BioMolExplorer está em execução.',size=12,color='#64748B') if self.page.web else
            ft.Container()],spacing=12,horizontal_alignment=ft.CrossAxisAlignment.STRETCH)

    async def export_project_dialog(self, project_id):
        token=self.token
        archive=await self.call(self.store.export_project,token,project_id)
        if self.token!=token:return
        await self.call(self.store.project,token,project_id)
        if self.page.web:
            if archive.stat().st_size>self.store.max_upload_bytes:
                self.notify('O pacote foi gerado em '+str(archive)+'. Para pacotes grandes, copie-o diretamente da pasta do projeto.')
                return
            data=await self.call(archive.read_bytes)
            if self.token!=token:return
            await self.picker.save_file(file_name=archive.name,src_bytes=data)
        else:
            destination=await self.picker.save_file(file_name=archive.name,dialog_title=self.tr('Exportar projeto'))
            if destination:await self.call(shutil.copyfile,archive,destination)
        self.notify('Projeto exportado com configuração, entradas, resultados e histórico. Contas e senhas ficam neste workspace.')

    async def pick_uploads(self, kind, project_id=None, archive=False):
        files=await self.picker.pick_files(allow_multiple=not archive)
        results=[]
        token=self.token
        for file in files or []:
            actual_kind='visualization' if kind!='other' and (file.name.lower().endswith(('.png','.jpg','.jpeg','.biomol-view.json'))) else kind
            if file.path and not self.page.web:
                results.append(Path(file.path) if archive else await self.call(self.store.import_local_file,token,project_id,file.path,actual_kind))
                continue
            if file.name in self.uploads:raise ValueError('Aguarde o envio anterior de '+file.name+'.')
            if file.size and file.size>self.store.max_upload_bytes:raise ValueError('Arquivo acima do limite de envio de 200 MB. Use a aplicação desktop para arquivos maiores.')
            ticket=uuid4().hex if archive else await self.call(self.store.prepare_upload,token,project_id,file.name,actual_kind)
            future=asyncio.get_running_loop().create_future()
            self.uploads[file.name]={'ticket':ticket,'project_id':project_id,'future':future,'archive':archive,'token':token}
            await self.picker.upload([ft.FilePickerUploadFile(name=file.name,id=file.id,upload_url=self.page.get_upload_url(ticket,600))])
            try:
                results.append(await asyncio.wait_for(future,timeout=600))
            finally:
                self.uploads.pop(file.name,None)
        return results

    async def import_project_dialog(self):
        token=self.token
        directory,folder=self.folder_controls()
        selected=ft.Text('Nenhum pacote selecionado',size=12)
        path=None
        async def pick(e):
            async def action():
                nonlocal path
                if self.token!=token:return
                selected.value='Recebendo e validando o pacote…';self.page.update()
                files=await self.pick_uploads('other',archive=True)
                if self.token!=token:return
                if files:path=files[0];selected.value='Pacote selecionado: '+Path(path).name;self.page.update()
            await self.guard(action)
        async def run(e):
            async def action():
                if self.token!=token:return
                if not path:raise ValueError('Selecione um projeto .bme.zip exportado pela plataforma.')
                if not directory.value or not directory.value.strip():raise ValueError('Informe a pasta de destino.')
                run_button.disabled=True;run_button.content='Importando…';self.page.update()
                try:
                    project=await self.call(self.store.import_project,token,path,directory.value.strip())
                finally:
                    run_button.disabled=False;run_button.content='Importar';self.page.update()
                if self.token!=token:return
                close_dialog(self.page,dialog)
                if self.page.web:Path(path).unlink(missing_ok=True)
                await self.open_project(project['id'])
                self.notify('Projeto importado. Os resultados concluídos serão reaproveitados quando seus dados e parâmetros forem compatíveis.')
            await self.guard(action)
        run_button=ft.Button('Importar',on_click=run)
        dialog=ft.AlertDialog(modal=True,title=ft.Text('Importar projeto',size=24),
            content=ft.Container(width=560,content=ft.Column([ft.Text('Migre seus experimentos sem repetir as etapas concluídas.'),
                ft.TextButton('Selecionar .bme.zip',icon=ft.Icons.UPLOAD_FILE,on_click=pick),selected,folder],spacing=20,tight=True)),
            actions=[ft.TextButton('Cancelar',on_click=lambda e:close_dialog(self.page,dialog)),run_button])
        self.page.show_dialog(dialog)

    async def history_dialog(self,project_id):
        token=self.token
        project=await self.call(self.store.project,token,project_id)
        events=await self.call(self.store.history,token,project_id)
        rows=[]
        async def confirm(event):
            async def restore(e):
                async def action():
                    restored=await self.call(self.store.rollback,token,project_id,event['id'])
                    close_dialog(self.page,confirmation);close_dialog(self.page,dialog)
                    if self.current and self.current['id']==project_id:await self.open_project(project_id)
                    else:await self.show_workspace()
                    self.notify('Versão restaurada. O histórico mantém o registro desta restauração.')
                await self.guard(action)
            confirmation=ft.AlertDialog(modal=True,title=ft.Text('Restaurar versão anterior?',size=24),
                content=ft.Container(width=540,content=ft.Column([ft.Text(event['summary']),
                    ft.Text('O projeto voltará ao estado anterior a esta alteração, incluindo configuração, arquivos, resultados e permissões. Mudanças posteriores sairão da versão atual e continuarão preservadas no histórico.')],spacing=18,tight=True)),
                actions=[ft.TextButton('Cancelar',on_click=lambda e:close_dialog(self.page,confirmation)),ft.Button('Restaurar versão',on_click=restore)])
            self.page.show_dialog(confirmation)
        for event in events:
            async def rollback(e,event=event):await self.guard(lambda:confirm(event))
            rows.append(ft.DataRow(cells=[ft.DataCell(ft.Text(datetime.fromtimestamp(event['created']).astimezone().strftime('%d/%m/%Y %H:%M:%S %Z'),size=12)),
                ft.DataCell(ft.Column([verbatim(ft.Text(event['actor_name'],size=13)),verbatim(ft.Text(event['actor_email'],size=11))],spacing=2,tight=True)),
                ft.DataCell(ft.Container(width=300,content=ft.Text(event['summary'],size=13))),
                ft.DataCell(ft.TextButton('Restaurar antes',on_click=rollback,disabled=not event['has_before'] or project['role']!='owner'))]))
        table=ft.DataTable(columns=[ft.DataColumn(ft.Text(v)) for v in ('Data e hora','Usuário','Alteração','Versão anterior')],rows=rows,column_spacing=28,data_row_min_height=58,data_row_max_height=90)
        dialog=ft.AlertDialog(title=ft.Text('Histórico · '+project['name'],size=24),content=ft.Container(width=1040,height=500,
            content=ft.Column([ft.Text('Mais recente primeiro. Restaurar antes desfaz a alteração selecionada e todas as posteriores.',size=12),
                ft.Row([table],scroll=ft.ScrollMode.AUTO)] if rows else [ft.Text('Nenhuma alteração registrada ainda.')],scroll=ft.ScrollMode.AUTO,spacing=16)),
            actions=[ft.TextButton('Fechar',on_click=lambda e:close_dialog(self.page,dialog))])
        self.page.show_dialog(dialog)

    def schedule_save(self):
        task=getattr(self,'_save_task',None)
        if task and not task.done():task.cancel()
        async def deferred():
            await asyncio.sleep(.6)
            if self.current and self.dirty and not getattr(self,'editing_stage',False):
                await self.guard(lambda:self.save_draft(silent=True))
        self._save_task=self.page.run_task(deferred)

    async def poll_project(self,project_id):
        token=self.token
        while self.token==token and self.current and self.current['id']==project_id:
            await asyncio.sleep(2)
            try:
                fresh=await self.call(self.store.project,token,project_id)
            except Exception:
                await self.guard(lambda:self.draw_project());return
            if self.token!=token or not self.current or self.current['id']!=project_id:return
            if fresh['role']!=self.current['role']:
                self.close_project_views();await self.open_project(project_id);return
            if not self.dirty and not getattr(self,'editing_stage',False) and fresh['updated']!=self.current['updated']:
                self.current.update(fresh);self.base_pipeline=copy.deepcopy(fresh['pipeline'])
                await self.draw_project()
