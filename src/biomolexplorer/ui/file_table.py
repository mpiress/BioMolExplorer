"""Reusable paginated file lists with actions aligned on the right."""
from pathlib import Path
import flet as ft
from .processing import save_file
from .removal import RemovalProgress
from .feedback import close_dialog
from .localization import verbatim


class FileTable:
    def __init__(self, ui, project_id, files, writable=False, remove=None, download=None, preview=None, extra_actions=None):
        self.ui,self.project_id,self.files=ui,project_id,files
        self.token=ui.token;self.active=True;self.offset=0;self.limit=25
        self.remove_action,self.download_action,self.preview_action=remove,download,preview
        self.extra_actions=extra_actions
        self.writable=writable
        self.table=ft.DataTable(columns=[ft.DataColumn(ft.Text(n)) for n in ('Arquivo','Tamanho','Ações')],rows=[],
            column_spacing=30,heading_row_color='#F1F5F9',data_row_min_height=72,data_row_max_height=92)
        origins=sorted({f.get('origin','') for f in files if f.get('origin')})
        self._origin_values=origins
        self.origin=ft.Dropdown(label='Origem dos resultados',value='',visible=len(origins)>1,
            options=[ft.DropdownOption(key='',text='Todas as origens')]+[verbatim(ft.DropdownOption(key=o,text=o),'text') for o in origins],
            enable_filter=True,enable_search=True,on_select=self.filter_origin)
        self.count=ft.Text(size=13,color='#64748B')
        self.previous=ft.TextButton('Página anterior',on_click=lambda e:self.move(-1))
        self.next=ft.TextButton('Próxima página',on_click=lambda e:self.move(1))
        self.page_size=ft.Dropdown(label='Itens por página',value='25',width=160,
            options=[ft.DropdownOption(key=str(n),text=str(n)) for n in (10,25,50,100)],on_select=self.resize)
        self.redraw()

    def valid(self):
        return self.active and self.ui.token==self.token and self.ui.current and self.ui.current['id']==self.project_id

    def filter_origin(self,e):
        self.offset=0;self.redraw();self.ui.page.update()

    def redraw(self):
        origins=sorted({f.get('origin','') for f in self.files if f.get('origin')})
        if origins!=self._origin_values:
            self.origin.options=[ft.DropdownOption(key='',text='Todas as origens')]+[
                verbatim(ft.DropdownOption(key=o,text=o),'text') for o in origins]
            self._origin_values=origins
        self.origin.visible=len(origins)>1
        if self.origin.value not in origins:self.origin.value=''
        visible=[f for f in self.files if not self.origin.value or f.get('origin')==self.origin.value]
        if self.offset>=len(visible):self.offset=max(0,((len(visible)-1)//self.limit)*self.limit)
        rows=[]
        for item in visible[self.offset:self.offset+self.limit]:
            async def download(e,item=item):
                async def action():
                    if not self.valid():return
                    if self.download_action:return await self.download_action(item)
                    data=await self.ui.call(self.ui.store.read_file,self.token,self.project_id,item['path'])
                    if self.valid():await save_file(self.ui,file_name=item['name'],src_bytes=data)
                await self.ui.guard(action)
            actions=[ft.IconButton(ft.Icons.DOWNLOAD,tooltip='Baixar '+item['name'],on_click=download)]
            if self.preview_action and (str(item.get('path','')).endswith('.biomol-view.json') or Path(item['name']).suffix.lower() in ('.png','.jpg','.jpeg')):
                async def preview(e,item=item):
                    if self.valid():await self.ui.guard(lambda:self.preview_action(item))
                actions.insert(0,ft.IconButton(ft.Icons.INSIGHTS,tooltip='Visualizar '+item['name'],on_click=preview))
            if self.remove_action:
                async def remove(e,item=item):
                    if self.valid():await self.confirm_remove(item)
                actions.append(ft.IconButton(ft.Icons.DELETE_OUTLINE,tooltip='Remover '+item['name'],on_click=remove,
                    disabled=not self.writable,icon_color='#B91C1C'))
            if self.extra_actions:actions.extend(self.extra_actions(item))
            name=ft.Column([verbatim(ft.Text(item.get('display_name',item['name']),size=13,max_lines=2,
                    overflow=ft.TextOverflow.ELLIPSIS,tooltip=item.get('provenance',item['name']))),
                *([verbatim(ft.Text(item['description'],size=11,color='#64748B',max_lines=1,overflow=ft.TextOverflow.ELLIPSIS,tooltip=item.get('provenance',item['description'])),'value')] if item.get('description') else [])],spacing=3)
            rows.append(ft.DataRow(cells=[ft.DataCell(ft.Container(name,width=420)),
                ft.DataCell(ft.Text(f"{item.get('size',0)/1024:.1f} KB",size=12)),ft.DataCell(ft.Row(actions,spacing=2))]))
        self.table.rows=rows
        self.count.value=f'{len(visible)} arquivos · Página {self.offset//self.limit+1}'
        self.previous.disabled=self.offset==0;self.next.disabled=self.offset+self.limit>=len(visible)

    def move(self,direction):
        if not self.valid():return
        self.offset=max(0,self.offset+direction*self.limit);self.redraw();self.ui.page.update()

    def resize(self,e):
        if not self.valid():return
        self.limit=int(self.page_size.value);self.offset=0;self.redraw();self.ui.page.update()

    async def confirm_remove(self,item):
        async def remove(e):
            async def action():
                if not self.valid():return
                await self.remove_action(item)
                if not self.valid():return
                self.files=[f for f in self.files if f is not item]
                self.redraw();close_dialog(self.ui.page,dialog);self.ui.page.update()
                self.ui.notify('Arquivo removido. A alteração foi registrada no histórico do projeto.')
            await progress.run(action)
        dialog=ft.AlertDialog(modal=True,title=ft.Text('Remover arquivo?',size=22),
            content=ft.Column([verbatim(ft.Text(item['name'])),ft.Text('O arquivo deixará de estar disponível. A próxima execução recalculará a etapa se necessário. O histórico permite restaurar o estado anterior.')],tight=True,spacing=16),
            actions=[ft.TextButton('Manter arquivo',on_click=lambda e:close_dialog(self.ui.page,dialog)),ft.TextButton('Remover arquivo',on_click=remove)])
        progress=RemovalProgress(self.ui,dialog,'Excluindo arquivo…')
        self.ui.page.show_dialog(dialog)

    def build(self):
        return ft.Column([self.origin,ft.Row([self.count,self.page_size],wrap=True,alignment=ft.MainAxisAlignment.SPACE_BETWEEN),
            ft.Row([self.table],scroll=ft.ScrollMode.AUTO),
            ft.Row([self.previous,self.next],wrap=True,alignment=ft.MainAxisAlignment.SPACE_BETWEEN)],spacing=16)
