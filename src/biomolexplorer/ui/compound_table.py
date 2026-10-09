"""Paginated retrieval results, local molecule previews and dataset curation."""
from pathlib import Path
import flet as ft
from biomolexplorer.compound_tables import CompoundTables
from biomolexplorer.visualizations import molecule_image
from .feedback import close_dialog
from .localization import verbatim


class CompoundTableViewer:
    def __init__(self,ui,project_id,run_id,stage_id,tables,inline=False):
        self.ui,self.project_id,self.run_id,self.stage_id=ui,project_id,run_id,stage_id
        self.token=ui.token
        self.inline=inline;self.active=True
        self.service=CompoundTables(ui.store)
        self.tables=tables;self.offset=0;self.limit=25;self.sequence=0;self.preview_version=0;self.data=None
        self.selector=ft.Dropdown(label='Tabela de compostos',value=tables[0]['path'],width=450,
            options=[ft.DropdownOption(key=t['path'],text=t['name']+(' · conjunto integrado' if t['integrated'] else '')) for t in tables],on_select=self.change_table)
        self.search=ft.TextField(label='Buscar composto ou SMILES',expand=True,on_submit=self.search_rows)
        self.count=ft.Text(size=13,color='#64748B')
        self.notice=ft.Text(size=12,color='#64748B')
        self.table=ft.DataTable(columns=[ft.DataColumn(ft.Text(name)) for name in ('Código do composto','SMILES','Visualização','Ações')],
            rows=[],data_row_min_height=68,data_row_max_height=80,column_spacing=24,heading_row_color='#F1F5F9')
        async def previous(e):await self.move(-1)
        async def next_page(e):await self.move(1)
        self.previous=ft.TextButton('Página anterior',on_click=previous)
        self.next=ft.TextButton('Próxima página',on_click=next_page)
        self.loading=ft.ProgressBar(visible=False)
        self.page_size=ft.Dropdown(label='Itens por página',value='25',width=160,
            options=[ft.DropdownOption(key=str(n),text=str(n)) for n in (10,25,50,100)],on_select=self.change_size)

    def valid(self):
        return self.active and self.ui.token==self.token and self.ui.current is not None and self.ui.current['id']==self.project_id and (self.inline or getattr(self.ui,'compound_viewer',None) is self)

    async def change_size(self,e):
        if not self.valid():return
        self.limit=int(self.page_size.value);self.offset=0
        await self.ui.guard(self.load)

    async def load(self):
        self.sequence+=1;sequence=self.sequence
        self.loading.visible=True
        self.ui.page.update()
        try:
            result=await self.ui.call(self.service.page,self.token,self.project_id,self.run_id,self.stage_id,
                self.selector.value,self.offset,self.limit,self.search.value or '')
        except Exception:
            if self.valid() and sequence==self.sequence:
                self.loading.visible=False
                self.ui.page.update()
            raise
        if not self.valid() or sequence!=self.sequence:
            return
        if self.offset>0 and not result['rows'] and result['matched']:
            self.offset=((result['matched']-1)//self.limit)*self.limit
            return await self.load()
        self.data=result;self.loading.visible=False
        label=self.ui.tr('composto encontrado' if result['matched']==1 else 'compostos encontrados') if hasattr(self.ui,'tr') else ('composto encontrado' if result['matched']==1 else 'compostos encontrados')
        self.count.value=f"{result['name']} · {result['matched']} {label} · {result['total']} na tabela · Página {self.offset//self.limit+1}"
        self.notice.value='A remoção afeta somente o arquivo selecionado. As análises que usam esse arquivo serão atualizadas na próxima execução.' if result['can_edit'] else 'Consulta disponível. Para remover compostos, é necessário ser editor e aguardar o pipeline terminar.'
        rows=[]
        for row in result['rows']:
            async def show_2d(e,row=row):await self.ui.guard(lambda:self.preview(row,False))
            async def show_3d(e,row=row):await self.ui.guard(lambda:self.preview(row,True))
            async def remove(e,row=row):await self.ui.guard(lambda:self.confirm_remove(row))
            smiles=verbatim(ft.Text(row['smiles'] or '—',size=12,max_lines=2,overflow=ft.TextOverflow.ELLIPSIS,selectable=True,tooltip=row['smiles']))
            identity=[verbatim(ft.Text(row['id'] or '—',size=13))]
            identity.extend(ft.IconButton(ft.Icons.OPEN_IN_NEW,tooltip='Abrir no '+link['provider'],
                url=ft.Url(link['url'],target=ft.UrlTarget.BLANK)) for link in row.get('links',[]))
            rows.append(ft.DataRow(cells=[ft.DataCell(ft.Row(identity,spacing=4)),
                ft.DataCell(ft.Container(verbatim(ft.Semantics(label=row['smiles'],exclude_semantics=True,content=smiles),'label'),width=300)),
                ft.DataCell(ft.Row([ft.TextButton('2D',on_click=show_2d,disabled=not bool(row['smiles'])),
                                    ft.TextButton('3D',on_click=show_3d,disabled=not bool(row['smiles']))],spacing=0)),
                ft.DataCell(ft.TextButton('Remover',icon=ft.Icons.DELETE_OUTLINE,on_click=remove,disabled=not result['can_edit'],style=ft.ButtonStyle(color='#B91C1C')))]))
        self.table.rows=rows
        self.previous.disabled=self.offset==0
        self.next.disabled=self.offset+self.limit>=result['matched']
        self.ui.page.update()

    async def change_table(self,e):
        if not self.valid():return
        self.offset=0;self.preview_version+=1
        await self.ui.guard(self.load)

    async def search_rows(self,e):
        if not self.valid():return
        self.offset=0
        self.preview_version+=1
        await self.ui.guard(self.load)

    async def move(self,direction):
        if not self.valid():return
        self.offset=max(0,self.offset+direction*self.limit)
        await self.ui.guard(self.load)

    async def preview(self,row,three_d):
        if not self.valid():return
        self.preview_version+=1;version=self.preview_version
        await self.ui.call(self.ui.store.project,self.token,self.project_id)
        if three_d:
            url=await self.ui.call(self.ui.compound_view_url,self.project_id,row['smiles'],row['id'],self.token)
            if not self.valid() or version!=self.preview_version:return
            await self.ui.call(self.ui.store.project,self.token,self.project_id)
            if not self.valid() or version!=self.preview_version:return
            await ft.UrlLauncher().launch_url(url,mode=ft.LaunchMode.EXTERNAL_APPLICATION,web_only_window_name='_blank')
            return
        model=await self.ui.call(molecule_image,row['smiles'])
        if not self.valid() or version!=self.preview_version:return
        await self.ui.call(self.ui.store.project,self.token,self.project_id)
        if not self.valid() or version!=self.preview_version:return
        width=max(280,min(620,(self.ui.page.width or 1440)-140))
        if model:
            content=ft.Image(src=model,width=width,height=300,fit=ft.BoxFit.CONTAIN)
        else:
            raise ValueError('A estrutura 2D está indisponível para este SMILES.')
        dialog=ft.AlertDialog(title=ft.Text(row['id']+' · '+('3D' if three_d else '2D'),size=22),
            content=ft.Container(width=width,content=ft.Column([content,
                ft.TextField(label='SMILES',value=row['smiles'],read_only=True,multiline=True,max_lines=3)],tight=True,spacing=14)),
            actions=[ft.TextButton('Fechar molécula',on_click=lambda e:close_dialog(self.ui.page,dialog))])
        self.ui.page.show_dialog(dialog)

    async def confirm_remove(self,row):
        if not self.valid():return
        self.preview_version+=1
        table=self.selector.value;expected_version=self.data['version']
        async def remove(e):
            async def action():
                if not self.valid():return
                await self.ui.call(self.service.remove,self.token,self.project_id,self.run_id,self.stage_id,table,row['index'],expected_version)
                if not self.valid():return
                close_dialog(self.ui.page,dialog)
                await self.load()
                self.ui.notify('Composto removido de '+self.data['name']+'.')
            await self.ui.guard(action)
        dialog=ft.AlertDialog(modal=True,title=ft.Text('Remover este composto?',size=22),
            content=ft.Column([verbatim(ft.Text(row['id'],weight=ft.FontWeight.W_600)),
                ft.Text('A remoção será aplicada apenas a '+self.data['name']+'. O registro original será preservado no histórico de alterações.')],tight=True,spacing=16),
            actions=[ft.TextButton('Manter composto',on_click=lambda e:close_dialog(self.ui.page,dialog)),ft.TextButton('Remover composto',on_click=remove)])
        self.ui.page.show_dialog(dialog)

    async def download(self,e):
        async def action():
            if not self.valid():return
            path,filename=self.selector.value,Path(self.selector.value).name
            data=await self.ui.call(self.ui.store.read_file,self.token,self.project_id,path)
            if self.valid():await self.ui.picker.save_file(file_name=filename,src_bytes=data)
        await self.ui.guard(action)

    def build(self):
        return ft.Column([ft.Row([self.selector,ft.TextButton('Baixar tabela',icon=ft.Icons.DOWNLOAD,on_click=self.download)],wrap=True),
            ft.Row([self.search,ft.TextButton('Buscar',on_click=self.search_rows)]),self.loading,
            ft.Row([self.count,self.page_size],wrap=True,alignment=ft.MainAxisAlignment.SPACE_BETWEEN),
            ft.Container(ft.Row([self.table],scroll=ft.ScrollMode.AUTO)),
            ft.Row([self.previous,self.next],wrap=True,alignment=ft.MainAxisAlignment.SPACE_BETWEEN),self.notice],spacing=14,expand=not self.inline)
