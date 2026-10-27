"""Paginated per-input chemical alerts with local molecular visualization."""
import flet as ft
from .feedback import close_dialog
from .localization import verbatim
from biomolexplorer.visualizations import molecule_image


class ADMETAlertTable:
    def __init__(self,results,filename):
        self.results=results;self.ui=results.ui;self.filename=filename
        self.active=True;self.offset=0;self.sequence=0
        self.search=ft.TextField(label='Pesquisar molécula ou alerta',expand=True,on_submit=self.search_changed)
        self.summary=ft.Text(size=12);self.table=ft.DataTable(columns=[ft.DataColumn(ft.Text('Molécula'))])
        self.previous=ft.TextButton('Anterior',on_click=self.previous_page)
        self.next=ft.TextButton('Próxima',on_click=self.next_page)
        self.source=ft.Text(size=12)
        self.dialog=ft.AlertDialog(modal=True,title=ft.Text('Alertas para revisão química'),
            content=ft.Container(width=max(260,min(1100,(self.ui.page.width or 1280)-100)),
                height=max(220,min(650,(getattr(self.ui.page,'height',None) or 900)-220)),content=ft.Column([
                    self.source,ft.Row([self.search,ft.IconButton(ft.Icons.SEARCH,on_click=self.search_changed)]),
                    ft.Text('true indica correspondência ao catálogo; false indica ausência de correspondência. Alertas não comprovam toxicidade.',size=12),
                    ft.Row([self.summary,self.previous,self.next],wrap=True),
                    ft.Column([ft.Row([self.table],scroll=ft.ScrollMode.AUTO)],scroll=ft.ScrollMode.AUTO,expand=True)],spacing=16)),
            actions=[ft.TextButton('Fechar alertas',on_click=self.close)],on_dismiss=self.dismissed)

    def valid(self):return self.active and self.results.valid()
    def dismissed(self,e):self.active=False;self.sequence+=1
    def close(self,e=None):self.dismissed(e);close_dialog(self.ui.page,self.dialog)

    async def open(self):
        await self.load()
        if self.valid():self.ui.page.show_dialog(self.dialog)

    async def search_changed(self,e):self.offset=0;await self.ui.guard(self.load)
    async def previous_page(self,e):await self.move(-1)
    async def next_page(self,e):await self.move(1)
    async def move(self,direction):self.offset=max(0,self.offset+direction*25);await self.ui.guard(self.load)

    async def load(self):
        self.sequence+=1;sequence=self.sequence;r=self.results
        data=await self.ui.call(r.service.admet_alert_page,r.token,r.project_id,r.run_id,r.stage['id'],
            self.filename,self.offset,25,self.search.value or '')
        if not self.valid() or sequence!=self.sequence:return
        self.source.value=data['label'];self.summary.value=f"{data['total']} moléculas · Página {self.offset//25+1}"
        self.table.columns=[ft.DataColumn(ft.Text('Código da molécula')),*[
            ft.DataColumn(ft.Text(name)) for name in data['catalogs']],ft.DataColumn(ft.Text('Alertas identificados')),ft.DataColumn(ft.Text('Estrutura'))]
        rows=[]
        for record in data['rows']:
            async def preview(e,record=record):
                async def action():
                    await self.ui.call(r.service.store.project,r.token,r.project_id)
                    if not self.valid():return
                    image=await self.ui.call(molecule_image,record['smiles'])
                    if not self.valid():return
                    if not image:raise ValueError('A estrutura 2D está indisponível para este SMILES.')
                    await self.ui.call(r.service.store.project,r.token,r.project_id)
                    if not self.valid():return
                    dialog=ft.AlertDialog(title=verbatim(ft.Text(record['id'])),
                        content=ft.Image(src=image,width=min(560,(self.ui.page.width or 1280)-100),height=320,fit=ft.BoxFit.CONTAIN),
                        actions=[ft.TextButton('Fechar molécula',on_click=lambda e:close_dialog(self.ui.page,dialog))])
                    self.ui.page.show_dialog(dialog)
                await self.ui.guard(action)
            rows.append(ft.DataRow(cells=[ft.DataCell(verbatim(ft.Text(record['id']))),*[
                ft.DataCell(verbatim(ft.Text('true' if record['flags'][name] else 'false',color='#B45309' if record['flags'][name] else '#15803D')))
                for name in data['catalogs']],ft.DataCell(verbatim(ft.Text(record['details'] or '—',width=340,max_lines=3,tooltip=record['details']))),
                ft.DataCell(ft.Semantics(label='Visualizar molécula',button=True,exclude_semantics=True,on_tap=preview,
                    content=ft.IconButton(ft.Icons.VISIBILITY_OUTLINED,tooltip='Visualizar molécula',on_click=preview,disabled=not record['smiles'])))]))
        self.table.rows=rows;self.previous.disabled=self.offset==0;self.next.disabled=self.offset+25>=data['total']
        self.ui.page.update()
