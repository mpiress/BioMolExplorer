"""Select, validate and remove typed files directly from the import block."""
import flet as ft
from biomolexplorer.input_validation import contract,validate_file
from biomolexplorer.import_inputs import file_kind
from .localization import verbatim


class ImportFiles:
    def __init__(self,form):
        from .app import KINDS
        self.form,self.ui=form,form.ui
        self.kinds=KINDS;self.entries={};self.type_controls={}
        params=form.stage['parameters']
        assets={a['id']:a for a in form.assets}
        for identifier in params.get('asset_ids',[]):
            asset=assets.get(identifier)
            self.entries[identifier]=dict(asset or {'id':identifier,'name':'Arquivo indisponível'},kind=file_kind(params,identifier))
        self.kind=ft.Dropdown(label='Tipo dos arquivos a adicionar',value=params.get('kind','compounds'),
            options=[ft.DropdownOption(key=k,text=v) for k,v in KINDS.items()],disabled=not form.writable)
        self.expected=ft.Text(contract(self.kind.value),size=12,color='#64748B')
        self.kind.on_select=lambda e:self.sync_contract()
        self.available=ft.Dropdown(label='Arquivo já enviado ao projeto',options=[],disabled=not form.writable,expand=True)
        self.table=ft.DataTable(columns=[ft.DataColumn(ft.Text(label)) for label in ('Tipo','Arquivo','Ações')],
            rows=[],column_spacing=24,data_row_min_height=72,data_row_max_height=88)
        self.count=ft.Text(size=12,color='#64748B')
        self.control=ft.Column([self.kind,self.expected,
            ft.TextButton('Selecionar arquivos no disco',icon=ft.Icons.UPLOAD_FILE,on_click=self.upload,disabled=not form.writable),
            ft.Row([self.available,ft.TextButton('Adicionar arquivo do projeto',on_click=self.add_available,disabled=not form.writable)],wrap=True),
            ft.TextButton('Adicionar todos os arquivos do projeto',on_click=self.add_all,disabled=not form.writable),
            ft.Row([self.table],scroll=ft.ScrollMode.AUTO),self.count],spacing=20,data='input')
        self.redraw()

    def sync_contract(self):
        self.expected.value=contract(self.kind.value);self.ui.page.update()

    def redraw(self):
        self.type_controls={};rows=[]
        for identifier,asset in self.entries.items():
            control=ft.Dropdown(value=asset['kind'],options=[ft.DropdownOption(key=k,text=v) for k,v in self.kinds.items()],
                disabled=not self.form.writable,width=285)
            async def change(e,identifier=identifier,control=control):
                previous=self.entries[identifier]['kind']
                async def action():
                    try:
                        await self.validate(identifier,control.value)
                    except Exception:
                        control.value=previous;self.ui.page.update();raise
                    self.entries[identifier]['kind']=control.value
                await self.ui.guard(action)
            control.on_select=change;self.type_controls[identifier]=control
            remove=ft.IconButton(icon=ft.Icons.REMOVE_CIRCLE_OUTLINE,tooltip='Remover da lista',disabled=not self.form.writable,
                on_click=lambda e,identifier=identifier:self.remove(identifier))
            rows.append(ft.DataRow(cells=[ft.DataCell(control),
                ft.DataCell(ft.Container(verbatim(ft.Text(asset['name'],tooltip=asset['name'])),width=360)),ft.DataCell(remove)]))
        self.table.rows=rows
        self.available.options=[verbatim(ft.DropdownOption(key=a['id'],text=a['name'])) for a in self.form.assets if a['id'] not in self.entries]
        if self.available.value not in {o.key for o in self.available.options}:self.available.value=None
        self.count.value=f'{len(self.entries)} arquivos selecionados'

    async def validate(self,identifier,kind):
        token,project_id=self.ui.token,self.ui.current['id']
        path=await self.ui.call(self.ui.store.asset_path,token,project_id,identifier)
        await self.ui.call(validate_file,path,kind)

    def remove(self,identifier):
        if not self.form.writable:return
        self.entries.pop(identifier,None);self.redraw();self.ui.page.update()

    async def add(self,identifiers,types=None):
        if not self.form.writable:return
        pending={}
        for identifier in identifiers:
            if identifier in self.entries:continue
            asset=next(a for a in self.form.assets if a['id']==identifier)
            kind=(types or {}).get(identifier,asset['kind'])
            await self.validate(identifier,kind)
            pending[identifier]=dict(asset,kind=kind)
        self.entries.update(pending)
        self.redraw();self.ui.page.update()

    async def upload(self,e):
        if not self.form.writable:return
        async def action():
            token,project_id=self.ui.token,self.ui.current['id'];kind=self.kind.value
            ids=await self.ui.pick_uploads(kind,project_id)
            if self.ui.token!=token or not self.ui.current or self.ui.current['id']!=project_id:return
            self.form.assets[:]=await self.ui.call(self.ui.store.assets,token,project_id)
            await self.add(ids)
        await self.ui.guard(action)

    async def add_available(self,e):
        async def action():
            if self.available.value:await self.add([self.available.value])
        await self.ui.guard(action)

    async def add_all(self,e):
        await self.ui.guard(lambda:self.add([a['id'] for a in self.form.assets]))

    def primary_kind(self):
        kinds={a['kind'] for a in self.entries.values()}
        return next(iter(kinds)) if len(kinds)==1 else 'other' if kinds else self.kind.value

    def asset_types(self):return {identifier:asset['kind'] for identifier,asset in self.entries.items()}
