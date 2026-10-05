"""Multiple typed sources and completed-result uploads in guided block forms."""
import flet as ft
from biomolexplorer.bindings import sources, pack
from biomolexplorer.flow import INPUTS, EXTERNAL_INPUTS, compatible, input_types
from biomolexplorer.input_validation import contract, output_kind
from .localization import verbatim


class InputEditor:
    def __init__(self,ui,stage,field,assets,writable):
        self.ui,self.stage,self.field,self.assets,self.writable=ui,stage,field,assets,writable
        self.rows=[]
        self.on_change=None
        self.box=ft.Column(spacing=24)
        self.value=stage['parameters'].get(field['name'])
        self.expected=input_types(stage,field['name']) or EXTERNAL_INPUTS.get(stage['operation'],{}).get(field['name'],{'other'})
        group=stage.get('bindings',{}).get(field['name'])
        for reference in sources(group) if group else [{}]:self.add(reference)
        if self.value and not group:self.rows[0][0].value='configured-path'
        self.control=ft.Column([ft.Text(field['label'],size=16,weight=ft.FontWeight.W_600),
            self.box,ft.Row([ft.TextButton('Adicionar entrada',icon=ft.Icons.ADD,on_click=self.append,disabled=not writable),
                ft.TextButton('Enviar meus arquivos',icon=ft.Icons.UPLOAD_FILE,on_click=self.upload,disabled=not writable)],wrap=True),
            ft.Text('Use a opção de processamento do bloco para manter os arquivos separados ou mesclar as entradas.',size=12,color='#64748B'),
            ft.Text('\n'.join(contract(k,'graphs' if stage['operation']=='graphs' else None)
                for k in sorted(self.expected) if k!='chembl'),size=12,color='#64748B')],spacing=18,data='input')

    def source_options(self):
        options=[ft.DropdownOption(key='',text='Conectar no canvas / selecionar arquivo')]
        options += [ft.DropdownOption(key='stage:'+s['id'],text=s['name']) for s in self.ui.current['pipeline'] if compatible(s,self.stage,self.field['name'])]
        options += [verbatim(ft.DropdownOption(key='asset:'+a['id'],text=a['name'])) for a in self.assets
            if a['kind'] in self.expected or (a['kind']=='other' and self.stage['operation']!='graphs')]
        if self.value:options.append(ft.DropdownOption(key='configured-path',text='Pasta configurada'))
        return options

    def selectors(self,value):
        names=getattr(self.ui,'artifact_choices',{}).get(value.split(':',1)[1],set()) if value and value.startswith('stage:') else set()
        if self.expected & {'compounds','fingerprints','similarity'}:names={n for n in names if n.endswith('.csv')}
        return [ft.DropdownOption(key='auto',text='Identificar automaticamente')]+[verbatim(ft.DropdownOption(key=n,text=n)) for n in sorted(names)]

    def add(self,reference=None):
        ref=reference or {}
        value='stage:'+ref['stage'] if 'stage' in ref else 'asset:'+ref['asset'] if 'asset' in ref else ''
        source=ft.Dropdown(label=self.field['label'],value=value,options=self.source_options(),disabled=not self.writable)
        selector=ft.Dropdown(label='Resultado usado nesta entrada',value=ref.get('selector','auto'),options=self.selectors(value),disabled=not self.writable)
        if selector.value not in {o.key for o in selector.options}:selector.options.append(ft.DropdownOption(key=selector.value,text=selector.value))
        def change(e):
            selector.options=self.selectors(source.value);selector.value='auto'
            self.changed()
        source.on_select=change
        selector.on_select=lambda e:self.changed()
        source.expand=selector.expand=True
        row=ft.Column([ft.Row([source]),ft.Row([selector])],spacing=20,horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
        item=(source,selector,row)
        def remove(e):self.rows.remove(item);self.box.controls.remove(row);self.changed()
        row.controls.append(ft.TextButton('Remover entrada',icon=ft.Icons.REMOVE_CIRCLE_OUTLINE,on_click=remove,disabled=not self.writable))
        self.rows.append(item);self.box.controls.append(row)

    def changed(self):
        if self.on_change:self.on_change()
        self.ui.page.update()

    def append(self,e):self.add();self.changed()

    async def upload(self,e):
        async def action():
            token,project_id=self.ui.token,self.ui.current['id']
            kind=next(iter(sorted(self.expected-{'chembl'})), 'other')
            ids=await self.ui.pick_uploads(kind,project_id)
            if self.ui.token!=token or not self.ui.current or self.ui.current['id']!=project_id:return
            fresh=await self.ui.call(self.ui.store.assets,token,project_id)
            self.assets[:]=fresh
            for source,_,_ in self.rows:source.options=self.source_options()
            for asset in ids:self.add({'asset':asset})
            self.changed()
        await self.ui.guard(action)

    def read(self):
        refs=[]
        for source,selector,_ in self.rows:
            if source.value and source.value!='configured-path':
                kind,identifier=source.value.split(':',1);refs.append({kind:identifier,'selector':selector.value or 'auto'})
        return pack(refs)

    def direct_path(self):
        return self.value if any(c.value=='configured-path' for c,_,_ in self.rows) else None


class ProvidedResults:
    def __init__(self,ui,stage,assets,writable):
        self.ui,self.stage,self.assets,self.writable=ui,stage,assets,writable
        self.kind=output_kind(stage['operation'])
        previous=stage.get('provided_results') or {}
        self.mode=ft.Dropdown(label='Como usar este bloco?',value='provided' if previous else 'execute',disabled=not writable,
            options=[ft.DropdownOption(key='execute',text='Executar a etapa com as entradas configuradas'),
                     ft.DropdownOption(key='provided',text='Usar meus resultados prontos · etapa concluída')])
        self.checks={}
        self.box=ft.Column(spacing=8)
        self.populate(previous.get('asset_ids',[]))
        self.details=ft.Column([ft.Text(contract(self.kind,stage['operation']),size=12,color='#64748B'),self.box,
            ft.TextButton('Enviar resultados prontos',icon=ft.Icons.UPLOAD_FILE,on_click=self.upload,disabled=not writable),
            ft.Text('Este bloco será considerado concluído. O pipeline usa estes arquivos nas etapas seguintes e não executa o cálculo deste bloco.',size=12)],spacing=18,visible=bool(previous))
        def change(e):self.details.visible=self.mode.value=='provided';ui.page.update()
        self.mode.on_select=change
        self.mode.expand=True
        self.control=ft.Column([ft.Row([self.mode]),self.details],spacing=20,data='input',horizontal_alignment=ft.CrossAxisAlignment.STRETCH)

    def populate(self,selected):
        self.checks={a['id']:verbatim(ft.Checkbox(label=a['name'],value=a['id'] in selected,disabled=not self.writable),'label')
            for a in self.assets if a['kind'] in (self.kind,'other','visualization')}
        self.box.controls=list(self.checks.values())

    async def upload(self,e):
        async def action():
            token,project_id=self.ui.token,self.ui.current['id']
            selected=[i for i,c in self.checks.items() if c.value]
            ids=await self.ui.pick_uploads(self.kind,project_id)
            if self.ui.token!=token or not self.ui.current or self.ui.current['id']!=project_id:return
            self.assets[:]=await self.ui.call(self.ui.store.assets,token,project_id)
            self.populate(selected+ids);self.ui.page.update()
        await self.ui.guard(action)

    def read(self):
        if self.mode.value!='provided':return None
        ids=[key for key,c in self.checks.items() if c.value]
        if not ids:raise ValueError('Selecione resultados prontos. Padrão esperado: '+contract(self.kind,self.stage['operation']))
        result={'kind':self.kind,'asset_ids':ids}
        if self.kind in ('structures','prepared_structures'):result['target']=self.stage['parameters'].get('target','MeuAlvo')
        return result
