"""Multiple typed sources and completed-result uploads in guided block forms."""
import flet as ft
from biomolexplorer.artifact_choices import logical_path, matches_selector
from biomolexplorer.bindings import sources, pack
from biomolexplorer.flow import INPUTS, EXTERNAL_INPUTS, compatible, input_types
from biomolexplorer.input_validation import contract, output_kind
from biomolexplorer.docking_inputs import target_input, prepared_receptor, ligand_input_file
from .localization import verbatim, stage_control


class InputEditor:
    def __init__(self,ui,stage,field,assets,writable):
        self.ui,self.stage,self.field,self.assets,self.writable=ui,stage,field,assets,writable
        self.is_target=target_input(stage,field['name'])
        self.is_ligand=stage['operation'] in ('docking_vina','docking_dock6','prepare_structures') and field['name']=='base_selected_mols'
        self.rows=[]
        self.compound_fields={}
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
            ft.Text('Selecione um complexo PDB ou reutilize um receptor preparado. Escolha os compostos ou resultados de docking em “Compostos selecionados”.'
                if stage['operation'] in ('docking_vina','docking_dock6') else '\n'.join(contract(k,'graphs' if stage['operation']=='graphs' else None)
                for k in sorted(self.expected) if k!='chembl'),size=12,color='#64748B')],spacing=18,data='input')

    def source_options(self):
        options=[ft.DropdownOption(key='',text='Conectar no canvas / selecionar arquivo')]
        options += [stage_control(ft.DropdownOption(key='stage:'+s['id'],text=s['name']),s,'text') for s in self.ui.current['pipeline'] if compatible(s,self.stage,self.field['name'])]
        options += [verbatim(ft.DropdownOption(key='asset:'+a['id'],text=a['name'])) for a in self.assets
            if (a['kind'] in self.expected or (a['kind']=='other' and self.stage['operation']!='graphs'))
            and (self.stage['operation'] in ('prepare_structures','docking_vina','docking_dock6') or not self.is_target or prepared_receptor(a['name']))
            and (not self.is_ligand or ligand_input_file(a['name']))
            and (self.stage['operation'] not in ('docking_vina','docking_dock6') or self.field['name']!='base_input_path'
                 or prepared_receptor(a['name']) or a['name'].endswith('.pdb') and a['name'].count('.')==1)]
        if self.value:options.append(ft.DropdownOption(key='configured-path',text='Pasta configurada'))
        return options

    def mixed_receptor_source(self,value):
        if self.stage['operation'] not in ('prepare_structures','docking_vina','docking_dock6') or self.field['name']!='base_input_path' or not value or not value.startswith('stage:'):return False
        from biomolexplorer.flow import output_types
        origin=next((s for s in self.ui.current['pipeline'] if s['id']==value.split(':',1)[1]),None)
        return bool(origin and origin['operation']=='import_results' and {'structures','prepared_structures'}<=output_types(origin))

    def prepared_source(self,value,selector=None):
        if self.stage['operation'] not in ('prepare_structures','docking_vina','docking_dock6') or self.field['name']!='base_input_path':return self.is_target
        if self.mixed_receptor_source(value) and selector and selector!='auto':return prepared_receptor(selector)
        if value and value.startswith('stage:'):
            from biomolexplorer.flow import output_types
            origin=next((s for s in self.ui.current['pipeline'] if s['id']==value.split(':',1)[1]),None)
            return bool(origin and 'prepared_structures' in output_types(origin))
        if value and value.startswith('asset:'):
            return any(a['id']==value.split(':',1)[1] and prepared_receptor(a['name']) for a in self.assets)
        return bool(value=='configured-path' and self.stage['parameters'].get('receptor_prepared',False))

    def selectors(self,value):
        is_target=self.prepared_source(value)
        names=getattr(self.ui,'artifact_choices',{}).get(value.split(':',1)[1],set()) if value and value.startswith('stage:') else set()
        if is_target and value and value.startswith('stage:') and value.split(':',1)[1] not in getattr(self.ui,'artifact_choices',{}):
            names={ref.get('selector','auto') for ref in sources(self.stage.get('bindings',{}).get(self.field['name'],{}))
                if ref.get('stage')==value.split(':',1)[1]}
        names={logical_path(n) for n in names}
        if self.stage['operation']=='retrieve_zinc':
            from pathlib import Path
            from biomolexplorer.zinc_retrieval import LIST_SUFFIXES
            names={n for n in names if Path(n).suffix.lower() in LIST_SUFFIXES}
        if self.is_ligand:names={n for n in names if ligand_input_file(n) and n.rsplit('/',1)[-1] not in ('pdb_codes.csv','centers.csv')}
        elif self.expected & {'compounds','fingerprints','similarity'}:names={n for n in names if n.endswith('.csv')}
        if self.mixed_receptor_source(value):
            names={n for n in names if prepared_receptor(n) or n.endswith('.pdb') and n.rsplit('/',1)[-1].count('.')==1}
        elif is_target:names={n for n in names if prepared_receptor(n)}
        if self.stage['operation'] in ('prepare_structures','docking_vina','docking_dock6') and self.field['name']=='base_input_path' and not is_target:
            names={n for n in names if n.endswith('.pdb') and n.rsplit('/',1)[-1].count('.')==1}
        return [ft.DropdownOption(key='auto',text='Identificar automaticamente')]+[verbatim(ft.DropdownOption(key=n,text=n)) for n in sorted(names)]

    def add(self,reference=None):
        ref=reference or {}
        value='stage:'+ref['stage'] if 'stage' in ref else 'asset:'+ref['asset'] if 'asset' in ref else ''
        source=ft.Dropdown(label=self.field['label'],value=value,options=self.source_options(),disabled=not self.writable)
        selector=ft.Dropdown(label='Receptor que deseja utilizar' if self.stage['operation'] in ('prepare_structures','docking_vina','docking_dock6') and self.field['name']=='base_input_path' else 'Resultado usado nesta entrada',value=logical_path(ref.get('selector','auto')),options=self.selectors(value),disabled=not self.writable)
        if selector.value not in {o.key for o in selector.options}:
            if self.is_target:selector.value='auto'
            else:selector.options.append(ft.DropdownOption(key=selector.value,text=selector.value))
        if self.is_target and source.value not in {o.key for o in source.options}:source.value=''
        def change(e):
            selector.options=self.selectors(source.value);selector.value='auto'
            self.changed()
        source.on_select=change
        selector.on_select=lambda e:self.changed()
        source.expand=selector.expand=True
        row=ft.Column([ft.Row([source]),ft.Row([selector])],spacing=20,horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
        if self.stage['operation'] in ('docking_vina','docking_dock6') and self.is_ligand:
            compound=ft.Dropdown(label='Composto para docking (opcional)',value=ref.get('compound_id',''),
                options=[],disabled=not self.writable,enable_filter=True,enable_search=True)
            self.compound_fields[id(source)]=compound
            def sync_compounds(reset=False):
                if reset:compound.value=''
                mapping=getattr(self.ui,'docking_compounds',{}).get(source.value,{})
                records={}
                for path,values in mapping.items():
                    if (selector.value=='auto' or source.value and source.value.startswith('asset:')
                            or matches_selector(path,selector.value)):
                        records.update(values)
                compound.options=[ft.DropdownOption(key='',text='Todos os compostos')]+[
                    verbatim(ft.DropdownOption(key=code,text=code)) for code in sorted(records)]
                if compound.value and compound.value not in records:
                    compound.options.append(verbatim(ft.DropdownOption(key=compound.value,text=compound.value)))
                compound.visible=bool(source.value and (source.value.startswith('asset:') or selector.value!='auto'))
                compound.disabled=not self.writable or not records
                compound.helper=None if records else 'Execute a origem para listar os compostos disponíveis.'
            def compound_source_changed(e):
                selector.options=self.selectors(source.value);selector.value='auto'
                sync_compounds(True);self.changed()
            def compound_file_changed(e):sync_compounds(True);self.changed()
            source.on_select=compound_source_changed
            selector.on_select=compound_file_changed
            compound.on_select=lambda e:self.changed()
            row.controls.append(compound)
            sync_compounds()
        item=(source,selector,row)
        def remove(e):self.rows.remove(item);self.compound_fields.pop(id(source),None);self.box.controls.remove(row);self.changed()
        row.controls.append(ft.TextButton('Remover entrada',icon=ft.Icons.REMOVE_CIRCLE_OUTLINE,on_click=remove,disabled=not self.writable))
        self.rows.append(item);self.box.controls.append(row)

    def changed(self):
        if self.on_change:self.on_change()
        self.ui.page.update()

    def append(self,e):self.add();self.changed()

    async def upload(self,e):
        async def action():
            token,project_id=self.ui.token,self.ui.current['id']
            kind='other' if self.stage['operation'] in ('docking_vina','docking_dock6') and self.field['name']=='base_input_path' else 'zinc_urls' if self.stage['operation']=='retrieve_zinc' else next(iter(sorted(self.expected-{'chembl'})), 'other')
            ids=await self.ui.pick_uploads(kind,project_id)
            if self.ui.token!=token or not self.ui.current or self.ui.current['id']!=project_id:return
            fresh=await self.ui.call(self.ui.store.assets,token,project_id)
            self.assets[:]=fresh
            if self.stage['operation']=='retrieve_pubchem' or self.is_ligand and self.stage['operation'] in ('docking_vina','docking_dock6'):
                from biomolexplorer.pubchem_retrieval import reference_choices
                attribute='pubchem_references' if self.stage['operation']=='retrieve_pubchem' else 'docking_compounds'
                if not hasattr(self.ui,attribute):setattr(self.ui,attribute,{})
                for asset in ids:
                    path=await self.ui.call(self.ui.store.asset_path,token,project_id,asset)
                    getattr(self.ui,attribute)['asset:'+asset]={str(path):await self.ui.call(reference_choices,[path])}
                if self.ui.token!=token or not self.ui.current or self.ui.current['id']!=project_id:return
            for source,_,_ in self.rows:source.options=self.source_options()
            allowed={o.key for o in self.source_options()}
            for asset in ids:
                if not self.is_target or 'asset:'+asset in allowed:self.add({'asset':asset})
            self.changed()
        await self.ui.guard(action)

    def read(self):
        if self.stage['operation'] in ('prepare_structures','docking_vina','docking_dock6') and self.field['name']=='base_input_path':
            modes=[self.prepared_source(source.value,selector.value) for source,selector,_ in self.rows if source.value]
            if any(modes) and not all(modes):
                raise ValueError('Use receptores brutos ou preparados em um mesmo bloco de preparação.')
        refs=[]
        for source,selector,_ in self.rows:
            if source.value and source.value!='configured-path':
                if self.prepared_source(source.value,selector.value) and (source.value not in {o.key for o in self.source_options()}
                        or selector.value not in {o.key for o in self.selectors(source.value)}):
                    raise ValueError('Selecione um receptor preparado para o alvo do docking.')
                kind,identifier=source.value.split(':',1)
                ref={kind:identifier,'selector':selector.value or 'auto'}
                compound=self.compound_fields.get(id(source))
                if compound and compound.value:ref['compound_id']=compound.value
                refs.append(ref)
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
