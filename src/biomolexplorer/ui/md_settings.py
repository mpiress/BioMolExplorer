"""Grouped AMBER settings and paged, exact receptor/ligand selection."""
import json
from pathlib import Path
import re
import flet as ft
from biomolexplorer.md_config import DEFAULTS,HELP,INTEGER,BOOLEANS,CHOICES
from biomolexplorer.md_inputs import MDInputs,pair_id
from .localization import verbatim
from .feedback import close_dialog
from .settings_layout import section
from .pair_table import PairTable

UI_HELP={key:value.replace('receptores','targets').replace('receptor','target') for key,value in HELP.items()}

LABELS=dict(zip(DEFAULTS,(
 'Protocolo de pontuação','Região móvel no AMBER Score','Convergência da minimização (RMS gradiente)',
 'Minimização antes do refinamento (ciclos)','Refinamento interno (passos)','Minimização após o refinamento (ciclos)',
 'Modelo de solvente implícito GB/SA','Corte de interações (Å)','Temperatura (K)',
 'Raio da região móvel (Å)','Expressão NAB do target','Expressão NAB do complexo',
 'Interromper o par sem preparação AMBER','Preservar as cargas da pose MOL2','Perfil de parâmetros do target',
)))


class MDSettings:
    def __init__(self,form):
        self.form=form;self.controls={};self.pairs=[];self.offset=0;self.limit=50
        self._input_cache={};self._annotation_cache={};self._match_cache=None;self._search_cache=None;self._pair_sources=None
        params=form.stage['parameters']
        self.selected=set(params.get('selected_pairs') or [])
        self.excluded=set(params.get('excluded_pairs') or [])
        self.summary=ft.Text();self.page_count=ft.Text(size=12)
        self.error=ft.Text(color='#B91C1C',visible=False)
        self.search=ft.TextField(label='Buscar target ou ligante',on_change=self.filter)
        self.target=ft.Dropdown(label='Alvo dos pares',value='all',options=[],on_select=self.filter,
            enable_filter=True,enable_search=True)
        self.pair_table=PairTable(('Alvo','PDB / target','Ligante','Cadeia','Ações'),
            form.ui.page,min_width=700,selectable=True,pin_actions=True)
        self.table=self.pair_table.table;self.table.column_spacing=24
        self.table.on_select_all=self.select_view
        self.box=self.pair_table.control
        self.previous=ft.TextButton('Anterior',on_click=lambda e:self.move(-1))
        self.next=ft.TextButton('Próxima',on_click=lambda e:self.move(1))
        self.empty=ft.Text('Os pares aparecerão quando as entradas PDB e MOL2 estiverem disponíveis.',size=12)
        from biomolexplorer.result_origins import ResultOrigins
        self.origins=ResultOrigins(list(getattr(form.ui,'artifact_items',{}).values()))
        self.pair_control=ft.Column([ft.Text('Ligantes e pares selecionados',size=16,weight=ft.FontWeight.W_600),
            ft.Text('Confirme as poses e os targets. Cada combinação selecionada é uma execução independente.',size=12),
            ft.ResponsiveRow([ft.Container(self.target,col={'xs':12,'md':4}),
                ft.Container(self.search,col={'xs':12,'md':8})],spacing=24,run_spacing=16),
            ft.Row([ft.TextButton('Selecionar pares desta visualização',on_click=self.all,disabled=not form.writable),
                ft.TextButton('Limpar seleção de pares',on_click=self.none,disabled=not form.writable)],wrap=True),
            self.summary,self.error,self.box,self.empty,ft.Row([self.previous,self.next,self.page_count],wrap=True)],
            spacing=18,data='input',horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
        for key,default in DEFAULTS.items():
            value=params.get(key,default)
            if key in BOOLEANS:
                control=ft.Switch(label=LABELS[key],value=value,disabled=not form.writable);getter=lambda c=control:bool(c.value)
            elif key in CHOICES:
                control=ft.Dropdown(label=LABELS[key],value=str(value),disabled=not form.writable,
                    options=[ft.DropdownOption(key=str(v),text={'dock6':'AMBER Score (DOCK6)',
                        'auto_curated':'Automático · ff94 / 32 PTMs FF_PTM (CC BY-NC 4.0)',
                        'native_ff94':'Nativo · ff94 e bibliotecas do usuário',
                        'ff03_ffptm':'ff03 + FF_PTM · perfil uniforme (CC BY-NC 4.0)',
                        'ligand':'Somente o ligante','everything':'Todo o complexo','nothing':'Somente energia (sem refinamento)',
                        'distance':'Resíduos próximos ao ligante','nab_atom_expression':'Expressões de átomos NAB'}.get(str(v),str(v))) for v in CHOICES[key]])
                getter=(lambda c=control:int(c.value)) if key=='amber_score_gb_model' else (lambda c=control:c.value)
            else:
                control=ft.TextField(label=LABELS[key],value=str(value),disabled=not form.writable)
                kind=int if key in INTEGER else float if isinstance(default,(float,int)) else str
                getter=lambda c=control,kind=kind:kind(c.value)
            control.tooltip=UI_HELP[key]
            form.readers[key]=getter;form.field_controls[key]=control;self.controls[key]=control
        form.readers['selected_pairs']=lambda:sorted(p['id'] for p in self.pairs if self.checked(p['id']))
        form.readers['excluded_pairs']=lambda:sorted(self.excluded & {p['id'] for p in self.pairs})
        def group(title,keys,expanded=False):
            cells=[]
            for key in keys:
                def help(e,key=key):
                    dialog=ft.AlertDialog(title=ft.Text(LABELS[key]),content=ft.Text(UI_HELP[key]),
                        actions=[ft.TextButton('Fechar',on_click=lambda e:close_dialog(form.ui.page,dialog))])
                    form.ui.page.show_dialog(dialog)
                self.controls[key].tooltip=UI_HELP[key]
                cells.append(self.controls[key])
            return section(title,cells,expanded=expanded)
        self.control=ft.Column([ft.Text('Amber Score: refinamento curto em solvente implícito GB/SA. Sem caixa de solvatação explícita, simulação convencional ou trajetória para RMSF.',size=13),
            group('Protocolo e preparação',['amber_score_type','residue_parameter_policy','use_existing_ligand_charges','amber_score_abort_on_unprepped_ligand'],True),
            group('Minimização e refinamento',['amber_score_before_md_minimization_cycles','amber_score_md_steps','amber_score_after_md_minimization_cycles','amber_score_minimization_rmsgrad','amber_score_temperature','amber_score_gb_model','amber_score_nonbonded_cutoff'],True),
            group('Região móvel',['amber_score_movable_region','amber_score_movable_distance_cutoff','amber_score_receptor_movable_atom_expr','amber_score_complex_movable_atom_expr'])],spacing=16)
        for editor in form.input_editors.values():editor.on_change=self.refresh
        self.restore_selected_inputs()
        self.refresh()

    def restore_selected_inputs(self):
        """Recover old automatic bindings only from explicitly saved pair IDs."""
        ui=self.form.ui
        if not self.selected or not getattr(ui,'store',None) or not ui.current.get('id'):return
        from biomolexplorer.bindings import sources
        if not any(ref.get('selector','auto')=='auto' for binding in self.form.stage.get('bindings',{}).values()
                   for ref in sources(binding)):return
        inputs=MDInputs(ui.store,ui.current['id'],{key:item.get('artifacts',[]) for key,item in getattr(ui,'artifact_items',{}).items()})
        from biomolexplorer.workspace import AccessDenied
        try:
            receptors=self.input_files(inputs,self.form.stage,'base_input_path')
            ligands=self.input_files(inputs,self.form.stage,'base_selected_mols')
            if len(receptors)*len(ligands)>100000:return
            paths={'base_input_path':set(),'base_selected_mols':set()}
            for receptor in receptors:
                for ligand in ligands:
                    key=pair_id(receptor,ligand)
                    if key in self.selected and key not in self.excluded:
                        paths['base_input_path'].add(str(receptor['path']))
                        paths['base_selected_mols'].add(str(ligand['path']))
            for field,editor in self.form.input_editors.items():
                if not editor.selection:continue
                rows=editor.selection.rows[field]
                for flag,ref in rows:
                    if getattr(flag,'physical_file',None) in paths.get(field,set()):flag.value=True
                editor.selection.pickers[field].redraw()
        except (ValueError,OSError,AccessDenied):
            # Missing historical files stay unselected; execution cannot guess.
            return

    def annotate(self,source):
        key=(source['identity'],str(source['path']))
        if key in self._annotation_cache:return self._annotation_cache[key]
        result=dict(source)
        identity=source['identity']
        origin=identity.split(':',2)[1] if identity.startswith('stage:') else None
        description=self.origins.describe(source['path'],origin) if origin else {}
        result['target']=description.get('source_label') if description.get('source_routes') else ''
        self._annotation_cache[key]=result
        return result

    def input_files(self,inputs,candidate,field):
        from biomolexplorer.result_index import file_stamp
        key=json.dumps([candidate.get('bindings',{}).get(field),candidate['parameters'].get(field)],sort_keys=True)
        cached=self._input_cache.get(field)
        if cached and cached[0]==key:
            try:
                if all(file_stamp(path)==stamp for path,stamp in cached[2]):return cached[1]
            except OSError:self._input_cache.pop(field,None)
        rows=inputs.files(candidate,field)
        paths={str(row['path']) for row in rows}
        direct=candidate['parameters'].get(field)
        if direct:paths.add(str(inputs.store.scoped_path(inputs.pid,direct)))
        stamps=[(path,file_stamp(path)) for path in sorted(paths)]
        self._input_cache[field]=(key,rows,stamps)
        return rows

    def refresh(self):
        previous=self.pairs
        self.pairs=[];self.error.visible=False;ui=self.form.ui
        if getattr(ui,'store',None) and ui.current.get('id'):
            candidate=dict(self.form.stage,parameters=dict(self.form.stage['parameters']),bindings=dict(self.form.stage.get('bindings',{})))
            for field,editor in self.form.input_editors.items():
                if editor.selection:
                    from biomolexplorer.bindings import pack
                    refs=pack(editor.selection.selected_references(editor.selection.rows[field],validate=False))
                else:refs=editor.read(validate=False)
                if refs:candidate.setdefault('bindings',{})[field]=refs
                else:candidate.setdefault('bindings',{}).pop(field,None)
                direct=editor.direct_path()
                if direct:candidate['parameters'][field]=direct
                else:candidate['parameters'].pop(field,None)
            inputs=MDInputs(ui.store,ui.current['id'],{key:item.get('artifacts',[]) for key,item in getattr(ui,'artifact_items',{}).items()})
            self._pose_inputs=inputs
            from biomolexplorer.workspace import AccessDenied
            try:
                receptors=[self.annotate(r) for r in self.input_files(inputs,candidate,'base_input_path')]
                ligands=self.input_files(inputs,candidate,'base_selected_mols')
                if len(receptors)*len(ligands)>100000:
                    raise ValueError('Selecione menos de 100000 combinações para a dinâmica.')
                source_key=tuple((field,self._input_cache[field][0],tuple(self._input_cache[field][2]))
                    for field in ('base_input_path','base_selected_mols'))
                if source_key==self._pair_sources:self.pairs=previous
                else:
                    self.pairs=[dict(id=pair_id(r,l),receptor=r,ligand=l) for r in receptors for l in ligands]
                    self._pair_sources=source_key
            except (ValueError,AccessDenied) as error:
                # A draft may temporarily reference unavailable outputs. Keep
                # editing usable; execution still performs strict validation.
                self.error.value=str(error);self.error.visible=True;self._pair_sources=None
        targets=sorted({p['receptor'].get('target','') for p in self.pairs}-{''})
        self.target.options=[ft.DropdownOption(key='all',text='Todos os alvos')]+[
            verbatim(ft.DropdownOption(key=name,text=name),'text') for name in targets]
        if self.target.value not in {'all',*targets}:self.target.value='all'
        self._match_cache=None
        if not self._search_cache or self._search_cache[:2]!=(id(self.pairs),len(self.pairs)):
            self._search_cache=(id(self.pairs),len(self.pairs),[self.search_text(pair) for pair in self.pairs])
        self.offset=0;self.redraw()

    def checked(self,key):return key not in self.excluded and key in self.selected

    @staticmethod
    def details(pair):
        receptor=Path(pair['receptor']['path']).name
        match=re.match(r'^([A-Za-z0-9]{4})_([A-Za-z0-9]+)\.',receptor)
        pdb,chain=(match[1],match[2]) if match else (pair['receptor']['name'],'—')
        ligand=pair['ligand'];name=ligand.get('compound',{}).get('molecule_chembl_id') or ligand['name']
        pose=next((ligand.get('compound',{}).get(k) for k in ('pose_file','conformer_file','prepared_mol2','prepared_pdbqt') if ligand.get('compound',{}).get(k)),ligand['path'])
        match=re.match(r'^([A-Za-z0-9]{4})_([A-Za-z0-9]+)_(-?\d+)([A-Za-z0-9])(?:\.|_)',Path(pose).name)
        residue='—'
        if match:name=match[2];residue=match[3];chain=match[4]
        return pdb,name,residue,chain

    def search_text(self,pair):
        return ' '.join((pair['receptor']['name'],pair['ligand']['name'],pair['receptor'].get('target',''),*self.details(pair))).casefold()

    def matches(self):
        term=(self.search.value or '').casefold()
        key=(id(self.pairs),len(self.pairs),term,self.target.value)
        if self._match_cache and self._match_cache[0]==key:return self._match_cache[1]
        if not self._search_cache or self._search_cache[:2]!=key[:2]:
            self._search_cache=(*key[:2],[self.search_text(pair) for pair in self.pairs])
        result=[pair for pair,text in zip(self.pairs,self._search_cache[2])
            if (self.target.value=='all' or pair['receptor'].get('target')==self.target.value) and term in text]
        self._match_cache=(key,result)
        return result

    def redraw(self):
        matches=self.matches()
        self.offset=min(self.offset,max(0,((len(matches)-1)//self.limit)*self.limit))
        self.summary.value=str(sum(self.checked(p['id']) for p in self.pairs))+' / '+str(len(self.pairs))+' pares selecionados'
        rendered=[]
        for pair in matches[self.offset:self.offset+self.limit]:
            def change(e,key=pair['id']):
                selected=getattr(e,'selected',None)
                if selected is None:selected=getattr(e,'data',False) in (True,'true')
                if selected:self.selected.add(key);self.excluded.discard(key)
                else:self.selected.discard(key)
                self.redraw();self.form.ui.page.update()
            actions=[]
            for kind in ('ligand',):
                source=pair[kind];path=source['path']
                if kind=='ligand' and source.get('compound'):
                    path=next((source['compound'].get(k) for k in ('pose_file','conformer_file','prepared_mol2','prepared_pdbqt') if source['compound'].get(k)),path)
                async def open(e,path=path,pair=pair):
                    async def preview():
                        selected=path
                        if pair['ligand'].get('compound',{}).get('target_scores'):
                            ligand=await self.form.ui.call(self._pose_inputs.targeted_ligand,pair['receptor'],pair['ligand'])
                            selected=ligand['compound']['conformer_file']
                        await self.form.ui.preview_artifact(self.form.ui.current['id'],selected)
                    await self.form.ui.guard(preview)
                actions.append(ft.IconButton(ft.Icons.VIEW_IN_AR,tooltip='Visualizar target' if kind=='receptor' else 'Visualizar ligante',on_click=open))
            def remove(e,key=pair['id']):
                if not self.form.writable:return
                self.selected.discard(key);self.excluded.add(key);self.redraw();self.form.ui.page.update()
            actions.append(ft.IconButton(ft.Icons.REMOVE_CIRCLE_OUTLINE,tooltip='Remover par',on_click=remove,
                disabled=not self.form.writable or not self.checked(pair['id'])))
            pdb,ligand_name,_,chain=self.details(pair)
            values=(pair['receptor'].get('target') or '—',pdb,ligand_name,chain)
            widths=(150,135,170,50)
            cells=[ft.DataCell(verbatim(ft.Text(value,width=width,max_lines=2,overflow=ft.TextOverflow.ELLIPSIS,
                tooltip=pair['receptor']['name'] if index==1 else pair['ligand']['name'] if index==2 else value),'value','tooltip'))
                for index,(value,width) in enumerate(zip(values,widths))]
            cells.append(ft.DataCell(ft.Row(actions,spacing=4)))
            checked=self.checked(pair['id'])
            rendered.append(ft.DataRow(cells=cells,selected=checked,color='#F0FDF4' if checked else None,
                on_select_change=change if self.form.writable else None))
        self.pair_table.set_rows(rendered)
        self.empty.visible=not matches
        self.previous.disabled=self.offset==0;self.next.disabled=self.offset+self.limit>=len(matches)
        self.page_count.value=f'Exibindo pares {self.offset+1 if matches else 0}–{min(self.offset+self.limit,len(matches))} de {len(matches)}.'

    def filter(self,e):self.offset=0;self.redraw();self.form.ui.page.update()
    def move(self,direction):self.offset=max(0,self.offset+direction*self.limit);self.redraw();self.form.ui.page.update()
    def all(self,e):
        if not self.form.writable:return
        keys={p['id'] for p in self.matches()};self.selected.update(keys);self.excluded.difference_update(keys);self.filter(e)
    def none(self,e):
        if not self.form.writable:return
        self.selected.clear();self.excluded.clear();self.filter(e)
    def select_view(self,e):
        if not self.form.writable:return
        selected=getattr(e,'selected',None)
        if selected is None:selected=getattr(e,'data',False) in (True,'true')
        if selected:self.all(e)
        else:
            self.selected.difference_update(p['id'] for p in self.matches());self.filter(e)
