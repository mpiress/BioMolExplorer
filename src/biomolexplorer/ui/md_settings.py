"""Grouped AMBER settings and paged, exact receptor/ligand selection."""
import copy
import flet as ft
from biomolexplorer.md_config import DEFAULTS,HELP,INTEGER,BOOLEANS,CHOICES
from biomolexplorer.md_inputs import MDInputs,pair_id
from .localization import verbatim
from .feedback import close_dialog
from .settings_layout import section

LABELS=dict(zip(DEFAULTS,(
 'Protocolo de dinâmica','Região móvel no AMBER Score','Convergência da minimização (RMS gradiente)',
 'Minimização antes da MD (ciclos)','Duração da dinâmica (passos)','Minimização após a MD (ciclos)',
 'Modelo de solvente implícito GB/SA','Corte de interações (Å)','Temperatura (K)',
 'Raio da região móvel (Å)','Expressão NAB do receptor','Expressão NAB do complexo',
 'Interromper o par sem preparação AMBER','Preservar as cargas da pose MOL2',
 'Passo de integração sander (ps)','Gravar a cada N passos','Aquecimento (passos)',
 'Equilibração (passos)','Semente da trajetória',
 'Máscara Amber da região móvel','Calcular RMSD, RMSF, ligações H e raio de giro')))


class MDSettings:
    def __init__(self,form):
        self.form=form;self.controls={};self.pairs=[];self.offset=0;self.limit=50
        params=form.stage['parameters']
        self.selected=set(params.get('selected_pairs') or [])
        self.excluded=set(params.get('excluded_pairs') or [])
        self.summary=ft.Text();self.box=ft.Column(spacing=6)
        self.search=ft.TextField(label='Buscar receptor ou ligante',on_change=self.filter)
        self.pair_control=ft.Column([ft.Text('Combinações receptor / ligante',size=16,weight=ft.FontWeight.W_600),
            ft.Text('Confirme as poses e os receptores. Cada combinação selecionada é uma execução independente.',size=12),self.search,
            ft.Row([ft.TextButton('Selecionar todos os pares',on_click=self.all,disabled=not form.writable),
                ft.TextButton('Limpar seleção de pares',on_click=self.none,disabled=not form.writable)],wrap=True),
            self.summary,self.box,ft.Row([ft.TextButton('Anterior',on_click=lambda e:self.move(-1)),
                ft.TextButton('Próxima',on_click=lambda e:self.move(1))])],spacing=12,data='input')
        for key,default in DEFAULTS.items():
            value=params.get(key,default)
            if key in BOOLEANS:
                control=ft.Switch(label=LABELS[key],value=value,disabled=not form.writable);getter=lambda c=control:bool(c.value)
            elif key in CHOICES:
                control=ft.Dropdown(label=LABELS[key],value=str(value),disabled=not form.writable,
                    options=[ft.DropdownOption(key=str(v),text={'dock6':'AMBER Score (DOCK6)','dock6_sander':'AMBER Score + trajetória sander'}.get(str(v),str(v))) for v in CHOICES[key]])
                getter=(lambda c=control:int(c.value)) if key=='amber_score_gb_model' else (lambda c=control:c.value)
            else:
                control=ft.TextField(label=LABELS[key],value=str(value),disabled=not form.writable)
                kind=int if key in INTEGER else float if isinstance(default,(float,int)) else str
                getter=lambda c=control,kind=kind:kind(c.value)
            control.tooltip=HELP[key]
            form.readers[key]=getter;form.field_controls[key]=control;self.controls[key]=control
        form.readers['selected_pairs']=lambda:sorted(self.selected)
        form.readers['excluded_pairs']=lambda:sorted(self.excluded)
        def group(title,keys,expanded=False):
            cells=[]
            for key in keys:
                def help(e,key=key):
                    dialog=ft.AlertDialog(title=ft.Text(LABELS[key]),content=ft.Text(HELP[key]),
                        actions=[ft.TextButton('Fechar',on_click=lambda e:close_dialog(form.ui.page,dialog))])
                    form.ui.page.show_dialog(dialog)
                self.controls[key].tooltip=HELP[key]
                cells.append(self.controls[key])
            return section(title,cells,expanded=expanded)
        self.control=ft.Column([ft.Text('O AMBER Score realiza refinamento curto em solvente implícito. Para trajetória e análises por quadro, selecione sander e configure AmberTools.',size=13),
            group('Protocolo e preparação',['amber_score_type','use_existing_ligand_charges','amber_score_abort_on_unprepped_ligand'],True),
            group('Minimização e dinâmica',['amber_score_before_md_minimization_cycles','amber_score_md_steps','amber_score_after_md_minimization_cycles','amber_score_minimization_rmsgrad','amber_score_temperature','amber_score_gb_model','amber_score_nonbonded_cutoff'],True),
            group('Região móvel',['amber_score_movable_region','amber_score_movable_distance_cutoff','amber_score_receptor_movable_atom_expr','amber_score_complex_movable_atom_expr']),
            group('Trajetória sander e análises',['md_time_step','md_heating_steps','md_equilibration_steps','trajectory_stride','random_seed','md_movable_mask','analyze_trajectory'])],spacing=16)
        for editor in form.input_editors.values():editor.on_change=self.refresh
        self.refresh()

    def refresh(self):
        self.pairs=[];ui=self.form.ui
        if getattr(ui,'store',None) and ui.current.get('id'):
            candidate=copy.deepcopy(self.form.stage)
            for field,editor in self.form.input_editors.items():
                refs=editor.read()
                if refs:candidate.setdefault('bindings',{})[field]=refs
                else:candidate.setdefault('bindings',{}).pop(field,None)
                direct=editor.direct_path()
                if direct:candidate['parameters'][field]=direct
            inputs=MDInputs(ui.store,ui.current['id'],{key:item.get('artifacts',[]) for key,item in getattr(ui,'artifact_items',{}).items()})
            receptors=inputs.files(candidate,'base_input_path');ligands=inputs.files(candidate,'base_selected_mols')
            if len(receptors)*len(ligands)<=100000:self.pairs=[dict(id=pair_id(r,l),receptor=r,ligand=l) for r in receptors for l in ligands]
        self.offset=0;self.redraw()

    def checked(self,key):return key not in self.excluded and key in self.selected
    def redraw(self):
        term=(self.search.value or '').lower();matches=[p for p in self.pairs if term in (p['receptor']['name']+' '+p['ligand']['name']).lower()]
        self.summary.value=str(sum(self.checked(p['id']) for p in self.pairs))+' / '+str(len(self.pairs))+' pares selecionados';self.box.controls=[]
        for pair in matches[self.offset:self.offset+self.limit]:
            def change(e,key=pair['id']):
                if e.control.value:self.selected.add(key);self.excluded.discard(key)
                else:self.selected.discard(key)
                self.redraw();self.form.ui.page.update()
            preview=[]
            for kind in ('receptor','ligand'):
                source=pair[kind]
                path=source['path']
                if kind=='ligand' and source.get('compound'):
                    path=next((source['compound'].get(k) for k in ('pose_file','conformer_file','prepared_mol2','prepared_pdbqt') if source['compound'].get(k)),path)
                async def open(e,path=path):
                    await self.form.ui.guard(lambda:self.form.ui.preview_artifact(self.form.ui.current['id'],path))
                preview.append(ft.IconButton(ft.Icons.VIEW_IN_AR,tooltip='Visualizar receptor' if kind=='receptor' else 'Visualizar ligante',on_click=open))
            checkbox=verbatim(ft.Checkbox(label=pair['receptor']['name']+' / '+pair['ligand']['name'],value=self.checked(pair['id']),disabled=not self.form.writable,on_change=change),'label')
            self.box.controls.append(ft.Row([checkbox,*preview],wrap=True))
        if not self.pairs:self.box.controls=[ft.Text('Os pares aparecerão quando as entradas PDB e MOL2 estiverem disponíveis.',size=12)]
    def filter(self,e):self.offset=0;self.redraw();self.form.ui.page.update()
    def move(self,direction):self.offset=max(0,min(max(0,len(self.pairs)-1),self.offset+direction*self.limit));self.redraw();self.form.ui.page.update()
    def all(self,e):self.selected={p['id'] for p in self.pairs};self.excluded.clear();self.filter(e)
    def none(self,e):self.selected=set();self.excluded.clear();self.filter(e)
