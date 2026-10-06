"""Structured stage configuration. Script generation stays behind typed controls."""
import copy
import json
import re
import flet as ft
from biomolexplorer.catalog import operation_fields, template_names, TITLES
from biomolexplorer.retrieval import MODES
from biomolexplorer.templates import RESOURCE_ROOT, validate_templates
from biomolexplorer.flow import compatible
from .localization import verbatim
from .template_labels import TEMPLATE_FIELDS

TEMPLATE_LABELS={'target':'Filtros do alvo molecular','bioactivity':'Filtros de atividade biológica','molecules':'Filtros das moléculas ChEMBL','similarmols':'Filtros dos similares ChEMBL','config':'Vina · opções do motor','prepare_complex':'Chimera · preparação do complexo','prepare_ligand':'Chimera · preparação do ligante','prepare_receptor':'Chimera · preparação do receptor','prepare_better_conform':'Chimera · conformação do ligante','prepare_md':'Chimera · preparação adicional','docking':'DOCK6 · docking','grid':'DOCK6 · grade','min':'DOCK6 · minimização','footprint':'DOCK6 · footprint','showbox':'DOCK6 · caixa','INSPH':'DOCK6 · esferas'}
LISTS={'organism':['Homo sapiens','Mus musculus','Rattus norvegicus'],
 'PolymerEntityTypeID':['Protein','DNA','RNA','NA-hybrid','Other'],
 'ExperimentalMethodID':['X-RAY DIFFRACTION','SOLUTION NMR','ELECTRON MICROSCOPY','ELECTRON CRYSTALLOGRAPHY','EPR','FIBER DIFFRACTION','FLUORESCENCE TRANSFER','INFRARED SPECTROSCOPY','NEUTRON DIFFRACTION','POWDER DIFFRACTION','SOLID-STATE NMR','SOLUTION SCATTERING','THEORETICAL MODEL']}
FILTER_LABELS={'organism':'Organismo','type__in':'Tipos de alvo','relationship_type':'Relação com o alvo','standard_type__in':'Medidas de atividade (Ki, IC50…)','molecule_type':'Tipo de molécula','max_value_ref':'Atividade máxima','standard_units':'Unidade de atividade','assay_type':'Tipo de ensaio','pchembl_value__isnull':'Ausência de pChEMBL (0 = exigir valor)','natural_product':'Produto natural (1 = sim, 0 = não)','similarity':'Similaridade mínima (%)','molecule_weight':'Massa molecular máxima'}
NUMERIC={'max_resolution','pH','threshold','pubchem_threshold','pubchem_max_records','max_records','radius','morgan_n_bits','exhaustiveness','num_modes','chunk_size','repulsion_weight','density','distance','plot_max_residues','max_targets','similarity_threshold'}

class GuidedForm:
    def __init__(self,ui,stage,assets,writable):
        self.ui,self.stage,self.assets,self.writable=ui,stage,assets,writable
        self.readers={}; self.binding_readers={}; self.direct_path_readers={}; self.template_readers={}
        self.input_editors={};self.field_controls={};self.conditional_cells={};self.filter_tiles={};self.filter_cells={}
        # Consensus uses its Vina input to derive the legacy base_input_path.
        # Do not keep a second, hidden dependency on the previously linked block.
        self.derived_consensus_input=(stage['operation']=='consensus'
            and not stage['parameters'].get('base_input_path')
            and stage.get('bindings',{}).get('base_input_path') in
                (None,stage.get('bindings',{}).get('base_vina_path')))
        self.title=ft.TextField(label='Nome do bloco',value=stage['name'],disabled=not writable,data='wide')
        self.enabled=ft.Switch(label='Incluir este bloco na execução',value=stage.get('enabled',True),disabled=not writable)
        self.controls=[self.title,self.enabled]
        if stage['operation']!='import_results':
            self.processing=ft.Dropdown(label='Como processar os arquivos selecionados?',
                value=stage.get('input_processing','individual'),disabled=not writable,options=[
                    ft.DropdownOption(key='individual',text='Processar individualmente'),
                    ft.DropdownOption(key='merge',text='Mesclar arquivos (merge)')])
            self.controls.append(self.processing)
        if stage['operation'] in ('retrieve_compounds','retrieve_structures','retrieve_zinc'):
            self.controls.append(ft.Text(
                'Ao concluir cada etapa, selecione no popup os arquivos que seguirão para o próximo bloco.',
                size=13,color='#64748B'))
        if stage['operation']=='retrieve_structures':
            self.controls.append(ft.Text('Escolha texto, IDs PDB, UniProt, EC ou filtros. Campos preenchidos são combinados; o nome da coleção só organiza a saída.',size=13))
        elif stage['operation']=='retrieve_compounds':
            self.controls.append(ft.Text('Busque alvos para obter bioatividades ou compostos diretamente por nome, IDs, similaridade ou subestrutura. Limpe filtros para ampliar a seleção.',size=13))
        if stage['operation']=='import_results':
            from .app import KINDS
            self.add('kind','Tipo de dados',stage['parameters'].get('kind','compounds'),list(KINDS))
            self.add('target','Nome do alvo',stage['parameters'].get('target','MeuAlvo'))
            checks={a['id']:verbatim(ft.Checkbox(label=a['name'],value=a['id'] in stage['parameters'].get('asset_ids',[]),disabled=not writable),'label') for a in assets}
            self.controls += [ft.Text('Selecione arquivos na aba Arquivos; depois volte ao pipeline.',size=12),*checks.values()]
            self.readers['asset_ids']=lambda:[k for k,c in checks.items() if c.value]
            return
        from .input_editor import ProvidedResults
        self.provided=ProvidedResults(ui,stage,assets,writable)
        self.controls.append(self.provided.control)
        if stage['operation']=='graphs':
            self.controls.append(ft.Text('Conecte um ou mais blocos Calcular similaridade. Use processamento individual para gerar um grafo por arquivo ou merge para reunir similaridades compatíveis. Opcionalmente, envie seus próprios CSVs de similaridade (source,target,value). Para entradas externas, uma tabela de códigos e SMILES permite visualizar moléculas, preservar nós isolados e calcular o fragmento comum.',size=13,color='#64748B',data='input'))
        if stage['operation']=='fingerprints':
            from biomolexplorer.fingerprint_selection import KINDS,LABELS
            selected=next((k for k in KINDS if stage['parameters'].get(k,k=='morgan')),'morgan')
            self.fingerprint_choice=ft.Dropdown(label='Tipo de fingerprint',value=selected,
                options=[ft.DropdownOption(key=k,text=LABELS[k]) for k in KINDS],disabled=not writable)
            self.controls.append(self.fingerprint_choice)
            for k in KINDS:self.readers[k]=lambda k=k:self.fingerprint_choice.value==k
        for field in operation_fields(stage['operation']):
            key=field['name']; value=stage['parameters'].get(key,field['default'])
            if key=='graph_inputs':continue
            if stage['operation']=='redocking' and key in ('pdb_codes','preparation_pairs','charge_type'):continue
            if stage['operation']=='fingerprints' and key in ('morgan','maccs','pharmacophore'):continue
            if key=='base_input_path' and self.derived_consensus_input:
                continue
            if field['path']:
                from .input_editor import InputEditor
                if stage['operation']=='graphs' and key=='base_input_path':
                    field=dict(field,label='Arquivo externo de compostos e SMILES · opcional')
                editor=InputEditor(ui,stage,field,assets,writable)
                self.input_editors[key]=editor
                self.controls.append(editor.control)
                self.direct_path_readers[key]=editor.direct_path
                self.binding_readers[key]=editor.read
            elif key=='chembl_filters': continue # Resource filters below also cover this parameter.
            elif key in ('pdb_code','pdb_codes'):
                self.records(key,field['label'],value)
            elif key=='sizeof_box':
                fields=[ft.TextField(label='Caixa '+axis+' (Å)',value=str(v),disabled=not writable,width=150) for axis,v in zip('XYZ',value or [20,20,20])]
                self.controls.append(ft.Row(fields,wrap=True,spacing=24,run_spacing=24)); self.readers[key]=lambda f=fields:[float(c.value) for c in f]
            elif key in LISTS or key=='files':
                self.list_field(key,field['label'],value,LISTS.get(key,[]))
            else: self.add(key,'Nome da coleção PDB (opcional)' if stage['operation']=='retrieve_structures' and key=='target' else field['label'],value,field['choices'])
        if stage['operation']=='redocking':
            from .redocking_pairs import RedockingPairs
            self.redocking_pairs=RedockingPairs(self)
            self.controls.append(self.redocking_pairs.control)
            self.input_editors['base_input_path'].on_change=self.redocking_pairs.refresh
            def redocking_inputs(e):
                from biomolexplorer.flow import input_types
                editor=self.input_editors['base_input_path']
                editor.stage=copy.deepcopy(stage)
                editor.stage['parameters']['prepare_complex']=bool(e.control.value)
                editor.expected=input_types(editor.stage,'base_input_path')
                for source,selector,_ in editor.rows:
                    source.options=editor.source_options()
                    if source.value not in {option.key for option in source.options}:
                        source.value='';selector.value='auto';selector.options=editor.selectors('')
                ui.page.update()
            self.field_controls['prepare_complex'].on_change=redocking_inputs
        if stage['operation']=='fingerprints':
            def choose(e=None):
                for key in ('radius','morgan_n_bits'):
                    self.field_controls[key].visible=self.fingerprint_choice.value=='morgan'
                    if key in self.conditional_cells:self.conditional_cells[key].visible=self.field_controls[key].visible
                if e is not None:ui.page.update()
            self.fingerprint_choice.on_select=choose;choose()
        if stage['operation']=='similarity':
            self.fingerprint_notice=ft.Text(size=12,color='#64748B')
            self.controls.append(self.fingerprint_notice)
            self.input_editors['base_input_path'].on_change=self.sync_fingerprint
            self.sync_fingerprint()
        if stage['operation']=='graphs':
            self.graph_notice=ft.Text(size=12,color='#64748B')
            self.controls.append(self.graph_notice)
            self.sync_graph_mode()
        for name in template_names(stage['operation']):
            if stage['operation']=='redocking' and name.startswith('chimera/'):continue
            self.template(name)
        if stage['operation']=='retrieve_compounds':
            self.field_controls['search_mode'].on_change=self.sync_retrieval
            self.field_controls['expand_chembl'].on_change=self.sync_retrieval
            self.sync_retrieval()
    def layout(self):
        """Group fields into roomy cards; expansion content has explicit spacing."""
        from .branding import BRAND_BLUE
        def fields(controls):
            cells=[]
            for control in controls:
                if isinstance(control,(ft.TextField,ft.Dropdown)):
                    control.width=None
                long_label=len(str(getattr(control,'label','')))>34
                col=12 if control.data=='wide' or long_label or not isinstance(control,(ft.TextField,ft.Dropdown)) else {'xs':12,'md':6}
                cells.append(ft.Container(content=control,col=col,padding=ft.Padding.symmetric(vertical=6),visible=control.visible))
                # Hidden conditional fields also release their space in the form.
                conditional=('radius','morgan_n_bits') if self.stage['operation']=='fingerprints' else ('max_targets','expand_chembl','similarity_threshold') if self.stage['operation']=='retrieve_compounds' else ()
                if control in [self.field_controls.get(k) for k in conditional]:
                    key=next(k for k in conditional if self.field_controls[k] is control)
                    self.conditional_cells[key]=cells[-1]
            return ft.ResponsiveRow(cells,spacing=24,run_spacing=20)
        def card(title,description,controls):
            return ft.Container(padding=24,bgcolor='#FFFFFF',border_radius=16,border=ft.Border.all(1,'#E2E8F0'),content=ft.Column([
                ft.Text(title,size=17,weight=ft.FontWeight.W_600,color='#172B4D'),
                ft.Text(description,size=12,color='#64748B'),ft.Container(height=4),fields(controls)],spacing=12))
        sections=[card('Identificação','Nomeie a etapa e escolha se ela participa da execução.',self.controls[:2])]
        inputs=[c for c in self.controls[2:] if c.data=='input']
        if self.stage['operation']=='graphs':
            editors=[self.input_editors[k].control for k in ('similarity_path','base_input_path')]
            inputs=[c for c in inputs if c not in editors]+editors
        parameters=[c for c in self.controls[2:] if c.data!='input' and not isinstance(c,ft.ExpansionTile)]
        if inputs: sections.append(card('Entradas','Escolha a origem e o resultado que esta etapa deve receber.',inputs))
        if parameters: sections.append(card('Parâmetros','Ajuste as opções da etapa. Os filtros específicos estão abaixo.',parameters))
        for control in self.controls[2:]:
            if not isinstance(control,ft.ExpansionTile): continue
            control.controls=[ft.Container(padding=ft.Padding.only(left=24,right=24,top=20,bottom=28),content=fields(control.controls))]
            control.bgcolor=control.collapsed_bgcolor='#FFFFFF'
            control.tile_padding=24
            control.text_color=control.collapsed_text_color='#172B4D'
            control.icon_color=control.collapsed_icon_color=BRAND_BLUE
            cell=ft.Container(control,border_radius=16,border=ft.Border.all(1,'#E2E8F0'),visible=control.visible)
            for key,tile in self.filter_tiles.items():
                if tile is control:self.filter_cells[key]=cell
            sections.append(cell)
        return ft.Column(sections,spacing=24,scroll=ft.ScrollMode.AUTO,horizontal_alignment=ft.CrossAxisAlignment.STRETCH)

    def add(self,key,label,value,choices=None,into=None):
        controls=self.controls if into is None else into
        if isinstance(value,bool):
            control=ft.Switch(label=label,value=value,disabled=not self.writable); getter=lambda:bool(control.value)
        elif choices:
            from .app import KINDS
            control=ft.Dropdown(label=label,value=value,options=[ft.DropdownOption(key=str(v),text=MODES.get(v,str(v)) if key=='search_mode' else KINDS.get(v,str(v)) if key=='kind' else str(v)) for v in choices],disabled=not self.writable); getter=lambda:control.value
        else:
            control=ft.TextField(label=label,value='' if value is None else str(value),disabled=not self.writable)
            if key=='dock6_app_path': control.value=str(self.ui.service.dock6_path or ''); control.read_only=True; control.helper='Configuração da instalação feita pelo administrador.'
            integers={'threshold','pubchem_threshold','pubchem_max_records','max_records','morgan_n_bits','num_modes','exhaustiveness','chunk_size','plot_max_residues','max_targets','similarity_threshold'}
            if key=='radius' and self.stage['operation']=='fingerprints': integers.add('radius')
            kind=int if key in integers else float if key in NUMERIC else int if type(value)is int else float if type(value)is float else str
            getter=lambda:None if not control.value else kind(control.value)
        if key in ('search_term','target'): control.data='wide'
        controls.append(control); self.readers[key]=getter;self.field_controls[key]=control
    def sync_retrieval(self,e=None):
        mode=self.field_controls['search_mode'].value or 'target'
        direct=mode in ('molecule_id','molecule_name','similarity','substructure')
        for key,visible in (('max_targets',not direct),('expand_chembl',not direct),('similarity_threshold',mode=='similarity')):
            self.field_controls[key].visible=visible
            if key in self.conditional_cells:self.conditional_cells[key].visible=visible
        self.field_controls['search_term'].helper={
            'target':'Ex.: CHEMBL220, P00533 ou acetylcholinesterase',
            'target_id':'Ex.: CHEMBL220, CHEMBL240', 'uniprot':'Ex.: P00533, P22303',
            'molecule_id':'Ex.: CHEMBL25, CHEMBL50',
            'similarity':'SMILES ou um único ID ChEMBL de composto',
            'substructure':'SMILES do fragmento que os compostos devem conter',
        }.get(mode,'Digite o nome ou texto da consulta.')
        for group,tile in self.filter_tiles.items():
            visible=group=='molecules' or (not direct and (group!='similarmols' or self.field_controls['expand_chembl'].value))
            tile.visible=visible
            if group in self.filter_cells:self.filter_cells[group].visible=visible
        if e:self.ui.page.update()

    def sync_fingerprint(self):
        from biomolexplorer.fingerprint_selection import generated_kind,LABELS
        editor=self.input_editors['base_input_path']
        candidate=copy.deepcopy(self.stage)
        candidate['bindings']['base_input_path']=editor.read()
        if not editor.direct_path():candidate['parameters'].pop('base_input_path',None)
        control=self.field_controls['fingerprint']
        try:
            kind,custom=generated_kind(candidate,self.ui.current['pipeline'])
            if kind:control.value=kind
            control.disabled=not self.writable or kind is not None or not custom
            self.fingerprint_notice.value=('Identificado pela entrada: '+LABELS[kind]+'.' if kind else
                'Escolha o tipo de fingerprint contido nos seus arquivos.' if custom else 'Conecte uma origem de fingerprints.')
        except ValueError as exc:
            control.disabled=True;self.fingerprint_notice.value=str(exc)
    def sync_graph_mode(self):
        self.graph_notice.value='A métrica e o limiar são definidos em Calcular similaridade. Este bloco preserva as relações de cada entrada e calcula o MCC e seu fragmento comum.'
    def list_field(self,key,label,value,suggestions):
        values=list(value or []) if isinstance(value,list) else [value] if value else []
        available=list(dict.fromkeys(suggestions+values)); box=ft.Column(); checks={}
        def add_item(item,selected):
            if item in checks: checks[item].value=True; return
            checks[item]=ft.Checkbox(label=item,value=selected,disabled=not self.writable); box.controls.append(checks[item])
        for item in available: add_item(item,item in values)
        entry=ft.TextField(label='Adicionar item a '+label,disabled=not self.writable)
        def append(e):
            if entry.value and entry.value.strip(): add_item(entry.value.strip(),True); entry.value=''; self.ui.page.update()
        self.controls.append(ft.ExpansionTile(title=ft.Text(label),controls=[box,entry,ft.TextButton('Adicionar item',on_click=append,disabled=not self.writable)]))
        self.readers[key]=lambda:[k for k,c in checks.items() if c.value] or None
    def records(self,key,label,value):
        multi=key=='pdb_codes' or self.stage['operation']=='docking_vina'
        rows=[]; box=ft.Column()
        automatic=ft.Checkbox(label='Detectar complexos pelos metadados dos arquivos',value=value is None,disabled=not self.writable or self.stage['operation']=='docking_dock6')
        if self.stage['operation']=='docking_dock6': automatic.value=False
        def add(record=None):
            record=record or ['', '', '', 'A']; fields=[ft.TextField(label=l,value=str(v),width=w,disabled=not self.writable) for l,v,w in zip(['Código PDB','Ligante','Número do resíduo','Cadeia'],record[:4],[150,120,150,90])]
            fields.append(ft.TextField(label='Resolução (Å, opcional)',value=str(record[4]) if len(record)>4 and record[4] is not None else '',width=170,disabled=not self.writable))
            row=ft.Row(fields,wrap=True,spacing=20,run_spacing=24)
            def remove(e): rows.remove(fields); box.controls.remove(row); self.ui.page.update()
            row.controls.append(ft.IconButton(icon=ft.Icons.DELETE_OUTLINE,on_click=remove,disabled=not self.writable)); rows.append(fields); box.controls.append(row)
        records=(value if isinstance(value[0],(list,tuple)) else [value]) if value else []
        for record in records: add(record)
        def append(e): add(); automatic.value=False; self.ui.page.update()
        self.controls.append(ft.ExpansionTile(title=ft.Text(label),controls=[automatic,box,ft.TextButton('Adicionar complexo',on_click=append,disabled=not self.writable)]))
        def read():
            if automatic.value: return None
            parsed=[[f[0].value.strip(),f[1].value.strip(),int(f[2].value),f[3].value.strip()]+([float(f[4].value)] if f[4].value else []) for f in rows]
            return parsed if multi else parsed[0] if parsed else None
        self.readers[key]=read
    def template(self,name):
        source=self.stage.get('templates',{}).get(name,(RESOURCE_ROOT/name).read_text()); controls=[]
        if name.endswith('.json'):
            group={'target':'target','bioactivity':'bioactivity','molecules':'molecules','similarmols':'similars'}[name.split('/')[-1][:-5]]
            data=copy.deepcopy(self.stage['parameters'].get('chembl_filters',{}).get(group,json.loads(source))); getters={}
            if group=='target':data.setdefault('organism','')
            if group=='bioactivity':
                data.setdefault('standard_type__in',[]);data.setdefault('assay_type','')
            for key,value in data.items():
                if key=='standard_type__in':
                    from .activity_measures import ActivityMeasures
                    picker=ActivityMeasures(self.ui.page,value or [],self.writable)
                    field=picker.control;getter=picker.values
                elif key=='organism':
                    organisms=list(dict.fromkeys(LISTS['organism']+['Escherichia coli','Saccharomyces cerevisiae','Danio rerio']+([value] if value else [])))
                    field=ft.Dropdown(label='Organismo',value=value or '',hint_text='Ex.: Homo sapiens',
                        helper_text='Nome científico do organismo do alvo. Selecione ou digite; vazio aceita qualquer organismo.',
                        editable=True,enable_filter=True,enable_search=True,
                        options=[ft.DropdownOption(key='',text='Qualquer')]+[ft.DropdownOption(key=v,text=v) for v in organisms],disabled=not self.writable)
                    typed=[None]
                    field.on_text_change=lambda e,t=typed:t.__setitem__(0,e.control.text or '')
                    field.on_select=lambda e,t=typed:t.__setitem__(0,None)
                    def getter(c=field,t=typed):
                        text=(t[0] if t[0] is not None else c.value or '').strip()
                        return '' if text==c.options[0].text else text
                elif key=='assay_type':
                    choices=[('B','Ligação'),('F','Funcional'),('A','ADMET'),('T','Toxicidade'),('P','Físico-químico'),('U','Não atribuído')]
                    field=ft.Dropdown(label='Tipo de ensaio',value=value or '',helper_text='Classificação do ensaio na ChEMBL.',
                        options=[ft.DropdownOption(key='',text='Qualquer')]+[ft.DropdownOption(key=k,text=k+' - '+label) for k,label in choices],disabled=not self.writable)
                    getter=lambda c=field:c.value or None
                elif isinstance(value,list) and key=='type__in':
                    suggested=['SINGLE PROTEIN','PROTEIN FAMILY','PROTEIN COMPLEX','CELL-LINE','TISSUE','ORGANISM'] if key=='type__in' else ['Ki','IC50','EC50','Kd']
                    checks=[ft.Checkbox(label=v,value=v in value,disabled=not self.writable,data=v) for v in dict.fromkeys(suggested+value)]
                    field=ft.Column([ft.Text(FILTER_LABELS.get(key,key)),ft.Row(checks,wrap=True,spacing=20,run_spacing=16)])
                    getter=lambda checks=checks:[c.data for c in checks if c.value]
                elif key in ('natural_product','pchembl_value__isnull') and value in (None,0,1):
                    field=ft.Dropdown(label='Produto natural' if key=='natural_product' else 'Exigir pChEMBL',value='' if value is None else str(value),options=[ft.DropdownOption(key='',text='Qualquer')]+[ft.DropdownOption(key=str(k),text=('Sim' if k==1 else 'Não') if key=='natural_product' else ('Sim' if k==0 else 'Não')) for k in (0,1)],disabled=not self.writable)
                    getter=lambda c=field:int(c.value) if c.value else None
                elif isinstance(value,list):
                    field=ft.TextField(label=FILTER_LABELS.get(key,key)+' (itens separados por vírgula)',value=', '.join(map(str,value)),disabled=not self.writable)
                    getter=lambda c=field:[v.strip() for v in c.value.split(',') if v.strip()]
                elif isinstance(value,(dict,tuple)): continue
                else:
                    field=ft.Switch(label=FILTER_LABELS.get(key,key),value=value,disabled=not self.writable) if isinstance(value,bool) else ft.TextField(label=FILTER_LABELS.get(key,key),value='' if value is None else str(value),disabled=not self.writable)
                    kind=int if key=='similarity' else float if key in ('molecule_weight','max_value_ref') else type(value) if value is not None else str
                    getter=lambda c=field,k=kind:k(c.value) if c.value not in (None,'') else None
                controls.append(field); getters[key]=getter
            def read(data=data,getters=getters):
                result=copy.deepcopy(data)
                for k,g in getters.items():
                    value=g()
                    if value is None or value=='' or value==[]: result.pop(k,None)
                    else: result[k]=value
                return json.dumps(result,indent=2,ensure_ascii=False)
            self.template_readers[name]=read
        elif name.startswith('chimera/'):
            lines=source.splitlines(keepends=True); toggles={}
            for command,label in [('addh','Adicionar hidrogênios'),('minimize','Minimizar energia'),('delete','Remover solvente / seleções do protocolo')]:
                indices=[i for i,line in enumerate(lines) if line.strip().startswith(command+' ') or line.strip()==command]
                if indices and all('{' not in lines[i] for i in indices):
                    c=ft.Switch(label=label,value=True,disabled=not self.writable); controls.append(c); toggles[command]=(c,indices)
            charges={}
            for i,line in enumerate(lines):
                match=re.search(r'\bmethod (gas|am1)\b',line)
                if match:
                    c=ft.Dropdown(label='Método de cargas do receptor',value=match.group(1),options=[ft.DropdownOption(key=v,text=v) for v in ('gas','am1')],disabled=not self.writable)
                    controls.append(c); charges[i]=c
            def read_chimera(lines=lines,t=toggles,charges=charges):
                return ''.join(re.sub(r'\bmethod (gas|am1)\b','method '+charges[i].value,line) if i in charges else line for i,line in enumerate(lines) if not any(i in ids and not c.value for c,ids in t.values()))
            self.template_readers[name]=read_chimera
            # Preserve markers and structural commands when producing the protocol.
        else:
            lines=source.splitlines(keepends=True); readers={}
            for i,line in enumerate(lines):
                stripped=line.strip()
                if not stripped or stripped.startswith('#') or '{' in line: continue
                if re.match(r'^(verbose|verbosity)\b',stripped,re.I):
                    lines[i]=re.sub(r'^(verbose|verbosity)(\s*=\s*|\s+).*',r'\1\g<2>0',stripped,flags=re.I)+'\n'
                    continue
                parts=re.split(r'\s*=\s*|\s+',stripped,maxsplit=1)
                positional=name in ('dock6/showbox.template','dock6/INSPH.template')
                if positional:
                    try: float(stripped)
                    except ValueError: continue
                    key=('Margem da caixa (Å)' if name.endswith('showbox.template') and i==1 else {3:'Raio mínimo (Å)',4:'Raio máximo (Å)',5:'Raio da esfera (Å)'}.get(i,'Parâmetro '+str(i+1)))
                    value=stripped
                else:
                    if len(parts)!=2: continue
                    key,value=parts
                if value.lower() in ('yes','no'):
                    c=ft.Switch(label=TEMPLATE_FIELDS.get(key,key.replace('_', ' ')),value=value.lower()=='yes',disabled=not self.writable); getter=lambda c=c:'yes' if c.value else 'no'
                else:
                    c=ft.TextField(label=TEMPLATE_FIELDS.get(key,key.replace('_', ' ')),value=value,disabled=not self.writable)
                    try: float(value); numeric=True
                    except ValueError: numeric=False
                    def getter(c=c,numeric=numeric):
                        if numeric: float(c.value)
                        if '\n' in c.value or '\r' in c.value: raise ValueError('Use um valor por configuração.')
                        return c.value
                controls.append(c); readers[i]=('' if positional else key,getter,'' if positional else ' = ' if '=' in line else ' ')
            def read(lines=lines,readers=readers):
                result=lines[:]
                for i,(key,getter,sep) in readers.items(): result[i]=key+sep+getter()+'\n'
                return ''.join(result)
            self.template_readers[name]=read
        if controls:
            tile=ft.ExpansionTile(title=ft.Text(TEMPLATE_LABELS.get(name.split('/')[-1].split('.')[0],name)),controls=controls)
            self.controls.append(tile)
            if name.startswith('crawlers/'):
                self.filter_tiles[name.split('/')[-1].split('.')[0]]=tile
    def read(self):
        result=copy.deepcopy(self.stage); result['name']=self.title.value or TITLES[result['operation']][0]; result['enabled']=bool(self.enabled.value)
        if 'process_all' in result:result['process_all']=False
        if hasattr(self,'processing'):result['input_processing']=self.processing.value
        for key,getter in self.readers.items():
            if result['operation']=='fingerprints' and key in ('radius','morgan_n_bits') and self.fingerprint_choice.value!='morgan':
                result['parameters'].pop(key,None);continue
            value=getter()
            if value is None: result['parameters'].pop(key,None)
            else: result['parameters'][key]=value
        for key,getter in self.binding_readers.items():
            direct_path=self.direct_path_readers[key]()
            if direct_path is not None:
                result['parameters'][key]=direct_path
                result['bindings'].pop(key,None)
                continue
            value=getter(); result['parameters'].pop(key,None)
            if value: result['bindings'][key]=value
            else: result['bindings'].pop(key,None)
        if self.derived_consensus_input:
            result['parameters'].pop('base_input_path',None)
            result['bindings'].pop('base_input_path',None)
        for name,getter in self.template_readers.items():
            value=getter()
            original=(RESOURCE_ROOT/name).read_text()
            different=json.loads(value)!={k:v for k,v in json.loads(original).items() if v is not None and v!='' and v!=[]} if name.endswith('.json') else value!=original
            if different: result.setdefault('templates',{})[name]=value
            else: result.get('templates',{}).pop(name,None)
        # Resource-based guided filters replace the corresponding parameter filters.
        if any(n.startswith('crawlers/') for n in self.template_readers): result['parameters'].pop('chembl_filters',None)
        result['parameters'].pop('verbose',None)
        if hasattr(self,'provided'):
            provided=self.provided.read()
            if provided:result['provided_results']=provided
            else:result.pop('provided_results',None)
        if hasattr(self,'redocking_pairs') and not result.get('provided_results'):
            records,settings=self.redocking_pairs.read()
            result['parameters']['pdb_codes']=records
            result['parameters']['preparation_pairs']=settings
            result['parameters'].pop('charge_type',None)
        validate_templates(result.get('templates',{}))
        if result['operation']=='similarity':
            from biomolexplorer.fingerprint_selection import generated_kind
            from biomolexplorer.bindings import sources
            refs=sources(result['bindings']['base_input_path']) if result['bindings'].get('base_input_path') else []
            if result.get('input_processing')!='individual' or len(refs)<=1:
                kind,_=generated_kind(result,self.ui.current['pipeline'])
                if kind:result['parameters']['fingerprint']=kind
        return result
