"""Explicit per-run file choices at pipeline stage boundaries."""
import copy
from types import SimpleNamespace
from pathlib import Path

import flet as ft

from biomolexplorer.artifacts import records as artifact_records, reference as artifact_reference
from .artifact_picker import ArtifactPicker
from biomolexplorer.artifact_choices import matches_selector, selector_for
from biomolexplorer.bindings import sources, pack
from biomolexplorer.catalog import field_label
from biomolexplorer.docking_inputs import target_input, prepared_receptor, ligand_input_file
from biomolexplorer.flow import INPUTS, EXTERNAL_INPUTS, input_types
from biomolexplorer.input_validation import columns
from .localization import verbatim, stage_control


def compatible_file(path, kinds):
    path = Path(path)
    required = {'compounds': {'canonical_smiles'}, 'chembl': {'canonical_smiles'},
                'fingerprints': {'fingerprint', 'molecule_chembl_id'},
                'similarity': {'source', 'target', 'value'}}
    for kind in kinds:
        if kind in required and path.suffix.lower() == '.csv' and required[kind] <= columns(path):
            return True
        if kind in ('vina','dock6') and path.suffix=='.csv' and {'engine','score','conformer_file'}<=columns(path):return True
        if kind == 'structures' and path.suffix.lower() == '.pdb' and path.name.count('.') == 1:
            return True
        if kind == 'prepared_structures' and prepared_receptor(path):
            return True
        if kind == 'vina' and path.suffix.lower() == '.pdbqt':
            return True
        if kind == 'dock6' and path.name.endswith('_scored.mol2'):
            return True
        if kind in ('other','zinc_urls'):
            from biomolexplorer.zinc_retrieval import LIST_SUFFIXES
            if path.suffix.lower() in LIST_SUFFIXES:return True
    return False


def file_selector(path, files):
    """Use the shortest suffix that identifies exactly one output."""
    return selector_for(path,files)


class FileSelection:
    def __init__(self, run, pending, assets=(), state=None, file_actions=None, compound_references=None, page=None, writable=True):
        self.stage = pending['configuration']
        from biomolexplorer.result_origins import ResultOrigins
        origin_labels=ResultOrigins(run['stages'],run.get('id'))
        self.rows = {}
        self.pickers = {}
        self.extras = {}
        self.reference_records={}
        self.compound_fields={}
        bindings=self.stage.get('bindings',{})
        self.derived_consensus=(self.stage['operation']=='consensus' and 'base_vina_path' in bindings
            and bindings.get('base_input_path')==bindings['base_vina_path'])
        self.processing=ft.Dropdown(label='Como processar os arquivos selecionados?',
            value=(state or {}).get('input_processing',self.stage.get('input_processing','individual')),options=[
                ft.DropdownOption(key='individual',text='Processar individualmente'),
                ft.DropdownOption(key='merge',text='Mesclar arquivos (merge)')])
        asset_names={a['id']:a['name'] for a in assets}
        completed = {s['id']: s for s in run['stages'] if s['status'] == 'succeeded'}
        controls = [ft.Text('Selecione os arquivos que serão usados em “' + pending['name'] + '”. '
                            'A execução continuará após sua confirmação.'),
                    self.processing,
                    ft.Text('Individual: cada arquivo gera resultados separados. Merge: as entradas são reunidas. '
                            'Com várias entradas, o modo individual processa cada combinação de arquivos. '
                            'Metadados necessários acompanham as estruturas escolhidas.',
                            size=12, color='#64748B')]
        if self.stage['operation'] in ('similarity','graphs'):
            self.follow_inputs=ft.Switch(label='Acompanhar novos conjuntos da origem',
                value=(state or {}).get('follow_inputs',self.stage.get('follow_inputs',True)),disabled=not writable)
            controls.append(self.follow_inputs)
        for field, group in self.stage.get('bindings', {}).items():
            if self.derived_consensus and field=='base_input_path':continue
            kinds = input_types(self.stage,field) or EXTERNAL_INPUTS.get(self.stage['operation'], {}).get(field, {'other'})
            if self.stage['operation']=='consensus' and field=='base_input_path':kinds={'vina'}
            rows = []
            controls.append(ft.Text(field_label(self.stage['operation'],field), weight=ft.FontWeight.W_600))
            seen = set()
            selected = (state or {}).get('selected',{}).get(field)
            configured = sources(group) if selected is None else selected
            for ref in sources(group):
                if 'asset' in ref:
                    if 'asset:'+ref['asset'] in seen:continue
                    seen.add('asset:'+ref['asset'])
                    if target_input(self.stage,field) and not prepared_receptor(asset_names.get(ref['asset'],'')):
                        continue
                    check = SimpleNamespace(label=asset_names.get(ref['asset'],'Arquivo enviado'),on_change=None,available=ref['asset'] in asset_names,
                        value=any(r.get('asset')==ref['asset'] for r in configured))
                    verbatim(check,'label')
                    rows.append((check, copy.deepcopy(ref)))

                    self.add_docking_choice(field,check,ref,(compound_references or {}).get('asset:'+ref['asset'],{}),controls,page)
                    continue
                if ref['stage'] in seen:
                    continue
                seen.add(ref['stage'])
                item = completed.get(ref['stage'])
                if item is None:continue
                files = list(dict.fromkeys(Path(p) for p in item['artifacts']))
                catalog={r['path']:r for r in artifact_records(item,run.get('id'))}
                matching = [p for p in files if compatible_file(p, kinds)]
                if self.stage['operation'] in ('prepare_structures','docking_vina','docking_dock6') and field=='base_input_path' and any(prepared_receptor(p) for p in files):
                    matching=[p for p in matching if prepared_receptor(p)]
                if self.stage['operation'] in ('docking_vina','docking_dock6') and field=='base_selected_mols':
                    matching=[p for p in files if ligand_input_file(p) and
                        (p.suffix!='.csv' or {'molecule_chembl_id','canonical_smiles'}<=columns(p))]
                if item.get('operation') in ('docking_vina','docking_dock6'):
                    tables=[p for p in matching if p.name=='docking_results.csv']
                    if tables:
                        from biomolexplorer.docking_data import result_tables
                        matching=result_tables(tables)
                if item.get('operation')=='retrieve_pubchem' and not item.get('configuration',{}).get('provided_results'):
                    matching=[p for p in matching if p.name=='compounds.csv' and p.parent.parent.name=='compounds']
                if 'prepared_structures' in kinds and any(p.parent.name == 'Prepared' for p in matching):
                    matching = [p for p in matching if p.parent.name == 'Prepared']
                for path in sorted(matching):
                    selector = catalog[str(path)]['selector']
                    # Preserve explicit choices from the popup or block settings.
                    # An automatic connection still starts with no files checked.
                    selectors = [r['selector'] for r in configured if r.get('stage')==ref['stage']
                        and r.get('selector','auto')!='auto']
                    chosen = any(matches_selector(path,choice)
                        and sum(matches_selector(p,choice) for p in files)==1 for choice in selectors)
                    record=catalog[str(path)]
                    chosen=any(r.get('artifact_id')==record['artifact_id'] or (not r.get('artifact_id') and chosen and matches_selector(path,r.get('selector','auto'))) for r in configured if r.get('stage')==item['id'])
                    description=origin_labels.describe(path,item['id'],run.get('id'))
                    check = SimpleNamespace(label=selector,display_label=description['origin']+' — '+selector,
                        value=chosen,on_change=None,available=path.is_file())
                    reference=artifact_reference(record)
                    previous=next((r for r in configured if r.get('stage')==item['id'] and
                        (r.get('artifact_id')==record['artifact_id'] if r.get('artifact_id') else matches_selector(path,r.get('selector','auto')))),{})
                    if previous.get('compound_id'):reference['compound_id']=previous['compound_id']
                    rows.append((check,reference))
                    if self.stage['operation']=='retrieve_pubchem':
                        from biomolexplorer.pubchem_retrieval import reference_choices
                        self.reference_records[(item['id'],selector)]=(compound_references or {}).get(str(path),{}) if compound_references is not None else reference_choices([path])
                    actions=file_actions(item,path) if file_actions else []
                    self.extras[id(check)]=actions
                    records={}
                    if self.stage['operation'] in ('docking_vina','docking_dock6') and field=='base_selected_mols':
                        from biomolexplorer.pubchem_retrieval import reference_choices
                        records=(compound_references or {}).get(str(path),{}) if compound_references is not None else reference_choices([path])
                    self.add_docking_choice(field,check,reference,records,controls,page)
                if not matching:
                    controls.append(ft.Text('Esta origem não produziu arquivos compatíveis com esta entrada.', color='#B91C1C'))
            for saved in configured:
                if saved.get('selector','auto')=='auto':continue
                if not any(flag.value and all((True if saved.get('artifact_id') else matches_selector(ref.get('selector','auto'),v)) if k=='selector' else ref.get(k)==v for k,v in saved.items() if k!='compound_id') for flag,ref in rows):
                    rows.append((SimpleNamespace(label=saved.get('selector',asset_names.get(saved.get('asset'),'Arquivo indisponível')),value=True,on_change=None,available=False),copy.deepcopy(saved)))
            self.rows[field] = rows
            origins={}
            for ref in sources(group):
                if 'stage' in ref:
                    origin=next((s for s in run['stages'] if s['id']==ref['stage']),{})
                    origins['stage:'+ref['stage']]=origin.get('name',ref['stage'])
                elif 'asset' in ref:origins['asset:'+ref['asset']]=asset_names.get(ref['asset'],'Arquivo enviado')
            self.pickers[field]=ArtifactPicker(rows,origins,page,self.extras,writable=writable)
            controls.append(self.pickers[field].control)
        if self.stage['operation']=='retrieve_pubchem':
            from biomolexplorer.catalog import CHOICE_LABELS
            self.selection_mode=ft.Dropdown(label='Referências para buscar similares',
                value=(state or {}).get('selection_mode',self.stage['parameters'].get('selection_mode','all')),
                options=[ft.DropdownOption(key=k,text=v) for k,v in CHOICE_LABELS['selection_mode'].items()])
            self.compound_id=ft.Dropdown(label='Composto de referência',enable_filter=True,enable_search=True,
                value=(state or {}).get('compound_id',self.stage['parameters'].get('compound_id','')),options=[])
            def refresh(e=None):
                available={}
                for rows in self.rows.values():
                    for check,ref in rows:
                        if check.value:
                            available.update(self.reference_records.get((ref.get('stage'),ref.get('selector')),{}))
                            if 'asset' in ref:available.update((compound_references or {}).get('asset:'+ref['asset'],{}))
                self.compound_id.options=[verbatim(ft.DropdownOption(key=k,text=k+' · '+v),'text') for k,v in sorted(available.items())]
                if self.compound_id.value not in available:self.compound_id.value=''
                single=self.selection_mode.value=='single'
                self.compound_id.visible=single;self.compound_id.disabled=not single
                if single:self.processing.value='merge'
                self.processing.disabled=single
                if e is not None and page:page.update()
            self.refresh_references=refresh
            self.selection_mode.on_select=refresh
            for rows in self.rows.values():
                for check,_ in rows:check.on_change=refresh
            refresh();controls.extend([self.selection_mode,self.compound_id])
        self.control = ft.Column(controls, spacing=12, scroll=ft.ScrollMode.AUTO)

    def add_docking_choice(self,field,check,ref,records,controls,page):
        if self.stage['operation'] not in ('docking_vina','docking_dock6') or field!='base_selected_mols':return
        choice=ft.Dropdown(label='Composto para docking (opcional)',value=ref.get('compound_id',''),
            options=[ft.DropdownOption(key='',text='Todos os compostos')]+[
                verbatim(ft.DropdownOption(key=code,text=code)) for code in sorted(records)],
            enable_filter=True,enable_search=True,visible=bool(check.value),disabled=not records)
        if choice.value and choice.value not in records:
            choice.options.append(verbatim(ft.DropdownOption(key=choice.value,text=choice.value)))
        self.compound_fields[id(check)]=choice
        def changed(e):
            choice.visible=bool(check.value)
            if e is not None and page:page.update()
        check.on_change=changed
        self.extras.setdefault(id(check),[]).append(choice)

    def selected_references(self,rows,validate=True):
        refs=[]
        for check,ref in rows:
            if not check.value:continue
            if validate and not getattr(check,'available',True):raise ValueError('Arquivo selecionado indisponível. Remova ou selecione novamente a origem.')
            ref=copy.deepcopy(ref)
            choice=self.compound_fields.get(id(check))
            if choice:
                ref.pop('compound_id',None)
                if choice.value:ref['compound_id']=choice.value
            refs.append(ref)
        return refs

    def state(self):
        """Keep unchecked rows as well as choices while the user edits settings."""
        result={'input_processing':self.processing.value,
                'selected':{field:self.selected_references(rows,validate=False)
                            for field,rows in self.rows.items()}}
        if hasattr(self,'selection_mode'):result.update(selection_mode=self.selection_mode.value,compound_id=self.compound_id.value)
        if hasattr(self,'follow_inputs'):result['follow_inputs']=bool(self.follow_inputs.value)
        return result

    def read(self, require_selection=True):
        stage = copy.deepcopy(self.stage)
        if hasattr(self,'selection_mode'):
            stage['parameters']['selection_mode']=self.selection_mode.value
            stage['parameters']['compound_id']=self.compound_id.value or ''
            if require_selection and self.selection_mode.value=='single' and not self.compound_id.value:
                raise ValueError('Selecione um composto de referência para a busca PubChem.')
        stage['input_processing']=self.processing.value
        if hasattr(self,'follow_inputs'):stage['follow_inputs']=bool(self.follow_inputs.value)
        for field, rows in self.rows.items():
            refs = self.selected_references(rows)
            if not refs:
                if stage['operation']=='graphs' and field=='base_input_path':
                    stage['bindings'].pop(field,None)
                    continue
                if require_selection:
                    raise ValueError('Selecione ao menos um arquivo para ' + field_label(self.stage['operation'],field) + '.')
                # Keep the upstream origins available in the settings form when
                # configuration is opened before any files have been selected.
                refs = [{'stage':ref['stage'],'selector':'auto'}
                        for ref in sources(stage['bindings'][field]) if 'stage' in ref]
            group = pack(refs)
            if group:stage['bindings'][field] = group
            else:stage['bindings'].pop(field,None)
        if self.derived_consensus:
            if 'base_vina_path' in stage['bindings']:stage['bindings']['base_input_path']=copy.deepcopy(stage['bindings']['base_vina_path'])
            else:stage['bindings'].pop('base_input_path',None)
        return stage
