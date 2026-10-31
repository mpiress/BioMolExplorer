"""Explicit per-run file choices at pipeline stage boundaries."""
import copy
from types import SimpleNamespace
from pathlib import Path

import flet as ft

from biomolexplorer.artifacts import records as artifact_records, reference as artifact_reference
from .artifact_picker import ArtifactPicker
from biomolexplorer.artifact_choices import matches_selector, selector_for, produced_files
from biomolexplorer.bindings import sources, pack
from biomolexplorer.catalog import field_label
from biomolexplorer.docking_inputs import target_input, prepared_receptor, preferred_receptors, ligand_input_file
from biomolexplorer.flow import INPUTS, EXTERNAL_INPUTS, input_types
from biomolexplorer.input_validation import columns
from .localization import verbatim, stage_control


def compatible_file(path, kinds, headers=None):
    path = Path(path)
    def file_columns(path):
        return headers[str(path)] if headers is not None and str(path) in headers else columns(path)
    required = {'compounds': {'canonical_smiles'}, 'chembl': {'canonical_smiles'},
                'fingerprints': {'fingerprint', 'molecule_chembl_id'},
                'similarity': {'source', 'target', 'value'}}
    for kind in kinds:
        if kind=='scores' and path.suffix.lower()=='.csv' and {'molecule_chembl_id','vina','dock6'}<=file_columns(path):return True
        if kind=='receptor_pdb' and path.suffix.lower()=='.pdb':return True
        if kind=='molecular_poses' and path.suffix.lower() in ('.mol2','.pdbqt','.sdf') and ligand_input_file(path):return True
        if kind=='dynamics' and path.suffix.lower() in ('.nc','.prmtop','.rst7','.json','.csv'):return True
        if kind in required and path.suffix.lower() == '.csv' and required[kind] <= file_columns(path):
            return True
        if kind in ('vina','dock6') and path.suffix=='.csv' and {'engine','score','conformer_file'}<=file_columns(path):return True
        if kind == 'structures' and path.suffix.lower() in ('.pdb','.cif','.mmcif') and path.name.count('.') == 1:
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


def selection_dialog_size(page,redocking=False):
    width=page.width or 1160;height=page.height or 900
    reserve=260 if redocking and width<700 else 220 if redocking else 250
    return dict(width=max(220,min(1180 if redocking else 860,width-100)),
                height=max(140,min(760 if redocking else 580,height-reserve)))


def selection_file_state(paths):
    """Cheap freshness check; never hashes or mutates scientific outputs."""
    result=[]
    for value in paths:
        path=Path(value)
        try:
            stat=path.stat();stamp=(stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns)
        except OSError:stamp=None
        result.append((str(path),stamp))
    return result


def prepare_selection_data(store,token,run,pending,assets):
    """Prepare only data in a worker thread; no dialogs, controls or writes."""
    operation=pending['configuration']['operation'];pid=run['project_id']
    paths=list(dict.fromkeys(p for item in run['stages'] for p in item.get('artifacts',[])))
    asset_paths={asset['id']:store.asset_path(token,pid,asset['id']) for asset in assets
                 if operation in ('retrieve_pubchem','docking_vina','docking_dock6') and asset['kind'] in ('compounds','other')}
    paths+=list(asset_paths.values());before=selection_file_state(paths)
    data=dict(assets=assets,paths=paths,files_before=before,records={},metadata={},reviews={},compound_references=None,outputs={})
    from biomolexplorer.result_origins import ResultOrigins
    data['origins']=ResultOrigins(run['stages'],run.get('id'))
    data['headers']={str(p):columns(p) for p in paths if Path(p).suffix.lower()=='.csv'}
    for item in run['stages']:
        if item['status']!='succeeded':continue
        data['records'][item['id']]=artifact_records(item,run.get('id'))
        data['outputs'][item['id']]=produced_files(item)
        for record in data['records'][item['id']]:
            data['origins'].describe(record['path'],item['id'],run.get('id'))
        if operation=='redocking':
            from biomolexplorer.redocking_config import metadata_records
            data['metadata'][item['id']]=metadata_records(item.get('artifacts',[]))
    if operation=='redocking':
        from biomolexplorer.result_files import ResultFiles
        data['reviews']=ResultFiles(store).pdb_review_statuses(token,pid,run['id'])
    if operation in ('retrieve_pubchem','docking_vina','docking_dock6'):
        from biomolexplorer.pubchem_retrieval import reference_choices
        from biomolexplorer.docking_inputs import compound_identifiers
        read_references=reference_choices if operation=='retrieve_pubchem' else compound_identifiers
        connected={r.get('stage') for r in sources(pending['configuration'].get('bindings',{}).get('base_selected_mols',{}))}
        references={};data['compound_references']=references
        for item in run['stages']:
            if item['status']!='succeeded' or operation!='retrieve_pubchem' and item['id'] not in connected:continue
            for filename in data['outputs'].get(item['id'],[]):
                path=store.scoped_path(pid,filename)
                if path.suffix.lower()=='.csv':references[str(path)]=read_references([path])
        for identifier,path in asset_paths.items():references['asset:'+identifier]=read_references([path])
    data['file_state']=selection_file_state(paths)
    return data


class FileSelection:
    def __init__(self, run, pending, assets=(), state=None, file_actions=None, compound_references=None, page=None, writable=True, prepared=None):
        self.prepared=prepared or {}
        self.stage = pending['configuration']
        from biomolexplorer.result_origins import ResultOrigins
        origin_labels=self.prepared.get('origins') or ResultOrigins(run['stages'],run.get('id'))
        self.rows = {}
        self.pickers = {}
        self.extras = {}
        self.reference_records={}
        self.pdb_reviews={}
        self.redocking_metadata={}
        self.compound_fields={}
        bindings=self.stage.get('bindings',{})
        self.derived_consensus=(self.stage['operation']=='consensus' and 'base_vina_path' in bindings
            and bindings.get('base_input_path')==bindings['base_vina_path'])
        self.processing=ft.Dropdown(label='Como processar os arquivos selecionados?',
            disabled=not writable,value=(state or {}).get('input_processing',self.stage.get('input_processing','individual')),options=[
                ft.DropdownOption(key='individual',text='Processar individualmente'),
                ft.DropdownOption(key='merge',text='Mesclar arquivos (merge)'),
                ft.DropdownOption(key='both',text='Individual + conjunto completo (ALL)')])
        if self.stage['operation'] in ('consensus','dynamics_dock6'):
            self.processing.value='merge'
            self.processing.disabled=True
            self.processing.options=[ft.DropdownOption(key='merge',text='Reunir as entradas selecionadas')]
        asset_names={a['id']:a['name'] for a in assets}
        completed = {s['id']: s for s in run['stages'] if s['status'] == 'succeeded'}
        controls = [ft.Text('Selecione os arquivos que serão usados em “' + pending['name'] + '”. '
                            'A execução continuará após sua confirmação.'),
                    self.processing,
                    ft.Text('Individual: cada arquivo gera resultados separados. Merge: as entradas são reunidas. '
                            'Com várias entradas, o modo individual processa cada combinação de arquivos. '
                            'Individual + ALL acrescenta a análise do conjunto completo e reaproveita lotes já calculados. '
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
                    if file_actions and self.stage['operation']=='dynamics_dock6':self.extras[id(check)]=file_actions({'asset_id':ref['asset']},Path(asset_names.get(ref['asset'],'')))

                    self.add_docking_choice(field,check,ref,(compound_references or {}).get('asset:'+ref['asset'],{}),controls,page)
                    continue
                if ref['stage'] in seen:
                    continue
                seen.add(ref['stage'])
                item = completed.get(ref['stage'])
                if item is None:continue
                files = list(dict.fromkeys(Path(p) for p in item['artifacts']))
                if self.stage['operation']=='redocking':
                    from biomolexplorer.redocking_config import metadata_records
                    self.redocking_metadata[item['id']]=(self.prepared['metadata'][item['id']]
                        if item['id'] in self.prepared.get('metadata',{}) else metadata_records(files))
                records=(self.prepared['records'][item['id']] if item['id'] in self.prepared.get('records',{})
                         else artifact_records(item,run.get('id')))
                catalog={r['path']:r for r in records}
                outputs=(self.prepared['outputs'][item['id']] if item['id'] in self.prepared.get('outputs',{}) else produced_files(item))
                matching = [p for p in outputs if compatible_file(p,kinds,self.prepared.get('headers'))]
                if item.get('operation')=='redocking':
                    from biomolexplorer.artifact_choices import redocking_receptors
                    allowed=set(redocking_receptors(files)) if kinds & {'structures','prepared_structures','receptor_pdb'} else set()
                    matching=[p for p in matching if p in allowed]
                if self.stage['operation']=='dynamics_dock6':
                    if field=='base_input_path':
                        from biomolexplorer.md_inputs import receptor_input_files
                        matching=receptor_input_files(matching)
                    else:
                        from biomolexplorer.molecular_roles import small_molecule_file
                        matching=[p for p in matching if small_molecule_file(p)]
                if item.get('operation') in ('docking_vina','docking_dock6','consensus') and not item.get('configuration',{}).get('provided_results') and not kinds<={'structures','receptor_pdb','prepared_structures'}:
                    from biomolexplorer.artifact_choices import consolidated_docking_tables
                    allowed=set(consolidated_docking_tables(files,item['operation']))
                    matching=[p for p in matching if p in allowed]
                if item.get('operation') in ('retrieve_compounds','expand_similar_compounds') and not item.get('configuration',{}).get('provided_results'):
                    from biomolexplorer.artifact_choices import choices
                    allowed=choices(item)
                    if allowed:matching=[p for p in matching if any(matches_selector(p,s) for s in allowed)]
                if self.stage['operation'] in ('prepare_structures','docking_vina','docking_dock6') and field=='base_input_path' and any(prepared_receptor(p) for p in files):
                    matching=preferred_receptors(p for p in matching if prepared_receptor(p))
                if self.stage['operation'] in ('docking_vina','docking_dock6') and field=='base_selected_mols':
                    matching=[p for p in matching if ligand_input_file(p) and
                        (p.suffix!='.csv' or {'molecule_chembl_id','canonical_smiles'}<=columns(p))]
                if item.get('operation') in ('docking_vina','docking_dock6'):
                    tables=[p for p in matching if p.name=='docking_results.csv']
                    if tables:
                        from biomolexplorer.docking_data import result_tables
                        matching=result_tables(tables)
                if item.get('operation')=='retrieve_pubchem' and not item.get('configuration',{}).get('provided_results'):
                    matching=[p for p in matching if (p.name=='compounds.csv' or '__pubchem__compounds' in p.stem) and p.parent.parent.name=='compounds']
                if 'prepared_structures' in kinds and any(p.parent.name == 'Prepared' for p in matching):
                    matching = [p for p in matching if p.parent.name == 'Prepared']
                if item.get('operation')=='graphs' and not item.get('configuration',{}).get('provided_results'):
                    from biomolexplorer.graph_exports import molecular_files
                    matching=molecular_files(matching)
                if item.get('operation')=='admet':
                    allowed=set(produced_files(item))
                    matching=[p for p in matching if p in allowed]
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
                    scope=''
                    if item.get('operation')=='graphs':
                        from biomolexplorer.graph_exports import subset_label
                        scope=subset_label(path,getattr(page,'translator',None))
                    from biomolexplorer.artifact_naming import human_name
                    title=scope or human_name(record)
                    if item.get('operation')=='admet' and not title.startswith('ADMET ·'):title='ADMET · '+title
                    count=record.get('record_count')
                    display=title+' · '+description['origin']+(f' · {count:,} registros' if count is not None else '')
                    check = SimpleNamespace(label=selector,display_label=display,physical_file=str(path),
                        value=chosen,on_change=None,available=path.is_file(),
                        target_label=record.get('context',{}).get('target') or record.get('context',{}).get('search_term') or
                            item.get('configuration',{}).get('parameters',{}).get('target') or item.get('name',item['id']),
                        origin_label=(record.get('block_id') or item['id'][:8])+' · '+item.get('name',item['id']))
                    reference=artifact_reference(record)
                    previous=next((r for r in configured if r.get('stage')==item['id'] and
                        (r.get('artifact_id')==record['artifact_id'] if r.get('artifact_id') else matches_selector(path,r.get('selector','auto')))),{})
                    if previous.get('compound_id'):reference['compound_id']=previous['compound_id']
                    rows.append((check,reference))
                    if self.stage['operation']=='retrieve_pubchem':
                        from biomolexplorer.pubchem_retrieval import reference_choices
                        self.reference_records[(item['id'],selector)]=(compound_references or {}).get(str(path),{}) if compound_references is not None else reference_choices([path])
                    if file_actions and self.stage['operation']=='redocking' and field=='base_input_path':
                        actions=lambda item=item,path=path:file_actions(item,path)
                    else:actions=file_actions(item,path) if file_actions else []
                    self.extras[id(check)]=actions
                    records={}
                    if self.stage['operation'] in ('docking_vina','docking_dock6') and field=='base_selected_mols':
                        from biomolexplorer.docking_inputs import compound_identifiers
                        records=(compound_references or {}).get(str(path),{}) if compound_references is not None else compound_identifiers([path])
                    self.add_docking_choice(field,check,reference,records,controls,page)
                if not matching:
                    controls.append(ft.Text('Esta origem não produziu arquivos compatíveis com esta entrada.', color='#B91C1C'))
            for saved in configured:
                if saved.get('selector','auto')=='auto':continue
                if not any(flag.value and all((True if saved.get('artifact_id') else matches_selector(ref.get('selector','auto'),v)) if k=='selector' else ref.get(k)==v for k,v in saved.items() if k!='compound_id') for flag,ref in rows):
                    flag=SimpleNamespace(label=saved.get('selector',asset_names.get(saved.get('asset'),'Arquivo indisponível')),value=True,on_change=None,available=False)
                    producer=completed.get(saved.get('stage'),{})
                    if producer.get('operation') in ('docking_vina','docking_dock6','consensus'):
                        # Keep a previously confirmed pose usable without offering
                        # intermediate poses again in the primary source catalog.
                        from biomolexplorer.artifacts import validate_reference
                        try:
                            path=Path(validate_reference(saved,producer.get('artifacts',[])))
                            if path.is_file() and compatible_file(path,kinds) and (self.stage['operation']!='dynamics_dock6' or field=='base_input_path' or ligand_input_file(path)):
                                flag.available=True;flag.legacy_selection=True;flag.physical_file=str(path)
                                self.extras[id(flag)]=file_actions(producer,path) if file_actions else []
                        except (OSError,ValueError,TypeError):pass
                    rows.append((flag,copy.deepcopy(saved)))
            self.rows[field] = rows
            origins={}
            for ref in sources(group):
                if 'stage' in ref:
                    origin=next((s for s in run['stages'] if s['id']==ref['stage']),{})
                    origins['stage:'+ref['stage']]=origin.get('name',ref['stage'])
                elif 'asset' in ref:origins['asset:'+ref['asset']]=asset_names.get(ref['asset'],'Arquivo enviado')
            picker_class=ArtifactPicker
            if self.stage['operation']=='dynamics_dock6':
                from .molecular_files import MolecularFilePicker
                picker_class=MolecularFilePicker
            if self.stage['operation']=='redocking' and field=='base_input_path':
                from .redocking_files import RedockingFilePicker
                picker_class=RedockingFilePicker
            self.pickers[field]=picker_class(rows,origins,page,self.extras,writable=writable)
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
        self.control = ft.Column(controls, spacing=12, scroll=ft.ScrollMode.AUTO,
            horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
        if self.stage['operation']=='redocking':
            self.processing.content_padding=ft.Padding.symmetric(horizontal=16,vertical=18)
            self.control.spacing=20
            self.control.horizontal_alignment=ft.CrossAxisAlignment.STRETCH
        for field,view in (state or {}).get('picker_views',{}).items():
            picker=self.pickers.get(field)
            if picker and hasattr(picker,'restore_view'):picker.restore_view(view)

    def apply_pdb_reviews(self,statuses):
        self.pdb_reviews.update(copy.deepcopy(statuses))
        for picker in self.pickers.values():
            if hasattr(picker,'apply_reviews'):picker.apply_reviews(statuses)

    def add_docking_choice(self,field,check,ref,records,controls,page):
        if self.stage['operation'] not in ('docking_vina','docking_dock6') or field!='base_selected_mols':return
        from .compound_choice import CompoundChoice
        choice=ft.Dropdown(label='Composto para docking (opcional)',value=ref.get('compound_id',''),
            enable_filter=True,enable_search=True,visible=bool(check.value),disabled=not records)
        picker=CompoundChoice(choice,page);picker.set_records(records)
        picker.control.visible=bool(check.value)
        self.compound_fields[id(check)]=choice
        def changed(e):
            choice.visible=picker.control.visible=bool(check.value)
            if e is not None and page:page.update()
        check.on_change=changed
        self.extras.setdefault(id(check),[]).append(picker.control)

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
        views={field:picker.view_state() for field,picker in self.pickers.items() if hasattr(picker,'view_state')}
        if views:result['picker_views']=views
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
        if stage['operation']=='redocking':
            from biomolexplorer.redocking_config import pair_key
            existing=stage['parameters'].get('pdb_codes') or []
            records=[]
            for check,ref in self.rows.get('base_input_path',[]):
                if not check.value:continue
                filename=getattr(check,'physical_file',ref.get('selector',check.label))
                code=Path(filename).name.split('.',1)[0].upper()
                review=self.pdb_reviews.get(filename)
                if review is not None:
                    chosen=[[review.get('pdb_id',code),r['LIGAND'],int(r['RESNUM']),r['CHAIN']]
                            + ([float(r.get('RESOLUTION') or review.get('resolution'))] if r.get('RESOLUTION') or review.get('resolution') else [])
                            for r in review.get('ligands',[]) if r.get('LIGAND') and r.get('RESNUM') and r.get('CHAIN')]
                else:
                    chosen=[r for r in existing if str(r[0]).upper()==code or f'{r[0]}_{r[3]}'.upper()==code]
                    if not chosen:
                        chosen=[r for r in self.redocking_metadata.get(ref.get('stage'),[])
                                if str(r[0]).upper()==code or f'{r[0]}_{r[3]}'.upper()==code]
                for record in chosen:
                    if pair_key(record) not in {pair_key(r) for r in records}:records.append(copy.deepcopy(record))
            stage['parameters']['pdb_codes']=records
            keys={pair_key(r) for r in records}
            stage['parameters']['preparation_pairs']={k:v for k,v in
                (stage['parameters'].get('preparation_pairs') or {}).items() if k in keys}
        if self.derived_consensus:
            if 'base_vina_path' in stage['bindings']:stage['bindings']['base_input_path']=copy.deepcopy(stage['bindings']['base_vina_path'])
            else:stage['bindings'].pop('base_input_path',None)
        return stage
