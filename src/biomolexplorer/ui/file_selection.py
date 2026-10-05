"""Explicit per-run file choices at pipeline stage boundaries."""
import copy
from pathlib import Path

import flet as ft

from biomolexplorer.bindings import sources, pack
from biomolexplorer.catalog import LABELS
from biomolexplorer.flow import INPUTS, EXTERNAL_INPUTS, input_types
from biomolexplorer.input_validation import columns
from .localization import verbatim


def compatible_file(path, kinds):
    path = Path(path)
    required = {'compounds': {'canonical_smiles'}, 'chembl': {'canonical_smiles'},
                'fingerprints': {'fingerprint', 'molecule_chembl_id'},
                'similarity': {'source', 'target', 'value'}}
    for kind in kinds:
        if kind in required and path.suffix.lower() == '.csv' and required[kind] <= columns(path):
            return True
        if kind == 'structures' and path.suffix.lower() == '.pdb' and path.name.count('.') == 1:
            return True
        if kind == 'prepared_structures' and path.name.endswith('.dockprep.pdbqt'):
            return True
        if kind == 'vina' and path.suffix.lower() == '.pdbqt':
            return True
        if kind == 'dock6' and path.name.endswith('_scored.mol2'):
            return True
        if kind == 'other' and path.suffix.lower() == '.txt':
            return True
    return False


def file_selector(path, files):
    """Use the shortest suffix that identifies exactly one output."""
    for length in range(1, len(path.parts)):
        suffix = '/'.join(path.parts[-length:])
        if sum(p.as_posix().endswith('/' + suffix) for p in files) == 1:
            return suffix
    raise ValueError('Os resultados contêm caminhos de arquivo duplicados.')


class FileSelection:
    def __init__(self, run, pending, assets=(), state=None):
        self.stage = pending['configuration']
        self.rows = {}
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
        for field, group in self.stage.get('bindings', {}).items():
            if self.derived_consensus and field=='base_input_path':continue
            kinds = input_types(self.stage,field) or EXTERNAL_INPUTS.get(self.stage['operation'], {}).get(field, {'other'})
            if self.stage['operation']=='consensus' and field=='base_input_path':kinds={'vina'}
            rows = []
            controls.append(ft.Text(LABELS.get(field, field), weight=ft.FontWeight.W_600))
            seen = set()
            selected = (state or {}).get('selected',{}).get(field)
            configured = sources(group) if selected is None else selected
            for ref in sources(group):
                if 'asset' in ref:
                    check = ft.Checkbox(label=asset_names.get(ref['asset'],'Arquivo enviado'),
                        value=any(r.get('asset')==ref['asset'] for r in configured))
                    verbatim(check,'label')
                    rows.append((check, copy.deepcopy(ref)))
                    controls.append(check)
                    continue
                if ref['stage'] in seen:
                    continue
                seen.add(ref['stage'])
                item = completed[ref['stage']]
                files = list(dict.fromkeys(Path(p) for p in item['artifacts']))
                controls.append(verbatim(ft.Text(item['name'], size=13)))
                matching = [p for p in files if compatible_file(p, kinds)]
                if 'prepared_structures' in kinds and any(p.parent.name == 'Prepared' for p in matching):
                    matching = [p for p in matching if p.parent.name == 'Prepared']
                for path in sorted(matching):
                    selector = file_selector(path, files)
                    # Preserve explicit choices from the popup or block settings.
                    # An automatic connection still starts with no files checked.
                    selectors = [r['selector'] for r in configured if r.get('stage')==ref['stage']
                        and r.get('selector','auto')!='auto']
                    chosen = any(path.as_posix().endswith('/'+choice)
                        and sum(p.as_posix().endswith('/'+choice) for p in files)==1 for choice in selectors)
                    check = verbatim(ft.Checkbox(label=selector, value=chosen),'label')
                    rows.append((check, {'stage': item['id'], 'selector': selector}))
                    controls.append(check)
                if not matching:
                    controls.append(ft.Text('Esta origem não produziu arquivos compatíveis com esta entrada.', color='#B91C1C'))
            self.rows[field] = rows
        self.control = ft.Column(controls, spacing=12, scroll=ft.ScrollMode.AUTO)

    def state(self):
        """Keep unchecked rows as well as choices while the user edits settings."""
        return {'input_processing':self.processing.value,
                'selected':{field:[copy.deepcopy(ref) for check,ref in rows if check.value]
                            for field,rows in self.rows.items()}}

    def read(self, require_selection=True):
        stage = copy.deepcopy(self.stage)
        stage['input_processing']=self.processing.value
        for field, rows in self.rows.items():
            refs = [ref for check, ref in rows if check.value]
            if not refs:
                if require_selection:
                    raise ValueError('Selecione ao menos um arquivo para ' + LABELS.get(field, field) + '.')
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
