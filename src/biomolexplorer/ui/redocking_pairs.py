"""Explicit redocking selection and preparation beside input sources."""
import copy
import re
import flet as ft
from biomolexplorer.redocking_config import DEFAULTS, pair_key, validate_pairs


class RedockingPairs:
    def __init__(self, form):
        self.form = form
        self.rows = []
        self.box = ft.Column(spacing=20)
        self.summary = ft.Text()
        self.available = ft.Dropdown(label='PDB / ligante / resíduo / cadeia disponíveis',
            disabled=not form.writable, enable_filter=True, enable_search=True)
        self.control = ft.Column([
            ft.Text('Pares para redocking', size=16, weight=ft.FontWeight.W_600),
            ft.Text('Selecione e valide cada par. Configure a cadeia, os cofatores e a preparação abaixo. A resolução é preservada dos metadados.', size=12),
            self.available,
            ft.TextButton('Selecionar par e configurar preparação', on_click=self.choose, disabled=not form.writable),
            ft.TextButton('Informar par manualmente', on_click=lambda e: self.add(), disabled=not form.writable),
            self.summary, self.box], spacing=18, data='input')
        for record in form.stage['parameters'].get('pdb_codes') or []:
            self.add(record, update=False)
        self.refresh()

    def candidates(self):
        editor = self.form.input_editors['base_input_path']
        result = []
        for source, selector, _ in editor.rows:
            for record in getattr(self.form.ui, 'redocking_records', {}).get(source.value, []):
                # A selected structure restricts the tuples; auto includes all metadata.
                if selector.value not in (None, 'auto') and selector.value.rsplit('/', 1)[-1].split('.')[0].lower() not in (record[0].lower(), f'{record[0]}_{record[3]}'.lower()):
                    continue
                if record not in result: result.append(record)
        return result

    def refresh(self):
        self.records = self.candidates()
        self.available.options = [ft.DropdownOption(key=str(i), text=f'{r[0]} / {r[1]} / {r[2]} / {r[3]}') for i, r in enumerate(self.records)]
        self.available.value = None
        self.summarize()

    def summarize(self):
        self.summary.value = 'Pares selecionados: ' + ('; '.join(
            f'{fields[0].value} / {fields[1].value} / {fields[2].value} / {fields[3].value}'
            for fields, _, _, _, _ in self.rows) if self.rows else 'nenhum')

    def choose(self, e):
        if self.available.value is not None:
            self.add(self.records[int(self.available.value)])

    def add(self, record=None, update=True):
        record = list(record or ['', '', '', 'A'])
        saved = self.form.stage['parameters'].get('preparation_pairs') or {}
        config = copy.deepcopy(saved.get(pair_key(record), {}))
        fields = [ft.TextField(label=label, value=str(value), disabled=not self.form.writable, width=170,
            on_change=lambda e: self.changed()) for label, value in zip(
                ('Código PDB', 'Ligante', 'Número do resíduo', 'Cadeia'), record[:4])]
        cofactors = ft.TextField(label='Cofatores do receptor (códigos separados por vírgula)', value=', '.join(config.get('cofactors', [])), disabled=not self.form.writable)
        has_cofactors = ft.Switch(label='Considerar cofatores como parte do receptor', value=bool(config.get('cofactors')), disabled=not self.form.writable)
        cofactors.visible = has_cofactors.value
        def toggle(e):
            cofactors.visible = has_cofactors.value
            self.form.ui.page.update()
        has_cofactors.on_change = toggle
        settings = {}
        preparation = []
        for role, title in (('receptor', 'Preparação do receptor'), ('ligand', 'Preparação e conformação do ligante')):
            values = dict(DEFAULTS, **config.get(role, {}))
            # Preserve the previous global charge method on initial migration.
            if role not in config:
                values['charge_type'] = self.form.stage['parameters'].get('charge_type', 'gas') if role=='ligand' else 'gas'
                templates=self.form.stage.get('templates',{})
                legacy=next((templates[n] for n in (('chimera/prepare_ligand.template','chimera/prepare_better_conform.template') if role=='ligand' else ('chimera/prepare_receptor.template',)) if n in templates),None)
                if legacy is not None:
                    commands=legacy.splitlines()
                    values['add_hydrogens']='addh' in commands
                    values['minimize']=any(line.startswith('minimize') for line in commands)
                    match=re.search(r'method (gas|am1)',legacy)
                    if match:values['charge_type']=match.group(1)
                complex_source=templates.get('chimera/prepare_complex.template')
                if complex_source is not None:
                    values['remove_solvent']='delete solvent' in complex_source.splitlines()
                    values['remove_hydrogens']='delete element.H' in complex_source.splitlines()
            options = {key: ft.Switch(label=label, value=values[key], disabled=not self.form.writable) for key, label in (
                ('remove_solvent', 'Remover solvente'), ('remove_hydrogens', 'Remover hidrogênios existentes'),
                ('add_hydrogens', 'Adicionar hidrogênios'), ('minimize', 'Minimizar energia'))}
            options['charge_type'] = ft.Dropdown(label='Método de cargas', value=values['charge_type'],
                options=[ft.DropdownOption(key=v, text=v) for v in ('gas', 'am1')], disabled=not self.form.writable)
            settings[role] = options
            preparation.append(ft.Column([ft.Text(title, weight=ft.FontWeight.W_600), *options.values()], spacing=12))
        panel = ft.Container(content=ft.Column([ft.Row(fields, wrap=True), has_cofactors, cofactors,
            *preparation], spacing=16), padding=20, border=ft.Border.all(1, '#E2E8F0'), border_radius=12)
        entry = (fields, record, (has_cofactors, cofactors), settings, panel)
        def remove(e):
            self.rows.remove(entry); self.box.controls.remove(panel); self.changed()
        panel.content.controls.append(ft.TextButton('Remover par', on_click=remove, disabled=not self.form.writable))
        self.rows.append(entry); self.box.controls.append(panel)
        self.summarize()
        if update: self.form.ui.page.update()

    def changed(self):
        self.summarize(); self.form.ui.page.update()

    def read(self):
        records = []; settings = {}
        for fields, original, (has_cofactors, cofactors), options, _ in self.rows:
            try: record = [fields[0].value.strip(), fields[1].value.strip().upper(), int(fields[2].value), fields[3].value.strip()]
            except (ValueError, TypeError): raise ValueError('Informe o número do resíduo do ligante para cada par selecionado.')
            # Only the selected metadata tuple supplies resolution.
            match = next((r for r in self.candidates() if r[:3] == record[:3] and r[3] == record[3]), None)
            if match is None and original[:4] == record[:4]: match = original
            if self.candidates() and match is None:
                raise ValueError('O par selecionado não corresponde ao PDB, ligante, resíduo e cadeia dos dados de entrada.')
            if match and len(match) > 4: record.append(match[4])
            config = {'cofactors': [c.strip().upper() for c in cofactors.value.split(',') if c.strip()] if has_cofactors.value else [],
                **{role: {key: c.value for key, c in controls.items()} for role, controls in options.items()}}
            if has_cofactors.value and not config['cofactors']: raise ValueError('Informe quais cofatores devem ser mantidos no receptor.')
            records.append(record); settings[pair_key(record)] = config
        validate_pairs(records, settings)
        return records, settings
