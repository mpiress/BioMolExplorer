"""Explicit redocking selection and preparation beside input sources."""
import copy
import flet as ft
from biomolexplorer.redocking_config import pair_key, validate_pairs


class RedockingPairs:
    def __init__(self, form):
        self.form = form
        self.rows = []
        from .preparation_settings import PreparationSettings
        saved=form.stage['parameters'].get('preparation_pairs') or {}
        seed=next(iter(saved.values()),{})
        self.preparation=PreparationSettings(form,seed)
        general_keys=('receptor','ligand','preparation_backend','ambertools_path','ambertools_env','vina_charge_policy')
        differing=any(any(config.get(key)!=seed.get(key) for key in general_keys) for config in saved.values())
        self.migration_notice=ft.Text('Existem opções de preparo antigas diferentes entre pares. As opções gerais exibidas partem do primeiro par e serão aplicadas a todos ao confirmar. Revise-as antes de continuar.',
            size=12,color='#92400E',visible=differing)

        # Cofactors and structural identities stay attached to each curated pair.
        self.preparation.has_cofactors.visible=False
        self.preparation.cofactors.visible=False
        self.preparation.groups['receptor'].title=ft.Text('Preparação de todos os receptores')
        self.preparation.groups['ligand'].title=ft.Text('Preparação de todos os ligantes')

        self.preparation_panels={}
        self.box = ft.Column(spacing=20)
        self.summary = ft.Text()
        self.available = ft.Dropdown(label='PDB / ligante / resíduo / cadeia disponíveis',
            disabled=not form.writable, enable_filter=True, enable_search=True)
        self.select_all_button = ft.TextButton('Selecionar todos os pares disponíveis',
            on_click=self.choose_all, disabled=not form.writable)
        self.control = ft.Column([
            ft.Text('Pares para redocking', size=16, weight=ft.FontWeight.W_600),
            ft.Text('Os pares das estruturas selecionadas são preservados da curadoria. Revise as identidades abaixo; a preparação geral será aplicada a todo o conjunto.', size=12),
            ft.Text('Para PDBs recuperados com triagem, esta lista contém os candidatos priorizados. Use Ligantes nos resultados da recuperação para revisar ou incluir alternativas.',size=12),
            self.available,
            ft.TextButton('Selecionar par e configurar preparação', on_click=self.choose, disabled=not form.writable),
            self.select_all_button,
            ft.TextButton('Informar par manualmente', on_click=lambda e: self.add(), disabled=not form.writable),
            self.summary, self.box, self.migration_notice, self.preparation.control], spacing=18, data='input')
        for record in form.stage['parameters'].get('pdb_codes') or []:
            self.add(record, update=False)
        self.refresh()

    def candidates(self):
        editor = self.form.input_editors['base_input_path']
        if editor.selection is not None:
            return editor.selection.read(require_selection=False)['parameters'].get('pdb_codes') or []
        result = []
        for source, selector, _ in editor.rows:
            for record in getattr(self.form.ui, 'redocking_records', {}).get(source.value, []):
                if selector.value not in (None, 'auto') and selector.value.rsplit('/', 1)[-1].split('.')[0].lower() not in (record[0].lower(), f'{record[0]}_{record[3]}'.lower()):
                    continue
                if record not in result:result.append(record)
        return result

    def refresh(self):
        self.records = self.candidates()
        editor=self.form.input_editors['base_input_path']
        explicit=editor.selection is not None or any(selector.value not in (None,'auto') for _,selector,_ in editor.rows)
        signature=tuple(pair_key(r) for r in self.records)
        if explicit and signature!=getattr(self,'_selection_signature',None):
            previous=getattr(self,'_selection_signature',())
            for entry in list(self.rows):
                if previous and pair_key([f.value for f in entry[0]]) in previous and pair_key([f.value for f in entry[0]]) not in signature:
                    self.rows.remove(entry);self.box.controls.remove(entry[-1])
                    self.preparation_panels.pop(id(entry[-1]),None)
            selected={pair_key([f.value for f in entry[0]]) for entry in self.rows}
            for record in self.records:
                if pair_key(record) not in selected:self.add(record,update=False)
        self._selection_signature=signature
        for control in self.control.controls[3:7]:control.visible=not (explicit and bool(self.records))

        self.available.options = [ft.DropdownOption(key=str(i), text=f'{r[0]} / {r[1]} / {r[2]} / {r[3]}') for i, r in enumerate(self.records)]
        self.available.value = None
        self.select_all_button.disabled=not self.form.writable or not self.records
        self.sync_preparation()
        self.summarize()

    def summarize(self):
        self.summary.value = 'Pares selecionados: ' + ('; '.join(
            f'{fields[0].value} / {fields[1].value} / {fields[2].value} / {fields[3].value}'
            for fields, _, _, _, _ in self.rows) if self.rows else 'nenhum')

    def choose(self, e):
        if self.available.value is not None:
            self.add(self.records[int(self.available.value)])

    def choose_all(self, e):
        if not self.form.writable:return
        selected=set();receptors={}
        for fields, _, (has_cofactors,cofactors), options, _ in self.rows:
            record=[field.value.strip() for field in fields]
            record[1]=record[1].upper()
            try:record[2]=int(record[2])
            except ValueError:pass  # Keep unfinished manual rows for the user.
            selected.add(pair_key(record))
            receptors.setdefault((record[0],record[3]),{
                'cofactors':[c.strip().upper() for c in cofactors.value.split(',') if c.strip()] if has_cofactors.value else [],
                'receptor':{key:control.value for key,control in options['receptor'].items()}})
        for record in self.candidates():
            key=pair_key(record)
            if key in selected:continue
            self.add(record,update=False,config=receptors.get((record[0],record[3])))
            selected.add(key)
        self.sync_preparation();self.summarize();self.form.ui.page.update()

    def add(self, record=None, update=True, config=None):
        record = list(record or ['', '', '', 'A'])
        saved = self.form.stage['parameters'].get('preparation_pairs') or {}
        config = copy.deepcopy(saved.get(pair_key(record), {}) if config is None else config)
        fields = [ft.TextField(label=label, value=str(value), disabled=not self.form.writable,
            on_change=lambda e: self.changed()) for label, value in zip(
                ('Código PDB', 'Ligante', 'Número do resíduo', 'Cadeia'), record[:4])]
        has_cofactors=ft.Switch(label='Manter cofatores neste receptor',value=bool(config.get('cofactors')),disabled=not self.form.writable)
        cofactors=ft.TextField(label='Códigos dos cofatores',value=', '.join(config.get('cofactors',[])))
        settings=self.preparation.settings
        has_cofactors.on_change=lambda e:(self.sync_preparation(),self.form.ui.page.update())
        from .settings_layout import grid
        panel=ft.Container(content=ft.Column([grid(fields),has_cofactors,cofactors],spacing=16),
            padding=20,border=ft.Border.all(1,'#E2E8F0'),border_radius=12)
        entry = (fields, record, (has_cofactors, cofactors), settings, panel)
        self.preparation_panels[id(panel)]=config
        def remove(e):
            self.preparation_panels.pop(id(panel),None)
            self.rows.remove(entry); self.box.controls.remove(panel); self.changed()
        panel.content.controls.append(ft.TextButton('Remover par', on_click=remove, disabled=not self.form.writable))
        self.rows.append(entry); self.box.controls.append(panel)
        if update:
            self.sync_preparation();self.summarize();self.form.ui.page.update()

    def changed(self):
        self.summarize(); self.form.ui.page.update()

    def sync_preparation(self):
        inactive=not self.form.writable or not self.form.field_controls['prepare_complex'].value
        self.preparation.sync(inactive)
        for _, _, (has_cofactors,cofactors), settings, panel in self.rows:
            has_cofactors.disabled=inactive
            cofactors.disabled=inactive or not has_cofactors.value
            for options in settings.values():
                for control in options.values():control.disabled=inactive

    def read(self):
        records = []; settings = {}
        candidates={pair_key(r):r for r in self.candidates()}
        general=self.preparation.read()
        for fields, original, (has_cofactors, cofactors), options, panel in self.rows:
            try: record = [fields[0].value.strip(), fields[1].value.strip().upper(), int(fields[2].value), fields[3].value.strip()]
            except (ValueError, TypeError): raise ValueError('Informe o número do resíduo do ligante para cada par selecionado.')
            # Only the selected metadata tuple supplies resolution.
            match = candidates.get(pair_key(record))
            if match is None and original[:4] == record[:4]: match = original
            if candidates and match is None:
                raise ValueError('O par selecionado não corresponde ao PDB, ligante, resíduo e cadeia dos dados de entrada.')
            if match and len(match) > 4: record.append(match[4])
            config = copy.deepcopy(self.preparation_panels[id(panel)])
            for key in ('receptor','ligand','preparation_backend','ambertools_path','ambertools_env','vina_charge_policy'):
                config.pop(key,None)
                if key in general:config[key]=copy.deepcopy(general[key])
            config['cofactors']=[c.strip().upper() for c in cofactors.value.split(',') if c.strip()] if has_cofactors.value else []
            if has_cofactors.value and not config['cofactors']: raise ValueError('Informe quais cofatores devem ser mantidos no receptor.')
            records.append(record); settings[pair_key(record)] = config
        # An empty draft is valid while the user is assembling the pipeline.
        # PipelineService.resume()/execution performs the mandatory check.
        if records:validate_pairs(records, settings)
        return records, settings
