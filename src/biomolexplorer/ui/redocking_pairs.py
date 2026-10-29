"""Explicit redocking selection and preparation beside input sources."""
import copy
from types import SimpleNamespace
import flet as ft
from biomolexplorer.redocking_config import pair_key, validate_pairs
from .localization import verbatim
from .feedback import close_dialog


class RedockingPairs:
    def __init__(self, form):
        self.form = form
        self.rows = []
        self.table=ft.DataTable(columns=[ft.DataColumn(ft.Text(t)) for t in ('PDB','Ligante','Número do resíduo','Cadeia','Ações')],rows=[],column_spacing=36,heading_row_color='#F1F5F9',data_row_min_height=56,data_row_max_height=64,divider_thickness=.5)
        self.box = ft.Row([self.table],scroll=ft.ScrollMode.AUTO,on_size_change=self.resize)
        self.preparation=form.preparation_settings
        self.summary = ft.Text()
        self.available = ft.Dropdown(label='PDB / ligante / resíduo / cadeia disponíveis',
            disabled=not form.writable, enable_filter=True, enable_search=True)
        self.select_all_button = ft.TextButton('Selecionar todos os pares disponíveis',
            on_click=self.choose_all, disabled=not form.writable)
        self.control = ft.Column([
            ft.Text('Pares para redocking', size=16, weight=ft.FontWeight.W_600),
            ft.Text('Selecione e valide cada par. A cadeia identifica o par; a preparação comum está nas seções de receptor e ligante. A resolução é preservada dos metadados.', size=12),
            ft.Text('Para PDBs recuperados com triagem, esta lista contém os candidatos priorizados. Use Ligantes nos resultados da recuperação para revisar ou incluir alternativas.',size=12),
            self.available,
            ft.TextButton('Adicionar par selecionado', on_click=self.choose, disabled=not form.writable),
            self.select_all_button,
            ft.TextButton('Informar par manualmente', on_click=lambda e:self.open_editor(), disabled=not form.writable),
            self.summary, self.box], spacing=18, data='pairs')
        for record in form.stage['parameters'].get('pdb_codes') or []:
            self.add(record, update=False)
        self.refresh()

    def resize(self,e):
        # Fill the available card width; preserve readable columns and horizontal
        # scrolling on small viewports. Avoid updates when layout is unchanged.
        width=max(680,round(e.width))
        if self.table.width!=width:
            self.table.width=width
            self.form.ui.page.update()

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
        self.summary.value = f'Pares selecionados: {len(self.rows)}'

    def choose(self, e):
        if self.available.value is not None:
            self.add(self.records[int(self.available.value)])

    def choose_all(self, e):
        if not self.form.writable:return
        selected={pair_key([field.value.strip() for field in fields]) for fields,_,_,_,_ in self.rows}
        for record in self.candidates():
            key=pair_key(record)
            if key in selected:continue
            self.add(record,update=False)
            selected.add(key)
        self.sync_preparation();self.summarize();self.form.ui.page.update()

    def add(self, record=None, update=True, config=None):
        record = list(record or ['', '', '', 'A'])
        fields=[SimpleNamespace(value=str(value)) for value in record[:4]]
        has_cofactors,cofactors=self.preparation.has_cofactors,self.preparation.cofactors
        settings=self.preparation.settings
        panel=ft.DataRow(cells=[ft.DataCell(verbatim(ft.Text(field.value))) for field in fields])
        entry=(fields,record,(has_cofactors,cofactors),settings,panel)
        def remove(e):
            self.rows.remove(entry);self.table.rows.remove(panel);self.changed()
        panel.cells.append(ft.DataCell(ft.Row([
            ft.IconButton(ft.Icon(ft.Icons.EDIT_OUTLINED,semantics_label='Editar par'),on_click=lambda e:self.open_editor(entry),disabled=not self.form.writable,tooltip='Editar par'),
            ft.IconButton(ft.Icon(ft.Icons.DELETE_OUTLINE,semantics_label='Remover par'),on_click=remove,disabled=not self.form.writable,tooltip='Remover par')],spacing=8)))
        self.rows.append(entry);self.table.rows.append(panel)
        if update:
            self.sync_preparation();self.summarize();self.form.ui.page.update()

    def changed(self):
        for fields,_,_,_,row in self.rows:
            for field,cell in zip(fields,row.cells[:4]):cell.content.value=field.value
        self.summarize(); self.form.ui.page.update()

    def open_editor(self,entry=None):
        if not self.form.writable:return
        original=[field.value for field in entry[0]] if entry else ['', '', '', 'A']
        controls=[ft.TextField(label=label,value=value) for label,value in zip(
            ('Código PDB','Ligante','Número do resíduo','Cadeia'),original)]
        error=ft.Text(color='#B91C1C',visible=False)
        def close(e=None):close_dialog(self.form.ui.page,dialog)
        def save(e):
            try:
                try:record=[controls[0].value.strip(),controls[1].value.strip().upper(),int(controls[2].value),controls[3].value.strip()]
                except (TypeError,ValueError):raise ValueError('Informe o número do resíduo do ligante para cada par selecionado.') from None
                others=[[f.value.strip() for f in fields] for fields,_,_,_,_ in self.rows if entry is None or fields is not entry[0]]
                validate_pairs(others+[record],{})
                candidates=self.candidates()
                match=next((r for r in candidates if list(r[:4])==record),None)
                if candidates and match is None:raise ValueError('O par selecionado não corresponde ao PDB, ligante, resíduo e cadeia dos dados de entrada.')
                if entry:
                    for field,value in zip(entry[0],record):field.value=str(value)
                    # Resolution belongs to the validated tuple, never to an edited number.
                    entry[1][:]=list(match or (entry[1] if pair_key(entry[1])==pair_key(record) else record))
                else:self.add(list(match or record),update=False)
                self.changed();close()
            except ValueError as exc:
                error.value=str(exc);error.visible=True;self.form.ui.page.update()
        dialog=ft.AlertDialog(modal=True,title=ft.Text('Editar par' if entry else 'Adicionar par'),
            content=ft.Container(width=min(620,max(260,(getattr(self.form.ui.page,'width',None) or 900)-120)),
                content=ft.Column([ft.ResponsiveRow([ft.Container(c,col={'xs':12,'sm':6},padding=ft.Padding.symmetric(vertical=8)) for c in controls],spacing=24,run_spacing=20),error],tight=True,spacing=16,scroll=ft.ScrollMode.AUTO)),
            actions=[ft.TextButton('Cancelar',on_click=close),ft.TextButton('Salvar par',on_click=save)])
        self.form.ui.page.show_dialog(dialog)

    def sync_preparation(self):
        self.preparation.sync(inactive=not self.form.field_controls['prepare_complex'].value)

    def read(self):
        records = []; settings = {}
        shared=self.preparation.read()
        for fields, original, (has_cofactors, cofactors), options, _ in self.rows:
            try: record = [fields[0].value.strip(), fields[1].value.strip().upper(), int(fields[2].value), fields[3].value.strip()]
            except (ValueError, TypeError): raise ValueError('Informe o número do resíduo do ligante para cada par selecionado.')
            # Only the selected metadata tuple supplies resolution.
            match = candidates.get(pair_key(record))
            if match is None and original[:4] == record[:4]: match = original
            if candidates and match is None:
                raise ValueError('O par selecionado não corresponde ao PDB, ligante, resíduo e cadeia dos dados de entrada.')
            if match and len(match) > 4: record.append(match[4])
            config=copy.deepcopy(shared)
            records.append(record); settings[pair_key(record)] = config
        # An empty draft is valid while the user is assembling the pipeline.
        # PipelineService.resume()/execution performs the mandatory check.
        if records:validate_pairs(records, settings)
        return records, settings
