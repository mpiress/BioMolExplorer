"""Targets and compounds of the pose validation block.

Each section follows the same layout: the origin of the selection, a combo
with every valid item of that origin, a search field with autocomplete beside
it, buttons to include the chosen item or the whole list, and the table of the
included items.  Targets are the prepared receptors of the connected redocking
or preparation block, listed by PDB code; compounds are the molecules that
produced the connected docking or consensus results.

The catalog is read from the artifact registry and the published score tables,
so opening the block never loads execution JSON nor walks the pose tree.  Every
change patches only the controls of its own section; the dialog never repaints
the page.
"""
import asyncio
import copy
from pathlib import Path
import flet as ft

from biomolexplorer.pose_validation_config import HELP, LABELS
from .localization import verbatim
from .feedback import close_dialog
from .selected_items_table import SelectedItemsTable

DEBOUNCE = .18
SUGGESTIONS = 50
SUGGESTION_ROWS = 6
SUGGESTION_HEIGHT = 44
COMBO_LIMIT = 1000
PAGE_SIZE = 25
ALL_ORIGINS = '__all__'
# Vina first: its scene also classifies the interaction of each residue.
ENGINE_ORDER = ('vina','dock6')
MUTED = '#64748B'
WARNING = '#B45309'


def action_chip(label, tooltip, on_click, unavailable=None):
    """Compact labelled action, in the style of the interface badges."""
    enabled = on_click is not None
    tooltip = tooltip if enabled else (unavailable or tooltip)
    return ft.Container(width=38,height=30,border_radius=8,ink=enabled,tooltip=tooltip,
        bgcolor='#EFF6FF' if enabled else '#F1F5F9',alignment=ft.Alignment.CENTER,
        border=ft.Border.all(1,'#C7D7FE' if enabled else '#E2E8F0'),on_click=on_click,
        content=verbatim(ft.Text(label,size=12,weight=ft.FontWeight.W_600,
                                 color='#2552E8' if enabled else MUTED)))


def card_decimals(properties, digits=3):
    """Card values with at most three decimals; other values are untouched.

    Providers and RDKit return full double precision, which the generic
    property panels print verbatim: this keeps every number of the card on the
    same scale without changing how other screens present their data.
    """
    result = {}
    for key, value in (properties or {}).items():
        if isinstance(value,float):
            result[key] = round(value,digits)
        elif isinstance(value,str) and not isinstance(value,bool):
            try:
                number = float(value.strip())
            except (TypeError,ValueError):
                result[key] = value
            else:
                result[key] = value if value.strip().lstrip('-+').isdigit() else round(number,digits)
        else:
            result[key] = value
    return result


def rounded(value, digits=2):
    """Scientific metadata may be absent or non-numeric; never break the table."""
    try:
        return format(round(float(value),digits),'g')
    except (TypeError,ValueError):
        return str(value or '')


class PosePanel:
    """Origin, combo, search, include buttons and selection table of one section."""

    TEXT = {
        'targets':dict(origin='Origem dos targets',combo='PDBs preparados disponíveis',
                       search='Pesquisar código PDB',hint='Ex.: 4M0E',
                       noun='target',plural='targets',none='Nenhum PDB preparado disponível na origem.',
                       empty='Nenhum target incluído: todos os PDBs preparados da origem serão validados.',
                       connect='Conecte um bloco de redocking ou preparação em Targets para listar os PDBs preparados.'),
        'compounds':dict(origin='Origem dos compostos',combo='Compostos disponíveis',
                         search='Pesquisar composto',hint='Ex.: CHEMBL92663',
                         noun='composto',plural='compostos',none='Nenhum composto com resultado na origem.',
                         empty='Nenhum composto incluído: todos os compostos com resultados serão validados.',
                         connect='Conecte um bloco Vina, DOCK6 ou Consenso em Compounds para listar os compostos.'),
    }

    def __init__(self, settings, kind, field, columns, weights, actions_width):
        self.settings, self.kind, self.field = settings, kind, field
        self.text = self.TEXT[kind]
        self.entries = []
        self.offset = 0
        self.focus = None
        self._task = None
        writable = settings.form.writable
        self.editor = settings.form.input_editors.get(field)
        self.initial_origins = self.origin_values()
        self.origin = ft.Dropdown(label=self.text['origin'],options=[],on_select=self.origin_changed,
            disabled=not writable,expand=True,menu_height=320)
        self.combo = ft.Dropdown(label=self.text['combo'],options=[],enable_filter=True,enable_search=True,
            menu_height=320,expand=True,disabled=not writable,hint_text='Selecione um item da lista')
        self.include_one = ft.FilledTonalButton('Incluir selecionado',icon=ft.Icons.ADD,
            on_click=self.include_selected,disabled=not writable)
        self.include_list = ft.TextButton('Incluir todos',icon=ft.Icons.PLAYLIST_ADD,
            on_click=self.include_all,disabled=not writable)
        self.search = ft.TextField(label=self.text['search'],hint_text=self.text['hint'],
            prefix_icon=ft.Icons.SEARCH,on_change=self.changed,on_submit=self.include_search,
            disabled=not writable,expand=True)
        self.search_add = ft.IconButton(ft.Icons.ADD_CIRCLE_OUTLINE,tooltip='Incluir item pesquisado',
            on_click=self.include_search,disabled=not writable)
        # A zero-height ListView inside an unbounded dialog can expand to an
        # enormous viewport in Flet Web. Suggestions are capped below, so a
        # content-sized Column is sufficient and never owns a scroll view.
        self.suggestions = ft.Column(spacing=0,tight=True)
        self.suggestion_count = ft.Text(size=11,color=MUTED)
        self.suggestion_box = ft.Column([self.suggestions,self.suggestion_count],
                                        spacing=4,tight=True,visible=False)
        self.available = ft.Text(size=12,color=MUTED)
        self.feedback = ft.Text(size=12,color=WARNING,visible=False)
        self.pair_table = SelectedItemsTable(columns,weights,actions_width)
        self.table = self.pair_table.table
        self.table_box = self.pair_table.control
        self.empty = ft.Text(self.text['empty'],size=12,color=MUTED)
        self.selected_count = ft.Text(size=12,color=MUTED)
        self.clear_button = ft.TextButton('Remover todos',icon=ft.Icons.CLEAR_ALL,on_click=self.clear,
            disabled=not writable)
        self.previous = ft.TextButton('Anterior',on_click=lambda e:self.move(-1))
        self.next = ft.TextButton('Próxima',on_click=lambda e:self.move(1))
        self.pager = ft.Row([self.previous,self.next],wrap=True,visible=False)
        self.control = ft.Column([
            self.origin,
            self.combo,
            ft.Row([self.include_one,self.include_list],wrap=True,spacing=8),
            ft.Row([self.search,self.search_add],spacing=4,
                   vertical_alignment=ft.CrossAxisAlignment.CENTER),
            self.suggestion_box,
            self.available,self.feedback,
            ft.Row([ft.Text('Itens incluídos',size=14,weight=ft.FontWeight.W_600),
                    self.selected_count,self.clear_button],wrap=True,spacing=12,
                   vertical_alignment=ft.CrossAxisAlignment.CENTER),
            self.empty,self.table_box,self.pager],
            spacing=12,data=kind,horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
        self.draw_origin()

    # ---------------------------------------------------------------- origin
    def origin_values(self):
        if self.editor is None:
            return []
        return [source.value for source,_,_ in self.editor.rows
                if source.value and source.value != 'configured-path']

    def origin_options(self):
        if self.editor is None:
            return []
        options = [option for option in self.editor.source_options() if option.key]
        if len(self.initial_origins) > 1:
            names = {option.key:option.text for option in options}
            label = ' + '.join(names.get(value,value) for value in self.initial_origins)
            options.insert(0,verbatim(ft.DropdownOption(key=ALL_ORIGINS,text='Todas as origens conectadas: '+label),'text'))
        return options

    def draw_origin(self):
        self.origin.options = self.origin_options()
        current = self.origin_values()
        if len(current) > 1:
            self.origin.value = ALL_ORIGINS if ALL_ORIGINS in {o.key for o in self.origin.options} else current[0]
        else:
            self.origin.value = current[0] if current else None
        self.origin.helper_text = None if current else self.text['connect']

    def set_origins(self, values):
        """Point the block input at these origins, as the canvas link would."""
        editor = self.editor
        for source,_,row in list(editor.rows):
            editor.entity_fields.pop(id(source),None); editor.compound_fields.pop(id(source),None)
            editor.reference_metadata.pop(id(source),None)
            if row in editor.box.controls:
                editor.box.controls.remove(row)
        editor.rows.clear()
        for value in values:
            kind,identifier = value.split(':',1)
            editor.add({kind:identifier})
        if not values:
            editor.add({})
        editor.selection_signature = None

    async def origin_changed(self, event=None):
        if not self.settings.form.writable:
            return
        value = self.origin.value
        values = list(self.initial_origins) if value == ALL_ORIGINS else ([value] if value else [])
        if values == self.origin_values():
            return
        self.set_origins(values)
        await self.settings.reload(prune=True)

    # ------------------------------------------------------------- selection
    @property
    def selected(self):
        return self.settings.selected_targets if self.kind == 'targets' else self.settings.selected_compounds

    def identifier(self, entry):
        return entry['receptor_id'] if self.kind == 'targets' else entry['compound_id']

    def known(self):
        return {self.identifier(entry):entry for entry in self.entries}

    def include(self, identifiers):
        if not self.settings.form.writable:
            return 0
        known = self.known(); added = 0; first = None
        for identifier in identifiers:
            if identifier in known and identifier not in self.selected:
                self.selected.append(identifier); added += 1
                first = first or identifier
        if added:
            for identifier in identifiers:
                self.settings.drop_exclusions(**({'receptor':identifier} if self.kind == 'targets'
                                                 else {'compound':identifier}))
            self.combo.value = None; self.search.value = ''
            # Show the page where the first included item landed.
            self.focus = first
            self.message('')
            self.settings.selection_changed(self.kind)
        elif identifiers:
            missing = [value for value in identifiers if value not in known]
            if not missing:
                self.message('O item escolhido já está incluído.')
            elif self.kind == 'compounds' and self.settings.selected_targets:
                # Consensus keeps one target per compound: a compound may have
                # results only for targets that are not included.
                self.message(f'{missing[0]} não tem resultado nos targets incluídos nesta origem.')
            else:
                self.message(f'{missing[0]} não está disponível na origem selecionada.')
        return added

    def include_selected(self, event=None):
        if not self.combo.value:
            self.message(f'Escolha um {self.text["noun"]} na lista antes de incluir.'); return
        self.include([self.combo.value])

    def include_all(self, event=None):
        candidates = [self.identifier(entry) for entry in self.visible()]
        if not candidates:
            self.message('Não há itens disponíveis para incluir.'); return
        self.include(candidates)

    def include_search(self, event=None):
        term = (self.search.value or '').strip().casefold()
        if not term:
            self.message('Digite um código para pesquisar.'); return
        exact = [entry for entry in self.visible() if term in self.keys(entry)]
        found = exact or self.matches()
        if len(found) == 1 or exact:
            self.include([self.identifier(entry) for entry in found]); return
        self.message(f'Nenhum {self.text["noun"]} disponível corresponde a “{self.search.value.strip()}”.'
                     if not found else
                     f'{len(found)} correspondências: escolha uma sugestão ou refine a pesquisa.')

    def remove(self, identifier, receptor=None):
        """Remove this row only: one target, or one ligand-receptor complex."""
        if not self.settings.form.writable:
            return
        self.focus = None
        if self.kind == 'compounds' and receptor:
            self.settings.exclude_pair(identifier,receptor)
        else:
            self.selected[:] = [value for value in self.selected if value != identifier]
            if self.kind == 'targets':
                self.settings.drop_exclusions(receptor=identifier)
            else:
                self.settings.drop_exclusions(compound=identifier)
        self.settings.selection_changed(self.kind)

    def clear(self, event=None):
        if not self.settings.form.writable or not self.selected:
            return
        for identifier in list(self.selected):
            self.settings.drop_exclusions(**({'receptor':identifier} if self.kind == 'targets'
                                             else {'compound':identifier}))
        self.selected.clear(); self.offset = 0; self.focus = None
        self.settings.selection_changed(self.kind)

    def move(self, step):
        offset = self.offset+step*PAGE_SIZE
        if 0 <= offset < len(self.rows()):
            self.offset = offset; self.draw_table(); self.settings.patch(self.control)

    def message(self, text):
        self.feedback.value = text; self.feedback.visible = bool(text)
        self.settings.patch(self.feedback)

    # ---------------------------------------------------------------- search
    def keys(self, entry):
        if self.kind == 'targets':
            return {entry['receptor_id'].casefold(),entry['pdb_code'].casefold()}
        return {entry['compound_id'].casefold()}

    def haystack(self, entry):
        if self.kind == 'targets':
            return ' '.join([entry['receptor_id'],entry.get('label',''),str(entry.get('ligand',''))]).casefold()
        return ' '.join([entry['compound_id'],*entry.get('receptors',[])]).casefold()

    def visible(self):
        """Items not yet included that match the current search."""
        term = (self.search.value or '').strip().casefold()
        chosen = set(self.selected)
        found = [entry for entry in self.entries if self.identifier(entry) not in chosen
                 and (not term or term in self.haystack(entry))]
        if term:
            # Codes that start with the typed text come first.
            found.sort(key=lambda entry:not any(key.startswith(term) for key in self.keys(entry)))
        return found

    def matches(self):
        return self.visible() if (self.search.value or '').strip() else []

    def changed(self, event=None):
        if self._task is not None and not self._task.done():
            self._task.cancel()
        page = getattr(self.settings.form.ui,'page',None)
        if page is None or not hasattr(page,'run_task'):
            self.filter(); return
        self._task = page.run_task(self._deferred)

    async def _deferred(self):
        try:
            await asyncio.sleep(DEBOUNCE)
        except asyncio.CancelledError:
            return
        self.filter()

    def filter(self, patch=True):
        """Autocomplete: suggestions and the combo follow every typed character."""
        if self.feedback.visible:
            self.message('')
        self.draw_choices()
        if patch:
            self.settings.patch(self.combo,self.include_list,self.suggestion_box,self.available)

    # ------------------------------------------------------------------ draw
    def option_text(self, entry):
        if self.kind == 'targets':
            return ' · '.join(filter(None,[entry['pdb_code'],('cadeia '+entry['chain']) if entry.get('chain') else '',
                ('ligante '+entry['ligand']) if entry.get('ligand') else '',
                ('RMSD '+rounded(entry['rmsd'],3)+' Å') if entry.get('rmsd') else '']))
        targets = entry.get('receptors',[])
        return ' · '.join(filter(None,[entry['compound_id'],
            ('1 target' if len(targets) == 1 else f'{len(targets)} targets') if targets else '',
            '/'.join(entry.get('engines',[]))]))

    def suggestion(self, entry):
        identifier = self.identifier(entry)
        title = entry['pdb_code'] if self.kind == 'targets' else identifier
        detail = self.option_text(entry).split(' · ',1)[1] if ' · ' in self.option_text(entry) else ''
        return ft.Container(height=SUGGESTION_HEIGHT,border_radius=8,ink=True,
            padding=ft.Padding.symmetric(horizontal=10),on_click=lambda e,value=identifier:self.include([value]),
            tooltip='Incluir '+identifier,
            content=ft.Row([ft.Icon(ft.Icons.ADD,size=16,color='#2552E8'),
                verbatim(ft.Text(title,size=13,weight=ft.FontWeight.W_500)),
                verbatim(ft.Text(detail,size=11,color=MUTED,expand=True,max_lines=1,overflow=ft.TextOverflow.ELLIPSIS))],
                spacing=8,vertical_alignment=ft.CrossAxisAlignment.CENTER))

    def draw_choices(self):
        visible = self.visible()
        self.combo.options = [verbatim(ft.DropdownOption(key=self.identifier(entry),text=self.option_text(entry)),'text')
                              for entry in visible[:COMBO_LIMIT]]
        if self.combo.value not in {option.key for option in self.combo.options}:
            self.combo.value = None
        self.combo.disabled = not self.settings.form.writable or not visible
        self.include_one.disabled = self.combo.disabled
        self.include_list.content = f'Incluir todos ({len(visible)})' if visible else 'Incluir todos'
        self.include_list.disabled = self.combo.disabled
        matches = self.matches()
        shown = matches[:SUGGESTION_ROWS]
        self.suggestions.controls = [self.suggestion(entry) for entry in shown]
        term = (self.search.value or '').strip()
        if term and not matches:
            self.suggestion_count.value = f'Nenhum {self.text["noun"]} disponível corresponde à pesquisa.'
        elif len(matches) > len(shown):
            self.suggestion_count.value = (f'Mostrando {len(shown)} de {len(matches)} sugestões; '
                                           'continue digitando para refinar.')
        else:
            self.suggestion_count.value = '1 sugestão' if len(matches) == 1 else f'{len(matches)} sugestões'
        self.suggestion_box.visible = bool(term)
        if not self.entries:
            self.available.value = self.text['connect'] if not self.origin_values() else self.text['none']
        else:
            hidden = len(visible)-len(self.combo.options)
            self.available.value = (f'{len(self.entries)} {self.text["plural"]} na origem · '
                                    f'{len(visible)} disponíveis para incluir'
                                    +(f' · a lista mostra {len(self.combo.options)}; pesquise para ver os demais'
                                      if hidden > 0 else ''))

    def rows(self):
        """Table rows: one per target, and one per ligand-receptor complex.

        A compound is validated once for every included target, so the table
        is already the list of complexes the execution will analyse.
        """
        known = self.known()
        if self.kind == 'targets':
            return [(identifier,known.get(identifier,{}),None,None) for identifier in self.selected]
        receptors = list(self.settings.selected_targets)
        rows = []
        for identifier in self.selected:
            entry = known.get(identifier,{})
            # Without a target selection every target of the compound is validated.
            chosen = receptors or list(entry.get('receptors') or [None])
            for receptor in chosen:
                if self.settings.is_excluded(identifier,receptor):
                    continue
                rows.append((identifier,entry,receptor,self.settings.pair_summary(identifier,receptor)))
        return rows

    def row_cells(self, identifier, entry, receptor=None, summary=None):
        if self.kind == 'targets':
            return [entry.get('pdb_code') or identifier,entry.get('chain') or '—',entry.get('ligand') or '—',
                    (rounded(entry.get('resolution'),2)+' Å') if entry.get('resolution') else '—',
                    (rounded(entry.get('rmsd'),3)+' Å') if entry.get('rmsd') else '—']
        summary = summary or {}
        score = summary.get('best_score')
        return [identifier,receptor or '—','/'.join(sorted(summary.get('engines') or ())) or '—',
                rounded(score,3) if score is not None else '—',str(summary.get('poses',0))]

    def draw_table(self):
        known = self.known()
        entries = self.rows()
        if self.focus is not None:
            index = next((position for position,row in enumerate(entries) if row[0] == self.focus),None)
            self.offset = (index or 0)//PAGE_SIZE*PAGE_SIZE
            self.focus = None
        if self.offset >= len(entries):
            self.offset = max(0,(len(entries)-1)//PAGE_SIZE*PAGE_SIZE)
        rows = []
        writable = self.settings.form.writable
        for identifier,entry,receptor,summary in entries[self.offset:self.offset+PAGE_SIZE]:
            if self.kind == 'targets':
                view = action_chip('3D','Visualizar target em 3D',
                    self.settings.view_structure(entry['file']) if entry.get('file') else None)
            else:
                # The complex, as the docking results present it; only when no
                # receptor is available does the pose alone open.
                opener = self.settings.view_complex(identifier,receptor,summary)
                if opener is None:
                    # Without a published receptor, open the pose of the engine
                    # that would be shown — never another engine's pose.
                    structure = self.settings.engine_pose(summary) or entry.get('prepared_file')
                    opener = self.settings.view_structure(structure) if structure else None
                view = action_chip('3D','Visualizar complexo receptor-ligante em 3D',opener,
                    unavailable='Nenhuma pose deste complexo está disponível nos resultados do pipeline.')
                card = action_chip('2D','Abrir card 2D do composto',
                                   self.settings.open_card(identifier,receptor))
            label = ('Remover este complexo' if self.kind == 'compounds' and receptor
                     else 'Remover '+self.text['noun'])
            remove = ft.IconButton(ft.Icon(ft.Icons.DELETE_OUTLINE,semantics_label=label),
                tooltip=label,disabled=not writable,icon_size=20,
                on_click=lambda e,value=identifier,target=receptor:self.remove(value,target))
            incomplete = identifier not in known or (self.kind == 'compounds' and not (summary or {}).get('poses'))
            rows.append(ft.DataRow(color='#FEF3C7' if incomplete else None,cells=[
                *[ft.DataCell(verbatim(ft.Text(str(value),size=13,max_lines=2,tooltip=str(value),
                                               overflow=ft.TextOverflow.ELLIPSIS)))
                  for value in self.row_cells(identifier,entry,receptor,summary)],
                ft.DataCell(ft.Row([view,card,remove] if self.kind == 'compounds' else [view,remove],spacing=4))]))
        self.pair_table.set_rows(rows)
        self.table_box.visible = bool(rows)
        self.empty.visible = not rows
        missing = sum(1 for identifier in self.selected if identifier not in known)
        empty_pairs = sum(1 for row in entries if self.kind == 'compounds' and not (row[3] or {}).get('poses'))
        self.selected_count.value = (f'{len(self.selected)} incluídos' if self.selected else '')+(
            f' · {len(entries)} complexos ligante-receptor' if self.kind == 'compounds' and entries else '')+(
            f' · {missing} sem resultado na origem atual' if missing and self.entries else '')+(
            f' · {empty_pairs} sem pose no target' if empty_pairs else '')
        self.clear_button.visible = bool(self.selected)
        self.pager.visible = len(entries) > PAGE_SIZE
        self.previous.disabled = self.offset == 0
        self.next.disabled = self.offset+PAGE_SIZE >= len(entries)

    def draw(self):
        self.draw_origin(); self.draw_choices(); self.draw_table()


class PoseSettings:
    def __init__(self, form):
        self.form = form
        parameters = form.stage['parameters']
        self.selected_targets = [str(v) for v in (parameters.get('selected_receptors') or [])]
        self.selected_compounds = [str(v) for v in (parameters.get('selected_ligands') or [])]
        # Complexes dropped by their own row, without losing the compound.
        self.excluded_pairs = [str(v) for v in (parameters.get('excluded_pairs') or [])]
        self.records = []
        self.origin_records = []
        self.pair_index = {}
        self.error = ft.Text(color='#B91C1C',visible=False,size=12)
        self.plan = ft.Container(ft.Text(size=12,color='#1E3A8A'),bgcolor='#EFF6FF',border_radius=10,
                                 padding=ft.Padding.symmetric(horizontal=12,vertical=10),visible=False)
        # The column shares keep every section inside the dialog width; only
        # the actions column is fixed, sized for its buttons.
        self.targets = PosePanel(self,'targets','base_input_path',
            ('PDB','Cadeia','Ligante de referência','Resolução','RMSD','Ações'),(3,2,4,3,3),100)
        self.compounds = PosePanel(self,'compounds','base_selected_mols',
            ('Composto','Target','Motores','Melhor score','Poses','Ações'),(5,3,3,3,2),146)
        self.targets_control = ft.Column([
            ft.Text('Escolha os PDBs preparados que o redocking usou para obter o RMSD. '
                    'A lista e a pesquisa usam o código PDB.',size=12),
            self.targets.control,self.error],spacing=12,data='targets',
            horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
        self.compounds_control = ft.Column([
            ft.Text('Escolha as moléculas que geraram os resultados de docking. Serão validadas as poses '
                    'dos compostos incluídos nos targets incluídos.',size=12),
            self.compounds.control,self.plan],spacing=12,data='compounds',
            horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
        form.readers['selected_receptors'] = lambda: list(self.selected_targets) or None
        form.readers['selected_ligands'] = lambda: list(self.selected_compounds) or None
        form.readers['excluded_pairs'] = lambda: self.visible_exclusions() or None
        for key in LABELS:
            control = form.field_controls.get(key)
            if control is not None and key in HELP:
                control.tooltip = HELP[key]
                if isinstance(control,ft.TextField):
                    control.helper = HELP[key]
        for editor in form.input_editors.values():
            # The origin fields of the sections drive these inputs; the generic
            # source rows would list the producers' output files instead.
            editor.control.visible = False
            editor.on_change = self.refresh
        self.refresh()

    # --------------------------------------------------------------- catalog
    def candidate_stage(self):
        """The stage as currently edited, so the catalog follows the connections."""
        stage = copy.deepcopy(self.form.stage)
        for field,editor in self.form.input_editors.items():
            references = editor.read(validate=False)
            if references:
                stage.setdefault('bindings',{})[field] = references
            else:
                stage.setdefault('bindings',{}).pop(field,None)
            direct = editor.direct_path()
            if direct:
                stage['parameters'][field] = direct
            else:
                stage['parameters'].pop(field,None)
        stage['parameters'].update(selected_receptors=list(self.selected_targets),
                                   selected_ligands=list(self.selected_compounds),
                                   excluded_pairs=list(self.excluded_pairs))
        return stage

    def inputs(self):
        from biomolexplorer.pose_inputs import PoseInputs
        resolver = getattr(self,'_resolver',None)
        if resolver is None:
            resolver = self._resolver = PoseInputs(self.form.ui.store,self.form.ui.current['id'],
                pipeline=(self.form.ui.current or {}).get('pipeline'))
        return resolver

    def load(self, prune=False):
        """Read the catalog of the current origins; safe off the UI thread."""
        ui = self.form.ui
        self.error.visible = False
        self.records = []; self.origin_records = []; self.pair_index = {}
        self.targets.entries = []; self.compounds.entries = []
        if not (getattr(ui,'store',None) and (ui.current or {}).get('id')):
            return
        from biomolexplorer.workspace import AccessDenied
        try:
            stage = self.candidate_stage()
            resolver = self.inputs()
            self.targets.entries = resolver.targets(stage)
            if prune and self.targets.entries:
                available = {entry['receptor_id'] for entry in self.targets.entries}
                self.selected_targets[:] = [value for value in self.selected_targets if value in available]
            # The lists follow the connected origin; the complexes of the table
            # follow the same resolution the execution uses.
            self.origin_records = resolver.records(stage)
            self.records = resolver.combined(stage)
            self.pair_index = self.build_pairs(self.records)
            self.compounds.entries = resolver.compounds(stage,self.selected_targets)
            if prune and self.compounds.entries:
                available = {entry['compound_id'] for entry in self.compounds.entries}
                self.selected_compounds[:] = [value for value in self.selected_compounds if value in available]
        except (ValueError,OSError,AccessDenied) as error:
            # A draft may reference results that are not available yet;
            # execution still validates the selection strictly.
            self.error.value = str(error); self.error.visible = True

    @staticmethod
    def build_pairs(records):
        """Engines, poses, best score and a pose file of each ligand-receptor complex."""
        index = {}
        for record in records:
            compound = str(record.get('compound_id') or record.get('ligand_id') or '').strip()
            if not compound:
                continue
            entry = index.setdefault((compound,str(record.get('receptor_id') or '')),
                                     {'engines':set(),'poses':0,'best_score':None,'pose':None,'by_engine':{}})
            entry['poses'] += 1
            engine = str(record.get('engine') or '')
            if engine:
                entry['engines'].add(engine)
            if entry['pose'] is None and record.get('pose') and Path(record['pose']).is_file():
                entry['pose'] = record['pose']
            try:
                score = float(record.get('score'))
            except (TypeError,ValueError):
                score = None
            if score is not None and (entry['best_score'] is None or score < entry['best_score']):
                entry['best_score'] = score
            # The best pose of each engine: the 3D view opens the complex, and
            # Vina comes first because its scene also classifies the residues.
            best = entry['by_engine'].setdefault(engine or 'score',
                {'pose':None,'receptor':None,'reference':None,'score':None,'smiles':None,'poses':0})
            best['poses'] += 1
            if best['pose'] is None or (score is not None and (best['score'] is None or score < best['score'])):
                best.update(pose=record.get('pose'),receptor=record.get('receptor'),
                            reference=record.get('reference'),score=score,
                            smiles=record.get('canonical_smiles'))
        return index

    def pair_summary(self, compound, receptor, index=None):
        """Summary of one complex; a PDB code without chain matches its chains."""
        if receptor is None:
            return None
        index = self.pair_index if index is None else index
        exact = index.get((compound,receptor))
        if exact is not None:
            return exact
        from biomolexplorer import pose_catalog
        found = [value for (code,identity),value in index.items()
                 if code == compound and pose_catalog.matches_target([receptor],identity)]
        if not found:
            return None
        merged = {'engines':set(),'poses':0,'best_score':None,'pose':None,'by_engine':{}}
        for value in found:
            merged['engines'] |= value['engines']; merged['poses'] += value['poses']
            merged['pose'] = merged['pose'] or value['pose']
            for engine,chosen in value['by_engine'].items():
                current = merged['by_engine'].get(engine)
                if current is None or (chosen['score'] is not None
                                       and (current['score'] is None or chosen['score'] < current['score'])):
                    merged['by_engine'][engine] = chosen
            if value['best_score'] is not None and (merged['best_score'] is None
                                                    or value['best_score'] < merged['best_score']):
                merged['best_score'] = value['best_score']
        return merged

    # ------------------------------------------------------------ exclusions
    @staticmethod
    def pair_key(compound, receptor):
        from biomolexplorer.pose_catalog import complex_key
        return complex_key(compound,receptor)

    def is_excluded(self, compound, receptor):
        return receptor is not None and self.pair_key(compound,receptor) in set(self.excluded_pairs)

    def exclude_pair(self, compound, receptor):
        """Drop one complex; a compound without any complex left is dropped too."""
        key = self.pair_key(compound,receptor)
        if key not in self.excluded_pairs:
            self.excluded_pairs.append(key)
        targets = self.selected_targets or [
            entry['receptor_id'] for entry in self.targets.entries]
        if all(self.is_excluded(compound,target) for target in targets):
            self.selected_compounds[:] = [value for value in self.selected_compounds if value != compound]
            self.drop_exclusions(compound=compound)

    def drop_exclusions(self, compound=None, receptor=None):
        """Forget exclusions of a compound or target that left the selection."""
        def keep(key):
            identifier,_,target = key.rpartition(':')
            return not (compound is not None and identifier == compound
                        or receptor is not None and target == receptor)
        self.excluded_pairs[:] = [key for key in self.excluded_pairs if keep(key)]

    def visible_exclusions(self):
        """Only exclusions of complexes the current selection still offers."""
        targets = set(self.selected_targets)
        compounds = set(self.selected_compounds)
        kept = []
        for key in self.excluded_pairs:
            compound,_,target = key.rpartition(':')
            if compound in compounds and (not targets or target in targets):
                kept.append(key)
        return kept

    def engine_rows(self, compound, receptor):
        """One row per docking engine of the pipeline for this complex.

        The complexes come from the same resolution the execution uses, where
        the blocks that produced each pose are authoritative.  Every engine
        always has its row: one without a result for the complex is reported
        with no values instead of being hidden.
        """
        summary = self.pair_summary(compound,receptor) or {}
        rows = []
        for engine in ENGINE_ORDER:
            chosen = (summary.get('by_engine') or {}).get(engine)
            rows.append({'engine':engine,'receptor':receptor,
                         'score':chosen['score'] if chosen else None,
                         'poses':chosen['poses'] if chosen else None,
                         'pose':chosen.get('pose') if chosen else None,
                         'summary':summary if chosen else None,'available':bool(chosen)})
        return rows

    def refresh(self):
        self.load(); self.draw()

    async def reload(self, prune=False):
        """Reload after an origin change without blocking the dialog."""
        ui = self.form.ui
        call = getattr(ui,'call',None)
        if call is not None:
            await call(self.load,prune)
        else:
            self.load(prune)
        self.draw()

    def selection_changed(self, kind):
        if kind == 'targets':
            # Compounds only make sense for the chosen targets.
            from biomolexplorer import pose_catalog
            self.compounds.entries = pose_catalog.compounds(self.origin_records,self.selected_targets)
            keep = {entry['compound_id'] for entry in self.compounds.entries}
            if self.selected_targets and self.targets.entries:
                self.selected_compounds[:] = [code for code in self.selected_compounds if code in keep]
        self.draw()

    def draw(self):
        self.targets.draw(); self.compounds.draw()
        text = self.plan_text()
        self.plan.content.value = text; self.plan.visible = bool(text)
        self.patch(self.targets_control,self.compounds_control)

    def plan_text(self):
        from biomolexplorer import pose_catalog
        by_id = {target['receptor_id']:target for target in self.targets.entries}
        try:
            selected = pose_catalog.selected_poses(self.records,by_id,self.selected_targets,
                                                   self.selected_compounds,self.visible_exclusions())
        except (OSError,ValueError):
            return ''
        if not selected:
            return ('Nenhuma pose corresponde à seleção atual. Revise os targets e compostos.'
                    if self.records else '')
        engines = sorted({str(record.get('engine') or '—') for record in selected})
        removed = len(self.visible_exclusions())
        return (f'{len(selected)} poses serão validadas · motores: '+', '.join(engines)+
                f' · receptores: {len({r["receptor_id"] for r in selected})}'+
                (f' · {removed} complexos removidos' if removed else ''))

    # --------------------------------------------------------------- viewers
    def view_structure(self, path):
        async def action():
            await self.form.ui.preview_artifact(self.form.ui.current['id'],path)

        async def handler(event):
            await self.form.ui.guard(action)
        return handler

    @staticmethod
    def engine_order(summary, engine=None):
        """Engines to try: only the named one, or Vina first and then the rest.

        A named engine is never replaced by another: a row of the card offers
        the pose of its own engine or nothing at all.
        """
        by_engine = (summary or {}).get('by_engine') or {}
        if engine is not None:
            return [engine] if engine in by_engine else []
        order = [name for name in ENGINE_ORDER if name in by_engine]
        return order+[name for name in sorted(by_engine) if name not in order]

    def engine_pose(self, summary, engine=None):
        """Published pose file of that engine, when the receptor is missing."""
        by_engine = (summary or {}).get('by_engine') or {}
        for name in self.engine_order(summary,engine):
            pose = by_engine[name].get('pose')
            if pose and Path(pose).is_file():
                return pose
        return None

    def complex_scene(self, compound, receptor, summary, engine=None):
        """Receptor, pose and reference of one complex, Vina first.

        The scene is the one the docking results open: the viewer computes the
        residue contacts and, for a Vina pose, the interaction of each residue.
        Each engine docks and minimizes its own geometry, so naming an engine
        opens that engine's pose and receptor and never another one's.
        """
        from biomolexplorer.pose_catalog import _preference
        by_engine = (summary or {}).get('by_engine') or {}
        order = self.engine_order(summary,engine)
        target = next((entry for entry in self.targets.entries
                       if entry['receptor_id'] == receptor),None) or {}
        prepared = target.get('file')
        for engine in order:
            chosen = by_engine[engine]
            pose = chosen.get('pose')
            if not pose or not Path(pose).is_file():
                continue
            structure = chosen.get('receptor')
            if not structure or not Path(structure).is_file() or (
                    prepared and Path(prepared).is_file() and _preference(prepared) < _preference(structure)):
                structure = prepared
            if not structure or not Path(structure).is_file():
                continue
            layers = [{'role':'receptor','path':str(structure)},{'role':'pose','path':str(pose)}]
            reference = chosen.get('reference')
            if reference and Path(reference).is_file():
                layers.append({'role':'reference','path':str(reference)})
            return dict(layers=layers,engine=engine,smiles=chosen.get('smiles'),
                        reference_kind='prepared' if len(layers) > 2 else None,
                        name=' · '.join(filter(None,[compound,receptor,engine.upper()])))
        return None

    def view_complex(self, compound, receptor, summary, engine=None):
        scene = self.complex_scene(compound,receptor,summary,engine)
        if scene is None:
            return None

        async def action():
            await self.form.ui.preview_complex(self.form.ui.current['id'],scene['layers'],scene['name'],
                                               ligand_smiles=scene.get('smiles'),
                                               reference_kind=scene.get('reference_kind'))

        async def handler(event):
            await self.form.ui.guard(action)
        return handler

    def compound_records(self, identifier):
        return [record for record in self.records
                if str(record.get('compound_id') or record.get('ligand_id') or '') == identifier]

    def card_data(self, identifier, receptor=None):
        """Identity, properties and docking results of one compound (worker thread).

        ``results`` holds the complexes of the target the card was opened for,
        and ``others`` the remaining targets of the compound: the card never
        presents another target's complex as the one of its own row.
        """
        from biomolexplorer.compound_properties import compound_properties
        records = self.compound_records(identifier)
        entry = next((e for e in self.compounds.entries if e['compound_id'] == identifier),{})
        smiles = entry.get('smiles') or next((r.get('canonical_smiles') for r in records if r.get('canonical_smiles')),'')
        properties = {}
        structure = {}
        if smiles:
            try:
                properties = card_decimals(compound_properties(
                    {'molecule_chembl_id':identifier,'canonical_smiles':smiles},calculate=True))
            except (ImportError,ValueError,RuntimeError):
                properties = {}
            try:
                # The 2D depiction and the Bemis-Murcko core of this molecule.
                from biomolexplorer.network_scaffolds import individual
                structure = individual(smiles)
            except (ImportError,ValueError,RuntimeError):
                structure = {}
        if receptor:
            results = self.engine_rows(identifier,receptor)
        else:
            results = [{'engine':engine,'receptor':target,'score':chosen['score'],
                        'poses':chosen['poses'],'pose':chosen['pose'],'summary':summary,'available':True}
                       for (compound,target),summary in sorted(self.pair_index.items())
                       if compound == identifier
                       for engine,chosen in sorted(summary['by_engine'].items())]
        pose = next((item['pose'] for item in results
                     if item['pose'] and Path(item['pose']).is_file()),None)
        return dict(entry=entry,records=records,smiles=smiles,properties=properties,structure=structure,
                    results=results,pose=pose)

    def panel(self, title, *controls, subtitle=None):
        """Card section: one title, an optional note and its own white surface."""
        return ft.Container(padding=14,bgcolor='#FFFFFF',border_radius=12,border=ft.Border.all(1,'#E2E8F0'),
            content=ft.Column([ft.Text(title,size=13,weight=ft.FontWeight.W_600,color='#172B4D'),
                *([ft.Text(subtitle,size=11,color=MUTED)] if subtitle else []),*controls],
                spacing=10,horizontal_alignment=ft.CrossAxisAlignment.STRETCH))

    def result_table(self, identifier, results, receptor):
        """Score, poses and complex of every engine and target of the compound.

        Each row opens its own engine's pose: DOCK6 refines the geometry Vina
        found, so the two poses of a complex are different structures.
        """
        best = min((item['score'] for item in results if item['score'] is not None),default=None)
        table = SelectedItemsTable(('Motor','Target','Score','Poses','3D'),(3,4,3,2),52)
        rows = []
        for item in results:
            available = item.get('available',True)
            score = item['score']
            values = [item['engine'] or '—',
                      (item['receptor'] or '—') if available else '—',
                      rounded(score,3) if score is not None else '—',
                      str(item['poses']) if item['poses'] is not None else '—']
            opener = self.view_complex(identifier,item['receptor'],item['summary'],item['engine'])
            cells = [ft.DataCell(verbatim(ft.Text(value,size=12,max_lines=1,tooltip=value,
                color=MUTED if not available else None,overflow=ft.TextOverflow.ELLIPSIS,
                weight=ft.FontWeight.W_600 if available and score is not None and score == best else None)))
                for value in values]
            cells.append(ft.DataCell(action_chip('3D','Visualizar complexo receptor-ligante em 3D',opener,
                unavailable='Este motor não tem resultado para o complexo nos blocos de docking do pipeline.')))
            # An engine without results in the pipeline keeps its row, in grey.
            rows.append(ft.DataRow(cells=cells,color=None if available else '#F1F5F9'))
        table.set_rows(rows)
        return table.control

    def open_card(self, identifier, receptor=None):
        async def action():
            from .molecular_cards import badge, image_panel, metric, property_sections, quick_metrics
            from .settings_layout import dialog_size, fit_dialog
            data = await self.form.ui.call(self.card_data,identifier,receptor)
            entry,properties,structure = data['entry'],data['properties'],data['structure']
            targets = entry.get('receptors') or []

            def shown(key):
                value = properties.get(key)
                return '—' if value in (None,'') else value
            badges = [badge('ID',identifier)]
            if receptor:
                # The complex this card was opened for, never another target.
                badges.append(badge('Target',receptor))
            if properties.get('Molecular_Formula'):
                badges.append(badge('Fórmula',properties['Molecular_Formula']))
            badges.append(badge('Targets',len(targets) or '—'))
            badges.append(badge('Poses',entry.get('poses','—')))
            for engine in entry.get('engines',[]):
                badges.append(badge('Motor',engine))
            # The best docking score of this complex, not the aggregate of the
            # connected origin, which may be a consensus relevance.
            scores = [item['score'] for item in data['results'] if item['score'] is not None]
            best = min(scores) if scores else (None if receptor else entry.get('best_score'))
            if best is not None:
                badges.append(badge('Melhor score',rounded(best,3)))
            header = ft.Column([ft.Row([ft.Icon(ft.Icons.SCIENCE_OUTLINED,color='#2552E8',size=22),
                verbatim(ft.Text(identifier,size=20,weight=ft.FontWeight.W_700,selectable=True,expand=True))],spacing=8),
                *([verbatim(ft.Text(str(properties['pref_name']),size=13,color='#475569',selectable=True))]
                  if properties.get('pref_name') else []),
                ft.Row(badges,wrap=True,spacing=6,run_spacing=6)],spacing=8)
            drawing = self.panel('Estrutura 2D',image_panel(structure.get('image'),height=240,
                notice=structure.get('notice') or 'Representação 2D indisponível: a origem não forneceu SMILES.'),
                subtitle='Desenho gerado com RDKit a partir do SMILES dos resultados.')
            core = self.panel('Núcleo Bemis–Murcko',
                ft.Image(src=structure['scaffold'],height=130,fit=ft.BoxFit.CONTAIN),
                subtitle=f'{structure.get("atoms",0)} átomos no núcleo de anéis e conectores.'
                ) if structure.get('scaffold') else None
            geometry = [metric('Doadores de H',shown('HBD'),col={'xs':6,'md':4}),
                        metric('Aceptores de H',shown('HBA'),col={'xs':6,'md':4}),
                        metric('Ligações rotacionáveis',shown('RotatableBonds'),col={'xs':6,'md':4}),
                        metric('Átomos pesados',shown('HeavyAtomCount'),col={'xs':6,'md':4}),
                        metric('Carga formal',shown('FormalCharge'),col={'xs':6,'md':4}),
                        metric('Targets com resultado',len(targets) or '—',col={'xs':6,'md':4})]
            chemistry = self.panel('Propriedades físico-químicas',quick_metrics(properties,col={'xs':12,'md':4}),
                                   ft.ResponsiveRow(geometry,spacing=8,run_spacing=8))
            results = self.panel(
                f'Resultados de docking em {receptor}' if receptor else 'Resultados de docking',
                self.result_table(identifier,data['results'],receptor),
                subtitle='Melhor score de cada motor neste complexo. Cada motor abre a sua própria pose; '
                         'um motor sem resultado no pipeline fica sem valores.'
                         if receptor else 'Melhor score de cada motor e target. '
                                          'Cada motor abre a sua própria pose.') if data['results'] else None
            # Two rows of two panels: the drawing beside the properties, and the
            # Bemis-Murcko core beside the docking results.
            columns = []
            for first,second in ((drawing,chemistry),(core,results)):
                cells = []
                for panel,width,other in ((first,5,second),(second,7,first)):
                    if panel is not None:
                        cells.append(ft.Container(panel,col={'xs':12,'md':12 if other is None else width}))
                if cells:
                    columns.append(ft.ResponsiveRow(cells,spacing=12,run_spacing=12))
            sections = property_sections(properties,{}) if properties else []
            content = [header,*columns]
            if sections:
                content.append(self.panel('Identidade e demais atributos',*sections,
                    subtitle='Cada valor traz a origem: provedor do resultado ou cálculo RDKit.'))
            # The complex of the row, from whichever block published it.
            summary = next((item['summary'] for item in data['results'] if item.get('available')),
                           None) or self.pair_summary(identifier,receptor)
            scene = self.complex_scene(identifier,receptor,summary)
            fallback = self.engine_pose(summary) or data['pose']
            actions = [ft.TextButton('Fechar',on_click=lambda e:close_dialog(self.form.ui.page,dialog))]
            if scene is not None:
                actions.insert(0,ft.TextButton('Abrir complexo em 3D',icon=ft.Icons.VIEW_IN_AR,
                                               on_click=self.view_complex(identifier,receptor,summary)))
            elif fallback:
                actions.insert(0,ft.TextButton('Abrir pose em 3D',icon=ft.Icons.VIEW_IN_AR,
                                               on_click=self.view_structure(fallback)))
            dialog = ft.AlertDialog(title=ft.Text('Detalhes do composto',size=20),
                content=ft.Container(width=dialog_size(self.form.ui.page)['width'],bgcolor='#F8FAFC',
                    padding=ft.Padding.all(4),content=ft.Column(content,spacing=14,tight=True,
                    scroll=ft.ScrollMode.AUTO,horizontal_alignment=ft.CrossAxisAlignment.STRETCH)),actions=actions)
            fit_dialog(dialog,self.form.ui.page)
            self.form.ui.page.show_dialog(dialog)

        async def handler(event):
            await self.form.ui.guard(action)
        return handler

    # ----------------------------------------------------------------- patch
    def patch(self, *controls):
        """Patch only these controls; the dialog never repaints the page."""
        for control in controls:
            update = getattr(control,'update',None)
            if callable(update):
                try:
                    update()
                except (AssertionError,AttributeError,RuntimeError):
                    # Not mounted yet (first build): the initial render carries the state.
                    continue
