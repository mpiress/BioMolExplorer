"""Session-owned localization at the presentation boundary.

Translate labels without changing persisted pipelines, dropdown keys, editable
values, molecular data or backend requests. No process-global language state.
"""
import re
from functools import lru_cache

import flet as ft
import flet.canvas as canvas
from biomolexplorer.i18n import ROOT,dictionary,locales,language_info

CATALOG = ROOT / 'en.json'
PATTERNS = CATALOG.with_name('patterns.json')
LANGUAGES = tuple(locales())

# Translate system-message payloads while keeping filenames, IDs and user titles intact.
SYSTEM_ARGUMENTS = {
    'Conecte {0}': {0},
    'Anéis com anéis: {0}': {0},
    'Anéis completos: {0}': {0},
    'Tempo máximo: {0} s': {0},
    '{0}Resultado reutilizado do pipeline; os critérios originais estão preservados.': {0},
    '{0} · {1} Todas as métricas ignoram escores como distância. Centralidade não indica atividade biológica.': {1},
    '{0} · Cor: {1} · Espessura: magnitude do escore disponível.': {0,1},
    'Cinza: indisponível · Contorno azul: MCC · {0}': {0},
    '{0} · {1}: {2} · Comunidade: {3}': {2,3},
    'Comunidade {0} · {1}': {1},
    'Subconjunto filtrado · {0}': {0},
    '{0} · {1}: {2} / {3} valores · Pearson r: {4} · associação exploratória, sem inferência causal.': {4},
    'O bloco “{0}” precisa de ajustes: {1}': {1},
    'Arquivo {0} ({1}): {2}': {2},
    'Falha no envio: {0}': {0},
    'Valor inválido para {0}': {0},
    '{0} deve indicar uma pasta do projeto.': {0},
    'A pasta de {0} não está disponível. Escolha uma entrada do projeto.': {0},
    'Conecte um bloco de origem ou selecione seus arquivos em: {0}': {0},
    'Selecione ao menos um arquivo para {0}.': {0},
    'Selecione ao menos um arquivo para {0}': {0},
    'A entrada {0} do bloco “{1}” está conectada a dados incompatíveis.': {0},
    '{0} (itens separados por vírgula)': {0},
    'Adicionar item a {0}': {0},
    '{0} - {1}': {1},
}



@lru_cache(maxsize=16)
def catalogs(language='en'):
    language_info(language)
    english=dictionary('en');target=dictionary(language)
    phrases={source:(source if language=='pt' else translated) for source,translated in english.items()}
    phrases.update(target)
    for source,translated in english.items():phrases.setdefault(translated,phrases[source])
    canonical=dictionary('en','patterns');localized=dictionary(language,'patterns')
    mapping={}
    for source,translated in canonical.items():
        result=localized.get(source,localized.get(translated,source if language=='pt' else translated))
        if result!=source:mapping[source]=(result,source)
        if result!=translated:mapping[translated]=(result,source)
    for source,translated in localized.items():
        if source not in mapping and source!=translated:mapping[source]=(translated,source)
    patterns = []
    for source, (translated,canonical_source) in mapping.items():
        parts = re.split(r'(\{\d+\})', source)
        groups=set();expression=''
        for part in parts:
            if re.fullmatch(r'\{\d+\}',part):
                name='p'+part[1:-1]
                expression+=f'(?P={name})' if name in groups else f'(?P<{name}>.*?)'
                groups.add(name)
            else:expression+=re.escape(part)
        # Prefer specific sentences to short prefixes such as "File {0}".
        specificity = sum(len(part) for part in parts if not re.fullmatch(r'\{\d+\}', part))
        patterns.append((specificity, re.compile(expression, re.DOTALL), translated, canonical_source))
    return phrases, sorted(patterns, key=lambda value: -value[0])


def verbatim(control, *attributes):
    """Protect user content and scientific text from presentation translation."""
    control._biomol_verbatim = set(attributes) if attributes else True
    return control


@lru_cache(maxsize=1)
def default_stage_names():
    from biomolexplorer.catalog import TITLES
    return {dictionary(code).get(title,title) for title,_,_ in TITLES.values() for code in locales()} | {title for title,_,_ in TITLES.values()} | {'Recuperar compostos','Retrieve compounds','Preparar meus complexos','Prepare my complexes'}


def stage_control(control,stage,attribute='value'):
    """Localize catalog stage names and preserve names entered by the user."""
    if stage.get('name') in default_stage_names():
        control._biomol_verbatim=set()
        return control
    return verbatim(control,attribute)


class Translator:
    def __init__(self, language='pt'):
        # Canvas updates revisit the same labels. Bound the session cache so
        # dynamic progress messages cannot accumulate indefinitely.
        self._translate=lru_cache(maxsize=4096)(self._translate)
        self.set_language(language)

    def set_language(self, language):
        language_info(language)
        self.language = language
        self._translate.cache_clear()

    def __call__(self, value):
        if not isinstance(value, str) or not value:
            return value
        return self._translate(value)

    def _translate(self, value):
        if value in ('Preparar meus complexos','Prepare my complexes'):
            value='Preparar para docking'
        phrases, patterns = catalogs(self.language)
        if value in phrases:
            return phrases[value]
        # Progress/error messages can contain several independently translated
        # lines, including messages emitted by a separate scientific worker.
        if '\n' in value:
            return '\n'.join(self(line) for line in value.split('\n'))
        for prefix in ('Etapa atual: ', 'Padrão esperado: ', 'Current stage: ', 'Expected format: '):
            if value.startswith(prefix):
                remainder=value[len(prefix):]
                return phrases.get(prefix,prefix)+(self(remainder) if prefix in ('Padrão esperado: ','Expected format: ') or remainder in default_stage_names() else remainder)
        for _, pattern, translated, canonical in patterns:
            match = pattern.fullmatch(value)
            if match:
                # Interpolation data is preserved exactly; it can be a filename,
                # a project name, an identifier or a scientific value.
                system_slots=SYSTEM_ARGUMENTS.get(canonical,set())
                if 'Padrão esperado: ' in canonical:
                    suffix=canonical.split('Padrão esperado: ',1)[1]
                    system_slots=system_slots | {int(i) for i in re.findall(r'\{(\d+)\}',suffix)}
                def substitute(marker):
                    index=int(marker[1]);payload=match.group('p'+str(index))
                    stage_patterns=('bloco “','de ajustes:','Configurar {','Resultados de {','Conectar saída de {','Escolha uma entrada compatível para {','Etapa: {','Entrada · {')
                    default_name=payload in default_stage_names() and any(part in canonical for part in stage_patterns)
                    return self(payload) if index in system_slots or default_name else payload
                return re.sub(r'\{(\d+)\}', substitute, translated)
        # Validation errors prepend a user filename to another known system
        # message. Keep the prefix intact and translate only the error payload.
        if ': ' in value:
            prefix,_,payload=value.partition(': ')
            translated=self(payload)
            if translated!=payload:return prefix+': '+translated
        return value


class LocalizedPage:
    """Page facade that localizes newly built and subsequently updated controls."""
    _children = ('controls', 'content', 'title', 'subtitle', 'actions', 'leading',
        'trailing', 'label', 'error', 'helper', 'rows', 'cells', 'columns', 'options', 'spans',
        'shapes', 'tabs', 'tab_bar', 'body', 'items', 'badge', 'icon', 'menu',
        'prefix', 'suffix', 'prefix_icon', 'suffix_icon', 'secondary_label',
        'bottom_axis', 'top_axis', 'left_axis', 'right_axis', 'labels', 'titles', 'segments', 'tooltip',
        'action', 'destinations', 'leading_icon', 'trailing_icon', 'error_content')
    _labels = ('label', 'secondary_label', 'hint_text', 'error', 'helper', 'helper_text', 'tooltip', 'semantics_label', 'message',
        'on_tap_hint_text', 'on_long_press_hint_text')

    def __init__(self, page, language='pt'):
        object.__setattr__(self, '_page', page)
        object.__setattr__(self, 'translator', Translator(language))
        object.__setattr__(self, '_dialogs', [])
        self.set_language(language)

    def __getattr__(self, name):
        return getattr(self._page, name)

    def __setattr__(self, name, value):
        setattr(self._page, name, value)

    def set_language(self, language):
        self.translator.set_language(language)
        self._page.locale_configuration = ft.LocaleConfiguration(
            supported_locales=[ft.Locale(info['language_code'],info.get('country_code')) for info in locales().values()],
            current_locale=ft.Locale(language_info(language)['language_code'],language_info(language).get('country_code')))

    def localize(self, control, seen=None):
        if control is None or isinstance(control, (str, int, float, bool)):
            return
        seen = set() if seen is None else seen
        if id(control) in seen:
            return
        seen.add(id(control))
        if isinstance(control, (list, tuple)):
            for child in control:
                self.localize(child, seen)
            return
        if not isinstance(control, (ft.Control, canvas.Shape, ft.TextSpan, ft.Tooltip)):
            return
        # Flet owns immutable controls created by its component renderer.
        # Translating their runtime snapshots would abort the entire update,
        # including unrelated graphs and charts. Localize our mutable controls.
        if hasattr(control, '_frozen'):
            return
        protected = getattr(control, '_biomol_verbatim', set())
        saved = getattr(control, '_biomol_messages', {})
        attributes = list(self._labels)
        if isinstance(control, (ft.Text, canvas.Text)):
            attributes.append('value')
        if isinstance(control,ft.Semantics):
            attributes.extend(('value','increased_value','decreased_value'))
        if isinstance(control, (ft.TextSpan, ft.DropdownOption)):
            attributes.append('text')
        if isinstance(control, (ft.Button, ft.TextButton,ft.PopupMenuItem)):
            attributes.append('content')
        for name in attributes:
            if protected is True or name in protected:
                continue
            value = getattr(control, name, None)
            if not isinstance(value, str):
                continue
            original, rendered = saved.get(name, (value, None))
            if value != rendered:
                original = value
            translated = self.translator(original)
            saved[name] = (original, translated)
            if translated != value:
                setattr(control, name, translated)
        control._biomol_messages = saved
        for name in self._children:
            self.localize(getattr(control, name, None), seen)

    def add(self, *controls):
        self.localize(controls)
        return self._page.add(*controls)

    def show_dialog(self, dialog):
        from .settings_layout import fit_dialog
        fit_dialog(dialog,self._page)
        self.localize(dialog)
        result = self._page.show_dialog(dialog)
        if not any(item is dialog for item in self._dialogs):
            self._dialogs.append(dialog)
        return result

    def pop_dialog(self):
        dialog = self._page.pop_dialog()
        self._dialogs[:] = [item for item in self._dialogs if item is not dialog]
        return dialog

    def update(self, *controls):
        from .settings_layout import fit_dialog
        for dialog in self._dialogs:fit_dialog(dialog,self._page)
        seen = set()
        self.localize(controls or getattr(self._page, 'controls', []), seen)
        self.localize(getattr(self._page,'overlay',[]),seen)
        self.localize(self._dialogs, seen)
        result = self._page.update(*controls)
        self._dialogs[:] = [item for item in self._dialogs if item.open]
        return result
