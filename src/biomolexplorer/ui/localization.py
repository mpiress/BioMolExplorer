"""Session-owned localization at the presentation boundary.

Translate labels without changing persisted pipelines, dropdown keys, editable
values, molecular data or backend requests. No process-global language state.
"""
import json
import re
from functools import lru_cache
from pathlib import Path

import flet as ft
import flet.canvas as canvas

CATALOG = Path(__file__).resolve().parents[1] / 'resources' / 'i18n' / 'en.json'
PATTERNS = CATALOG.with_name('patterns.json')
LANGUAGES = ('en', 'pt')


@lru_cache(maxsize=1)
def catalogs():
    phrases = json.loads(CATALOG.read_text(encoding='utf-8'))
    patterns = []
    for source, translated in json.loads(PATTERNS.read_text(encoding='utf-8')).items():
        parts = re.split(r'(\{\d+\})', source)
        expression = ''.join('(.*?)' if re.fullmatch(r'\{\d+\}', part) else re.escape(part) for part in parts)
        # Prefer specific sentences to short prefixes such as "File {0}".
        specificity = sum(len(part) for part in parts if not re.fullmatch(r'\{\d+\}', part))
        patterns.append((specificity, re.compile(expression, re.DOTALL), translated))
    return phrases, sorted(patterns, key=lambda value: -value[0])


def verbatim(control, *attributes):
    """Protect user content and scientific text from presentation translation."""
    control._biomol_verbatim = set(attributes) if attributes else True
    return control


class Translator:
    def __init__(self, language='pt'):
        self.set_language(language)

    def set_language(self, language):
        if language not in LANGUAGES:
            raise ValueError('Language must be en or pt.')
        self.language = language

    def __call__(self, value):
        if self.language == 'pt' or not isinstance(value, str) or not value:
            return value
        phrases, patterns = catalogs()
        if value in phrases:
            return phrases[value]
        if '\nPadrão esperado: ' in value:
            message,_,expected=value.partition('\nPadrão esperado: ')
            return self(message)+'\nExpected format: '+self(expected)
        # Progress/error messages can contain several independently translated
        # lines, including messages emitted by a separate scientific worker.
        if '\n' in value:
            return '\n'.join(self(line) for line in value.split('\n'))
        for prefix in ('Etapa atual: ', 'Padrão esperado: '):
            if value.startswith(prefix):
                remainder=value[len(prefix):]
                return phrases.get(prefix,prefix)+(self(remainder) if prefix=='Padrão esperado: ' else remainder)
        for _, pattern, translated in patterns:
            match = pattern.fullmatch(value)
            if match:
                # Interpolation data is preserved exactly; it can be a filename,
                # a project name, an identifier or a scientific value.
                return re.sub(r'\{(\d+)\}', lambda m: match.group(int(m[1])+1), translated)
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
        'prefix', 'suffix', 'prefix_icon', 'suffix_icon')
    _labels = ('label', 'hint_text', 'error', 'helper', 'helper_text', 'tooltip', 'semantics_label', 'message')

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
            supported_locales=[ft.Locale('en', 'US'), ft.Locale('pt', 'BR')],
            current_locale=ft.Locale('en', 'US') if language == 'en' else ft.Locale('pt', 'BR'))

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
        if not isinstance(control, (ft.Control, canvas.Shape, ft.TextSpan)):
            return
        protected = getattr(control, '_biomol_verbatim', set())
        saved = getattr(control, '_biomol_messages', {})
        attributes = list(self._labels)
        if isinstance(control, (ft.Text, canvas.Text)):
            attributes.append('value')
        if isinstance(control, (ft.TextSpan, ft.DropdownOption)):
            attributes.append('text')
        if isinstance(control, (ft.Button, ft.TextButton)):
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
        seen = set()
        self.localize(controls or getattr(self._page, 'controls', []), seen)
        self.localize(self._dialogs, seen)
        result = self._page.update(*controls)
        self._dialogs[:] = [item for item in self._dialogs if item.open]
        return result
