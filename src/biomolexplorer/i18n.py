"""Explicit runtime messages shared by workers and session-owned presenters."""
import json
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).with_name('resources') / 'i18n'


@lru_cache(maxsize=2)
def runtime_catalog(language):
    if language not in ('pt', 'en'):
        raise ValueError('Language must be en or pt.')
    return json.loads((ROOT / f'runtime_{language}.json').read_text(encoding='utf-8'))


class RuntimeMessages:
    def __init__(self, language):
        self.language = language
        self.catalog = runtime_catalog(language)

    def __call__(self, key, **arguments):
        # Missing keys fail explicitly instead of silently changing language.
        return self.catalog[key].format(**arguments)
