"""Explicit runtime messages shared by workers and session-owned presenters."""
import json
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).with_name('resources') / 'i18n'


@lru_cache(maxsize=1)
def locales():
    """Locale metadata and dictionaries are packaged together; no UI branches."""
    return json.loads((ROOT/'locales.json').read_text(encoding='utf-8'))


def language_info(language):
    if language not in locales():
        raise ValueError('Unknown language: '+str(language))
    return locales()[language]


@lru_cache(maxsize=32)
def dictionary(language,kind='phrases'):
    language_info(language)
    filename=language+'.json' if kind=='phrases' else ('patterns.json' if kind=='patterns' and language=='en' else kind+'_'+language+'.json')
    path=ROOT/filename
    return json.loads(path.read_text(encoding='utf-8')) if path.is_file() else {}


@lru_cache(maxsize=16)
def localized_dictionary(language,kind):
    # Incomplete new locales fall back to English, never implicitly Portuguese.
    return {**dictionary('en',kind),**dictionary(language,kind)}


def runtime_catalog(language):
    return localized_dictionary(language,'runtime')


class RuntimeMessages:
    def __init__(self, language):
        self.language = language
        self.catalog = runtime_catalog(language)

    def __call__(self, key, **arguments):
        # Missing keys fail explicitly instead of silently changing language.
        return self.catalog[key].format(**arguments)
