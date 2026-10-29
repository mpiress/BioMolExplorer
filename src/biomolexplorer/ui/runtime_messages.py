"""Render worker feedback in the viewing session, never in a global language."""
import re
from biomolexplorer.i18n import RuntimeMessages
from .localization import Translator


def session_language(ui):
    return getattr(getattr(ui.page, 'translator', None), 'language', getattr(ui, 'language', 'pt'))


def stage_log(value, language):
    """Localize application messages while retaining paths and native tool output.

    Raw downloaded logs remain unchanged for scientific diagnostics. The worker's
    event IDs select session messages even for logs saved in a different language.
    """
    messages, translate = RuntimeMessages(language), Translator(language)
    lines = []
    for line in value.splitlines(keepends=True):
        event = re.match(r'^(.*?\bevent=([\w.]+) \[[^\n]*?\] )(.*?)(\n?)$', line)
        if event:
            prefix, name, payload, newline = event.groups()
            key = 'log.event.' + name
            if key in messages.catalog:
                error = payload.partition(': ')[2]
                rendered = messages(key, python=payload.partition('python=')[2], error=translate(error))
            else:
                rendered = translate(payload)
            lines.append(prefix + rendered + newline)
        else:
            lines.append(translate(line.rstrip('\n')) + ('\n' if line.endswith('\n') else ''))
    return ''.join(lines)
