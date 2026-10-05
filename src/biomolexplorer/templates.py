"""Validate per-stage scientific configuration overlays, never arbitrary code."""
import json
import re
import shutil
from pathlib import Path
from string import Formatter

RESOURCE_ROOT = Path(__file__).parent / 'resources'
CHIMERA_COMMANDS = {'open','delete','select','write','close','addh','addcharge','minimize'}


def validate_templates(templates):
    if not isinstance(templates, dict):
        raise ValueError('As configurações de scripts devem ser um objeto JSON.')
    for name, text in templates.items():
        reference = (RESOURCE_ROOT / name).resolve()
        if not reference.is_relative_to(RESOURCE_ROOT.resolve()) or not reference.is_file():
            raise ValueError('Template desconhecido.')
        if not isinstance(text,str) or len(text) > 100000:
            raise ValueError('Template inválido ou muito grande.')
        if name.endswith('.json'):
            if not isinstance(json.loads(text),dict):
                raise ValueError('Filtros devem ser objetos JSON.')
            continue
        original = reference.read_text()
        expected = {field for _,field,_,_ in Formatter().parse(original) if field}
        actual = {field for _,field,_,_ in Formatter().parse(text) if field}
        if actual != expected or any(not re.fullmatch('[a-zA-Z_][a-zA-Z_0-9]*', field) for field in actual):
            raise ValueError('Mantenha os marcadores de entrada e saída do template: ' + ', '.join(sorted(expected)))
        if any(c in text for c in (';', '`', '$', '\x00', '\\')) or '..' in text or '~' in text:
            raise ValueError('Comandos externos e caminhos fora do projeto não são permitidos nos templates.')
        # Static absolute paths are disallowed; placeholders are resolved to scoped inputs.
        if re.search(r'(?<![\w{])/(?:\w|\.)', text):
            raise ValueError('Use os marcadores do template para os caminhos de arquivos.')
        if name.startswith('chimera/'):
            if re.search(r'\b(?:https?|ftp|file):|(?:^|\s)(?:python|system|run)\b',text,re.I):
                raise ValueError('Use apenas os arquivos e comandos científicos do projeto.')
            for line in text.splitlines():
                if line.strip() and not line.lstrip().startswith('#') and line.split()[0].lower() not in CHIMERA_COMMANDS:
                    raise ValueError('O editor aceita os comandos científicos do template, não scripts executáveis arbitrários.')


def materialize_templates(destination, templates):
    validate_templates(templates)
    destination = Path(destination)
    shutil.copytree(RESOURCE_ROOT, destination)
    for name,text in templates.items():
        (destination / name).write_text(text)
    return destination
