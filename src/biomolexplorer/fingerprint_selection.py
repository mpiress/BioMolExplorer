"""Infer fingerprint algorithms from generated inputs, including older projects."""
from pathlib import Path
from .bindings import sources

KINDS = ('morgan', 'maccs', 'pharmacophore')
LABELS = {'morgan': 'Morgan', 'maccs': 'MACCS', 'pharmacophore': 'Farmacóforo'}


def filename_kind(filename):
    name = Path(filename).name.lower()
    return next((kind for kind in KINDS if name.startswith(kind + '_')), None)


def generated_kind(stage, pipeline, artifacts=None):
    """Return (algorithm, has_custom_inputs). Reject incompatible generated sources."""
    lookup = {s['id']: s for s in pipeline}
    kinds = set()
    binding=stage.get('bindings', {}).get('base_input_path')
    refs = sources(binding) if binding else []
    custom = bool(stage.get('parameters', {}).get('base_input_path'))
    for ref in refs:
        producer = lookup.get(ref.get('stage'))
        if not producer or producer['operation'] != 'fingerprints' or producer.get('provided_results'):
            custom = True
            continue
        selector = ref.get('selector', 'auto')
        kind = filename_kind(selector) if selector != 'auto' else None
        if kind is None and artifacts is not None:
            files = sorted(p for p in artifacts.get(producer['id'], []) if filename_kind(p))
            if files:
                kind = filename_kind(files[0])
        if kind is None:
            enabled = [k for k in KINDS if producer['parameters'].get(k, k == 'morgan')]
            if len(enabled) == 1:
                kind = enabled[0]
            elif enabled:
                raise ValueError('Esta origem contém mais de um tipo de fingerprint. Selecione um arquivo específico na entrada.')
        if kind:
            kinds.add(kind)
    if len(kinds) > 1:
        raise ValueError('Use entradas com o mesmo tipo de fingerprint. Morgan, MACCS e farmacóforo não podem ser combinados.')
    return next(iter(kinds), None), custom
