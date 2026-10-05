"""Three-way merging preserves independent edits by project collaborators."""
from copy import deepcopy

_MISSING = object()


def merge(base, local, remote, label):
    if local == remote or remote == base:
        return deepcopy(local) if local is not _MISSING else _MISSING
    if local == base:
        return deepcopy(remote) if remote is not _MISSING else _MISSING
    if all(isinstance(value, dict) for value in (base, local, remote)):
        result = {}
        for key in base.keys() | local.keys() | remote.keys():
            value = merge(base.get(key, _MISSING), local.get(key, _MISSING), remote.get(key, _MISSING), label + '/' + key)
            if value is not _MISSING:
                result[key] = value
        return result
    raise ValueError(f'Outro colaborador alterou o mesmo campo ({label}). Sua edição foi mantida; reabra o projeto e confira a versão compartilhada antes de reaplicar a alteração.')


def merge_pipeline(base, local, remote):
    indexes = [{s['id']: s for s in stages} for stages in (base, local, remote)]
    result = merge(*indexes, 'bloco')
    order = list(dict.fromkeys([s['id'] for s in remote] + [s['id'] for s in local]))
    return [result[i] for i in order if i in result]
