"""Shared retrieval criteria and collection naming, independent of the UI."""
import hashlib
import re

MODES = {
    'target': 'Alvo · identificação automática',
    'target_name': 'Alvo · nome parcial',
    'target_id': 'Alvo · IDs ChEMBL',
    'uniprot': 'Alvo · acesso UniProt',
    'target_text': 'Alvo · texto livre',
    'molecule_id': 'Composto · IDs ChEMBL',
    'molecule_name': 'Composto · nome parcial',
    'similarity': 'Composto · similaridade por SMILES ou ID',
    'substructure': 'Composto · subestrutura SMILES',
}
UNIPROT = re.compile(r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\Z')


def terms(value):
    if value is None:
        return []
    if isinstance(value, str):
        value = re.split(r'[,;\s]+', value.strip())
    if not isinstance(value, (list, tuple)) or any(not isinstance(v, str) for v in value):
        raise ValueError('Informe identificadores como lista ou texto separado por vírgulas.')
    return list(dict.fromkeys(v.strip().upper() for v in value if v.strip()))


def identifiers(value, kind):
    values = terms(value)
    pattern = re.compile(r'[1-9][A-Z0-9]{3}\Z') if kind == 'pdb' else UNIPROT if kind == 'uniprot' else re.compile(r'CHEMBL[0-9]+\Z') if kind == 'chembl' else re.compile(r'[A-Z0-9]{1,5}\Z')
    if not values or any(not pattern.fullmatch(v) for v in values):
        raise ValueError(f'Identificador {kind} inválido. Confira os códigos informados.')
    return values


def collection_name(value):
    """Keep existing simple collection paths; hash structural query syntax."""
    value = value.strip().replace(' ', '')
    if re.fullmatch(r'[\w.-]{1,100}', value) and value not in ('.', '..'):
        return value
    return 'consulta_'+hashlib.sha256(value.encode()).hexdigest()[:16]


def target_criteria(term, mode='target'):
    term = term.strip()
    if mode == 'target':
        if re.fullmatch(r'CHEMBL[0-9]+', term.upper()):
            mode = 'target_id'
        elif UNIPROT.fullmatch(term.upper()):
            mode = 'uniprot'
        else:
            mode = 'target_name'
    if mode == 'target_id':
        return {'target_chembl_id__in': identifiers(term, 'chembl')}
    if mode == 'uniprot':
        return {'target_components__accession__in': identifiers(term, 'uniprot')}
    if mode == 'target_name':
        return {'pref_name__icontains': term}
    if mode == 'target_text':
        return {'_search': term}
    raise ValueError('Modo de busca de alvo inválido.')
