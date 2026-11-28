"""Validation contracts per artifact type: invalid content is never published.

`artifact_contracts` declares what each type *is*; this module decides whether a
file actually is one.  The contract runs at the publication boundary and again
before reuse, so content that fails it never enters the store and never comes
out of it as a cache hit.

Each contract is as deep as its type justifies and no deeper.  A score table is
checked by its declared columns, because reading fifteen thousand rows to
publish them twice would cost more than it protects.  A prepared engine input
is checked by the validators that already exist for it — native Vina aborts the
process instead of raising on a malformed PDBQT, which is why
``ligand_cache.validate_prepared_ligand`` exists and why it belongs here too.
"""
from pathlib import Path

# How much of a molecular file the structural contracts examine.  A pose file
# can hold thousands of models; the first one answers whether the file is a
# molecular artifact at all.
SCANNED_LINES=4096


class InvalidContent(ValueError):
    """Content that does not satisfy the contract of its artifact type."""


def _present(path):
    path=Path(path)
    if not path.is_file() or not path.stat().st_size:
        raise InvalidContent('O conteúdo está ausente ou vazio: '+path.name)
    return path


def _header(path):
    import csv
    with path.open(encoding='utf-8-sig',newline='') as stream:
        for row in csv.reader(stream):
            return [column.strip() for column in row]
    return []


def _table(path,artifact_type):
    """A tabular artifact carries the columns its type declares."""
    from .artifact_contracts import TYPES
    columns=set(_header(path))
    if not columns:
        raise InvalidContent('A tabela não contém cabeçalho: '+path.name)
    declared=getattr(TYPES.get(artifact_type),'schema',None) or ()
    missing=[column for column in declared if column not in columns]
    if missing:
        raise InvalidContent('A tabela não declara as colunas do seu tipo ('+', '.join(missing)+'): '+path.name)
    return path


def _coordinates(line):
    import math
    try:
        values=[float(line[start:end]) for start,end in ((30,38),(38,46),(46,54))]
    except ValueError:
        raise InvalidContent('Coordenadas atômicas inválidas.') from None
    if not all(math.isfinite(value) for value in values):
        raise InvalidContent('Coordenadas atômicas não finitas.')
    return values


def _structure(path):
    """A molecular artifact has at least one atom with finite coordinates."""
    suffix=path.suffix.lower()
    if suffix in ('.cif','.mmcif'):
        with path.open(encoding='utf-8',errors='replace') as stream:
            for index,line in enumerate(stream):
                if line.startswith('_atom_site.'):return path
                if index>=SCANNED_LINES:break
        raise InvalidContent('O arquivo mmCIF não declara posições atômicas: '+path.name)
    if suffix=='.mol2':
        with path.open(encoding='utf-8',errors='replace') as stream:
            for index,line in enumerate(stream):
                if line.startswith('@<TRIPOS>ATOM'):return path
                if index>=SCANNED_LINES:break
        raise InvalidContent('O arquivo MOL2 não contém a seção de átomos: '+path.name)
    with path.open(encoding='utf-8',errors='replace') as stream:
        for index,line in enumerate(stream):
            if line.startswith(('ATOM  ','HETATM')):
                _coordinates(line);return path
            if index>=SCANNED_LINES:break
    raise InvalidContent('O arquivo não contém registros atômicos: '+path.name)


def _pose(path):
    """The first model of a pose file is a molecular artifact."""
    return _structure(path)


def _prepared(path):
    """A prepared engine input goes through the validator the engine needs."""
    _structure(path)
    suffix=path.suffix.lower()
    if suffix in ('.pdbqt','.mol2'):
        from .ligand_cache import validate_prepared_ligand
        try:
            validate_prepared_ligand(path,suffix.lstrip('.'))
        except InvalidContent:
            raise
        except ValueError as error:
            raise InvalidContent(str(error)) from None
    return path


def _download_list(path):
    from .zinc_retrieval import read_download_list
    try:
        read_download_list(path)
    except ValueError as error:
        raise InvalidContent(str(error)) from None
    return path


CONTRACTS={'present':_present,'table':_table,'structure':_structure,
           'pose':_pose,'prepared':_prepared,'download_list':_download_list}


def validate(artifact_type,path):
    """Apply the contract of ``artifact_type`` to ``path``; raise to refuse it."""
    from .artifact_contracts import CONTENT_CONTRACTS
    path=_present(path)
    name=CONTENT_CONTRACTS.get(artifact_type,'present')
    contract=CONTRACTS.get(name)
    if contract is None:
        raise InvalidContent('Contrato de conteúdo desconhecido: '+str(name))
    if contract is _table:
        return _table(path,artifact_type)
    return contract(path)


def validator(artifact_type):
    """The contract as a one-argument callable, for the content store."""
    return lambda path:validate(artifact_type,path)
