"""Typed scientific ports shared by the editor, validator and executor.

The contract describes scientific values.  File suffixes remain an adapter
concern and are intentionally absent from the public port definitions.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class ArtifactType:
    name: str
    schema_version: int = 1
    entity_kind: str | None = None
    required_metadata: tuple[str, ...] = ()
    schema: tuple[str, ...] = ()


@dataclass(frozen=True)
class PortContract:
    name: str
    types: tuple[str, ...]
    cardinality: str = 'many'
    required: bool | str = False
    schema_version: int = 1
    required_metadata: tuple[str, ...] = ()


TYPES = {
    'compounds': ArtifactType('compounds', entity_kind='compound', required_metadata=('identifier',),schema=('molecule_chembl_id','canonical_smiles')),
    'chembl': ArtifactType('chembl', entity_kind='compound', required_metadata=('identifier',),schema=('molecule_chembl_id','canonical_smiles')),
    'fingerprints': ArtifactType('fingerprints', entity_kind='compound', required_metadata=('algorithm',),schema=('molecule_chembl_id','fingerprint')),
    'similarity': ArtifactType('similarity', entity_kind='molecular_pair', required_metadata=('metric',),schema=('source','target','value')),
    'structures': ArtifactType('structures', entity_kind='structure', required_metadata=('pdb_id',)),
    'prepared_structures': ArtifactType('prepared_structures', entity_kind='receptor', required_metadata=('preparation',)),
    'receptor_pdb': ArtifactType('receptor_pdb', entity_kind='receptor'),
    'molecular_poses': ArtifactType('molecular_poses', entity_kind='pose', required_metadata=('receptor_id','ligand_id')),
    'vina': ArtifactType('vina', entity_kind='pose', required_metadata=('engine',)),
    'dock6': ArtifactType('dock6', entity_kind='pose', required_metadata=('engine',)),
    'scores': ArtifactType('scores', entity_kind='docking_pair'),
    'dynamics': ArtifactType('dynamics', entity_kind='trajectory'),
    'zinc_urls': ArtifactType('zinc_urls', entity_kind='compound'),
    'other': ArtifactType('other'),
}


OUTPUT_PORTS = {
    'retrieve_pubchem': {'compounds': ('compounds',)},
    'retrieve_compounds': {'compounds': ('compounds','chembl')},
    'expand_similar_compounds': {'compounds': ('compounds',)},
    'retrieve_zinc': {'compounds': ('compounds',)},
    'retrieve_structures': {'receptor_structures': ('structures','receptor_pdb')},
    'prepare_structures': {'prepared_receptors': ('prepared_structures','receptor_pdb'), 'compounds': ('compounds',), 'poses': ('molecular_poses',)},
    'redocking': {'structures': ('structures','prepared_structures','receptor_pdb'), 'redocking_results': ('molecular_poses',)},
    'admet': {'filtered_compounds': ('compounds',)},
    'fingerprints': {'fingerprints': ('fingerprints',), 'compounds': ('compounds',)},
    'similarity': {'similarity': ('similarity',)},
    'graphs': {'selected_compounds': ('compounds',)},
    'docking_vina': {'docking_results': ('vina','molecular_poses'), 'compounds': ('compounds',), 'receptors': ('receptor_pdb',)},
    'docking_dock6': {'docking_results': ('dock6','molecular_poses'), 'compounds': ('compounds',), 'receptors': ('receptor_pdb',)},
    'consensus': {'consensus_scores': ('scores','compounds')},
    'dynamics_dock6': {'dynamics': ('dynamics',), 'poses': ('molecular_poses',), 'receptors': ('receptor_pdb',)},
}


INPUT_PORTS = {
    'dynamics_dock6': {'base_input_path': ('receptor_pdb',), 'base_selected_mols': ('compounds','molecular_poses')},
    'retrieve_pubchem': {'base_input_path': ('compounds',)},
    'expand_similar_compounds': {'base_input_path': ('chembl',)},
    'retrieve_zinc': {'base_input_path': ('zinc_urls',)},
    'admet': {'base_input_path': ('compounds',)},
    'fingerprints': {'base_input_path': ('compounds',)},
    'similarity': {'base_input_path': ('fingerprints',)},
    'graphs': {'similarity_path': ('similarity',)},
    'prepare_structures': {'base_input_path': ('structures','prepared_structures'), 'base_selected_mols': ('compounds',)},
    'redocking': {'base_input_path': ('structures','prepared_structures',)},
    'docking_vina': {'base_input_path': ('structures','prepared_structures'), 'base_selected_mols': ('compounds','vina','dock6')},
    'docking_dock6': {'base_input_path': ('structures','prepared_structures'), 'base_selected_mols': ('compounds','vina','dock6')},
    'consensus': {'base_vina_path': ('vina',), 'base_dock6_path': ('dock6',)},
}

INPUT_RULES = {
    'dynamics_dock6': {'base_input_path':('one_or_many',True),'base_selected_mols':('one_or_many',True)},
    'retrieve_pubchem': {'base_input_path':('one_or_many','reference_source=compounds')},
    'expand_similar_compounds': {'base_input_path':('one_or_many',True)},
    'retrieve_zinc': {'base_input_path':('one_or_many',True)},
    'admet': {'base_input_path':('one_or_many',True)},
    'fingerprints': {'base_input_path':('one_or_many',True)},
    'similarity': {'base_input_path':('one_or_many',True)},
    'graphs': {'similarity_path':('one_or_many',True)},
    'prepare_structures': {'base_input_path':('one_or_many',True),'base_selected_mols':('one_or_many',True)},
    'redocking': {'base_input_path':('one_or_many',True)},
    'docking_vina': {'base_input_path':('one_or_many',True),'base_selected_mols':('one_or_many',True)},
    'docking_dock6': {'base_input_path':('one_or_many',True),'base_selected_mols':('one_or_many',True)},
    'consensus': {'base_vina_path':('one_or_many',True),'base_dock6_path':('one_or_many',True)},
}

OUTPUT_CARDINALITY = {
    # A logical port may publish one occurrence per input batch.  The
    # consumer's processing mode decides whether those occurrences remain
    # individual or are merged; the scientific contract must accept both.
    operation:{port:'one_or_many' for port in ports}
    for operation,ports in OUTPUT_PORTS.items()
}

ROLE_PORTS = {
    'retrieve_pubchem': {'compounds':'compounds','mol':'compounds','sim':'compounds','merged':'compounds'},
    'retrieve_compounds': {'compounds':'compounds','mol':'compounds','sim':'compounds','merged':'compounds'},
    'expand_similar_compounds': {'compounds':'compounds','mol':'compounds','sim':'compounds','merged':'compounds'},
    'retrieve_zinc': {'compounds':'compounds','mol':'compounds'},
    'retrieve_structures': {'structure':'receptor_structures'},
    'prepare_structures': {'receptor':'prepared_receptors','complexes':'prepared_receptors','centers':'prepared_receptors',
                           'ligand':'poses','compounds':'compounds'},
    'redocking': {'structure':'structures','complexes':'structures','receptor':'structures','pose':'redocking_results',
                  'scores':'redocking_results','rmsd':'redocking_results','ligand':'redocking_results'},
    'admet': {'compounds':'filtered_compounds','selected':'filtered_compounds','rejected':'filtered_compounds',
              'evaluated':'filtered_compounds','bbb-pos':'filtered_compounds','bbb-neg':'filtered_compounds',
              'hia-pos':'filtered_compounds','hia-neg':'filtered_compounds'},
    'fingerprints': {'fingerprints':'fingerprints','morgan':'fingerprints','maccs':'fingerprints',
                     'pharmacophore':'fingerprints','compounds':'compounds'},
    'graphs': {'lcc':'selected_compounds','component':'selected_compounds','compounds':'selected_compounds'},
    'docking_vina': {'scores':'docking_results','pose':'docking_results','compounds':'compounds','receptor':'receptors'},
    'docking_dock6': {'scores':'docking_results','pose':'docking_results','compounds':'compounds','receptor':'receptors'},
    'consensus': {'ranking':'consensus_scores','scores':'consensus_scores'},
    'dynamics_dock6': {'trajectory':'dynamics','topology':'dynamics','coordinates':'dynamics','rmsd':'dynamics',
                       'rmsf':'dynamics','gyration':'dynamics','energies':'dynamics','protocol':'dynamics',
                       'pose':'poses','ligand':'poses','receptor':'receptors'},
}


def output_types(operation):
    return {kind for kinds in OUTPUT_PORTS.get(operation,{}).values() for kind in kinds}


def output_port(operation, accepted_types):
    """Return the single compatible semantic port, rejecting ambiguity."""
    matches=[name for name,kinds in OUTPUT_PORTS.get(operation,{}).items() if set(kinds)&set(accepted_types)]
    return matches[0] if len(matches)==1 else None


def port_types(operation, port):
    return set(OUTPUT_PORTS.get(operation,{}).get(port,()))


def infer_output_port(operation, role, path='',legacy=False):
    ports=OUTPUT_PORTS.get(operation,{})
    if len(ports)==1:return next(iter(ports),None)
    role=str(role or '').lower();name=str(path).lower()
    direct=ROLE_PORTS.get(operation,{}).get(role)
    if direct:return direct
    for prefix,port in ROLE_PORTS.get(operation,{}).items():
        if role.startswith(prefix+'-'):return port
    if not legacy:return None
    hints={
        'fingerprints': ('fingerprint',), 'similarity': ('similarity','edges'),
        'prepared_receptors': ('receptor','structure','complex','center'),
        'structures': ('receptor','structure','complex','center'),
        'redocking_results': ('redocking','pose','score'),
        'docking_results': ('vina','dock6','pose','score','ranking'),
        'receptors': ('receptor','structure'), 'poses': ('pose','ligand'),
        'dynamics': ('trajectory','topology','rmsd','rmsf','gyration','energy','protocol'),
        'consensus_scores': ('score','ranking','consensus'),
        'filtered_compounds': ('compound','selected','rejected','evaluated','bbb','hia','rule'),
        'selected_compounds': ('compound','component','lcc','mol'),
        'compounds': ('compound','mol','sim','merged'),
    }
    value=role+' '+name
    matches=[port for port in ports if any(token in value for token in hints.get(port,()))]
    return matches[0] if len(matches)==1 else None


def operation_contract(operation):
    def input_contract(name,kinds):
        cardinality,required=INPUT_RULES.get(operation,{}).get(name,('one_or_many',False))
        metadata=tuple(sorted({value for kind in kinds for value in TYPES[kind].required_metadata}))
        return PortContract(name,tuple(kinds),cardinality,required,required_metadata=metadata).__dict__
    def output_contract(name,kinds):
        metadata=tuple(sorted({value for kind in kinds for value in TYPES[kind].required_metadata}))
        return PortContract(name,tuple(kinds),OUTPUT_CARDINALITY[operation][name],False,required_metadata=metadata).__dict__
    return {
        'operation': operation,
        'schema_version': 1,
        'inputs': {name: input_contract(name,kinds) for name,kinds in INPUT_PORTS.get(operation,{}).items()},
        'outputs': {name: output_contract(name,kinds) for name,kinds in OUTPUT_PORTS.get(operation,{}).items()},
    }

