"""Lazy operation registry: add providers/analyses without coupling to Flet."""
from importlib import import_module
from pathlib import Path
import shutil
import math

from .contracts import OperationResult, OperationSpec
from .paths import resolve_path

OPERATIONS = {
    'retrieve_compounds': OperationSpec('wrappers.crawlers', 'retrieve_compounds',
        ('search_term',), ('include_pubchem', 'pubchem_threshold', 'pubchem_max_records', 'chembl_filters')),
    'expand_similar_compounds': OperationSpec('wrappers.crawlers', 'expand_similar_compounds',
        ('search_term',), ('threshold', 'max_records', 'base_input_path', 'input_file')),
    'retrieve_structures': OperationSpec('wrappers.crawlers', 'load_pdb', ('target',),
        ('pdb_ec', 'organism', 'PolymerEntityTypeID', 'ExperimentalMethodID', 'max_resolution', 'must_have_ligand'),
        {'PolymerEntityTypeID': 'crawlers.complex:PolymerEntityType',
         'ExperimentalMethodID': 'crawlers.complex:ExperimentalMethod'}),
    'retrieve_zinc': OperationSpec('wrappers.crawlers', 'load_zinc', ('filename', 'base_input_path'), ('verbose',)),
    'prepare_structures': OperationSpec('wrappers.redocking', 'prepare_structures', ('base_input_path', 'target'),
        ('pdb_codes', 'pH', 'charge_type')),
    'admet': OperationSpec('wrappers.admet', 'ADMETWrapper', ('base_input_path',), ('input_file', 'verbose')),
    'fingerprints': OperationSpec('wrappers.molecular_analyzer', 'generate_fingerprints',
        ('base_input_path',), ('morgan_n_bits', 'radius', 'files', 'morgan', 'maccs', 'pharmacophore', 'chunk_size')),
    'similarity': OperationSpec('wrappers.molecular_analyzer', 'compute_similarity',
        ('base_input_path',), ('metric', 'fingerprint', 'filename', 'threshold', 'approximate'),
        {'metric': 'kernel.descriptors:similarityFunctions', 'fingerprint': 'kernel.descriptors:fingerprints'}),
    'graphs': OperationSpec('wrappers.molecular_analyzer', 'analyze_graphs',
        (), ('similarity_path','base_input_path',
             'mcs_timeout','mcs_ring_matches_ring_only','mcs_complete_rings_only','graph_inputs')),
    'redocking': OperationSpec('wrappers.redocking', 'perform_redocking', ('base_input_path', 'target'),
        ('pdb_codes', 'pH', 'sizeof_box', 'exhaustiveness', 'num_modes', 'prepare_complex', 'charge_type')),
    'docking_vina': OperationSpec('wrappers.docking', 'perform_docking_vina',
        ('base_input_path', 'target', 'base_selected_mols', 'mol_filename'),
        ('pdb_code', 'pH', 'sizeof_box', 'exhaustiveness', 'num_modes')),
    'docking_dock6': OperationSpec('wrappers.docking', 'perform_docking_dock6',
        ('base_input_path', 'target', 'base_selected_mols', 'dock6_app_path', 'charge_type', 'mol_filename', 'pdb_code', 'base_vina_path'),
        ('density', 'radius', 'distance', 'conformer_search_type', 'plot_max_residues')),
    'consensus': OperationSpec('wrappers.docking', 'generate_consensus', ('base_input_path', 'target'),
        ('repulsion_weight', 'base_vina_path', 'base_dock6_path')),
}


def validate_operation(operation, parameters):
    if operation not in OPERATIONS:
        raise ValueError(f'Unknown operation: {operation}')
    if not isinstance(parameters, dict):
        raise ValueError('Parameters must be a dictionary')
    spec = OPERATIONS[operation]
    unknown = parameters.keys() - set(spec.required + spec.optional)
    missing = set(spec.required) - parameters.keys()
    if unknown or missing:
        raise ValueError(f'Invalid parameters: unknown={sorted(unknown)}, missing={sorted(missing)}')
    for name in spec.required:
        if parameters[name] is None or parameters[name] == '':
            raise ValueError(f'{name} cannot be empty')
    for name in ('search_term', 'target'):
        if name in parameters:
            value = parameters[name]
            if not isinstance(value, str) or not value.strip() or any(c in value for c in '/\\\x00') or value.strip() in ('.', '..'):
                raise ValueError(f'{name} must be a non-empty name without path separators')
    for name in ('pubchem_threshold', 'threshold'):
        if name in parameters and parameters[name] is not None:
            value = parameters[name]
            minimum = 0 if operation in ('similarity','graphs') else 1
            if type(value) is not int or not minimum <= value <= 100:
                raise ValueError(f'{name} must be an integer between {minimum} and 100')
    if operation == 'fingerprints' and 'radius' in parameters and (type(parameters['radius']) is not int or parameters['radius'] < 0):
        raise ValueError('radius must be a non-negative integer')
    for name in ('density', 'distance', 'max_resolution') + (('radius',) if operation == 'docking_dock6' else ()):
        if name in parameters:
            value = parameters[name]
            if name == 'max_resolution' and value is None:
                continue
            if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
                raise ValueError(f'{name} must be a positive finite number')
    if 'chembl_filters' in parameters:
        filters = parameters['chembl_filters']
        if not isinstance(filters, dict) or filters.keys() - {'target', 'bioactivity', 'molecules', 'similars'}:
            raise ValueError('chembl_filters must contain only known ChEMBL stages')
        if not all(isinstance(value, dict) for value in filters.values()):
            raise ValueError('Each ChEMBL stage filter must be a dictionary')
    for name in ('max_records', 'pubchem_max_records', 'morgan_n_bits', 'num_modes', 'exhaustiveness', 'chunk_size','mcs_timeout'):
        if name in parameters and (type(parameters[name]) is not int or parameters[name] < 1):
            raise ValueError(f'{name} must be a positive integer')
    if 'pH' in parameters:
        value = parameters['pH']
        if type(value) not in (int,float) or not math.isfinite(value) or not 0 <= value <= 14:
            raise ValueError('pH must be a finite number between 0 and 14')
    if 'repulsion_weight' in parameters:
        value=parameters['repulsion_weight']
        if type(value) not in (int,float) or not math.isfinite(value) or value<0:
            raise ValueError('repulsion_weight must be a non-negative finite number')
    if 'sizeof_box' in parameters:
        value = parameters['sizeof_box']
        if not isinstance(value,list) or len(value)!=3 or any(type(v) not in (int,float) or not math.isfinite(v) or v<=0 for v in value):
            raise ValueError('sizeof_box must contain three positive finite dimensions')
    for name in ('include_pubchem', 'verbose', 'must_have_ligand', 'prepare_complex', 'morgan', 'maccs', 'pharmacophore', 'approximate',
                 'mcs_ring_matches_ring_only','mcs_complete_rings_only'):
        if name in parameters and type(parameters[name]) is not bool:
            raise ValueError(f'{name} must be boolean')
    if operation=='graphs' and 'graph_inputs' in parameters:
        entries=parameters['graph_inputs']
        if not isinstance(entries,list) or not entries or len(entries)>1000:
            raise ValueError('Selecione de 1 a 1000 entradas de grafos.')
        for entry in entries:
            if (not isinstance(entry,dict) or entry.get('kind') != 'similarity'
                    or not isinstance(entry.get('file'),str) or not isinstance(entry.get('compound_files'),list)
                    or any(not isinstance(p,str) for p in entry['compound_files'])):
                raise ValueError('Entrada de grafo inválida. Use apenas similaridade pronta: CSV com source,target,value; tabelas de compostos são opcionais.')
            if any(k in entry and not isinstance(entry[k],str) for k in ('label','metric','fingerprint')):
                raise ValueError('Nome, métrica e tipo da entrada de grafo devem ser textos.')
    return spec


def execute_operation(operation, parameters, output_path):
    from .molecule_quality import CURRENT_REPORT, QualityReport, prepare_inputs, sanitize_csv, table_kind
    validate_operation(operation,parameters)
    output = Path(output_path).resolve();output.mkdir(parents=True,exist_ok=True)
    report = QualityReport(operation);token = CURRENT_REPORT.set(report)
    try:
        clean = prepare_inputs(operation,parameters,output,report)
        result = _execute_operation(operation,clean,output)
        for artifact in result.artifacts if operation!='graphs' else []:
            if Path(artifact).suffix=='.csv' and table_kind(artifact):sanitize_csv(artifact,artifact,report=report)
        path = report.write(output)
        details = dict(result.details,excluded_records=len(report.records))
        if path:details['exclusion_report']=str(path)
        return OperationResult(operation,result.artifacts+([str(path)] if path and str(path) not in result.artifacts else []),details)
    finally:
        # Also retain exclusions when a separate provider/engine error occurs.
        if report.records:report.write(output)
        CURRENT_REPORT.reset(token)


def _execute_operation(operation, parameters, output_path):
    spec = validate_operation(operation, parameters)
    kwargs = dict(parameters)
    output = Path(output_path).resolve()
    output.mkdir(parents=True, exist_ok=True)
    for name in ('base_input_path', 'base_selected_mols', 'dock6_app_path', 'base_vina_path', 'base_dock6_path', 'similarity_path','fingerprints_path'):
        if kwargs.get(name) is not None:
            kwargs[name] = str(resolve_path(kwargs[name]))
    for name, reference in spec.enums.items():
        if name not in kwargs:
            continue
        module, enum_name = reference.split(':')
        enum = getattr(import_module(module), enum_name)
        def convert(value):
            return enum[value] if value in enum.__members__ else enum(value)
        kwargs[name] = [convert(v) for v in kwargs[name]] if isinstance(kwargs[name], list) else convert(kwargs[name])
    if operation in ('redocking', 'prepare_structures'):
        # The legacy redocking routine filters/deletes structures in its input
        # directory. Work on a snapshot, preserving the user's retrieval data.
        target = kwargs['target'].replace(' ', '')
        snapshot = output / 'structures'
        shutil.copytree(Path(kwargs['base_input_path']) / target, snapshot / target)
        kwargs['base_input_path'] = str(snapshot)
    if operation == 'fingerprints':
        # The legacy wrapper puts outputs beside its input; explicitly isolate them.
        from wrappers.molecular_analyzer import generate_fingerprints
        result = generate_fingerprints(**kwargs, base_output_path=str(output))
    else:
        kwargs['base_output_path'] = str(output)
        result = getattr(import_module(spec.module), spec.function)(**kwargs)
        if operation == 'admet':
            result = result.run_pipeline()
    artifacts = [str(p) for p in sorted(output.rglob('*')) if p.is_file()
                 and 'cache' not in p.parts and '.quality-inputs' not in p.parts
                 and p.suffix in ('.csv', '.png', '.json', '.pdb', '.pdbqt', '.mol2')]
    details = {'rows': len(result)} if hasattr(result, 'columns') else {}
    return OperationResult(operation, artifacts, details)
