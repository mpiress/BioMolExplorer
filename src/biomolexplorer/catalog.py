"""UI metadata derived from operation signatures without importing scientific code."""
import ast
import copy
from functools import lru_cache
from pathlib import Path
from .operations import OPERATIONS
from .paths import SOURCE_ROOT

TITLES = {
 'import_results': ('Importar meus arquivos', 'Use compostos, PDBs ou resultados preparados por você.', 'Entradas'),
 'retrieve_compounds': ('Recuperar compostos', 'Busque moléculas na ChEMBL e amplie a seleção com PubChem.', 'Recuperação'),
 'expand_similar_compounds': ('Expandir similares', 'Encontre novos compostos a partir de downloads ChEMBL existentes.', 'Recuperação'),
 'retrieve_structures': ('Recuperar estruturas PDB', 'Selecione estruturas de proteínas e seus ligantes.', 'Recuperação'),
 'retrieve_zinc': ('Recuperar ZINC', 'Baixe uma seleção usando um arquivo de endereços ZINC.', 'Recuperação'),
 'prepare_structures': ('Preparar meus complexos', 'Prepare PDBs próprios com Chimera e Open Babel, sem recuperar PDBs.', 'Docking'),
 'admet': ('Avaliar ADMET', 'Calcule propriedades e aplique os filtros moleculares disponíveis.', 'Análise'),
 'fingerprints': ('Gerar fingerprints', 'Descreva as moléculas para comparações de similaridade.', 'Análise'),
 'similarity': ('Calcular similaridade', 'Construa relações entre moléculas usando seus fingerprints.', 'Análise'),
 'graphs': ('Filtrar por grafos', 'Explore relações e selecione o maior componente conectado.', 'Análise'),
 'redocking': ('Executar redocking', 'Avalie o protocolo de docking com ligantes de referência.', 'Docking'),
 'docking_vina': ('Docking com Vina', 'Execute a seleção de compostos contra receptores preparados.', 'Docking'),
 'docking_dock6': ('Docking com DOCK6', 'Refine conformações e obtenha scores e footprints.', 'Docking'),
 'consensus': ('Consenso de docking', 'Combine resultados Vina e DOCK6.', 'Docking'),
}
LABELS = {
 'search_term': 'Alvo molecular (nome ou ID ChEMBL)', 'target': 'Nome do alvo / pasta das estruturas',
 'include_pubchem': 'Ampliar a seleção com PubChem', 'pubchem_threshold': 'Similaridade mínima PubChem (%)',
 'pubchem_max_records': 'Máximo de similares por referência', 'threshold': 'Limiar de similaridade (%)',
 'max_records': 'Máximo de registros', 'base_input_path': 'Dados de entrada',
 'base_selected_mols': 'Compostos selecionados', 'base_vina_path': 'Resultados Vina',
 'base_dock6_path': 'Resultados DOCK6',
 'similarity_path': 'Relações de similaridade', 'dock6_app_path': 'Instalação DOCK6',
 'mcs_timeout': 'Tempo máximo da busca do fragmento (s)',
 'mcs_ring_matches_ring_only': 'Comparar átomos de anéis somente com anéis',
 'mcs_complete_rings_only': 'Exigir anéis completos no fragmento',
 'pdb_ec': 'Número EC', 'organism': 'Organismos', 'max_resolution': 'Resolução máxima (Å)',
 'must_have_ligand': 'Exigir ligante', 'input_file': 'Arquivo de compostos (opcional)',
 'mol_filename': 'Nome da tabela de compostos (sem .csv)', 'filename': 'Nome do arquivo',
 'metric': 'Métrica de similaridade', 'fingerprint': 'Tipo de fingerprint', 'radius': 'Raio',
 'morgan_n_bits': 'Número de bits Morgan', 'approximate': 'Busca aproximada (mais rápida)',
 'pH': 'pH de preparação', 'sizeof_box': 'Dimensões da caixa [x, y, z]', 'exhaustiveness': 'Esforço de busca',
 'num_modes': 'Número de poses', 'prepare_complex': 'Preparar complexos antes do redocking',
 'charge_type': 'Método de atribuição de cargas', 'pdb_codes': 'Complexos [PDB, ligante, resíduo, cadeia, resolução]',
 'pdb_code': 'Complexo de referência', 'conformer_search_type': 'Busca de conformações',
 'chunk_size': 'Moléculas por bloco', 'chembl_filters': 'Filtros ChEMBL por etapa',
 'repulsion_weight': 'Peso de repulsão no consenso', 'density': 'Densidade da superfície',
 'distance': 'Distância de seleção', 'plot_max_residues': 'Máximo de resíduos na figura',
 'PolymerEntityTypeID': 'Tipos de polímero', 'ExperimentalMethodID': 'Métodos experimentais',
 'files': 'Arquivos de compostos', 'morgan': 'Morgan', 'maccs': 'MACCS', 'pharmacophore': 'Farmacóforo',
 'verbose': 'Detalhar o processamento',
}
PATH_FIELDS = {'base_input_path','base_selected_mols','base_vina_path','base_dock6_path','similarity_path'}
ENUMS = {
 'metric': ['Tanimoto','Dice','Cosine','Sokal','Russel','RogotGoldberg','AllBit','Kulczynski','McConnaughey','Asymmetric','BraunBlanquet'],
 'fingerprint': ['morgan','maccs','pharmacophore'], 'charge_type': ['gas','am1'], 'conformer_search_type': ['flex','rigid'],
}
ENUM_DEFAULTS = {'metric': 'Tanimoto', 'fingerprint': 'morgan'}
PRESETS = {
 'Compostos → ADMET': ['retrieve_compounds', 'admet'],
 'Seleção por grafos → ADMET': ['retrieve_compounds', 'fingerprints', 'similarity', 'graphs', 'admet'],
 'Meus compostos → ADMET': ['import_results', 'admet'],
 'Meus PDBs → preparação': ['import_results', 'prepare_structures'],
 'Redocking → Vina': ['retrieve_structures', 'redocking', 'retrieve_compounds', 'docking_vina'],
 'Exploração completa → consenso': ['retrieve_structures', 'redocking', 'retrieve_compounds', 'fingerprints', 'similarity', 'graphs', 'admet', 'docking_vina', 'docking_dock6', 'consensus'],
}


@lru_cache(maxsize=64)
def operation_fields(operation):
    spec = OPERATIONS[operation]
    tree = ast.parse((SOURCE_ROOT / Path(*spec.module.split('.')).with_suffix('.py')).read_text())
    node = next(n for n in tree.body if isinstance(n, (ast.FunctionDef,ast.ClassDef)) and n.name == spec.function)
    if isinstance(node, ast.ClassDef):
        node = next(n for n in node.body if isinstance(n, ast.FunctionDef) and n.name == '__init__')
    defaults = dict(zip([a.arg for a in node.args.args][-len(node.args.defaults):], node.args.defaults)) if node.args.defaults else {}
    fields = []
    for key in spec.required + spec.optional:
        default = ENUM_DEFAULTS.get(key)
        if key in defaults:
            try:
                default = ast.literal_eval(defaults[key])
            except ValueError:
                default = ENUM_DEFAULTS.get(key)
        if key == 'include_pubchem':
            default = True
        fields.append({'name': key, 'label': LABELS.get(key,key), 'required': key in spec.required,
                       'default': default, 'path': key in PATH_FIELDS, 'choices': ENUMS.get(key)})
    return fields


def template_names(operation):
    names = []
    if operation == 'retrieve_compounds':
        names += [f'crawlers/{name}.json' for name in ('target','bioactivity','molecules','similarmols')]
    if operation in ('prepare_structures','redocking','docking_vina','docking_dock6'):
        names += [f'chimera/{p.name}' for p in sorted((SOURCE_ROOT / 'biomolexplorer/resources/chimera').glob('*.template'))]
    if operation in ('redocking','docking_vina'):
        names += ['vina/config.template']
    if operation == 'docking_dock6':
        names += [f'dock6/{p.name}' for p in sorted((SOURCE_ROOT / 'biomolexplorer/resources/dock6').glob('*.template'))]
    return names


def new_stage(operation):
    from uuid import uuid4
    if operation == 'import_results':
        parameters = {'kind': 'compounds', 'target': 'MeuAlvo', 'asset_ids': []}
    else:
        parameters = {field['name']: copy.deepcopy(field['default']) for field in operation_fields(operation)
                      if field['default'] is not None and not field['path']}
        if operation == 'fingerprints':
            parameters.update(morgan=True, maccs=False, pharmacophore=False)
        if 'search_term' in OPERATIONS[operation].required:
            parameters['search_term'] = 'CHEMBL220'
        if 'target' in OPERATIONS[operation].required:
            parameters['target'] = 'MeuAlvo'
        if 'mol_filename' in OPERATIONS[operation].required:
            parameters['mol_filename'] = 'molecules'
    stage = {'id': uuid4().hex, 'operation': operation, 'name': TITLES[operation][0],
            'enabled': True, 'parameters': parameters, 'bindings': {}, 'depends_on': [], 'templates': {}}
    if operation in ('retrieve_compounds', 'retrieve_structures', 'retrieve_zinc'):
        stage['process_all'] = False
    return stage
