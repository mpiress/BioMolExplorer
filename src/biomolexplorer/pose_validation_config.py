"""Parameter contract of the pose validation block, shared by UI, CLI and runner."""
import math
import re

# Tolerance defaults mirror the wrapper signature; the UI reads these labels.
NUMERIC = {'vdw_clash_scale':(0.3,1.0),'delta_energy_limit':(0.0,100000.0),
           'bond_tolerance':(0.0,2.0),'angle_tolerance':(0.0,180.0),'planarity_tolerance':(0.0,5.0)}
INTEGER = {'minimization_steps':(1,20000)}
BOOLEANS = ('strict_filter',)
IDENTIFIER = re.compile(r'^[A-Za-z0-9_.+:-]{1,200}$')
# One ligand-receptor complex removed from the selection: compound:receptor.
PAIR = re.compile(r'^[A-Za-z0-9_.+-]{1,100}:[A-Za-z0-9_.+-]{1,100}$')
MAX_SELECTED = 5000
MAX_POSE_INPUTS = 20000
LABELS = {
    'strict_filter':'Entregar apenas poses válidas aos próximos blocos',
    'vdw_clash_scale':'Fração mínima da soma dos raios de van der Waals',
    'delta_energy_limit':'Energia interna máxima acima do confôrmero relaxado (kcal/mol)',
    'bond_tolerance':'Desvio relativo máximo nos comprimentos de ligação',
    'angle_tolerance':'Desvio máximo nos ângulos de ligação (°)',
    'planarity_tolerance':'Desvio máximo de planaridade (Å)',
    'minimization_steps':'Iterações máximas da minimização (MMFF94 ou UFF)',
}
HELP = {
    'strict_filter':'Ativado, somente as poses aprovadas seguem para Consenso ou Amber Score. '
                    'Desativado, todas as poses continuam disponíveis com os metadados da validação.',
    'vdw_clash_scale':'Choque estérico: um contato entre átomos pesados do ligante e do receptor ou cofatores '
                      'é choque quando a distância fica abaixo desta fração da soma dos raios de van der Waals. '
                      '0,75 preserva ligações de hidrogênio e contatos polares próximos.',
    'delta_energy_limit':'Diferença entre a energia da pose ancorada e a do mesmo ligante minimizado.',
    'bond_tolerance':'Comparação com a conformação relaxada: 0,25 aceita desvios de até 25% por ligação.',
    'angle_tolerance':'Comparação com a conformação relaxada, em graus.',
    'planarity_tolerance':'Distância máxima dos átomos de anéis aromáticos e ligações duplas ao plano ajustado.',
    'minimization_steps':'Limite de iterações da minimização usada como referência geométrica e energética.',
}


def _identifiers(value, field):
    if value is None:
        return
    if not isinstance(value,list) or len(value) > MAX_SELECTED:
        raise ValueError('Selecione até '+str(MAX_SELECTED)+' identificadores em '+field+'.')
    for item in value:
        if not isinstance(item,str) or not IDENTIFIER.fullmatch(item):
            raise ValueError('Identificador inválido em '+field+'.')


def _pairs(value, field):
    if value is None:
        return
    if not isinstance(value,list) or len(value) > MAX_SELECTED:
        raise ValueError('Selecione até '+str(MAX_SELECTED)+' complexos em '+field+'.')
    for item in value:
        if not isinstance(item,str) or not PAIR.fullmatch(item):
            raise ValueError('Complexo inválido em '+field+'.')


def validate_parameters(parameters):
    for field in BOOLEANS:
        if field in parameters and type(parameters[field]) is not bool:
            raise ValueError(LABELS[field]+' deve ser verdadeiro ou falso.')
    for field,(low,high) in NUMERIC.items():
        if parameters.get(field) is None:
            continue
        value = parameters[field]
        if type(value) not in (int,float) or not math.isfinite(value) or not low <= value <= high:
            raise ValueError(LABELS[field]+f' deve ser um número entre {low:g} e {high:g}.')
    for field,(low,high) in INTEGER.items():
        if parameters.get(field) is None:
            continue
        if type(parameters[field]) is not int or not low <= parameters[field] <= high:
            raise ValueError(LABELS[field]+f' deve ser um inteiro entre {low} e {high}.')
    _identifiers(parameters.get('selected_receptors'),'alvos selecionados')
    _identifiers(parameters.get('selected_ligands'),'ligantes selecionados')
    _pairs(parameters.get('excluded_pairs'),'complexos removidos')
    poses = parameters.get('pose_inputs')
    if poses is None or poses == 'pending-input':
        return
    if not isinstance(poses,list) or len(poses) > MAX_POSE_INPUTS:
        raise ValueError('Selecione até '+str(MAX_POSE_INPUTS)+' poses para validar em um bloco.')
    for entry in poses:
        if not isinstance(entry,dict) or not isinstance(entry.get('pose'),str) or not entry['pose']:
            raise ValueError('Cada pose a validar precisa de um arquivo de origem.')
        for key in ('receptor','reference'):
            if entry.get(key) is not None and not isinstance(entry[key],str):
                raise ValueError('Caminho inválido de receptor na validação de poses.')
        for key in ('pose_id','receptor_id','ligand_id','compound_id'):
            if entry.get(key) is not None and not isinstance(entry[key],str):
                raise ValueError('Identificador inválido na validação de poses.')
