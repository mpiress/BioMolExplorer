"""Validated native DOCK6 6.11 Amber Score settings."""
import math
import re

MANUAL='https://dock.docking.org/DOCK_6/dock6_manual.htm#AMBERScore'
DEFAULTS=dict(amber_score_type='dock6',amber_score_movable_region='ligand',
    amber_score_minimization_rmsgrad=.01,amber_score_before_md_minimization_cycles=100,
    amber_score_md_steps=3000,amber_score_after_md_minimization_cycles=100,
    amber_score_gb_model=5,amber_score_nonbonded_cutoff=18.,amber_score_temperature=300.,
    amber_score_movable_distance_cutoff=10.,amber_score_receptor_movable_atom_expr='::',
    amber_score_complex_movable_atom_expr='::',amber_score_abort_on_unprepped_ligand=True,
    use_existing_ligand_charges=False)
HELP={
 'amber_score_type':'Amber Score nativo do DOCK6 6.11: refinamento curto em GB/SA, sem simulação convencional ou trajetória externa.',
 'amber_score_movable_region':'Região móvel do AMBER Score: ligand (ligante), everything (todos), nothing (energia sem MD), distance (resíduos próximos ao sítio), nab_atom_expression (expressões NAB).',
 'amber_score_minimization_rmsgrad':'Critério de convergência RMS do gradiente na minimização; número finito positivo. Padrão DOCK6: 0.01.',
 'amber_score_before_md_minimization_cycles':'Ciclos de gradiente conjugado antes da dinâmica; inteiro maior ou igual a zero. Padrão: 100.',
 'amber_score_md_steps':'Passos do refinamento interno; inteiro maior ou igual a zero. O AMBER Score nativo usa 0.001 ps por passo. O padrão 3000 corresponde a um refinamento curto de 3 ps, não a uma simulação de equilíbrio validada.',
 'amber_score_after_md_minimization_cycles':'Ciclos de minimização após a dinâmica; inteiro maior ou igual a zero. Padrão: 100.',
 'amber_score_gb_model':'Solvente implícito: 1 (HCT), 2 (OBC), 5 (OBC modificado). A energia final do AMBER Score inclui o termo de área superficial LCPO.',
 'amber_score_nonbonded_cutoff':'Corte das interações não ligadas, em angstroms; número positivo. Padrão DOCK6: 18 Å.',
 'amber_score_temperature':'Temperatura em kelvin; número positivo. Padrão: 300 K.',
 'amber_score_movable_distance_cutoff':'Distância em angstroms para selecionar resíduos móveis no modo distance. As esferas do sítio são geradas das coordenadas do ligante selecionado.',
 'amber_score_receptor_movable_atom_expr':'Expressão NAB do receptor, usada em nab_atom_expression. Exemplo: ::C,CA,N,O. Não é uma máscara cpptraj.',
 'amber_score_complex_movable_atom_expr':'Expressão NAB do complexo. No complexo, o ligante vem antes do receptor na numeração de cadeias NAB.',
 'amber_score_abort_on_unprepped_ligand':'Interrompe a execução DOCK6 do par se faltar preparação AMBER. Os outros pares continuam e a falha permanece visível.',
 'use_existing_ligand_charges':'Reutiliza as cargas do MOL2 na preparação (-c do prepare_amber.pl). Desativado: calcula AM1-BCC com as ferramentas fornecidas pelo DOCK6.',
}
INTEGER={'amber_score_before_md_minimization_cycles','amber_score_md_steps','amber_score_after_md_minimization_cycles'}
BOOLEANS={'amber_score_abort_on_unprepped_ligand','use_existing_ligand_charges'}
CHOICES={'amber_score_type':['dock6'],
 'amber_score_movable_region':['ligand','everything','nothing','distance','nab_atom_expression'],
 'amber_score_gb_model':[1,2,5]}


def validate_settings(parameters):
    if parameters.get('amber_score_type','dock6') != 'dock6':
        raise ValueError('Only native DOCK6 6.11 Amber Score is supported. Reconfigure this legacy external-MD block explicitly; saved results are preserved.')
    legacy_defaults=dict(md_time_step=.001,trajectory_stride=10,md_heating_steps=1000,md_equilibration_steps=1000,random_seed=2026,md_movable_mask=':1',analyze_trajectory=True)
    for key,default in legacy_defaults.items():
        if key in parameters and parameters[key]!=default:
            raise ValueError(key+': conventional-MD settings are retired; explicitly reconfigure this block for native Amber Score.')
    settings={key:parameters.get(key,value) for key,value in DEFAULTS.items()}
    for key,choices in CHOICES.items():
        if settings[key] not in choices or type(settings[key]) is bool:raise ValueError(f'{key}: selecione uma opção válida.')
    for key in INTEGER:
        minimum=0
        if type(settings[key]) is not int or not minimum<=settings[key]<=2147483647:raise ValueError(f'{key}: informe um inteiro >= {minimum}.')
    for key in BOOLEANS:
        if type(settings[key]) is not bool:raise ValueError(f'{key}: use verdadeiro ou falso.')
    for key in ('amber_score_minimization_rmsgrad','amber_score_nonbonded_cutoff','amber_score_temperature','amber_score_movable_distance_cutoff'):
        value=settings[key]
        if type(value) not in (int,float) or not math.isfinite(value) or value<=0:raise ValueError(f'{key}: informe um número finito positivo.')
    for key in ('amber_score_receptor_movable_atom_expr','amber_score_complex_movable_atom_expr'):
        value=settings[key]
        if not isinstance(value,str) or any(c in value for c in ('\n','\r','\0',"'",'"')):raise ValueError(f'{key}: expressão inválida; use uma única linha sem aspas.')
    for key in ('selected_pairs','excluded_pairs'):
        rows=parameters.get(key) or []
        if not isinstance(rows,list) or len(rows)>100000 or any(not isinstance(v,str) or not re.fullmatch('[a-f0-9]{64}',v) for v in rows):raise ValueError(f'{key}: lista de pares inválida.')
    return settings


def execution_settings(settings):
    """Only options used by the chosen protocol determine per-pair reuse."""
    return {key:value for key,value in settings.items() if key in DEFAULTS}
