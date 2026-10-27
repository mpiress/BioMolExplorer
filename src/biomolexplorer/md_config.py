"""Validated DOCK6 AMBER-score and optional AmberTools trajectory settings."""
import math
import re

MANUAL='https://dock.docking.org/DOCK_6/dock6_manual.htm#AMBERScore'
DEFAULTS=dict(amber_score_type='dock6',amber_score_movable_region='ligand',
    amber_score_minimization_rmsgrad=.01,amber_score_before_md_minimization_cycles=100,
    amber_score_md_steps=3000,amber_score_after_md_minimization_cycles=100,
    amber_score_gb_model=5,amber_score_nonbonded_cutoff=18.,amber_score_temperature=300.,
    amber_score_movable_distance_cutoff=10.,amber_score_receptor_movable_atom_expr='::',
    amber_score_complex_movable_atom_expr='::',amber_score_abort_on_unprepped_ligand=True,
    use_existing_ligand_charges=False,md_time_step=.001,trajectory_stride=10,
    md_heating_steps=1000,md_equilibration_steps=1000,
    random_seed=2026,md_movable_mask=':1',analyze_trajectory=True)
HELP={
 'amber_score_type':'dock6: AMBER Score nativo, com minimização–MD–minimização e GB/SA. dock6_sander: também executa uma trajetória do complexo com sander; requer AmberTools. Esta opção da plataforma não é uma palavra-chave do DOCK6.',
 'amber_score_movable_region':'Região móvel do AMBER Score: ligand (ligante), everything (todos), nothing (energia sem MD), distance (resíduos próximos ao sítio), nab_atom_expression (expressões NAB).',
 'amber_score_minimization_rmsgrad':'Critério de convergência RMS do gradiente na minimização; número finito positivo. Padrão DOCK6: 0.01.',
 'amber_score_before_md_minimization_cycles':'Ciclos de gradiente conjugado antes da dinâmica; inteiro maior ou igual a zero. Padrão: 100.',
 'amber_score_md_steps':'Passos de dinâmica; inteiro positivo. O AMBER Score nativo usa 0.001 ps por passo. O padrão 3000 corresponde a um refinamento curto de 3 ps, não a uma simulação de equilíbrio validada.',
 'amber_score_after_md_minimization_cycles':'Ciclos de minimização após a dinâmica; inteiro maior ou igual a zero. Padrão: 100.',
 'amber_score_gb_model':'Solvente implícito: 1 (HCT), 2 (OBC), 5 (OBC modificado). A energia final do AMBER Score inclui o termo de área superficial LCPO.',
 'amber_score_nonbonded_cutoff':'Corte das interações não ligadas, em angstroms; número positivo. Padrão DOCK6: 18 Å.',
 'amber_score_temperature':'Temperatura em kelvin; número positivo. Padrão: 300 K.',
 'amber_score_movable_distance_cutoff':'Distância em angstroms para selecionar resíduos móveis no modo distance. As esferas do sítio são geradas das coordenadas do ligante selecionado.',
 'amber_score_receptor_movable_atom_expr':'Expressão NAB do receptor, usada em nab_atom_expression. Exemplo: ::C,CA,N,O. Não é uma máscara cpptraj.',
 'amber_score_complex_movable_atom_expr':'Expressão NAB do complexo. No complexo, o ligante vem antes do receptor na numeração de cadeias NAB.',
 'amber_score_abort_on_unprepped_ligand':'Interrompe a execução DOCK6 do par se faltar preparação AMBER. Os outros pares continuam e a falha permanece visível.',
 'use_existing_ligand_charges':'Reutiliza as cargas do MOL2 na preparação (-c do prepare_amber.pl). Desativado: calcula AM1-BCC com as ferramentas fornecidas pelo DOCK6.',
 'md_time_step':'Passo de integração sander, em ps. Use 0.001 ps como padrão; máximo 0.002 ps. Não altera o passo fixo do AMBER Score nativo.',
 'trajectory_stride':'Intervalo de gravação da trajetória sander, em passos; inteiro positivo, não maior que o total de passos.',
 'md_heating_steps':'Passos AmberTools/sander para aquecimento de 0 K até a temperatura selecionada.',
 'md_equilibration_steps':'Passos AmberTools/sander de equilíbrio na temperatura selecionada antes da produção.',
 'random_seed':'Semente inteira positiva para a trajetória sander. Faz parte da identidade do experimento e do cache.',
 'md_movable_mask':'Máscara Amber da região móvel na trajetória sander: :1 para o ligante de um resíduo; vazio para todos os átomos. Confirme o mapeamento no PDB AMBER preparado. É independente das expressões NAB.',
 'analyze_trajectory':'RMSD do backbone e do ligante após ajuste da proteína, RMSF C-alfa, ligações H da interface, raio de giro, SASA e contatos. Requer cpptraj e dois quadros. Não estima energia livre.',
}
INTEGER={'amber_score_before_md_minimization_cycles','amber_score_md_steps','amber_score_after_md_minimization_cycles',
         'md_heating_steps','md_equilibration_steps','trajectory_stride','random_seed'}
BOOLEANS={'amber_score_abort_on_unprepped_ligand','use_existing_ligand_charges','analyze_trajectory'}
CHOICES={'amber_score_type':['dock6','dock6_sander'],
 'amber_score_movable_region':['ligand','everything','nothing','distance','nab_atom_expression'],
 'amber_score_gb_model':[1,2,5]}


def validate_settings(parameters):
    settings={key:parameters.get(key,value) for key,value in DEFAULTS.items()}
    for key,choices in CHOICES.items():
        if settings[key] not in choices or type(settings[key]) is bool:raise ValueError(f'{key}: selecione uma opção válida.')
    for key in INTEGER:
        minimum=0 if key.endswith('minimization_cycles') else 1
        if type(settings[key]) is not int or settings[key]<minimum:raise ValueError(f'{key}: informe um inteiro >= {minimum}.')
    for key in BOOLEANS:
        if type(settings[key]) is not bool:raise ValueError(f'{key}: use verdadeiro ou falso.')
    for key in ('amber_score_minimization_rmsgrad','amber_score_nonbonded_cutoff','amber_score_temperature','amber_score_movable_distance_cutoff','md_time_step'):
        value=settings[key]
        if type(value) not in (int,float) or not math.isfinite(value) or value<=0:raise ValueError(f'{key}: informe um número finito positivo.')
    if settings['md_time_step']>.002:raise ValueError('md_time_step: o limite é 0.002 ps.')
    if settings['amber_score_type']=='dock6_sander':
        if settings['analyze_trajectory'] and settings['amber_score_md_steps']//settings['trajectory_stride']<2:
            raise ValueError('A análise requer ao menos dois quadros gravados; ajuste os passos ou trajectory_stride.')
        if settings['trajectory_stride']>settings['amber_score_md_steps']:raise ValueError('trajectory_stride não pode exceder amber_score_md_steps.')
        if settings['amber_score_movable_region']=='nothing':raise ValueError('Escolha uma região móvel para executar a trajetória sander.')
    for key in ('amber_score_receptor_movable_atom_expr','amber_score_complex_movable_atom_expr','md_movable_mask'):
        value=settings[key]
        if not isinstance(value,str) or any(c in value for c in ('\n','\r','\0',"'",'"')):raise ValueError(f'{key}: expressão inválida; use uma única linha sem aspas.')
    for key in ('selected_pairs','excluded_pairs'):
        rows=parameters.get(key) or []
        if not isinstance(rows,list) or len(rows)>100000 or any(not isinstance(v,str) or not re.fullmatch('[a-f0-9]{64}',v) for v in rows):raise ValueError(f'{key}: lista de pares inválida.')
    return settings


def execution_settings(settings):
    """Only options used by the chosen protocol determine per-pair reuse."""
    return {key:value for key,value in settings.items() if settings['amber_score_type']=='dock6_sander'
            or key.startswith('amber_score_') or key=='use_existing_ligand_charges'}
