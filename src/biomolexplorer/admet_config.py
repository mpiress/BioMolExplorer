"""Research criteria, independent of scientific/UI dependencies.

Profiles are editable screening hypotheses, not validated ADMET models.
"""
import copy
import math

PROFILES = {
    'exploratory': {'label':'Exploração ampla', 'description':'Descreva todos os compostos válidos; alertas não removem moléculas.'},
    'oral': {'label':'Candidatos orais', 'description':'Hipótese de absorção passiva favorável, Lipinski com até uma violação e Veber.',
             'hia_target':'HIA+', 'rule_set':'oral', 'max_rule_violations':1},
    'cns': {'label':'Acesso ao SNC', 'description':'Priorize BBB+ pela heurística local. Efluxo, ionização e exposição livre precisam de avaliação adicional.', 'bbb_target':'BBB+'},
    'peripheral': {'label':'Restrição periférica', 'description':'Priorize BBB− pela heurística local; isso não comprova ausência de exposição cerebral.', 'bbb_target':'BBB-'},
    'fragment': {'label':'Descoberta por fragmentos', 'description':'Use a regra de três com extensões de flexibilidade e polaridade; ajuste para sua biblioteca.', 'rule_set':'fragment', 'max_rule_violations':0},
    'natural_products': {'label':'Produtos naturais / além da regra de cinco', 'description':'Explore espaço químico amplo sem impor regras de fármacos orais. Ajuste os limites à hipótese.'},
    'custom': {'label':'Hipótese personalizada', 'description':'Combine barreiras, regras, alertas e intervalos definidos por você.'},
    'legacy': {'label':'Compatibilidade com ADMET anterior', 'description':'Reproduz a exclusão automática histórica. Os antigos padrões são alertas, não comprovação de toxicidade.'},
}
DESCRIPTORS = {
    'MW':'Massa molecular (g/mol)', 'TPSA':'TPSA (Å²)', 'WLOGP':'WLOGP (Wildman–Crippen)',
    'HBD':'Doadores de H', 'HBA':'Aceptores de H', 'RB':'Ligações rotacionáveis',
    'QED':'QED (0–1)', 'FractionCSP3':'Fração de carbonos sp³ (0–1)',
    'HeavyAtoms':'Átomos pesados', 'AromaticRings':'Anéis aromáticos', 'Rings':'Anéis',
    'FormalCharge':'Carga formal', 'ESOL_logS':'logS estimado (log10 mol/L)',
    'logKp':'logKp estimado (log10 cm/s)',
    'Lipinski_violations':'Violações de Lipinski', 'Veber_violations':'Violações de Veber',
    'Fragment_violations':'Violações da regra de três ampliada',
}
OPTIONS = {
    'admet_profile':list(PROFILES), 'admet_mode':['annotate','filter'],
    'prediction_backend':['rdkit','admet_ai'], 'bbb_target':['any','BBB+','BBB-'],
    'hia_target':['any','HIA+','HIA-'], 'rule_set':['none','lipinski','veber','oral','fragment'],
    'alert_policy':['annotate','exclude'],
}
PARAMETERS=set(OPTIONS)|{'max_rule_violations','alert_catalogs','property_limits'}
DEFAULTS=dict(admet_profile='legacy',admet_mode='filter',prediction_backend='rdkit',
              bbb_target='any',hia_target='any',rule_set='none',max_rule_violations=0,
              alert_policy='annotate',alert_catalogs=['PAINS','BRENK'],property_limits={})


def resolve_config(parameters):
    profile=parameters.get('admet_profile','legacy')
    if profile not in PROFILES:raise ValueError('Perfil ADMET inválido.')
    result=copy.deepcopy(DEFAULTS)
    result.update({k:v for k,v in PROFILES[profile].items() if k in PARAMETERS})
    result.update({k:v for k,v in parameters.items() if k in PARAMETERS and v is not None})
    for key,values in OPTIONS.items():
        if result[key] not in values:raise ValueError('Opção ADMET inválida: '+key)
    count=result['max_rule_violations']
    if type(count) is not int or not 0<=count<=6:raise ValueError('O número de violações permitidas deve ser um inteiro entre 0 e 6.')
    catalogs=result['alert_catalogs']
    if not isinstance(catalogs,list) or any(v not in ('PAINS','BRENK','NIH') for v in catalogs):
        raise ValueError('Escolha catálogos PAINS, BRENK ou NIH.')
    limits=result['property_limits']
    if not isinstance(limits,dict):raise ValueError('Os limites ADMET devem ser um mapa de intervalos.')
    for key,interval in limits.items():
        if key not in DESCRIPTORS and not (isinstance(key,str) and key.startswith('ML_') and len(key)>3 and result['prediction_backend']=='admet_ai'):
            raise ValueError('Propriedade ADMET desconhecida: '+str(key))
        if not isinstance(interval,(list,tuple)) or len(interval)!=2:
            raise ValueError('Cada propriedade ADMET deve ter mínimo e máximo.')
        if any(v is not None and (type(v) not in (int,float) or not math.isfinite(v)) for v in interval):
            raise ValueError('Os limites ADMET devem ser números finitos ou nulos.')
        lo,hi=interval
        if lo is not None and hi is not None and lo>hi:raise ValueError('O mínimo ADMET não pode exceder o máximo.')
    return result


def selection_reasons(properties,config):
    reasons=[]
    for key,target in (('BBB',config['bbb_target']),('HIA',config['hia_target'])):
        if target!='any' and properties.get(key)!=target:reasons.append(key+': requer '+target)
    rules={'none':[], 'lipinski':['Lipinski'], 'veber':['Veber'], 'oral':['Lipinski','Veber'], 'fragment':['Fragment']}[config['rule_set']]
    for rule in rules:
        # Oral combines the two independently: Veber has no allowed violations.
        allowed=0 if config['rule_set']=='oral' and rule=='Veber' else config['max_rule_violations']
        if properties[rule+'_violations']>allowed:reasons.append(rule+': violações acima de '+str(allowed))
    if config['alert_policy']=='exclude' and properties['StructuralAlert_count']:
        reasons.append('Alertas estruturais: '+properties['StructuralAlerts'])
    for key,(lo,hi) in config['property_limits'].items():
        value=properties.get(key)
        if type(value) not in (int,float) or not math.isfinite(value):reasons.append(key+': valor indisponível')
        elif (lo is not None and value<lo) or (hi is not None and value>hi):reasons.append(key+': fora do intervalo')
    return reasons
