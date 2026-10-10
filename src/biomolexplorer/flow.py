"""Visual workflow contracts, independent of UI controls and scientific imports."""
import copy
import math

from .catalog import LABELS, operation_fields, field_label
from .pipeline import validate_pipeline
from .bindings import sources, dependencies, pack

NODE_WIDTH = 290
CANVAS_WIDTH = 2200
CANVAS_HEIGHT = 40000
ROW_SPACING = 320
PORT_LABELS = {'base_input_path':'Dados de entrada','base_selected_mols':'Compostos',
               'similarity_path':'Similaridades',
               'base_vina_path':'Resultados Vina','base_dock6_path':'Resultados DOCK6'}
OUTPUTS = {
    'retrieve_pubchem': {'compounds'},
    'retrieve_compounds': {'compounds','chembl'}, 'expand_similar_compounds': {'compounds'},
    'retrieve_zinc': {'compounds'}, 'retrieve_structures': {'structures'},
    'prepare_structures': {'prepared_structures','compounds'}, 'redocking': {'structures','prepared_structures','redocking'},
    'admet': {'compounds'}, 'fingerprints': {'fingerprints','compounds'}, 'similarity': {'similarity'},
    'graphs': {'compounds'}, 'docking_vina': {'vina','compounds'}, 'docking_dock6': {'dock6','compounds'}, 'consensus': {'scores','compounds'},
}
INPUTS = {
    'retrieve_pubchem': {'base_input_path':{'compounds'}},
    'expand_similar_compounds': {'base_input_path':{'chembl'}},
    'retrieve_zinc': {'base_input_path':{'zinc_urls','other'}},
    'admet': {'base_input_path':{'compounds'}}, 'fingerprints': {'base_input_path':{'compounds'}},
    'similarity': {'base_input_path':{'fingerprints'}},
    'graphs': {'similarity_path':{'similarity'}},
    'prepare_structures': {'base_input_path':{'structures','prepared_structures'},'base_selected_mols':{'compounds'}},
    'redocking': {'base_input_path':{'structures','prepared_structures'}},
    'docking_vina': {'base_input_path':{'structures','prepared_structures'},'base_selected_mols':{'compounds','vina','dock6'}},
    'docking_dock6': {'base_input_path':{'structures','prepared_structures'},'base_selected_mols':{'compounds','vina','dock6'}},
    'consensus': {'base_vina_path':{'vina'},'base_dock6_path':{'dock6'}},
}
# These inputs are configured with uploaded files, never with canvas connections.
EXTERNAL_INPUTS = {'graphs': {'base_input_path': {'compounds'}}}


def input_types(stage, field):
    if stage['operation']=='redocking' and field=='base_input_path':
        return {'structures'} if stage['parameters'].get('prepare_complex',True) else {'prepared_structures'}
    return INPUTS.get(stage['operation'],{}).get(field,set())


def output_types(stage):
    if stage.get('provided_results'):
        return {stage['provided_results']['kind']}
    if stage['operation']=='import_results':
        params=stage['parameters']
        types=params.get('asset_types')
        return {types.get(asset,params.get('kind','other')) for asset in params.get('asset_ids',[])} if types else {params.get('kind','other')}
    return OUTPUTS[stage['operation']]


def input_ports(stage):
    operation=stage['operation']
    if operation=='import_results' or stage.get('provided_results'):
        return []
    if operation=='retrieve_pubchem' and stage['parameters'].get('reference_source')=='manual':return []
    fields={f['name']:f for f in operation_fields(operation)}
    return [{'field':field,'label':field_label(operation,field) if field=='base_input_path' else PORT_LABELS[field], 'types':input_types(stage,field),
             'required': fields[field]['required'] or operation in ('retrieve_pubchem','expand_similar_compounds','consensus') or operation=='prepare_structures' and field=='base_selected_mols'}
            for field,types in INPUTS.get(operation,{}).items()]


def compatible(source,target,field):
    if source['operation']=='prepare_structures' and target['operation'] in ('docking_vina','docking_dock6'):
        engine='vina' if target['operation']=='docking_vina' else 'dock6'
        if source['parameters'].get('docking_engines','both') not in ('both',engine):return False
    if target['operation']=='graphs' and source['operation']!='similarity':
        return False
    return source['id']!=target['id'] and bool(output_types(source) & input_types(target,field))


def connect(stages, source_id, target_id, field):
    by_id={s['id']:s for s in stages}
    if source_id not in by_id or target_id not in by_id:
        raise ValueError('O bloco da conexão não existe.')
    source,target=by_id[source_id],by_id[target_id]
    if not compatible(source,target,field):
        raise ValueError('Dados incompatíveis. Escolha uma saída do tipo esperado por esta entrada.')
    if not source.get('enabled',True):
        raise ValueError('Ative o bloco de origem antes de conectá-lo.')
    candidate=copy.deepcopy(stages)
    proposed=next(s for s in candidate if s['id']==target_id)
    old=proposed.setdefault('bindings',{}).get(field)
    proposed['bindings'][field]=pack((sources(old) if old else [])+[{'stage':source_id,'selector':'auto'}])
    if target['operation'] in ('prepare_structures','docking_vina','docking_dock6') and field=='base_input_path':
        proposed['parameters']['receptor_prepared']='prepared_structures' in output_types(source)
    if target['operation']=='consensus' and field=='base_vina_path':
        proposed['bindings']['base_input_path']=copy.deepcopy(proposed['bindings']['base_vina_path'])
    validate_pipeline(candidate)
    target['bindings']=proposed['bindings']
    target['parameters']=proposed['parameters']
    target['parameters'].pop(field,None)
    if field=='base_input_path' and 'target' in target['parameters'] and 'target' in source['parameters']:
        target['parameters']['target']=source['parameters']['target']
    if target['operation']=='expand_similar_compounds':
        target['parameters']['search_term']=source['parameters']['search_term']
    if target['operation']=='retrieve_pubchem':target['parameters']['reference_source']='compounds'


def disconnect(stage,field):
    stage.get('bindings',{}).pop(field,None)
    if stage['operation']=='consensus' and field=='base_vina_path':
        stage.get('bindings',{}).pop('base_input_path',None)


def remove(stages, stage_id):
    for stage in stages:
        for field,binding in list(stage.get('bindings',{}).items()):
            remaining=[b for b in sources(binding) if b.get('stage')!=stage_id]
            if remaining:
                stage['bindings'][field]=pack(remaining)
            else:
                disconnect(stage,field)
        stage['depends_on']=[s for s in stage.get('depends_on',[]) if s!=stage_id]
    stages[:]=[s for s in stages if s['id']!=stage_id]


def node_height(stage):
    return 168+len(input_ports(stage))*34


def position(stage,index=0):
    value=stage.get('position')
    if isinstance(value,dict) and all(type(value.get(k)) in (int,float) and math.isfinite(value[k]) for k in ('x','y')):
        return bounded_position(value['x'],value['y'],stage)
    return bounded_position(60+(index%5)*350,64+(index//5)*ROW_SPACING,stage)


def bounded_position(x,y,stage):
    return {'x':max(24,min(CANVAS_WIDTH-NODE_WIDTH-24,float(x))),
            'y':max(24,min(CANVAS_HEIGHT-node_height(stage)-24,float(y)))}


def available_position(stages,stage):
    """Find a free cell after deletion or movement instead of stacking new nodes."""
    for index in range(500):
        candidate=bounded_position(60+(index%5)*350,64+(index//5)*ROW_SPACING,stage)
        def overlaps(other,i):
            p=position(other,i)
            return abs(candidate['x']-p['x'])<NODE_WIDTH+20 and abs(candidate['y']-p['y'])<max(node_height(stage),node_height(other))+20
        if not any(overlaps(other,i) for i,other in enumerate(stages)):
            return candidate
    return bounded_position(60+(len(stages)%5)*350,64+(len(stages)//5)*ROW_SPACING,stage)


def arrange(stages):
    order=validate_pipeline(stages)
    levels,rows={},{}
    by_id={s['id']:s for s in stages}
    for stage_id in order:
        stage=by_id[stage_id]
        parents=dependencies(stage)
        level=max((levels[p]+1 for p in parents),default=0)
        levels[stage_id]=level
        row=rows.get(level,0)
        rows[level]=row+1
    band_y=64
    for band in range(max(levels.values(),default=0)//5+1):
        counters={}
        for stage in stages:
            level=levels[stage['id']]
            if level//5!=band: continue
            row=counters.get(level,0); counters[level]=row+1
            stage['position']=bounded_position(60+(level%5)*350,band_y+row*ROW_SPACING,stage)
        band_y+=max(counters.values(),default=1)*ROW_SPACING+80



def stage_issues(stage):
    issues=[]
    if not stage.get('enabled',True):
        return ['Bloco desativado']
    if stage.get('provided_results'):
        return []
    if stage['operation']=='import_results':
        return [] if stage['parameters'].get('asset_ids') else ['Selecione seus arquivos']
    if stage['operation']=='graphs' and not any(stage.get('bindings',{}).get(k) or stage['parameters'].get(k)
            for k in ('similarity_path','graph_inputs')):
        issues.append('Conecte Calcular similaridade ou envie um CSV de similaridade')
    for port in input_ports(stage):
        if port['required'] and not stage.get('bindings',{}).get(port['field']) and not stage['parameters'].get(port['field']):
            issues.append('Conecte '+port['label'].lower())
    for field in operation_fields(stage['operation']):
        if field['required'] and not field['path'] and stage['parameters'].get(field['name']) in (None,''):
            issues.append('Configure '+LABELS.get(field['name'],field['name']).lower())
    return issues
