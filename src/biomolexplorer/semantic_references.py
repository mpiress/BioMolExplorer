"""Validation and migration of logical artifact and entity references."""
import re
from dataclasses import dataclass,field

from .artifact_contracts import OUTPUT_PORTS, port_types

REFERENCE_KEYS={'stage','asset','output_port','selection','mode','run_id','artifact_id','selector','compound_id'}
SELECTION_KEYS={'entity_ids','pdb_id','chain_id','assembly_id','ligand_id','receptor_id','pose_id','community_id','properties'}


@dataclass(frozen=True)
class ArtifactReference:
    source_stage: str
    output_port: str
    mode: str = 'current'
    run_id: str | None = None
    artifact_id: str | None = None

    def to_binding(self):
        value={'stage':self.source_stage,'output_port':self.output_port,'mode':self.mode}
        if self.run_id:value['run_id']=self.run_id
        if self.artifact_id:value['artifact_id']=self.artifact_id
        return value


@dataclass(frozen=True)
class EntityReference:
    artifact: ArtifactReference
    entity_type: str
    entity_id: str
    qualifiers: dict=field(default_factory=dict)


@dataclass(frozen=True)
class SelectionReference:
    artifact: ArtifactReference
    selection: dict

    def to_binding(self):return dict(self.artifact.to_binding(),selection=dict(self.selection))


def validate_selection(value):
    if not isinstance(value,dict) or set(value)-SELECTION_KEYS:
        raise ValueError('Seleção científica inválida.')
    for key,item in value.items():
        if key=='properties':
            if not isinstance(item,dict) or any(not isinstance(k,str) or not isinstance(v,(str,int,float,bool,list,tuple,dict)) for k,v in item.items()):
                raise ValueError('Filtro de propriedades inválido.')
        elif key=='entity_ids':
            if not isinstance(item,list) or not item or any(not isinstance(v,str) or not v.strip() for v in item):
                raise ValueError('Selecione ao menos uma entidade válida.')
        elif not isinstance(item,(str,int)) or not str(item).strip():
            raise ValueError('Identificador de entidade inválido.')
    return value


def validate_reference(ref, producer=None):
    if not isinstance(ref,dict) or set(ref)-REFERENCE_KEYS or ('stage' in ref)==('asset' in ref):
        raise ValueError('Escolha uma origem científica válida.')
    identity=ref.get('stage',ref.get('asset'))
    if not isinstance(identity,str) or not re.fullmatch('[a-f0-9]{32}',identity):
        raise ValueError('Uma conexão contém um identificador inválido.')
    if ref.get('mode','current') not in ('current','pinned'):
        raise ValueError('Modo de referência inválido.')
    if ref.get('mode')=='pinned' and not (ref.get('run_id') and ref.get('artifact_id')):
        raise ValueError('Uma referência histórica deve fixar execução e artefato.')
    if 'output_port' in ref:
        if 'stage' not in ref or not isinstance(ref['output_port'],str):
            raise ValueError('Porta de saída inválida.')
        if producer and ref['output_port'] not in OUTPUT_PORTS.get(producer['operation'],{}):
            raise ValueError('A porta de saída não pertence ao bloco produtor.')
    if 'selection' in ref:validate_selection(ref['selection'])
    return ref


def migrate_reference(ref, producer=None, accepted_types=()):
    """Add an unambiguous semantic port without changing legacy selectors."""
    result=dict(ref)
    if producer and 'stage' in result and 'output_port' not in result:
        matches=[name for name in OUTPUT_PORTS.get(producer['operation'],{})
                 if port_types(producer['operation'],name)&set(accepted_types)]
        if len(matches)==1:result['output_port']=matches[0]
    if result.get('compound_id') and 'selection' not in result:
        result['selection']={'entity_ids':[result['compound_id']]}
    return result


def migrate_pipeline_references(stages):
    """Upgrade only unambiguous legacy bindings and report every ambiguity.

    The returned pipeline is a detached copy.  Selectors remain available to
    the legacy adapter, but a new semantic port is never guessed when two
    producer outputs satisfy the same consumer contract.
    """
    import copy
    from .artifact_contracts import INPUT_PORTS
    from .bindings import pack,sources
    migrated=copy.deepcopy(stages);producers={stage['id']:stage for stage in migrated}
    issues=[];changed=False
    for consumer in migrated:
        for field,binding in list(consumer.get('bindings',{}).items()):
            accepted=INPUT_PORTS.get(consumer.get('operation'),{}).get(field,())
            updated=[]
            for position,reference in enumerate(sources(binding)):
                producer=producers.get(reference.get('stage'))
                value=migrate_reference(reference,producer,accepted)
                if producer and 'stage' in reference and 'output_port' not in value:
                    candidates=[name for name in OUTPUT_PORTS.get(producer.get('operation'),{})
                                if port_types(producer.get('operation'),name)&set(accepted)]
                    issues.append({'consumer_stage_id':consumer['id'],'input_port':field,
                        'reference_index':position,'producer_stage_id':producer['id'],
                        'candidate_output_ports':candidates,'reason':'ambiguous' if len(candidates)>1 else 'incompatible'})
                changed|=value!=reference;updated.append(value)
            packed=pack(updated)
            if packed is not None:consumer['bindings'][field]=packed
    return migrated,issues,changed


def migrate_project_references(db,project_id,pipeline=None):
    """Persist an idempotent project migration and its explicit issues."""
    import json
    if pipeline is None:
        row=db.execute('SELECT pipeline FROM projects WHERE id=?',(project_id,)).fetchone()
        if row is None:return [],False
        pipeline=json.loads(row['pipeline'])
    migrated,issues,changed=migrate_pipeline_references(pipeline)
    db.execute('DELETE FROM artifact_migration_issues WHERE project_id=?',(project_id,))
    for issue in issues:
        db.execute('''INSERT INTO artifact_migration_issues(
            project_id,consumer_stage_id,input_port,reference_index,producer_stage_id,reason,details)
            VALUES (?,?,?,?,?,?,?)''',(project_id,issue['consumer_stage_id'],issue['input_port'],issue['reference_index'],
                issue['producer_stage_id'],issue['reason'],json.dumps(issue,sort_keys=True)))
    if changed:
        db.execute('UPDATE projects SET pipeline=? WHERE id=?',(json.dumps(migrated,allow_nan=False),project_id))
    return issues,changed


def parse_selection(text):
    """Parse compact UI filters such as ``pdb:4M0E, chain:A``."""
    values=[value.strip() for value in str(text or '').split(',') if value.strip()]
    if not values:return None
    aliases={'pdb':'pdb_id','chain':'chain_id','assembly':'assembly_id','ligand':'ligand_id',
             'receptor':'receptor_id','pose':'pose_id','community':'community_id'}
    result={};entities=[]
    for value in values:
        if ':' in value:
            key,item=value.split(':',1);field=aliases.get(key.strip().casefold())
            if field:
                if field in result and result[field]!=item.strip():
                    raise ValueError('Use apenas um valor por filtro científico nomeado.')
                result[field]=item.strip();continue
        entities.append(value)
    if entities:result['entity_ids']=entities
    return validate_selection(result)


def selection_text(selection):
    if not selection:return ''
    reverse={'pdb_id':'pdb','chain_id':'chain','assembly_id':'assembly','ligand_id':'ligand',
             'receptor_id':'receptor','pose_id':'pose','community_id':'community'}
    values=[f'{reverse[key]}:{value}' for key,value in selection.items() if key in reverse]
    values.extend(selection.get('entity_ids',[]))
    return ', '.join(values)

