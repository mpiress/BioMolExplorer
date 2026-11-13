"""Resolve target-specific consensus actions through authorized producer tables."""
import csv
import json
import math
import re
from functools import lru_cache
from pathlib import Path

from .workspace import AccessDenied


def target_details(row):
    raw=row.get('target_scores') or '[]'
    try:details=json.loads(raw) if isinstance(raw,str) else raw
    except (ValueError,TypeError):raise ValueError('Invalid consensus target metadata') from None
    if not isinstance(details,list):raise ValueError('Invalid consensus target metadata')
    result={}
    for item in details:
        if not isinstance(item,dict) or not isinstance(item.get('receptor_id'),str):
            raise ValueError('Invalid consensus target metadata')
        if item['receptor_id'] in result:raise ValueError('Duplicate consensus target metadata')
        result[item['receptor_id']]=item
    if not result:
        identity=row.get('receptor_id','')
        result[identity]={'receptor_id':identity}
    return result


@lru_cache(maxsize=4)
def _table_index(filename,stamp):
    from .input_validation import ALIASES
    from .result_index import file_stamp
    path=Path(filename);result={}
    with path.open(encoding='utf-8-sig',newline='') as stream:
        for index,original in enumerate(csv.DictReader(stream)):
            row={ALIASES.get(key,key):value for key,value in original.items()}
            code=row.get('molecule_chembl_id');target=row.get('receptor_id','')
            value=float(row['score'])
            if not code or not math.isfinite(value):raise ValueError('Invalid producer score table')
            pose=row.get('conformer_file') or ''
            if pose and not Path(pose).is_absolute():pose=str(path.parent/pose)
            candidate=(value,pose,index)
            key=(code,target)
            if key not in result or candidate<result[key]:result[key]=candidate
    if file_stamp(path)!=stamp:raise ValueError('Producer scores changed while reading')
    return result


def target_route(service,token,pid,rid,sid,table,index,version,target,column):
    """Return an independently authorized producer selection, or the saved pose."""
    from .result_index import file_stamp
    from .result_state import input_digest
    from .bindings import sources
    if type(index) is not int or index<0:raise AccessDenied('Conformação não autorizada.')
    path=service._path(token,pid,rid,sid,table)
    row=service._record(path,index,version,root=service.store.project_dir(pid))
    run,stage=service._stage(token,pid,rid,sid)
    if stage.get('operation')!='consensus' or not isinstance(target,str) or target not in target_details(row):
        raise AccessDenied('Alvo não autorizado para este resultado.')
    # A retained representative remains independently usable even if its
    # original producer was removed. Other targets never fall back to it.
    if target==row.get('receptor_id') and row.get(column):
        saved=service.store.scoped_path(pid,Path(row[column]) if Path(row[column]).is_absolute() else path.parent/row[column])
        if str(saved) in service._allowed(pid,stage) and saved.is_file():return None
    engine='vina' if column=='vina_pose' else 'dock6'
    port='base_'+engine+'_path'
    tables=stage.get('input_files',{}).get(port,[])
    if not tables and engine=='vina':tables=stage.get('input_files',{}).get('base_input_path',[])
    origins={ref.get('stage') for ref in sources(stage.get('configuration',{}).get('bindings',{}).get(port,{}))}
    if not origins and engine=='vina':
        origins={ref.get('stage') for ref in sources(stage.get('configuration',{}).get('bindings',{}).get('base_input_path',{}))}
    candidates=[]
    expected=target_details(row)[target].get(engine+'_raw')
    from .published_results import _input_hashes
    original_hashes=_input_hashes(stage)
    recorded=original_hashes.get(port) or (original_hashes.get('base_input_path') if engine=='vina' else None)
    for filename in dict.fromkeys(tables):
        source=service.store.scoped_path(pid,filename)
        if source.suffix.lower()!='.csv' or not source.is_file():continue
        producers=[item for item in run['stages'] if str(source) in item.get('artifacts',[]) and
                   (not origins or item['id'] in origins)]
        source_rid=rid;child=None
        if not producers:
            relative=source.relative_to(service.store.project_dir(pid)).parts
            if len(relative)<4 or relative[0]!='runs':continue
            source_rid,source_sid=relative[1:3]
            if origins and source_sid not in origins:continue
            cache=getattr(service,'_target_origins',{})
            child=cache.get((source_rid,source_sid))
            if child is None:
                child=type(service)(service.store)
                cache[(source_rid,source_sid)]=child
                service._target_origins=cache
            _,producer=child._stage(token,pid,source_rid,source_sid)
            if str(source) not in producer.get('artifacts',[]):continue
        elif len(producers)==1:producer=producers[0]
        else:continue
        if producer.get('operation')=='consensus':continue
        digest=input_digest(source)
        if recorded and digest not in recorded:continue
        manifest=producer.get('artifact_manifest',{})
        if not recorded and not manifest.get(str(source)):continue
        if manifest and manifest.get(str(source))!=digest:continue
        match=_table_index(str(source),file_stamp(source)).get((row['molecule_chembl_id'],target))
        if match is None or expected is not None and match[0]!=float(expected):continue
        candidates.append((match[0],match[1],str(source),match[2],digest,producer,source_rid,child))
    if not candidates:raise ValueError('Resultado deste alvo indisponível na origem.')
    value,pose,source,position,digest,producer,source_rid,child=min(candidates,key=lambda item:item[:4])
    if child is None:
        child=type(service)(service.store)
        key=list(service._stage_snapshot[0]);key[2]=producer['id']
        child._stage_snapshot=(tuple(key),run,producer)
    # Reuse this run's already-authorized read model; child._stage still checks
    # the live token and revision before using it.
    selection=dict(table=source,index=position,version=digest,column='conformer_file')
    actual=child._record(Path(source),position,digest,root=service.store.project_dir(pid))
    if actual.get('molecule_chembl_id')!=row['molecule_chembl_id'] or actual.get('receptor_id','')!=target:
        raise AccessDenied('Alvo não autorizado para este resultado.')
    fields=('footprint_file',) if column=='footprint_file' else ('conformer_file','receptor_file','reference_file')
    for field in fields:
        linked=actual.get(field)
        if not linked:continue
        linked=service.store.scoped_path(pid,Path(linked) if Path(linked).is_absolute() else Path(source).parent/linked)
        sha=producer.get('artifact_manifest',{}).get(str(linked))
        if sha and input_digest(linked)!=sha:raise ValueError('O arquivo de origem foi alterado.')
    return child,source_rid,producer['id'],selection


def target_label(service,pid,run,row,target):
    """Prefer recorded protein names, using local PDB headers as a fallback."""
    paths=[]
    pdb=target.split('_',1)[0]
    cache=getattr(service,'_target_label_paths',None)
    key=service._stage_snapshot[0]
    if cache is None or cache[0]!=key:
        by_name={}
        for item in run['stages']:
            for path in item.get('artifacts',[]):
                if str(path).lower().endswith('.pdb'):by_name.setdefault(Path(path).name.upper(),[]).append(path)
        service._target_label_paths=(key,by_name)
    paths.extend(service._target_label_paths[1].get(pdb.upper()+'.PDB',[]))
    for field in ('vina_receptor_file','dock6_receptor_file','receptor_file'):
        if target==row.get('receptor_id') and row.get(field) and Path(row[field]).is_absolute():paths.insert(0,row[field])
    text=' '.join(str(row.get(key) or '') for key in ('biological_target','target_name','receptor_name')) if target==row.get('receptor_id') else ''
    for value in dict.fromkeys(paths):
        path=service.store.scoped_path(pid,value)
        if not path.is_file():continue
        with path.open(errors='replace') as stream:text+=' '+stream.read(16384)
        label=_protein_label(text)
        if label:return label+' · '+target
    return (_protein_label(text)+' · '+target) if _protein_label(text) else target


def _protein_label(text):
    value=text.upper()
    for label,patterns in (
        ('BChE',('BUTYRYLCHOLINESTERASE',r'\bBCHE\b')),
        ('AChE',('ACETYLCHOLINESTERASE',r'\bACHE\b')),
        ('BACE1',('BETA-SECRETASE', 'BETA SECRETASE', r'\bBACE1\b'))):
        if any(re.search(pattern,value) for pattern in patterns):return label
    return None
