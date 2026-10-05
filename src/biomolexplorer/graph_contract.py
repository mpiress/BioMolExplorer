"""Upgrade saved graph configurations without changing stored experiment files."""
import copy
from .bindings import sources, pack


def normalize_legacy_graphs(pipeline):
    pipeline=copy.deepcopy(pipeline)
    if not isinstance(pipeline,list):return pipeline
    lookup={s['id']:s for s in pipeline if isinstance(s,dict) and isinstance(s.get('id'),str)}
    for stage in pipeline:
        if not isinstance(stage,dict) or stage.get('operation')!='graphs':continue
        if not isinstance(stage.get('parameters'),dict) or not isinstance(stage.get('bindings',{}),dict):continue
        for name in ('fingerprints_path','metric','fingerprint','threshold'):
            stage['parameters'].pop(name,None)
        stage.get('bindings',{}).pop('fingerprints_path',None)
        for field in ('base_input_path','similarity_path'):
            binding=stage.get('bindings',{}).get(field)
            if not isinstance(binding,dict) or not binding:continue
            references=sources(binding)
            if not isinstance(references,list):continue
            kept=[]
            for ref in references:
                if not isinstance(ref,dict) or ('stage' in ref and (not isinstance(ref['stage'],str) or ref['stage'] not in lookup)):
                    kept.append(ref)
                    continue
                if 'asset' in ref or (field=='similarity_path' and lookup.get(ref.get('stage'),{}).get('operation')=='similarity'):
                    kept.append(ref)
                elif field=='similarity_path':
                    producer=lookup.get(ref.get('stage'),{})
                    if producer.get('operation')=='import_results' and producer['parameters'].get('kind')=='similarity':
                        kept.extend({'asset':asset,'selector':'auto'} for asset in producer['parameters'].get('asset_ids',[]))
            if kept:stage['bindings'][field]=pack(kept)
            else:stage['bindings'].pop(field,None)
    return pipeline
