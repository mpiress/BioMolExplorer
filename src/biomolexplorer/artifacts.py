"""Artifact identities and lineage, independent of scientific formats and UI.

A slot identifies a producer and a qualified logical output; an occurrence also
identifies its producing run. Paths locate bytes and are never the identity.
Legacy runs are indexed lazily, without rewriting their scientific files.
"""
import hashlib
import json
from pathlib import Path
from .artifact_choices import selector_for, matches_selector, logical_path
from collections import Counter
from dataclasses import dataclass,field,asdict
from .bindings import sources


@dataclass
class ArtifactDescriptor:
    artifact_id: str
    block_id: str
    operation: str
    output: str
    physical_file: str
    source: dict=field(default_factory=dict)
    parent_artifact_ids: list=field(default_factory=list)
    lineage: list=field(default_factory=list)
    parameters: dict=field(default_factory=dict)
    record_count: int | None=None
    provenance: dict=field(default_factory=dict)

    def to_dict(self):return asdict(self)


def slot_id(stage_id, selector):
    return hashlib.sha256(json.dumps([stage_id, selector], separators=(',', ':')).encode()).hexdigest()


def locator(stage_id, path):
    """Canonical producer-relative path, independent of run/worker directories."""
    parts=Path(logical_path(path)).parts
    if stage_id in parts:
        parts=parts[len(parts)-1-list(reversed(parts)).index(stage_id)+1:]
        # A calculation-cache hit is materialized below ``cache/redocking-N``
        # and repeats part of the producer-relative path.  That implementation
        # detail must not change the identity selected by downstream docking.
        cache=next((i for i,value in enumerate(parts[:-1])
                    if value=='cache' and parts[i+1].startswith('redocking-')),None)
        if cache is not None:
            prefix=parts[:cache];tail=parts[cache+2:];overlap=0
            for size in range(1,min(len(prefix),len(tail))+1):
                if prefix[-size:]==tail[:size]:overlap=size
            parts=prefix+tail[overlap:]
        return '/'.join(parts)
    if 'artifacts' in parts:
        return '/'.join(parts[len(parts)-1-list(reversed(parts)).index('artifacts')+1:])
    # Legacy files outside a producer directory have no known relative root.
    # Preserve their qualified location rather than inventing an ambiguous ID.
    return logical_path(path)


def records(item, run_id=None, *, decorate=True):
    files=list(dict.fromkeys(str(p) for p in item.get('artifacts', [])))
    previous={r['path']:r for r in item.get('artifact_records', [])}
    config=item.get('configuration', item)
    parts={path:Path(logical_path(path)).parts for path in files}
    suffixes={path:['/'.join(value[-n:]) for n in range(1,len(value)+1)] for path,value in parts.items()}
    counts=Counter(suffix for options in suffixes.values() for suffix in options)
    result=[]
    for path in files:
        selector=next((suffix for suffix in suffixes[path] if counts[suffix]==1),None)
        if selector is None:raise ValueError('Os resultados contêm caminhos de arquivo duplicados.')
        from .artifact_naming import publication_metadata
        publication=publication_metadata(path)
        prior=previous.get(path, {})
        identity=prior.get('artifact_id') or publication.get('artifact_id') or slot_id(item['id'],locator(item['id'],path))
        batch=next((b for b in item.get('batches',[]) if path in b.get('artifacts',[])),{})
        origin=prior.get('run_id') or batch.get('reused_from_run') or item.get('reused_from_run') or run_id
        runtime=batch.get('processing_runtime') if batch else item.get('processing_runtime')
        result.append(dict(artifact_id=identity,
            block_id=prior.get('block_id') or config.get('block_id') or item.get('block_id'),
            output=prior.get('output') or publication.get('output'),
            original_name=prior.get('original_name') or publication.get('legacy_name') or Path(path).name,
            provenance=dict(prior.get('provenance',{}),legacy_name=publication.get('legacy_name',Path(path).name)),
            record_count=prior.get('record_count',publication.get('record_count')),
            naming=prior.get('naming') or publication.get('naming',{}),
            occurrence_id=prior.get('occurrence_id') or publication.get('naming',{}).get('identifier') or hashlib.sha256(json.dumps([identity, origin]).encode()).hexdigest(),
            stage_id=item['id'], stage_name=item.get('name', item['id']),
            operation=item.get('operation', config.get('operation', 'unknown')),
            content_sha256=item.get('artifact_manifest',{}).get(str(Path(path).resolve())),
            processing_version=batch.get('processing_version') if batch else item.get('processing_version'),
            producer_software=batch.get('producer_software') if batch else item.get('producer_software'),
            runtime_id=hashlib.sha256(json.dumps(runtime,sort_keys=True).encode()).hexdigest() if runtime else None,
            dependency_manifest=item.get('input_manifest',{}),
            run_id=origin, name=Path(path).name, selector=selector, locator=locator(item['id'],path), path=path,
            processing_scope='all' if batch.get('configuration',{}).get('_combined_input') else 'individual',
            batch_label=batch.get('label'),
            context=prior.get('context', {'target':config.get('parameters', {}).get('target'),
                     'search_term':config.get('parameters', {}).get('search_term')}),
            inputs=prior.get('inputs', {})))
    if decorate:
        from .artifact_naming import name_records
        name_records([(r,item) for r in result])
    return result


def reference(record):
    return {'stage':record['stage_id'], 'selector':record['selector'],
            'artifact_id':record['artifact_id']}


def validate_reference(ref, files):
    """Verify identity and availability before any consumer materializes inputs."""
    if 'stage' not in ref or ref.get('selector', 'auto') == 'auto':
        if ref.get('artifact_id'):raise ValueError('Um artefato identificado exige um arquivo específico.')
        return
    if ref.get('artifact_id'):
        from .artifact_naming import publication_metadata
        matching=[p for p in files if publication_metadata(p).get('artifact_id',slot_id(ref['stage'],locator(ref['stage'],p)))==ref['artifact_id']]
        if len(matching)!=1:
            from .admet_outputs import legacy_reference
            migrated=legacy_reference(ref,files)
            if migrated is not None:return migrated
            raise ValueError('A identidade do arquivo selecionado mudou. Revise a seleção da origem.')
    else:
        matching=[p for p in files if matches_selector(p, ref['selector'])]
        if not matching:
            from .admet_outputs import legacy_reference
            migrated=legacy_reference(ref,files)
            if migrated is not None:return migrated
    if len(matching)!=1 or not Path(matching[0]).is_file():
        raise ValueError('O arquivo escolhido não está nos resultados da etapa (ausente ou ambíguo). Revise a seleção da origem.')
    return matching[0]


def index_run(stages, run_id, assets=()):
    """Record exact file inputs, including per-batch lineage for individual mode."""
    assets_by_path={str(a['path']):a for a in assets}
    by_path={}
    for item in stages:
        item['artifact_records']=records(item, run_id,decorate=False)
        for record in item['artifact_records']:by_path.setdefault(record['path'], []).append(record)
    for item in stages:
        batches={path:batch for batch in item.get('batches',[]) for path in batch.get('artifacts',[])}
        for record in item['artifact_records']:
            batch=batches.get(record['path'])
            if not item.get('reused') and (batch or item).get('artifact_context'):
                record['context']=dict((batch or item)['artifact_context'])
            config=(batch or {}).get('configuration', item.get('configuration', {}))
            inputs=(batch or item).get('input_files', {})
            lineage={}
            for field, paths in inputs.items():
                refs=sources(config.get('bindings', {}).get(field, {}))
                origins={r['stage'] for r in refs if 'stage' in r}
                lineage[field]=[]
                for path in paths:
                    matches=[r for r in by_path.get(str(path), []) if r['stage_id'] in origins]
                    if matches:
                        lineage[field].extend({k:r[k] for k in ('artifact_id','occurrence_id','stage_id','stage_name','run_id','name','selector','path')} for r in matches)
                    else:
                        asset=assets_by_path.get(str(path),{})
                        lineage[field].append({'path':str(path), 'name':asset.get('name',Path(path).name),
                            'asset_id':asset.get('id'), 'asset_ids':[asset['id']] if asset else []})
            # Cached records retain the execution that originally produced them.
            if lineage and not item.get('reused'):record['inputs']=lineage
            parents=[r for rows in record.get('inputs',{}).values() for r in rows]
            record['parent_artifact_ids']=list(dict.fromkeys(r['artifact_id'] for r in parents if r.get('artifact_id')))
            record['lineage']=record['parent_artifact_ids']
            record['source']=({field:[dict(r) for r in rows] for field,rows in record.get('inputs',{}).items()}
                or {key:value for key,value in record['context'].items() if value is not None})

    from .artifact_naming import name_records
    name_records([(r,s) for s in stages for r in s['artifact_records']])


def annotate_files(item, files, run_id=None):
    """Attach provenance to display records without replacing original filenames."""
    catalog={r['path']:r for r in records(item,run_id)}
    for file in files:
        record=catalog.get(file['path'])
        if not record:continue
        inputs=sorted({r.get('stage_name',r.get('name',''))
            for rows in record['inputs'].values() for r in rows})
        context=record['context'].get('target') or record['context'].get('search_term')
        from .artifact_naming import human_name
        parts=[human_name(record),record.get('block_id') or record['stage_name']]
        if context:parts.append(str(context))
        if inputs:parts.append(' + '.join(inputs))
        file.update(artifact_id=record['artifact_id'],occurrence_id=record['occurrence_id'],
                    origin_run_id=record['run_id'],label=record['public_name'],public_name=record['public_name'],
                    display_name=record['public_name'],naming=record['naming'])
    return files
