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
from .bindings import sources


def slot_id(stage_id, selector):
    return hashlib.sha256(json.dumps([stage_id, selector], separators=(',', ':')).encode()).hexdigest()


def locator(stage_id, path):
    """Canonical producer-relative path, independent of run/worker directories."""
    parts=Path(logical_path(path)).parts
    if stage_id in parts:
        return '/'.join(parts[len(parts)-1-list(reversed(parts)).index(stage_id)+1:])
    if 'artifacts' in parts:
        return '/'.join(parts[len(parts)-1-list(reversed(parts)).index('artifacts')+1:])
    # Legacy files outside a producer directory have no known relative root.
    # Preserve their qualified location rather than inventing an ambiguous ID.
    return logical_path(path)


def records(item, run_id=None):
    files=list(dict.fromkeys(str(p) for p in item.get('artifacts', [])))
    previous={r['path']:r for r in item.get('artifact_records', [])}
    config=item.get('configuration', {})
    suffixes={path:['/'.join(Path(logical_path(path)).parts[-n:]) for n in range(1,len(Path(logical_path(path)).parts)+1)] for path in files}
    counts=Counter(suffix for options in suffixes.values() for suffix in options)
    result=[]
    for path in files:
        selector=next((suffix for suffix in suffixes[path] if counts[suffix]==1),None)
        if selector is None:raise ValueError('Os resultados contêm caminhos de arquivo duplicados.')
        identity=slot_id(item['id'], locator(item['id'], path))
        prior=previous.get(path, {})
        origin=prior.get('run_id') or item.get('reused_from_run') or run_id
        result.append(dict(artifact_id=identity,
            occurrence_id=hashlib.sha256(json.dumps([identity, origin]).encode()).hexdigest(),
            stage_id=item['id'], stage_name=item.get('name', item['id']),
            operation=item.get('operation', config.get('operation', 'unknown')),
            run_id=origin, name=Path(path).name, selector=selector, locator=locator(item['id'],path), path=path,
            context=prior.get('context', {'target':config.get('parameters', {}).get('target'),
                     'search_term':config.get('parameters', {}).get('search_term')}),
            inputs=prior.get('inputs', {})))
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
        matching=[p for p in files if slot_id(ref['stage'],locator(ref['stage'],p))==ref['artifact_id']]
        if len(matching)!=1:
            raise ValueError('A identidade do arquivo selecionado mudou. Revise a seleção da origem.')
    else:
        matching=[p for p in files if matches_selector(p, ref['selector'])]
    if len(matching)!=1 or not Path(matching[0]).is_file():
        raise ValueError('O arquivo escolhido não está nos resultados da etapa (ausente ou ambíguo). Revise a seleção da origem.')
    return matching[0]


def index_run(stages, run_id, assets=()):
    """Record exact file inputs, including per-batch lineage for individual mode."""
    assets_by_path={str(a['path']):a for a in assets}
    by_path={}
    for item in stages:
        item['artifact_records']=records(item, run_id)
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


def annotate_files(item, files, run_id=None):
    """Attach provenance to display records without replacing original filenames."""
    catalog={r['path']:r for r in records(item,run_id)}
    for file in files:
        record=catalog.get(file['path'])
        if not record:continue
        inputs=sorted({r.get('stage_name',r.get('name',''))+' · '+r.get('stage_id','')[:8]
            for rows in record['inputs'].values() for r in rows})
        context=record['context'].get('target') or record['context'].get('search_term')
        parts=[record['stage_name']+' · '+record['stage_id'][:8]]
        if context:parts.append(str(context))
        if inputs:parts.append(' + '.join(inputs))
        parts.append(record['selector'])
        file.update(artifact_id=record['artifact_id'],occurrence_id=record['occurrence_id'],
                    origin_run_id=record['run_id'],label=' — '.join(parts))
    return files
