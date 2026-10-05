"""Shared input references for validation, execution, caching and the canvas."""

def sources(binding):
    return binding.get('sources',[]) if 'sources' in binding else [binding]


def dependencies(stage):
    if stage.get('provided_results'):
        return set()
    return set(stage.get('depends_on',[])) | {source['stage']
        for binding in stage.get('bindings',{}).values() for source in sources(binding) if 'stage' in source}


def asset_references(stage):
    if stage.get('provided_results'):
        return set(stage['provided_results'].get('asset_ids',[]))
    return set(stage['parameters'].get('asset_ids',[])) | {source['asset']
        for binding in stage.get('bindings',{}).values() for source in sources(binding) if 'asset' in source}


def pack(references):
    unique=[]
    for reference in references:
        if reference not in unique:
            unique.append(reference)
    return unique[0] if len(unique)==1 else {'sources':unique} if unique else None


def input_labels(stage,assets=None,automatic=None):
    names={a['id']:a['name'] for a in assets or []}
    if stage.get('provided_results'):
        return [names.get(identifier,'Resultado fornecido') for identifier in stage['provided_results']['asset_ids']]
    labels=[];seen=[]
    for binding in stage.get('bindings',{}).values():
        for source in sources(binding):
            if source in seen:continue
            seen.append(source)
            if 'asset' in source:
                labels.append(names.get(source['asset'],'Arquivo próprio'))
            elif source.get('selector','auto')!='auto':
                labels.append(source['selector'])
            else:
                labels.append((automatic or {}).get(source['stage'],'Seleção automática'))
    return labels
