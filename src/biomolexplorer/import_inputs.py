"""Per-file import types with compatibility for legacy single-type blocks."""
KINDS={'compounds','structures','prepared_structures','fingerprints','similarity',
       'vina','dock6','scores','visualization','other','zinc_urls'}


def file_kind(params,asset_id):
    return params.get('asset_types',{}).get(asset_id,params.get('kind','other'))


def validate_types(params):
    types=params.get('asset_types',{})
    if (not isinstance(types,dict) or types and set(types)!=set(params.get('asset_ids',[]))
            or any(not isinstance(kind,str) or kind not in KINDS for kind in types.values())
            or not isinstance(params.get('kind'),str) or params.get('kind') not in KINDS):
        raise ValueError('Informe um tipo válido para cada arquivo selecionado.')


def validate_import_files(paths,params):
    from .input_validation import validate_bundle,validate_file
    grouped={}
    for asset_id,path in paths:
        grouped.setdefault(file_kind(params,asset_id),[]).append(path)
    for kind,files in grouped.items():
        if kind=='visualization':
            for path in files:validate_file(path,kind)
        else:validate_bundle(files,kind)
