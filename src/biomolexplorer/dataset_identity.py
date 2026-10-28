"""Dataset labels from persisted provenance, with conservative file-only fallback."""
import hashlib
from pathlib import Path
import re

TAGS = {'MOL','SIM','PUBCHEM','ZINC','MIX','INPUT'}


def source_label(step):
    target = str(step.get('target') or '').strip()
    if target and target not in ('MeuAlvo','Estruturas','PubChem','ZINC','Zinc'):return target
    name = str(step.get('stage_name') or '').strip()
    from .catalog import TITLES
    defaults = {value[0] for value in TITLES.values()}
    defaults |= {'Retrieve ChEMBL','Retrieve PDB structures','Retrieve structures','Retrieve PubChem','Retrieve ZINC'}
    if name and name not in defaults:return name
    return str(step.get('identifier') or '').strip()


def file_tag(path):
    path = Path(path)
    tokens = re.split(r'[^A-Za-z0-9]+',path.stem.upper())
    for tag in ('MIX','PUBCHEM','ZINC','SIM','MOL','INPUT'):
        if tag in tokens:return tag
    if any(p.casefold()=='pubchem' for p in path.parts):return 'PUBCHEM'
    if any(p.casefold()=='zinc' for p in path.parts):return 'ZINC'
    if path.stem.lower() in ('similars','similarmols','similar_molecules'):return 'SIM'
    if path.stem.lower() in ('molecules','mols'):return 'MOL'
    if 'ChEMBL' in path.parts:
        if 'similars' in path.parts:return 'SIM'
        if 'molecules' in path.parts:return 'MOL'
    return 'INPUT'


def admet_origin(paths, origins=None):
    tags=[];labels=[]
    for path in paths:
        routes = origins.describe(path).get('source_routes',[]) if origins else []
        if not routes:
            tags.append(file_tag(path));continue
        for route in routes:
            providers=[s for s in route if s.get('operation') in ('retrieve_pubchem','retrieve_zinc','retrieve_compounds','expand_similar_compounds')]
            provider=providers[-1] if providers else None
            if provider and provider['operation']=='retrieve_pubchem':tag='PUBCHEM'
            elif provider and provider['operation']=='retrieve_zinc':tag='ZINC'
            else:
                candidates=[file_tag(value) for s in reversed(route)
                    for value in (s.get('path',''),s.get('name',''))]
                tag=next((t for t in candidates if t!='INPUT'),file_tag(path))
                if tag=='INPUT' and provider and provider['operation'] in ('retrieve_compounds','expand_similar_compounds'):
                    tag='MIX'
            tags.append(tag)
            label=next((value for s in route if s.get('operation')!='external' and (value:=source_label(s))), '')
            if label and label not in labels:labels.append(label)
    unique=set(tags)
    return dict(tag=next(iter(unique)) if len(unique)==1 else 'MIX' if unique else 'INPUT',
                label=' + '.join(labels))


def output_stem(filename, origin=None):
    origin=origin or {}
    tag=origin.get('tag') or file_tag(filename)
    if tag not in TAGS:raise ValueError('Tag de conjunto ADMET inválida.')
    stem=Path(filename).stem
    label=str(origin.get('label') or '').strip()
    value=(label+'_' if label and label.casefold()!=stem.casefold() else '')+stem
    value=re.sub(r'_(MOL|SIM|PUBCHEM|ZINC|MIX|INPUT)$','',value,flags=re.I)
    normalized=re.sub(r'[^\w.+-]+','_',value,flags=re.UNICODE).strip('._') or 'dataset'
    if normalized!=value:normalized+='_'+hashlib.sha256(value.encode()).hexdigest()[:8]
    value=normalized
    if len(value)>120:value=value[:100]+'_'+hashlib.sha256(value.encode()).hexdigest()[:12]
    return value+'_'+tag


def redocking_targets(routes, pdb):
    """A merged collection must match the original PDB before choosing its target."""
    exact=[r for r in routes if any(Path(s.get('path','')).suffix.lower()=='.pdb'
        and Path(s.get('path','')).stem.upper()==pdb.upper() for s in r[:-1])]
    candidates=exact or routes
    labels=[]
    for route in candidates:
        sources=[s for s in route[:-1] if s.get('operation')=='retrieve_structures']
        for source in sources:
            value=source_label(source)
            if value and value not in labels:labels.append(value)
    # Metadata-only routes from several sources cannot safely identify a PDB.
    return sorted(labels,key=str.casefold) if exact or len(labels)<=1 else []
