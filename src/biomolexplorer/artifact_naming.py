"""Publication boundary for short output names, separate from engine work files."""
import csv
import hashlib
import json
from pathlib import Path
import re
from functools import lru_cache

OPERATIONS={'retrieve_compounds':'chembl','retrieve_pubchem':'pubchem',
    'expand_similar_compounds':'chembl-similar','retrieve_zinc':'zinc','admet':'admet',
    'fingerprints':'fingerprint','similarity':'similarity','graphs':'graph',
    'import_results':'import','prepare_structures':'preparation','redocking':'redocking',
    'docking_vina':'vina','docking_dock6':'dock6','consensus':'consensus','dynamics_dock6':'md'}
PATTERN=re.compile(r'^(B\d+)__([a-z0-9-]+)__(.+)$')


def assign_block_ids(stages, previous=()):
    """Allocate monotonically, preserving IDs through reordering and deletion."""
    known={s['id']:s['block_id'] for s in previous if s.get('block_id')}
    used=set(known.values())|{s['block_id'] for s in stages if s.get('block_id')}
    number=max([int(v[1:]) for v in used if re.fullmatch(r'B\d+',v)]+[0])+1
    assigned=set()
    for stage in stages:
        value=known.get(stage['id'],stage.get('block_id'))
        if value is not None and (not re.fullmatch(r'B\d{2,}',value) or value in assigned):
            raise ValueError('Cada bloco deve ter um block_id curto e único (B01, B02…).')
        if value is None:
            while f'B{number:02d}' in used:number+=1
            value=f'B{number:02d}';number+=1;used.add(value)
        stage['block_id']=value;assigned.add(value)
    return stages


def slug(value):
    return re.sub('[^a-z0-9-]+','-',str(value).lower()).strip('-')[:48] or 'result'


def output_role(operation,path,parameters):
    path=Path(path);name=path.stem.lower()
    canonical=PATTERN.match(path.stem)
    if canonical and canonical[2]==OPERATIONS.get(operation) and re.fullmatch('[a-z0-9-]+',canonical[3]):
        return canonical[3]
    if operation in ('retrieve_compounds','expand_similar_compounds','retrieve_pubchem'):
        if name.endswith('_mols'):return 'mol'
        if name.endswith('_sims'):return 'sim'
        if name.endswith('_full'):return 'merged'
    if operation=='fingerprints':
        return next((v for v in ('morgan','maccs','pharmacophore') if v in name), 'fingerprints')
    if operation=='similarity':return slug(parameters.get('metric','tanimoto'))
    if operation=='admet':
        if name.endswith('_alerts'):return 'alerts'
        rule=re.search(r'_rules_(lipinski-veber|lipinski|veber|fragment)$',name)
        if rule:return 'rules-'+rule[1]
        if name.endswith('_admet_report'):return 'report'
        for key,label in (('bbb+','bbb-pos'),('bbb-','bbb-neg'),('hia+','hia-pos'),('hia-','hia-neg'),
                          ('evaluated','evaluated'),('rejected','rejected'),('selected','selected')):
            if name.endswith('_'+key):return label
        return 'compounds'
    if operation=='graphs':
        if 'full_graph' in name:return 'full-graph'
        match=re.search(r'component_(\d+)',name)
        if match:return ('lcc-' if name.endswith('_mcc') else 'component-')+match[1]
        if path.parent.name=='Molecules':return 'lcc'
        return 'edges'
    roles={'retrieve_structures':{'.pdb':'structure','.cif':'structure','.mmcif':'structure'},
        'prepare_structures':{'.pdbqt':'receptor','.mol2':'ligand'},
        'docking_vina':{'.pdbqt':'pose'},'docking_dock6':{'.mol2':'pose'},
        'dynamics_dock6':{'.nc':'trajectory','.prmtop':'topology','.rst7':'coordinates'}}
    if path.name=='docking_results.csv':return 'scores'
    if path.name=='pdb_codes.csv':return 'rmsd' if operation=='redocking' else 'complexes'
    if path.name=='centers.csv':return 'centers'
    if operation=='consensus' and path.suffix=='.csv':return 'ranking'
    if operation=='retrieve_zinc' and path.suffix in ('.sdf','.smi','.mol2'):return 'mol'
    if path.suffix in roles.get(operation,{}):return roles[operation][path.suffix]
    if path.suffix=='.csv' and 'compound' in name:return 'compounds'
    return slug(name)


def canonical_filename(block_id,operation,output,extension):
    if not re.fullmatch(r'B\d{2,}',block_id):raise ValueError('block_id inválido.')
    extension=extension.lstrip('.')
    if not re.fullmatch(r'[a-zA-Z0-9]+(?:[.-][a-zA-Z0-9]+)*',extension):raise ValueError('Extensão de artefato inválida.')
    return f'{block_id}__{OPERATIONS.get(operation,slug(operation))}__{slug(output)}.{extension}'


def publication_metadata(path):
    """Legacy engine names remain metadata, never the public display label."""
    path=Path(path);manifest=path.parent/'.artifact-publication.json'
    stamp=_file_stamp(manifest)
    if stamp is None:return {}
    return dict(_publication_index(str(manifest),stamp).get(path.name,{}))


def _publish_legacy_analytical(stage,files,owned_root):
    """Rename materialized analytical tables once; references are never copied.

    Engine-coupled structure bundles keep their required internal filenames.
    Analytical CSVs and paired visualization files have independent public
    names. Structured references are updated without retaining duplicate data
    or filesystem aliases (project archives require regular files).
    """
    operation=stage['operation'];block=stage.get('block_id')
    analytical={'retrieve_compounds','retrieve_pubchem','expand_similar_compounds',
                'admet','fingerprints','similarity','graphs'}
    if not block or operation not in analytical or stage.get('provided_results'):return files
    owned_root=Path(owned_root).resolve();parameters=stage.get('parameters',{})
    def published(path):
        match=PATTERN.match(path.stem)
        return bool(match and match[1]==block and match[2]==OPERATIONS.get(operation)
            and re.fullmatch('[a-z0-9-]+',match[3]))
    mapping={};metadata={};occupied=set()
    for value in sorted(files):
        path=Path(value)
        if not path.is_relative_to(owned_root) or not path.is_file() or path.is_symlink():continue
        if path.suffix.lower()!='.csv' or published(path):continue
        if path.name in ('molecule_exclusions.csv','preparation_failures.csv'):continue
        role=output_role(operation,path,parameters)
        name=canonical_filename(block,operation,role,path.suffix)
        if (path.parent,name) in occupied:
            role+='-'+hashlib.sha256(path.relative_to(owned_root).as_posix().encode()).hexdigest()[:8]
            name=canonical_filename(block,operation,role,path.suffix)
        occupied.add((path.parent,name));destination=path.with_name(name)
        if destination.exists() and destination!=path:raise ValueError('Artifact publication collision: '+name)
        mapping[path]=destination
        for suffix in ('.biomol-view.json','.png'):
            companion=path.with_suffix(suffix)
            if str(companion) in files and companion.exists():mapping[companion]=destination.with_suffix(suffix)
        with path.open(encoding='utf-8-sig',newline='') as stream:
            reader=csv.reader(stream);next(reader,None);count=sum(1 for _ in reader)
        from .artifacts import slot_id,locator
        metadata[destination]=dict(output=role,record_count=count,legacy_name=path.name,
            artifact_id=slot_id(stage['id'],locator(stage['id'],destination)))
    # Reports and visualizations are public outputs too. Their roles derive
    # from the operation and file type, rather than an inherited input stem.
    for value in sorted(files):
        path=Path(value)
        if path in mapping or not path.is_relative_to(owned_root) or not path.is_file() or path.is_symlink():continue
        if published(path) or path.name in ('molecule_exclusions.json','preparation_failures.csv'):continue
        if path.suffix.lower() not in ('.json','.png','.svg','.html','.csv','.sdf','.mol2','.smi'):continue
        extension='biomol-view.json' if path.name.endswith('.biomol-view.json') else path.suffix.lstrip('.')
        role='view' if extension=='biomol-view.json' else 'report' if extension=='json' else 'figure' if extension in ('png','svg','html') else output_role(operation,path,parameters)
        destination=path.with_name(canonical_filename(block,operation,role,extension))
        if destination in mapping.values() or destination.exists() and destination!=path:
            role+='-'+hashlib.sha256(path.relative_to(owned_root).as_posix().encode()).hexdigest()[:8]
            destination=path.with_name(canonical_filename(block,operation,role,extension))
        mapping[path]=destination
    from .artifacts import slot_id,locator
    for source,destination in mapping.items():
        if destination not in metadata:
            metadata[destination]=dict(output='view' if destination.name.endswith('.biomol-view.json') else
                'figure' if destination.suffix in ('.png','.svg') else 'report',
                legacy_name=source.name,artifact_id=slot_id(stage['id'],locator(stage['id'],destination)))
        source.rename(destination)
    replacements={str(source):str(destination) for source,destination in mapping.items()}
    names={source.name:destination.name for source,destination in mapping.items()
        if sum(other.name==source.name for other in mapping)==1}
    def remap(value):
        if isinstance(value,str):
            if value in replacements:return replacements[value]
            if value in names:return names[value]
            # Views may store project-relative exports (Molecules/xxx.csv).
            matches=[destination for source,destination in mapping.items()
                if len(Path(value).parts)>1 and source.as_posix().endswith('/'+value)]
            if len(matches)==1:return str(Path(value).with_name(matches[0].name))
            return value
        if isinstance(value,list):return [remap(v) for v in value]
        if isinstance(value,dict):return {k:remap(v) for k,v in value.items()}
        return value
    from .storage import write_json
    for value in files:
        path=mapping.get(Path(value),Path(value))
        if not path.is_relative_to(owned_root):continue
        if path.suffix=='.json':
            try:
                original=json.loads(path.read_text());updated=remap(original)
                if original!=updated:write_json(updated,path)
            except (OSError,ValueError):pass
    for folder in {p.parent for p in metadata}:
        from .storage import write_json
        manifest=folder/'.artifact-publication.json'
        saved=json.loads(manifest.read_text()) if manifest.exists() else {}
        saved.update({p.name:v for p,v in metadata.items() if p.parent==folder})
        write_json(saved,manifest)
    return [str(mapping.get(Path(p),Path(p))) for p in files]


def legacy_human_name(record):
    output=record.get('output','result')
    legacy=record.get('provenance',{}).get('legacy_name','')
    if output=='compounds' and legacy.lower().endswith(('_mols.csv','_sims.csv','_full.csv')):
        output=output_role(record.get('operation'),legacy,{})
    if output=='compounds':
        if record.get('operation') in ('retrieve_compounds','retrieve_zinc'):output='mol'
        elif record.get('operation')=='retrieve_pubchem':output='sim'
    labels={'morgan':'Morgan fingerprints','maccs':'MACCS fingerprints','pharmacophore':'Pharmacophore fingerprints',
        'mol':'MOL · Moléculas','sim':'SIM · Similares','merged':'ALL · Conjunto completo',
        'compounds':'Compostos','bbb-pos':'BBB+','bbb-neg':'BBB−','hia-pos':'HIA+','hia-neg':'HIA−',
        'alerts':'Alertas químicos','rules-lipinski':'Regras · Lipinski','rules-veber':'Regras · Veber',
        'rules-lipinski-veber':'Regras · Lipinski + Veber','rules-fragment':'Regras · Fragmentos',
        'evaluated':'ADMET · Todos os avaliados','selected':'ADMET · Selecionados','rejected':'ADMET · Fora dos critérios',
        'structure':'PDB · Estrutura','receptor':'REC · Receptor preparado','ligand':'LIG · Ligante preparado',
        'pose':'POSE · Pose de docking','scores':'SCORE · Resultados de docking','ranking':'RANK · Consenso',
        'complexes':'Complexos e ligantes','rmsd':'RMSD · Redocking','centers':'Centros de docking','trajectory':'MD · Trajetória',
        'topology':'MD · Topologia','coordinates':'MD · Coordenadas','report':'Relatório','figure':'Figura','view':'Visualização',
        'tanimoto':'Similaridade Tanimoto','dice':'Similaridade Dice','cosine':'Similaridade cosseno',
        'lcc':'MCC / LCC','full-graph':'Grafo completo','edges':'Relações do grafo'}
    label=labels.get(output,record.get('name') or output.replace('-',' ').capitalize())
    if re.fullmatch(r'(lcc|component)-\d+',output):
        label=('MCC / LCC · ' if output.startswith('lcc-') else 'Componente · ')+output.rsplit('-',1)[1]
    if output in ('structure','receptor','ligand','pose','trajectory','topology','coordinates'):
        label+=' · '+record.get('name',Path(record.get('path','')).name)
    if record.get('operation')=='admet' and output=='compounds':label='Compostos processados'
    batch=record.get('batch_label')
    if batch and record.get('processing_scope')!='all':
        match=PATTERN.match(batch)
        label+=' · '+(labels.get(match[3],match[3]) if match else batch)
    return label+' · ALL' if record.get('processing_scope')=='all' else label


# Public names are presentation metadata. Native bundles and their logical slot
# identities retain their validated filenames, including legacy publications.
SCHEME_VERSION=1
PUBLIC_OPERATIONS={**OPERATIONS,'retrieve_structures':'PDB',
    'retrieve_compounds':'ChEMBL','expand_similar_compounds':'ChEMBL-Similar',
    'retrieve_pubchem':'PubChem','retrieve_zinc':'ZINC','admet':'ADMET',
    'fingerprints':'Fingerprints','similarity':'Similarity','graphs':'Graphs',
    'prepare_structures':'Preparation','redocking':'Redocking',
    'docking_vina':'Docking-Vina','docking_dock6':'Docking-DOCK6',
    'consensus':'Consensus','dynamics_dock6':'AmberScore-DOCK6','import_results':'Import'}
PUBLIC_PATTERN=re.compile(r'^TG\[([^\[\]/\\]+)\]-SRC\[([^\[\]/\\]+)\]-OP\[([^\[\]/\\]+)\]-(.+?)(?:-ID\[([a-f0-9]{8,64})\])?\.([A-Za-z0-9][A-Za-z0-9-]*(?:\.[A-Za-z0-9-]+)*)$')
GENERIC_TARGETS={'','meualvo','estruturas','structures','pubchem','zinc','chembl',
    'molecules','moleculas','compounds','fingerprints','similarity','admet','target','dados','results','na','none','unknown','undefined','input'}


def component(value,limit=32):
    """Safe, bounded labels; complete original labels remain in metadata."""
    import unicodedata
    original=str(value or '').strip().replace('−','-').replace('–','-')
    text=unicodedata.normalize('NFKD',original).encode('ascii','ignore').decode()
    text=re.sub(r'[^A-Za-z0-9+_.-]+','-',text).strip('.')
    if text and not re.search('[A-Za-z0-9]',text):text=''
    if not text:text='NA' if not original else 'U-'+hashlib.sha256(original.encode()).hexdigest()[:8]
    if len(text)>limit:text=text[:limit-9].rstrip('-')+'-'+hashlib.sha256(original.encode()).hexdigest()[:8]
    return text


def file_extension(filename):
    name=Path(filename).name
    if name.lower().endswith('.biomol-view.json'):return 'biomol-view.json'
    path=Path(name)
    if path.suffix.lower() in ('.gz','.bz2','.xz','.zst') and len(path.suffixes)>1:return ''.join(path.suffixes[-2:]).lstrip('.')
    return path.suffix.lstrip('.') or 'dat'


@lru_cache(maxsize=512)
def _publication_index(path,stamp):
    try:
        value=json.loads(Path(path).read_text())
        return value if isinstance(value,dict) else {}
    except (OSError,ValueError):return {}

def public_filename(target,source,operation,result_type,extension,identifier=None):
    extension=str(extension).lstrip('.')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9-]*(?:\.[A-Za-z0-9-]+)*',extension):raise ValueError('Extensão de artefato inválida.')
    tail='-ID['+identifier+']' if identifier else ''
    if identifier and not re.fullmatch(r'[a-f0-9]{8,64}',identifier):raise ValueError('Identificador de artefato inválido.')
    return f'TG[{component(target,24)}]-SRC[{component(source,24)}]-OP[{component(operation,24)}]-{component(result_type,40).replace('.', '-')}{tail}.{extension}'


def parse_filename(filename):
    """Parse public and Bxx compatibility names, never infer an artifact ID."""
    name=Path(filename).name;match=PUBLIC_PATTERN.fullmatch(name)
    if match:
        return dict(scheme=SCHEME_VERSION,target=match[1],source=match[2],operation=match[3],
            result_type=match[4],identifier=match[5],extension=match[6])
    legacy=PATTERN.match(Path(name).stem)
    if legacy:return dict(scheme=0,block_id=legacy[1],operation=legacy[2],result_type=legacy[3],extension=file_extension(name))
    return None


def _fields(path):
    # One shared header cache, invalidated on any normal file modification.
    return _cached_fields(str(path),_file_stamp(path))


def _file_stamp(path):
    try:
        info=Path(path).stat();return info.st_size,info.st_mtime_ns,info.st_ctime_ns
    except OSError:return None


@lru_cache(maxsize=2048)
def _cached_fields(path,stamp):
    if stamp is None or Path(path).suffix.lower()!='.csv':return frozenset()
    try:
        with Path(path).open(encoding='utf-8-sig',newline='') as stream:return frozenset(next(csv.reader(stream),[]))
    except (OSError,UnicodeError,csv.Error):return frozenset()


def result_type(record,config=None):
    """One content taxonomy for all producers, including old native bundles."""
    config=config or {};params=config.get('parameters',{})
    saved=record.get('naming') or {}
    if saved.get('scheme')==SCHEME_VERSION and saved.get('result_type') and saved.get('extension')==file_extension(record.get('name',record.get('path',''))):
        return saved['result_type']
    operation=record.get('operation',config.get('operation',''))
    original=record.get('original_name') or record.get('provenance',{}).get('legacy_name') or record.get('name') or Path(record.get('path','')).name
    parsed=parse_filename(original)
    if parsed and parsed['scheme']==SCHEME_VERSION:return parsed['result_type']
    path=Path(original)
    if path.suffix.lower() in ('.gz','.bz2','.xz','.zst'):path=path.with_suffix('')
    name=path.name.lower();suffix=path.suffix.lower()
    role=record.get('output') or output_role(operation,path,params)
    roles={'mol':'MOL','sim':'SIM','merged':'ALL','compounds':'MOL','alerts':'Alerts',
        'bbb-pos':'BBB+','bbb-neg':'BBB-','hia-pos':'HIA+','hia-neg':'HIA-',
        'evaluated':'Evaluated','selected':'Selected','rejected':'Rejected',
        'morgan':'Fingerprints-Morgan','maccs':'Fingerprints-MACCS','pharmacophore':'Fingerprints-Pharmacophore',
        'scores':'Scores','ranking':'Scores','centers':'Centers','trajectory':'Trajectory',
        'topology':'Topology','coordinates':'Coordinates','rmsd':'RMSD','complexes':'Pairs',
        'full-graph':'Graph-Full','lcc':'Molecules-LCC','edges':'Edges'}
    # Reports, grids and references have actual content types, not block names.
    if 'exclusions' in name:return 'Exclusions'
    if 'failures' in name:return 'Failures'
    if name in ('centers.csv','centers.json'):return 'Centers'
    if 'energies' in name or 'amber_score' in name:return 'AmberScore'
    if 'rmsf' in name:return 'RMSF'
    if 'rmsd' in name and suffix in ('.csv','.json','.png','.svg'):return 'RMSD'
    if 'gyration' in name:return 'RadiusOfGyration'
    if name=='protocol.json':return 'Protocol'
    if name.endswith('_admet_report.json'):return 'Methods'
    if name=='dynamics_results.json':return 'Pairs'
    if name=='analysis.json':return 'Analysis'
    if suffix in ('.log','.out','.err'):return 'Log'
    if suffix in ('.prmtop','.top'):return 'Topology'
    if suffix in ('.nc','.dcd','.xtc','.trr'):return 'Trajectory'
    if suffix in ('.rst7','.inpcrd','.rst'):return 'Coordinates'
    if suffix in ('.bmp','.nrg'):return 'Grid'
    if suffix=='.sph':return 'Spheres'
    if operation=='retrieve_zinc' and suffix in ('.sdf','.smi','.mol2'):return 'MOL'
    if suffix in ('.pdb','.cif','.mmcif','.pdbqt','.mol2','.sdf','.mol'):
        if operation in ('redocking','prepare_structures') and suffix in ('.pdb','.cif','.mmcif') and '.lig.' not in name and '.noh.' not in name:
            return 'Complex-'+component(path.stem.removesuffix('.complex'),24)
        if '.dockprep.' in name or '.noh.' in name:return 'Receptor-'+component(path.name.split('.',1)[0],24)
        if '.lig.' in name or '.ligand.' in name:
            parent=Path(record.get('path','')).parent.name
            prefix='ReferenceLigand-' if operation=='redocking' and parent=='Prepared' else 'Ligand-' if parent in ('Prepared','Molecules') else 'Pose-' if operation.startswith('docking_') or operation=='redocking' else 'Ligand-'
            return prefix+component(path.name.split('.',1)[0],28)
        if operation=='retrieve_structures':return 'Complex-'+component(path.stem,20)
        if operation in ('docking_vina','docking_dock6','redocking'):return 'Pose-'+component(path.stem.removesuffix('_scored'),28)
        return 'Ligand' if suffix in ('.sdf','.mol','.mol2') else 'Structure'
    if role.startswith('rules-'):return 'Rules-'+role[6:].replace('lipinski','Lipinski').replace('veber','Veber').replace('fragment','Fragment')
    if operation=='graphs':
        component_role=re.fullmatch(r'(component|lcc)-(\d+)',role)
        if component_role:return 'Molecules-C'+str(int(component_role[2]))+('-LCC' if component_role[1]=='lcc' else '')
    if operation=='retrieve_pubchem' and role in ('compounds','mol'):role='sim'
    if operation=='similarity' and suffix=='.csv':return 'Similarity-'+component(str(params.get('metric') or role).capitalize(),20)
    if suffix=='.csv':
        fields=_fields(record.get('path',''))
        if 'RMSD' in fields:return 'RMSD'
        if {'molecule_chembl_id','engine','score'}<=fields or {'molecule_chembl_id','vina','dock6'}<=fields:return 'Scores'
        if 'PDB_CODE' in fields and 'LIGAND' in fields:return 'Pairs'
        if operation=='dynamics_dock6' and ('amber_score' in fields or 'Amber_Score' in fields):return 'AmberScore'
        if role in roles:return roles[role]
        if 'fingerprint' in fields:return 'Fingerprints'
        if {'source','target','value'}<=fields:return 'Similarity'
        if {'molecule_chembl_id','canonical_smiles'}<=fields:return 'MOL'
        return 'Table'
    if suffix in ('.png','.svg','.pdf','.html') or name.endswith('.biomol-view.json'):
        base=roles.get(role)
        if base is None:
            # Paired visualizations retain the data category when recoverable.
            base_role=output_role(operation,path.with_suffix('.csv'),params)
            base=roles.get(base_role)
        prefix='View' if name.endswith('.biomol-view.json') else 'Figure'
        return prefix+('-'+base if base and base not in ('MOL','Table') else '')
    if suffix=='.json':return 'Report'
    if suffix in ('.txt','.in','.def','.template','.com'):return 'Parameters' if suffix!='.txt' else 'Text'
    return roles.get(role,'Data')


def _target_label(value):
    value=str(value or '').strip()
    value=re.sub(r'^(?:Target|Alvo)\s+|\s+(?:PDB|ChEMBL)$','',value,flags=re.I).strip()
    if component(value).casefold() in GENERIC_TARGETS:return None
    gene=re.search(r'(?<![A-Z0-9])(ACHE|BCHE|BACE1)(?![A-Z0-9])',value.upper())
    return gene[1] if gene else component(value,24) if value else None


def name_records(entries):
    """Enrich (record,stage) pairs in memory, preserving exact producer IDs.

    Immediate parents choose SRC; only their target metadata is propagated.
    Filenames are a compatibility fallback, never proof of identity/validity.
    """
    entries=list(entries);by_path={};by_id={}
    for record,stage in entries:
        by_path.setdefault(record['path'],[]).append((record,stage));by_id.setdefault(record.get('artifact_id'),[]).append((record,stage))
    memo={};visiting=set();kinds={}
    def kind_of(record,stage):
        key=(record.get('stage_id'),record['path'])
        if key not in kinds:kinds[key]=result_type(record,stage.get('configuration',stage))
        return kinds[key]
    def describe(record,stage):
        key=(record.get('stage_id'),record['path'])
        if key in memo:return memo[key]
        if key in visiting:return dict(targets=[],source='INPUT')
        visiting.add(key)
        config=stage.get('configuration',stage);op=record.get('operation',config.get('operation'));kind=kind_of(record,stage)
        batch=next((b for b in stage.get('batches',[]) if record['path'] in b.get('artifacts',[])),{})
        inputs=record.get('inputs') or {field:[dict(path=str(p)) for p in paths] for field,paths in (batch or stage).get('input_files',{}).items()}
        port='base_selected_mols' if op in ('prepare_structures','docking_vina','docking_dock6','dynamics_dock6') and not kind.startswith(('Receptor','Complex','Structure','Centers')) else 'base_vina_path' if op=='consensus' else 'similarity_path' if op=='graphs' else 'base_input_path'
        def parents(field):
            result=[]
            for ref in inputs.get(field,[]):
                choices=by_path.get(ref.get('path'),[]) or by_id.get(ref.get('artifact_id'),[])
                if ref.get('stage_id'):choices=[entry for entry in choices if entry[0].get('stage_id')==ref['stage_id']]
                if len(choices)==1:result.append(choices[0])
            return result
        immediate=parents(port)
        target_parents=parents('base_input_path') if op in ('docking_vina','docking_dock6','redocking','dynamics_dock6') and not (op.startswith('docking_') and kind.startswith('Ligand')) else immediate
        # Native per-pair protocol metadata disambiguates merged receptor sets.
        if op=='dynamics_dock6' and 'pairs' in Path(record['path']).parts:
            parts=Path(record['path']).parts;position=parts.index('pairs')
            pair_root=Path(*parts[:position+2]);protocol=pair_root/'protocol.json'
            details=_publication_index(str(protocol),_file_stamp(protocol)) if _file_stamp(protocol) else {}
            pair_inputs=details.get('inputs',{})
            receptor=pair_inputs.get('receptor',{}).get('source')
            ligand=pair_inputs.get('ligand',{}).get('source')
            if receptor:target_parents=[entry for entry in target_parents if entry[0]['path']==receptor]
            else:target_parents=[]
            if ligand:immediate=[entry for entry in immediate if entry[0]['path']==ligand]
        if op in ('redocking','docking_vina','docking_dock6','prepare_structures'):
            match=re.match(r'^([1-9][A-Za-z0-9]{3})(?:[_.-]|$)',record.get('original_name',record.get('name','')))
            if match:
                pdb=match[1].upper()
                exact=[entry for entry in target_parents if Path(entry[0]['path']).name.upper().startswith((pdb+'.',pdb+'_'))]
                if exact:target_parents=exact
        targets=sorted({t for r,s in target_parents for t in describe(r,s).get('targets',[])})
        context=record.get('context') or {};params=config.get('parameters',{})
        explicit=context.get('biological_target') or context.get('target_id') or (config.get('provided_results') or {}).get('target')
        if explicit and _target_label(explicit):targets=[_target_label(explicit)]
        biological_source=op=='retrieve_structures' or op=='retrieve_compounds' and (params.get('search_mode') or 'target').startswith(('target','uniprot'))
        if not targets and biological_source:
            candidate=_target_label(context.get('target') or params.get('target'))
            if not candidate:
                gene=re.search(r'(?<![A-Z0-9])(ACHE|BCHE|BACE1)(?![A-Z0-9])',str(record.get('stage_name',stage.get('name',''))).upper())
                if gene:candidate=gene[1]
            if not candidate:
                from .dataset_identity import source_label
                label=source_label(dict(target=context.get('target') or params.get('target'),stage_name=record.get('stage_name',stage.get('name','')),identifier=params.get('search_term','')))
                candidate=_target_label(label)
            if candidate:targets=[candidate]
        previous=record.get('naming',{})
        if not targets and previous.get('scheme')==SCHEME_VERSION:targets=previous.get('targets',[])
        source_kinds=sorted({kind_of(r,s) for r,s in immediate})
        # A structure bundle includes bookkeeping tables; those are supporting
        # metadata rather than additional biological input categories.
        categories=sorted({'REC' if k.startswith('Receptor') else 'PDB' if k.startswith(('Complex','Structure')) else 'LIG' if k.startswith(('Ligand','Pose','ReferenceLigand')) else k for k in source_kinds})
        if port=='base_input_path' and op in ('redocking','prepare_structures','docking_vina','docking_dock6','dynamics_dock6'):
            structures=[k for k in categories if k in ('PDB','REC')]
            if structures:categories=structures
        source=categories[0] if len(categories)==1 else 'MIX-'+hashlib.sha256(json.dumps(categories).encode()).hexdigest()[:8] if categories else None
        if source is None:
            from .dataset_identity import file_tag
            tags={file_tag(ref.get('path',ref.get('name',''))) for ref in inputs.get(port,[])}-{'INPUT'}
            source='+'.join(sorted(tags)) if tags else previous.get('source') if previous.get('scheme')==SCHEME_VERSION else None
        if not source:source={'retrieve_compounds':'ChEMBL','expand_similar_compounds':'ChEMBL','retrieve_structures':'PDB','retrieve_pubchem':'PubChem','retrieve_zinc':'ZINC'}.get(op,'INPUT')
        # Structure identifiers belong to TIPO; SRC is their immediate category.
        source='REC' if source.startswith('Receptor') else 'PDB' if source.startswith(('Complex','Structure')) else 'LIG' if source.startswith(('Ligand','Pose')) else source
        if len(source)>24:source='MIX-'+hashlib.sha256(source.encode()).hexdigest()[:8]
        operation=PUBLIC_OPERATIONS.get(op,component(op))
        target='+'.join(targets) if targets else 'NA'
        if len(target)>24:target='MULTI-'+hashlib.sha256(target.encode()).hexdigest()[:8]
        identity=record.get('occurrence_id') or record.get('artifact_id') or hashlib.sha256(record['path'].encode()).hexdigest()
        if not isinstance(identity,str) or not re.fullmatch('[a-f0-9]{8,64}',identity):identity=hashlib.sha256(str(identity).encode()).hexdigest()
        naming=dict(scheme=SCHEME_VERSION,target=target,targets=targets,source=source,operation=operation,
            result_type=kind,extension=file_extension(record.get('name',record['path'])),identifier=identity,
            source_artifact_ids=sorted({r['artifact_id'] for r,s in immediate if r.get('artifact_id')}),source_types=source_kinds)
        visiting.remove(key);memo[key]=naming;return naming
    identifiers={}
    for record,stage in entries:
        naming=dict(describe(record,stage));identifier=naming['identifier'];short=identifier[:8]
        identifiers.setdefault(short,set()).add(identifier)
        record['naming']=naming
    for record,stage in entries:
        naming=record['naming'];identifier=naming['identifier'];length=8
        while len({v[:length] for v in identifiers[identifier[:8]]})<len(identifiers[identifier[:8]]):length+=1
        record['public_name']=public_filename(naming['target'],naming['source'],naming['operation'],naming['result_type'],naming['extension'],identifier[:length])
        record['display_name']=record['public_name']
    return [r for r,s in entries]


def human_name(record):
    if record.get('public_name'):return record['public_name']
    item=dict(record);item.setdefault('path',item.get('name','artifact.dat'));item.setdefault('name',Path(item['path']).name)
    name_records([(item,dict(operation=item.get('operation'),parameters={}))])
    return item['public_name']


def publish(stage,files,owned_root,*,stages=None,run_id=None,input_files=None):
    """Publish contextual names for every block without changing bundle bytes.

    Existing analytical work names remain a private compatibility protocol.
    Public names are stable metadata aliases, including for native structures.
    """
    files=_publish_legacy_analytical(stage,files,owned_root)
    import copy
    from .artifacts import index_run,records
    view=copy.deepcopy(stages or [])
    current=dict(stage,configuration=copy.deepcopy(stage),artifacts=list(files),
        input_files=copy.deepcopy(input_files or {}),status='succeeded')
    view=[s for s in view if s['id']!=stage['id']]+[current]
    index_run(view,run_id)
    metadata={};root=Path(owned_root).resolve()
    for record in current['artifact_records']:
        path=Path(record['path'])
        if not path.is_relative_to(root) or not path.is_file() or path.is_symlink():continue
        old=publication_metadata(path)
        metadata[path]=dict(old,artifact_id=record['artifact_id'],output=record.get('output') or output_role(stage['operation'],path,stage.get('parameters',{})),
            legacy_name=old.get('legacy_name',path.name),public_name=record['public_name'],naming=record['naming'])
    from .storage import write_json
    for folder in {p.parent for p in metadata}:
        manifest=folder/'.artifact-publication.json'
        saved=json.loads(manifest.read_text()) if manifest.exists() else {}
        saved.update({p.name:data for p,data in metadata.items() if p.parent==folder})
        write_json(saved,manifest)
    return files


def download_name(context,*,result_type=None,extension=None,qualifier=None,operation=None):
    """One export convention for original files and generated dashboard views."""
    context=context or {};naming=dict(context.get('naming') or {})
    existing=context.get('public_name') or context.get('display_name')
    parsed=parse_filename(existing or '')
    if not naming and parsed and parsed['scheme']==SCHEME_VERSION:naming=parsed
    if not naming:
        record=dict(context);record.setdefault('operation',operation or 'export')
        record.setdefault('path',record.get('name','result.'+(extension or 'dat')))
        record.setdefault('name',Path(record['path']).name)
        name_records([(record,dict(operation=record['operation'],parameters={}))]);naming=record['naming'];existing=record['public_name']
    kind=result_type or naming['result_type'];ext=extension or naming['extension']
    if existing and kind==naming['result_type'] and ext==naming['extension'] and qualifier in (None,'all'):return existing
    identifier=hashlib.sha256(json.dumps([naming.get('identifier'),kind,ext,qualifier],sort_keys=True,default=str).encode()).hexdigest()[:8]
    return public_filename(naming.get('target','NA'),naming.get('source','INPUT'),naming.get('operation') or PUBLIC_OPERATIONS.get(operation,'Export'),kind,ext,identifier)



def aggregate_download_name(operation,result_type,extension,*,targets=(),source='INPUT',identifiers=(),qualifier=None):
    labels=sorted({label for target in targets if (label:=_target_label(target))})
    target='+'.join(labels) or 'NA'
    if len(target)>24:target='MULTI-'+hashlib.sha256(target.encode()).hexdigest()[:8]
    identity=hashlib.sha256(json.dumps([operation,result_type,extension,sorted(identifiers),qualifier],sort_keys=True).encode()).hexdigest()[:8]
    return public_filename(target,source,PUBLIC_OPERATIONS.get(operation,component(operation)),result_type,extension,identity)

def bundle_manifest(records,stages=()):
    """Portable public/private mapping; science bundle contents stay untouched."""
    return dict(scheme=SCHEME_VERSION,convention='TG[TARGET]-SRC[SOURCE]-OP[OPERATION]-TYPE-ID[IDENTIFIER].ext',
        producers={s['id']:s.get('configuration',s) for s in stages},
        artifacts=[{k:r.get(k) for k in ('artifact_id','occurrence_id','public_name','original_name','name','path','archive_name',
            'content_sha256','processing_version','producer_software','naming','parent_artifact_ids','inputs')} for r in records])


def native_name(path):
    """Compatibility hint for canonical exports; never a validity decision."""
    parsed=parse_filename(path)
    if not parsed or parsed['scheme']!=SCHEME_VERSION:return Path(path).name
    metadata=publication_metadata(path)
    if metadata.get('legacy_name') and not parse_filename(metadata['legacy_name']):return metadata['legacy_name']
    kind=parsed['result_type'];extension=parsed['extension']
    fixed={'Centers':'centers.csv','RMSD':'pdb_codes.csv','Pairs':'pdb_codes.csv','Scores':'docking_results.csv'}
    if kind in fixed and extension=='csv':return fixed[kind]
    if kind.startswith('Receptor-'):
        stem=kind[len('Receptor-'):]
        if re.fullmatch(r'[1-9][A-Za-z0-9]{3}_[A-Za-z0-9]+',stem):return stem+('.noH.' if extension.lower()=='pdb' else '.dockprep.')+extension
    if kind.startswith('Complex-') and re.fullmatch(r'[1-9][A-Za-z0-9]{3}',kind[8:]):return kind[8:]+'.'+extension
    for prefix in ('Pose-','Ligand-','ReferenceLigand-'):
        if kind.startswith(prefix):
            stem=kind[len(prefix):]
            return stem+('_scored.' if prefix=='Pose-' and 'DOCK6' in parsed['operation'] else '.lig.')+extension
    return Path(path).name
