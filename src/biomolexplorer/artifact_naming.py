"""Publication boundary for short output names, separate from engine work files."""
import csv
import hashlib
import json
from pathlib import Path
import re

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
    if not manifest.is_file():return {}
    try:return json.loads(manifest.read_text()).get(path.name,{})
    except (OSError,ValueError):return {}


def publish(stage,files,owned_root):
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


def human_name(record):
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
