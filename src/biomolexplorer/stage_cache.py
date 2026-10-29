"""Content-based reuse of project stage results, independent of canvas positions."""
import hashlib
import json
from pathlib import Path
from functools import lru_cache

from .catalog import PATH_FIELDS, template_names
from .paths import SOURCE_ROOT
from .templates import RESOURCE_ROOT
from .bindings import dependencies


def file_digest(path):
    with Path(path).open('rb') as stream:
        digest = hashlib.sha256()
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
        return digest.hexdigest()


def artifact_manifest(artifacts):
    return {str(Path(path).resolve()): file_digest(path) for path in artifacts}


def manifest_matches(manifest):
    if not manifest:
        return False
    try:
        return all(file_digest(path) == digest for path, digest in manifest.items())
    except OSError:
        return False


def implementation_digest():
    digest = hashlib.sha256()
    digest.update((SOURCE_ROOT/'biomolexplorer/preparation_failures.py').read_bytes())
    for name in ('docking_preparation.py','docking_data.py','dock6_validation.py','ligand_mol2.py','chimera_preparation.py','admet_outputs.py','molecular_surface.py','pdb_ligands.py','redocking_preparation.py','redocking_cache.py','redocking_dock6.py','redocking_config.py'):
        digest.update((SOURCE_ROOT/'biomolexplorer'/name).read_bytes())
    for folder in ('caad', 'crawlers', 'kernel', 'wrappers'):
        for path in sorted((SOURCE_ROOT / folder).rglob('*.py')):
            digest.update(str(path.relative_to(SOURCE_ROOT)).encode())
            digest.update(path.read_bytes())
    for name in ('operations.py','pipeline.py','input_validation.py','molecule_quality.py','fingerprint_selection.py','graph_inputs.py','docking_inputs.py','visualizations.py','admet_config.py','admet_predictors.py','admet_ai_runner.py','egg_plot.py','vina_session.py'):
        digest.update((SOURCE_ROOT / 'biomolexplorer' / name).read_bytes())
    return digest.hexdigest()


@lru_cache(maxsize=32)
def operation_digest(operation):
    """Version only the scientific implementation and its local import closure."""
    import ast
    from .operations import OPERATIONS
    spec=OPERATIONS.get(operation)
    # This module also contains RMSD/docking and graph code. Hashing its entire
    # import closure would invalidate fingerprints after an unrelated Vina fix.
    focused={
        'fingerprints':{'wrappers/molecular_analyzer.py':['generate_fingerprints'],
                        'kernel/descriptors.py':['Descriptors.__init__','Descriptors.get_fingerprints','Descriptors.compute_morgan_fingerprints',
                            'Descriptors.morgan_worker','Descriptors.compute_maccs_fingerprints','Descriptors.maccs_worker',
                            'Descriptors.compute_pharmacophore_fingerprints','Descriptors.pharmacophore_worker']},
        'similarity':{'wrappers/molecular_analyzer.py':['compute_similarity','read_filters'],
                      'kernel/descriptors.py':['MolSimilarity','similarityFunctions','fingerprints']}}
    if operation in focused and spec and spec.module=='wrappers.molecular_analyzer':
        digest=hashlib.sha256(b'focused-science-v1')
        for filename,names in focused[operation].items():
            tree=ast.parse((SOURCE_ROOT/filename).read_bytes())
            selected_nodes=[]
            for name in names:
                nodes=tree.body
                for part in name.split('.'):
                    selected=next(n for n in nodes if getattr(n,'name',None)==part)
                    nodes=getattr(selected,'body',[])
                digest.update(name.encode());digest.update(ast.dump(selected,include_attributes=False).encode())
                selected_nodes.append(selected)
            referenced={n.id for selected in selected_nodes for n in ast.walk(selected) if isinstance(n,ast.Name)}
            for node in tree.body:
                if isinstance(node,(ast.Import,ast.ImportFrom)) or isinstance(node,(ast.Assign,ast.AnnAssign)) and any(
                        isinstance(child,ast.Name) and child.id in referenced for child in ast.walk(node)):
                    digest.update(ast.dump(node,include_attributes=False).encode())
        helpers=['molecular_parsing.py','molecule_quality.py','input_validation.py','storage.py']
        if operation=='fingerprints':helpers+=['molecular_cache.py']
        for name in helpers:digest.update((SOURCE_ROOT/'biomolexplorer'/name).read_bytes())
        return digest.hexdigest()
    pending=[spec.module] if spec else []
    pending+=['biomolexplorer.input_validation','biomolexplorer.molecule_quality']
    infrastructure={'workspace','pipeline','jobs','stage_cache','result_sync','clear_results',
                    'project_state','project_folders','project_archive','project_merge','config','contracts','catalog','operations',
                    'bindings','paths','diagnostics','logging','progress','storage','artifacts','artifact_choices','artifact_naming',
                    'result_origins','result_files'}
    visited=set();digest=hashlib.sha256(b'operation-contract-v2')
    while pending:
        module=pending.pop()
        if module.startswith('biomolexplorer.ui.') or module.startswith('biomolexplorer.') and module.split('.')[1] in infrastructure:continue
        if module in visited:continue
        visited.add(module);path=SOURCE_ROOT.joinpath(*module.split('.')).with_suffix('.py')
        if not path.is_file():path=SOURCE_ROOT.joinpath(*module.split('.'),'__init__.py')
        if not path.is_file():continue
        data=path.read_bytes();digest.update(module.encode());digest.update(data)
        tree=ast.parse(data)
        for node in ast.walk(tree):
            names=[]
            if isinstance(node,ast.Import):names=[a.name for a in node.names]
            elif isinstance(node,ast.ImportFrom):
                base=node.module or ''
                if node.level:base='.'.join(module.split('.')[:-node.level])+('.'+base if base else '')
                names=[base]+[base+'.'+a.name for a in node.names]
            pending.extend(sorted(n for n in names if n.split('.')[0] in ('biomolexplorer','wrappers','kernel','caad','crawlers')))
    version=digest.hexdigest()
    # Audited routing-only changes may retain the previous scientific contract.
    # Later source changes produce another raw hash and cannot inherit approval.
    compatibility=RESOURCE_ROOT/'processing_compatibility.json'
    if compatibility.is_file():
        return json.loads(compatibility.read_text()).get('aliases',{}).get(version,version)
    return version


@lru_cache(maxsize=16)
def runtime_packages(python,include_models=False):
    """Read the scientific worker's versions, not the UI interpreter's packages."""
    import subprocess
    if not Path(python).is_file():return {'interpreter':python}
    code="import sys,json,importlib.metadata as m\nversions={}\nfor name in ['rdkit','numpy','pandas','scipy','networkx','biopython','torch','autogluon.tabular','admet-ai','chemprop','datasketch']:\n try: versions[name]=m.version(name)\n except m.PackageNotFoundError: pass\nmodels={}\n"
    if include_models:
        code+="import hashlib\nfrom pathlib import Path\ntry:\n distribution=m.distribution('admet-ai')\n for relative in distribution.files or []:\n  if Path(relative).suffix not in ('.pt','.ckpt','.pkl','.json','.py'):continue\n  path=distribution.locate_file(relative)\n  if not path.is_file():continue\n  digest=hashlib.sha256()\n  with path.open('rb') as stream:\n   for chunk in iter(lambda:stream.read(1048576),b''):digest.update(chunk)\n  models[str(relative)]=digest.hexdigest()\nexcept m.PackageNotFoundError:pass\n"
    code+="print(json.dumps(dict(python=sys.version,packages=versions,models=models)))"
    result=subprocess.run([python,'-c',code],capture_output=True,text=True,timeout=30,check=True)
    return json.loads(result.stdout)


def runtime_signature(operation,python,parameters):
    import os
    import shutil
    info=runtime_packages(python)
    groups={'fingerprints':('rdkit','numpy','pandas'),'similarity':('rdkit','numpy','pandas','scipy','datasketch'),
            'graphs':('rdkit','numpy','pandas','scipy','networkx'),
            'admet':('rdkit','numpy','pandas','scipy','torch','autogluon.tabular')}
    selected=groups.get(operation,('rdkit','numpy','pandas','biopython'))
    if operation=='admet' and parameters.get('prediction_backend','rdkit')!='admet_ai':selected=('rdkit','numpy','pandas','scipy')
    value=dict(python=info.get('python',python),packages={n:info.get('packages',{}).get(n) for n in selected},tools={})
    if operation=='admet' and parameters.get('prediction_backend')=='admet_ai':
        from .admet_predictors import isolated_python
        predictor=isolated_python()
        value['predictor']=runtime_packages(str(predictor or python),include_models=True)
    tools=[]
    if operation in ('prepare_structures','docking_vina','docking_dock6') or operation=='redocking' and parameters.get('prepare_complex',True):
        from .preparation_cache import tool_manifest
        value['tools'].update(tool_manifest(str(Path(python).parent)+os.pathsep+os.environ.get('PATH','')))
    if operation=='docking_vina' or operation=='redocking' and parameters.get('redocking_engine','vina')=='vina':tools+=['vina']
    if operation=='docking_dock6' or operation=='redocking' and parameters.get('redocking_engine')=='dock6':tools+=['dock6','grid','sphgen']
    search=str(Path(python).parent)+os.pathsep+os.environ.get('PATH','')
    for name in tools:
        configured=Path(parameters.get('dock6_app_path') or '')/'bin'/name
        path=str(configured) if name in ('dock6','grid','sphgen') and configured.is_file() else shutil.which(name,path=search)
        if path:value['tools'][name]=file_digest(path)
    return value


def path_manifest(path):
    path = Path(path)
    if path.is_file():
        return {path.name: file_digest(path)}
    ignored = {'.biomolexplorer', '.curation', '__pycache__', 'cache', 'logs', '.history', '.exports', '.input-cache', 'project.json'}
    return {p.relative_to(path).as_posix(): file_digest(p) for p in sorted(path.rglob('*'))
            if p.is_file() and not (set(p.relative_to(path).parts) & ignored)}


def stage_key(stage, parameters, inputs, software):
    config = {'version': 2, 'operation': stage['operation'], 'parameters': dict(parameters),
              'bindings': stage.get('bindings',{}), 'provided_results':stage.get('provided_results'),
              'inputs': dict(inputs), 'processing_version':operation_digest(stage['operation']),
              'runtime':runtime_signature(stage['operation'],software.split(':',1)[-1],parameters), 'templates': {},
              'input_processing':stage.get('input_processing','merge')}
    if stage.get('provided_results'):
        config['bindings']={}
    elif stage['operation']=='admet':
        # Schema 3 omitted independently selected rules/alerts. These outputs
        # must be regenerated once, never adopted as a complete schema 4 run.
        from .admet_outputs import SCHEMA_VERSION
        config['admet_output_schema']=SCHEMA_VERSION
    elif stage['operation']=='redocking':
        from .ligand_mol2 import PREPARATION_SCHEMA_VERSION
        config['redocking_preparation_schema']=PREPARATION_SCHEMA_VERSION
        if parameters.get('redocking_engine','vina')=='vina':
            unused=('dock6_app_path','density','radius','distance','conformer_search_type','max_orientations','num_scored_conformers','dock6_seed')
        else:unused=('sizeof_box','exhaustiveness','num_modes','dock6_app_path')
        for field in unused:config['parameters'].pop(field,None)
    if stage['operation']=='graphs' and parameters.get('graph_inputs'):
        # Source run paths must not change an otherwise identical experiment.
        config['parameters']['graph_inputs']=[dict({k:v for k,v in entry.items() if k!='label'},
            file={'name':Path(entry['file']).name,'sha256':file_digest(entry['file'])},
            compound_files=[{'name':Path(p).name,'sha256':file_digest(p)} for p in entry['compound_files']])
            for entry in parameters['graph_inputs']]
    # The installation directory is a machine setting. Scientific options and
    # immutable input/output hashes still determine whether an experiment fits.
    if stage['operation']=='docking_dock6':
        config['parameters'].pop('dock6_app_path',None)
    # Bound paths change with run IDs; input contents, selectors and inferred
    # filenames determine the scientific request instead of those output paths.
    for field in stage.get('bindings', {}):
        config['parameters'].pop(field, None)
    for field in PATH_FIELDS:
        if parameters.get(field) and field not in stage.get('bindings', {}):
            config['inputs']['path:' + field] = path_manifest(parameters[field])
            config['parameters'].pop(field,None)
    names=template_names(stage['operation'])
    if stage['operation']=='redocking':
        ignored='vina/' if parameters.get('redocking_engine','vina')=='dock6' else 'dock6/'
        names=[name for name in names if not name.startswith(ignored)]
    for name in names:
        config['templates'][name] = stage.get('templates', {}).get(name, (RESOURCE_ROOT / name).read_text())
    # Extra valid overlays can influence preparation protocols indirectly.
    config['templates'].update({name:value for name,value in stage.get('templates',{}).items()
                                if stage['operation']!='redocking' or not name.startswith(ignored)})
    encoded = json.dumps(config, sort_keys=True, ensure_ascii=False, allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def input_manifests(stage, results, selected=None):
    inputs = {}
    if selected:
        def identity(path):
            parts=Path(path).parts
            if 'artifacts' in parts:return '/'.join(parts[len(parts)-list(reversed(parts)).index('artifacts'):])
            return '/'.join(parts[-2:])
        return {'selected:'+field:sorted((identity(p),file_digest(p)) for p in paths)
                for field,paths in selected.items()}
    for identifier in sorted(dependencies(stage)):
        files = []
        for path in results[identifier]:
            parts = Path(path).parts
            # Keep the output's internal path while discarding its source run ID.
            start = len(parts) - 1 - list(reversed(parts)).index('artifacts') if 'artifacts' in parts else len(parts) - 2
            files.append(('/'.join(parts[start + 1:]), file_digest(path)))
        inputs['stage:' + identifier] = sorted(files)
    return inputs
