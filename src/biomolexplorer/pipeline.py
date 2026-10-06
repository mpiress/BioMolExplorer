"""Versioned DAG execution with project-scoped inputs, templates and permissions."""
import copy
import csv
import json
import os
import re
import shutil
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from itertools import product
from pathlib import Path
from uuid import uuid4

from .catalog import PATH_FIELDS, ENUMS, LABELS
from .config import AppConfig
from .jobs import JobManager, TERMINAL
from .operations import OPERATIONS, validate_operation
from .templates import validate_templates, materialize_templates
from .workspace import AccessDenied
from .diagnostics import get_logger, log_context, event, diagnose_exception
from .stage_cache import artifact_manifest, manifest_matches, implementation_digest, input_manifests, stage_key, path_manifest
from .bindings import sources as input_sources, dependencies as stage_dependencies, asset_references

ACTIVE_RUNS = ('queued', 'running', 'awaiting_input')


def execution_context(stages):
    return [dict(s['configuration'],_execution_batches=s.get('batches',[])) for s in stages]


def execution_modes(stages):
    """Every stage uses explicit choices, including legacy automatic pipelines."""
    return {identifier:'curated' for identifier in validate_pipeline(stages)}


def execution_stage(stage, mode):
    effective = copy.deepcopy(stage)
    effective.pop('_process_all_inputs',None)
    return effective


def validate_pipeline(stages):
    if not isinstance(stages,list) or len(stages) > 100:
        raise ValueError('O pipeline deve conter até 100 etapas.')
    ids = [stage.get('id') for stage in stages if isinstance(stage,dict)]
    if len(ids) != len(stages) or any(not isinstance(i,str) or not re.fullmatch('[a-f0-9]{32}',i) for i in ids) or len(ids) != len(set(ids)):
        raise ValueError('As etapas precisam de identificadores únicos.')
    dependencies = {}
    for stage in stages:
        if stage.get('operation') not in OPERATIONS and stage.get('operation') != 'import_results':
            raise ValueError('Etapa desconhecida.')
        if not isinstance(stage.get('parameters'),dict) or not isinstance(stage.get('bindings',{}),dict):
            raise ValueError('Parâmetros e conexões devem ser objetos JSON.')
        if type(stage.get('enabled',True)) is not bool:
            raise ValueError('Estado da etapa inválido.')
        if stage.get('input_processing','merge') not in ('merge','individual'):
            raise ValueError('Escolha mesclar arquivos ou processar individualmente.')
        if 'process_all' in stage and (type(stage['process_all']) is not bool or
                stage['operation'] not in ('retrieve_compounds','retrieve_structures','retrieve_zinc')):
            raise ValueError('A opção de processar tudo deve ser marcada ou desmarcada em um bloco de recuperação.')
        if not isinstance(stage.get('name'),str) or not stage['name'].strip() or len(stage['name']) > 180:
            raise ValueError('Dê um nome válido a cada etapa.')
        validate_templates(stage.get('templates',{}))
        provided=stage.get('provided_results')
        if provided is not None:
            from .input_validation import output_kind
            if (not isinstance(provided,dict) or set(provided)-{'kind','asset_ids','target'}
                    or provided.get('kind')!=output_kind(stage['operation'])
                    or not isinstance(provided.get('asset_ids'),list) or not provided['asset_ids']
                    or any(not isinstance(i,str) or not re.fullmatch('[a-f0-9]{32}',i) for i in provided['asset_ids'])):
                raise ValueError('Selecione resultados prontos compatíveis com o tipo deste bloco.')
        depends_on = stage.get('depends_on',[])
        if not isinstance(depends_on,list) or any(not isinstance(source,str) for source in depends_on):
            raise ValueError('As dependências de uma etapa devem ser uma lista de identificadores.')
        dependency = set(depends_on)
        spec = OPERATIONS.get(stage['operation'])
        accepted_inputs = PATH_FIELDS & set(spec.required + spec.optional) if spec else set()
        for field, group in stage.get('bindings',{}).items():
            if field not in accepted_inputs or not isinstance(group,dict):
                raise ValueError('Conexão de entrada inválida.')
            if 'sources' in group and (set(group)!={'sources'} or not isinstance(group['sources'],list) or not group['sources']):
                raise ValueError('Selecione ao menos um arquivo por entrada.')
            for binding in input_sources(group):
                if not isinstance(binding,dict) or set(binding) - {'stage','asset','selector'} or ('stage' in binding) == ('asset' in binding):
                    raise ValueError('Escolha uma etapa ou um arquivo como entrada.')
                reference = binding.get('stage',binding.get('asset'))
                if not isinstance(reference,str) or not re.fullmatch('[a-f0-9]{32}',reference):
                    raise ValueError('Uma conexão contém um identificador inválido.')
                selector = binding.get('selector','auto')
                if (not isinstance(selector,str) or not selector or '\\' in selector or
                        any(c in selector for c in '\n\r\x00') or Path(selector).is_absolute() or '..' in Path(selector).parts):
                    raise ValueError('Escolha um nome de arquivo válido para a conexão.')
                if 'stage' in binding:
                    dependency.add(binding['stage'])
        if provided:
            dependency=set()
        if dependency - set(ids) or stage['id'] in dependency:
            raise ValueError('Uma conexão aponta para uma etapa inexistente ou para ela mesma.')
        dependencies[stage['id']] = dependency
    ordered = []
    pending = dict(dependencies)
    while pending:
        ready = [stage_id for stage_id,deps in pending.items() if deps.issubset(ordered)]
        if not ready:
            raise ValueError('Há um ciclo nas conexões do pipeline.')
        for stage_id in ready:
            ordered.append(stage_id)
            del pending[stage_id]
    # Share the visual editor's data contracts with programmatic submissions.
    # The import stays local because flow.connect() also validates the DAG.
    from .flow import compatible
    by_id = {stage['id']:stage for stage in stages}
    for stage in stages:
        if stage.get('provided_results'):
            continue
        for field,group in stage.get('bindings',{}).items():
          for binding in input_sources(group):
            if 'stage' not in binding:continue
            port = 'base_vina_path' if stage['operation'] == 'consensus' and field == 'base_input_path' else field
            source = by_id[binding['stage']]
            if source['operation'] == 'import_results' and not isinstance(source['parameters'].get('kind','other'),str):
                raise ValueError('Tipo de importação inválido.')
            if not compatible(source,stage,port):
                raise ValueError(f'A entrada {LABELS.get(field,field)} do bloco “{stage["name"]}” está conectada a dados incompatíveis.')
    return ordered


def _csv_columns(path):
    try:
        with path.open(encoding='utf-8-sig', newline='') as stream:
            return set(next(csv.reader(stream),[]))
    except (OSError,UnicodeError):
        return set()


def select_input(artifacts, field, selector='auto'):
    files = [Path(p) for p in artifacts]
    if selector != 'auto':
        candidates = [p for p in files if p.name == selector or p.as_posix().endswith('/' + selector)]
        if not candidates:
            raise ValueError('O arquivo escolhido não está nos resultados da etapa.')
        if len(candidates)>1:
            raise ValueError('O nome escolhido identifica mais de um resultado. Selecione um caminho de arquivo específico.')
        candidates.sort(key=lambda p:(not (p.name=='compounds.csv' and p.parent.parent.name=='compounds'),len(p.parts),p.as_posix()))
        return candidates[0].parent
    if field in ('base_selected_mols', 'base_input_path'):
        compounds = [p for p in files if p.suffix == '.csv' and {'canonical_smiles','molecule_chembl_id'}.issubset(_csv_columns(p))]
        if compounds:
            # Prefer consolidated / selected datasets, then unfiltered ADMET results.
            compounds.sort(key=lambda p: (p.name != 'compounds.csv',p.name != 'molecules.csv', bool(re.search(r'_(BBB|HIA)',p.stem)),len(p.parts)))
            return compounds[0].parent
        fingerprints = [p for p in files if p.suffix == '.csv' and {'fingerprint','molecule_chembl_id'}.issubset(_csv_columns(p))]
        if fingerprints:
            return fingerprints[0].parent
        metadata = [p for p in files if p.name == 'pdb_codes.csv']
        if metadata:
            metadata.sort(key=lambda p: (not (p.parent / 'Prepared').is_dir(),len(p.parts)))
            return metadata[0].parent.parent
        pdbs = [p for p in files if p.suffix in ('.pdb','.pdbqt')]
        if pdbs:
            return pdbs[0].parent.parent
    if field == 'similarity_path':
        matching = [p for p in files if p.suffix == '.csv' and {'source','target','value'}.issubset(_csv_columns(p))]
        if matching:
            return matching[0].parent
    if field == 'base_vina_path':
        matching = [p for p in files if p.suffix == '.pdbqt']
        if matching:
            matching.sort(key=lambda p: ('Vina' not in p.parts,len(p.parts),p.name))
            return matching[0].parent
    if field == 'base_dock6_path':
        matching = [p for p in files if p.name.endswith('_scored.mol2')]
        if matching:
            return matching[0].parent
        matching = [p for p in files if p.suffix == '.pdbqt']
        if matching:
            return matching[0].parent
    raise ValueError('Não foi possível identificar esta entrada. Escolha um arquivo ou uma etapa compatível.')


class PipelineService:
    def __init__(self, store, worker_python=None, max_runs=2, cpu_workers=2, job_timeout=86400, dock6_path=None):
        self.store = store
        config = AppConfig(workspace=store.root,worker_python=worker_python,cpu_workers=cpu_workers,job_timeout=job_timeout)
        self.worker_python = config.worker_python
        self.cpu_workers = config.cpu_workers
        self.job_timeout = config.job_timeout
        self.dock6_path = Path(dock6_path).resolve() if dock6_path else None
        self._executor = ThreadPoolExecutor(max_workers=max_runs, thread_name_prefix='biomol-pipeline')
        self._lock = threading.RLock()
        self._active = {}
        self._closed = False
        import fcntl
        self._owner = (store.root / 'pipeline.lock').open('a+')
        try:
            fcntl.flock(self._owner, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            self._owner.close()
            self._executor.shutdown(wait=False)
            raise RuntimeError('Outra aplicação já está usando este workspace.') from None
        with store.connect() as db:
            interrupted_at = time.time()
            for run in db.execute("SELECT id,stages FROM runs WHERE status IN ('queued','running')").fetchall():
                stages = json.loads(run['stages'])
                for stage in stages:
                    if stage['status'] == 'running':
                        stage.update(status='interrupted', error='A aplicação foi interrompida durante esta etapa.',
                                     finished_at=interrupted_at)
                    elif stage['status'] == 'queued':
                        stage.update(status='skipped', error='A execução foi interrompida antes de iniciar esta etapa.')
                db.execute("UPDATE runs SET status='interrupted',stages=?,error='A aplicação foi interrompida.',updated=? WHERE id=?",
                           (json.dumps(stages),interrupted_at,run['id']))

    def submit(self, token, project_id, stage_ids=None, reuse_results=True):
        if type(reuse_results) is not bool:
            raise ValueError('A opção de reaproveitamento deve ser verdadeira ou falsa.')
        project = self.store.project(token, project_id, 'editor')
        user = self.store.user(token)
        stages = project['pipeline']
        order = validate_pipeline(stages)
        by_id = {stage['id']:stage for stage in stages}
        if not stages:
            raise ValueError('Adicione uma etapa antes de executar.')
        if stage_ids is not None:
            if not isinstance(stage_ids,(list,tuple,set)) or not stage_ids or any(not isinstance(i,str) for i in stage_ids):
                raise ValueError('Selecione pelo menos uma etapa para executar.')
            requested = set(stage_ids)
            if not requested.issubset(by_id):
                raise ValueError('Etapa inexistente.')
            # Execute selected stages and all required ancestors.
            while True:
                expanded = set(requested)
                for stage_id in requested:
                    stage = by_id[stage_id]
                    expanded.update(stage_dependencies(stage))
                if expanded == requested:
                    break
                requested = expanded
            stages = [by_id[i] for i in order if i in requested]
        else:
            stages = [by_id[i] for i in order]
        if not any(stage.get('enabled',True) for stage in stages):
            raise ValueError('Ative pelo menos um bloco antes de executar.')
        for stage in stages:
            if not stage.get('enabled',True):
                continue
            sources = stage_dependencies(stage)
            if any(not by_id[source].get('enabled',True) for source in sources):
                raise ValueError('Uma etapa ativa depende de uma etapa desativada.')
            try:
                self._validate_parameters(project_id, stage, partial=True)
            except (ValueError,FileNotFoundError) as exc:
                raise ValueError(f'O bloco “{stage["name"]}” precisa de ajustes: {exc}') from None
            if stage['operation'] == 'import_results' or stage.get('provided_results'):
                names = set()
                params=stage.get('provided_results') or stage['parameters']
                paths=[]
                for asset in params.get('asset_ids',[]):
                    path = self.store.asset_path(token,project_id,asset)
                    paths.append(path)
                    if not path.is_file():
                        raise ValueError(f'O bloco “{stage["name"]}” contém um arquivo indisponível. Envie-o novamente.')
                    if path.name in names:
                        raise ValueError(f'O bloco “{stage["name"]}” contém arquivos com o mesmo nome. Renomeie um deles antes de importar.')
                    names.add(path.name)
                from .input_validation import validate_bundle
                validate_bundle(paths,params['kind'],stage['operation'] if stage.get('provided_results') else None)
            for group in ({} if stage.get('provided_results') else stage.get('bindings',{})).values():
              for binding in input_sources(group):
                if 'asset' in binding:
                    if not self.store.asset_path(token,project_id,binding['asset']).is_file():
                        raise ValueError(f'A entrada do bloco “{stage["name"]}” está indisponível. Envie o arquivo novamente.')
        run_id = uuid4().hex
        snapshot = [{'id':stage['id'],'name':stage['name'],'operation':stage['operation'],
                     'status':'queued' if stage.get('enabled',True) else 'skipped','configuration':stage,'artifacts':[],
                     'input_mode':'curated', 'reuse_results':reuse_results, 'requires_curation':
                         bool(stage_dependencies(stage))}
                    for stage in stages]
        with self._lock:
            if self._closed:
                raise RuntimeError('A aplicação está sendo encerrada.')
            if len(self._active) >= 20:
                raise ValueError('A fila de execução está cheia. Tente novamente após a conclusão de uma execução.')
            with self.store.connect() as db:
                db.execute('BEGIN IMMEDIATE')
                # Project deletion and collaboration changes can happen while
                # validating files. Recheck before registering this execution.
                self.store._require_user(user['id'],project_id,'editor')
                if db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')", (project_id,)).fetchone():
                    raise ValueError('Este projeto já tem uma execução ativa.')
                db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,NULL)',
                           (run_id,project_id,user['id'],'queued',json.dumps(snapshot),time.time(),time.time()))
            event = threading.Event()
            self._active[run_id] = {'event': event, 'manager': None, 'job_id': None, 'reuse_results': reuse_results}
            try:
                self._executor.submit(self._run, run_id, project_id, user['id'], snapshot)
            except Exception:
                self._update(run_id,'failed',snapshot,'Não foi possível iniciar a execução.')
                self._active.pop(run_id,None)
                raise
        return self.store.get_run(token,run_id)

    def _validate_parameters(self, project_id, stage, partial=False, cache_only=False):
        if stage.get('provided_results'):
            return
        params = stage['parameters']
        operation = stage['operation']
        for field,value in params.items():
            if field in PATH_FIELDS | {'dock6_app_path'} and value is not None and not isinstance(value,str):
                raise ValueError(LABELS.get(field,field) + ' deve indicar uma pasta do projeto.')
        if operation == 'import_results':
            if set(params) - {'kind','asset_ids','target'} or not params.get('asset_ids'):
                raise ValueError('Selecione arquivos para a etapa de importação.')
            if not isinstance(params['asset_ids'],list) or any(not isinstance(asset,str) or not re.fullmatch('[a-f0-9]{32}',asset) for asset in params['asset_ids']):
                raise ValueError('Selecione arquivos válidos para a etapa de importação.')
            if params.get('kind') not in ('compounds','structures','prepared_structures','fingerprints','similarity','vina','dock6','scores','other'):
                raise ValueError('Tipo de importação inválido.')
        else:
            supplied = dict(params)
            for field in stage.get('bindings',{}):
                supplied[field] = 'pending-input'
            if operation == 'consensus' and not supplied.get('base_input_path') and supplied.get('base_vina_path'):
                supplied['base_input_path'] = supplied['base_vina_path']
            if partial:
                required = set(OPERATIONS[operation].required)
                if operation == 'expand_similar_compounds':
                    required.add('base_input_path')
                fallback = None
                if params.get('base_input_path'):
                    candidate = Path(params['base_input_path'])
                    fallback = self.store.scoped_path(project_id,candidate if candidate.is_absolute() else self.store.project_dir(project_id) / candidate)
                if operation == 'graphs' and not any(supplied.get(k) for k in ('similarity_path','graph_inputs')):
                    raise ValueError('Conecte um bloco Calcular similaridade ou forneça um CSV externo com source,target,value para gerar os grafos.')
                if operation == 'consensus':
                    for field,folder in (('base_vina_path','Vina'),('base_dock6_path','Dock6')):
                        if not (fallback and (fallback / folder).is_dir()):
                            required.add(field)
                missing = {field for field in required if not supplied.get(field)}
                if missing:
                    raise ValueError('Conecte um bloco de origem ou selecione seus arquivos em: ' +
                                     ', '.join(LABELS.get(field,field) for field in sorted(missing)))
            # Connected inputs are curated after retrieval; execution validates the selected pairs.
            validate_operation(operation,supplied,defer_redocking_selection=partial and bool(stage.get('bindings',{})))
        for field,value in params.items():
            if field == 'target' and operation == 'import_results' and (not isinstance(value,str) or not re.fullmatch('[A-Za-z0-9_-]{1,80}',value)):
                raise ValueError('O nome da pasta do alvo deve conter apenas letras, números, _ e -.')
            if field=='target' and (not isinstance(value,str) or not re.fullmatch('[A-Za-z0-9 _-]{1,80}',value)):
                raise ValueError('Use letras, números, espaço, _ ou - no nome da pasta do alvo.')
            if field in PATH_FIELDS and value:
                path = Path(value)
                path = self.store.scoped_path(project_id,path if path.is_absolute() else self.store.project_dir(project_id) / path)
                if not path.is_dir():
                    raise ValueError('A pasta de ' + LABELS.get(field,field).lower() + ' não está disponível. Escolha uma entrada do projeto.')
            if field == 'dock6_app_path' and value:
                if partial or cache_only:
                    continue
                if self.dock6_path is None or Path(value).resolve() != self.dock6_path:
                    raise AccessDenied('Use a instalação DOCK6 configurada pelo administrador.')
            if field in ('filename','input_file','mol_filename') and value:
                if not isinstance(value,str) or Path(value).name != value or '\\' in value or any(c in value for c in '\n\r\x00'):
                    raise ValueError('O nome do arquivo deve pertencer à entrada selecionada.')
            if field == 'files' and value:
                if not isinstance(value,list) or not all(isinstance(v,str) and Path(v).name == v and '\\' not in v for v in value):
                    raise ValueError('Informe apenas nomes de arquivos da entrada selecionada.')
            if field in ENUMS and value is not None and value not in ENUMS[field]:
                raise ValueError('Valor inválido para ' + field)
            if field in ('pdb_code','pdb_codes') and value:
                if not isinstance(value,(list,tuple)):
                    raise ValueError('Informe os registros PDB como uma lista JSON.')
                records = value if isinstance(value[0],(list,tuple)) else [value]
                for record in records:
                    if len(record) < 4 or any(not re.fullmatch('[A-Za-z0-9_.+-]{1,40}',str(item)) for item in record[:4]):
                        raise ValueError('Use registros PDB com identificador, ligante, número de resíduo e cadeia.')

    def _update(self, run_id, status, stages, error=None):
        with self.store.connect() as db:
            db.execute('BEGIN IMMEDIATE')
            prior=db.execute('SELECT * FROM runs WHERE id=?',(run_id,)).fetchone()
            db.execute('UPDATE runs SET status=?,stages=?,updated=?,error=? WHERE id=?',
                       (status,json.dumps(stages),time.time(),error,run_id))
            if prior and prior['status'] in ACTIVE_RUNS and status not in ACTIVE_RUNS:
                from .project_state import read_snapshot,record
                last=db.execute('SELECT id FROM project_history WHERE project_id=? ORDER BY created DESC LIMIT 1',(prior['project_id'],)).fetchone()
                before=read_snapshot(self.store,prior['project_id'],last['id'],'after') if last else None
                record(self.store,db,prior['project_id'],prior['user_id'],'execution',f'Execução {run_id[:8]}: {status}',before)

    def _materialize_inputs(self, project_id, stage, field, refs, results, directory=None):
        from .input_validation import columns, merge_csv, validate_file
        from .project_state import digest
        import hashlib
        operation=stage['operation']
        kind=('similarity' if field=='similarity_path' else 'fingerprints' if operation=='similarity' else
              'vina' if field=='base_vina_path' or operation=='consensus' and field=='base_input_path' else 'dock6' if field=='base_dock6_path' else
              'compounds' if field=='base_selected_mols' or operation in ('admet','fingerprints','graphs','expand_similar_compounds') else
              'other' if operation=='retrieve_zinc' else 'prepared_structures' if operation in ('docking_vina','docking_dock6')
              or operation=='redocking' and not stage['parameters'].get('prepare_complex',True) else 'structures')
        required={'compounds':{'canonical_smiles'},'fingerprints':{'fingerprint','molecule_chembl_id'},'similarity':{'source','target','value'}}
        selected=[]
        for ref in refs:
            start=len(selected)
            if 'asset' in ref:
                with self.store.connect() as db:
                    asset=db.execute('SELECT path FROM assets WHERE project_id=? AND id=?',(project_id,ref['asset'])).fetchone()
                if asset is None:raise AccessDenied('Arquivo não autorizado.')
                selected.append(self.store.scoped_path(project_id,asset['path']))
                continue
            upstream=results.get(ref['stage'])
            if not upstream:raise ValueError('A origem desta entrada não possui resultados disponíveis.')
            selector=ref.get('selector','auto')
            if selector!='auto':
                folder=select_input(upstream,field,selector)
                selected.append(next(Path(p) for p in upstream if Path(p).parent==folder and Path(p).name==Path(selector).name))
            elif kind in required:
                matching=[Path(p) for p in upstream if Path(p).suffix.lower()=='.csv' and required[kind].issubset(columns(p))]
                matching.sort(key=lambda p:(p.name!='compounds.csv',p.name!='molecules.csv',bool(re.search(r'_(BBB|HIA)',p.stem)),len(p.parts)))
                if not matching:raise ValueError('A origem não possui um arquivo compatível. Selecione o resultado nas configurações do bloco.')
                selected.extend(matching if stage.get('_process_all_inputs') else matching[:1])
            else:
                allowed={'structures':{'.pdb'},'prepared_structures':{'.pdb','.pdbqt'},
                         'vina':{'.pdbqt'},'dock6':{'.mol2'},'other':{'.txt'}}[kind]
                candidates=[Path(p) for p in upstream if Path(p).suffix.lower() in allowed]
                if kind=='prepared_structures' and any(p.parent.name=='Prepared' for p in candidates):
                    candidates=[p for p in candidates if p.parent.name=='Prepared']
                if kind=='structures':
                    candidates=[p for p in candidates if p.name.count('.')==1]
                selected.extend(candidates)
            if kind in ('structures','prepared_structures'):
                chosen=selected[start:]
                files=[Path(p) for p in upstream]
                folders={p.parent for p in chosen}
                if kind=='prepared_structures':folders|={p.parent.parent for p in chosen if p.parent.name=='Prepared'}
                selected.extend(p for p in files if p.parent in folders and p.name in ('pdb_codes.csv','centers.csv'))
                if kind=='prepared_structures':
                    receptors={p.name.split('.',1)[0] for p in chosen if '.dockprep.' in p.name}
                    ligands=set()
                    for metadata in (p for p in files if p.parent in folders and p.name=='pdb_codes.csv'):
                        with metadata.open(encoding='utf-8-sig',newline='') as stream:
                            for record in csv.DictReader(stream):
                                if f"{record.get('PDB_CODE')}_{record.get('CHAIN')}" in receptors:
                                    ligands.update((f"{record['PDB_CODE']}_{record['LIGAND']}_{record['RESNUM']}{record['CHAIN']}",
                                        f"{record['PDB_CODE']}_{record['LIGAND']}_{record['RESNUM']}_{record['CHAIN']}"))
                    selected.extend(p for p in files if p.parent in folders and
                        p.suffix.lower() in ('.pdb','.pdbqt','.mol2') and
                        p.name.split('.',1)[0] in receptors|ligands)
                if not any(p.suffix.lower() in ('.pdb','.pdbqt') for p in chosen):
                    data_folders=folders|{folder/'Prepared' for folder in folders} if kind=='prepared_structures' else folders
                    selected.extend(p for p in files if p.parent in data_folders and p.suffix.lower() in ('.pdb','.pdbqt'))
        selected=list(dict.fromkeys(selected))
        if not selected:raise ValueError('Selecione ao menos um arquivo de entrada.')
        for file in selected:validate_file(file,kind,validate_rows=kind not in required)
        key=hashlib.sha256(json.dumps([(str(p),digest(p)) for p in selected]+[(operation,field,stage['parameters'].get('target'))],sort_keys=True).encode()).hexdigest()
        destination=(directory or self.store.project_dir(project_id)/'.input-cache'/key)/field
        destination.mkdir(parents=True,exist_ok=True)
        filename=None
        if kind in required:
            filename=selected[0].name if stage.get('_individual_input') and len(selected)==1 else 'selected_'+kind+'.csv'
            merged=destination/filename
            if not merged.exists():
                from .molecule_quality import QualityReport,merge_clean_csv
                quality=QualityReport(operation)
                merge_clean_csv(selected,merged,kind,quality)
                quality.write(destination)
        elif kind=='other':
            filename='selected_urls.txt'
            (destination/filename).write_text(''.join(p.read_text().rstrip()+'\n' for p in selected),encoding='utf-8')
        else:
            target=stage['parameters'].get('target','MeuAlvo').replace(' ','')
            data=destination/target if kind in ('structures','prepared_structures') else destination
            names={};metadata_sources={}
            for source in selected:
                file=(data/'Prepared' if kind=='prepared_structures' and source.name!='pdb_codes.csv' else data)/source.name
                if kind in ('structures','prepared_structures') and source.name in ('pdb_codes.csv','centers.csv'):
                    metadata_sources.setdefault(file,[]).append(source)
                    continue
                if file in names and digest(source)!=digest(names[file]):raise ValueError('Arquivos de entrada com o mesmo nome possuem conteúdos diferentes. Renomeie-os antes de combinar.')
                names[file]=source;file.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(source,file)
            for file,paths in metadata_sources.items():
                fields=[];records={};centers={}
                for source in paths:
                    with source.open(encoding='utf-8-sig',newline='') as stream:
                        reader=csv.DictReader(stream)
                        headers=reader.fieldnames or [];rows=list(reader)
                    fields.extend(f for f in headers if f not in fields)
                    if file.name=='centers.csv':
                        for field in headers:
                            values=[r[field] for r in rows]
                            if field in centers and centers[field]!=values:raise ValueError('Centros conflitantes para '+field+'.')
                            centers[field]=values
                    else:
                        for row in rows:
                            key=tuple(row[k] for k in ('PDB_CODE','LIGAND','RESNUM','CHAIN'))
                            old=records.setdefault(key,{})
                            if any(old.get(k) and v and old[k]!=v for k,v in row.items()):
                                raise ValueError('Metadados PDB conflitantes para '+str(key)+'.')
                            old.update({k:v for k,v in row.items() if v})
                file.parent.mkdir(parents=True,exist_ok=True)
                with file.open('w',encoding='utf-8',newline='') as stream:
                    writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
                    writer.writerows([{f:centers[f][i] for f in fields} for i in range(3)] if file.name=='centers.csv' else records.values())
            metadata=data/'pdb_codes.csv'
            if kind in ('structures','prepared_structures') and metadata.exists():
                with metadata.open(encoding='utf-8-sig',newline='') as stream:
                    reader=csv.DictReader(stream);fields=reader.fieldnames;rows=list(reader)
                filenames={p.name for p in selected}
                if kind=='structures':rows=[r for r in rows if r['PDB_CODE']+'.pdb' in filenames]
                else:rows=[r for r in rows if f'{r["PDB_CODE"]}_{r["CHAIN"]}.dockprep.pdbqt' in filenames]
                if not rows:raise ValueError('Os arquivos selecionados não possuem registros PDB compatíveis nos metadados.')
                with metadata.open('w',encoding='utf-8',newline='') as stream:
                    writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
                centers=data/'Prepared'/'centers.csv'
                if kind=='prepared_structures' and centers.exists():
                    with centers.open(encoding='utf-8-sig',newline='') as stream:
                        reader=csv.DictReader(stream)
                        keys={key for r in rows for key in (f'{r["PDB_CODE"]}_{r["LIGAND"]}_{r["RESNUM"]}{r["CHAIN"]}',
                            f'{r["PDB_CODE"]}_{r["LIGAND"]}_{r["RESNUM"]}_{r["CHAIN"]}')}
                        fields=[f for f in reader.fieldnames or [] if f in keys]
                        coordinates=list(reader)
                    if not fields:raise ValueError('Os centros selecionados não correspondem aos registros PDB.')
                    with centers.open('w',encoding='utf-8',newline='') as stream:
                        writer=csv.DictWriter(stream,fieldnames=fields,extrasaction='ignore')
                        writer.writeheader();writer.writerows(coordinates)
        return destination,filename,selected

    def _resolve(self, project_id, user_id, stage, results, input_directory=None, item=None, cache_only=False, pipeline=None):
        params = dict(stage['parameters'])
        if stage.get('_process_all_inputs') and 'base_input_path' in stage.get('bindings',{}):
            if stage['operation']=='prepare_structures':params.pop('pdb_codes',None)
            elif stage['operation']=='docking_vina':params.pop('pdb_code',None)
        if stage['operation']=='graphs':
            from .graph_inputs import GraphInputs
            if pipeline is None:
                with self.store.connect() as db:
                    row=db.execute('SELECT pipeline FROM projects WHERE id=?',(project_id,)).fetchone()
                pipeline=json.loads(row[0])
            params=GraphInputs(self.store,project_id,pipeline,results).resolve(stage,item)
            entries=params.get('graph_inputs',[])
            if stage.get('input_processing')=='merge' and len(entries)>1:
                if len({(e.get('metric','pronta'),e.get('fingerprint','morgan')) for e in entries})>1:
                    raise ValueError('Mescle apenas similaridades com a mesma métrica e fingerprint.')
                from .stage_cache import file_digest
                import hashlib
                paths=[Path(e['file']) for e in entries]
                key=hashlib.sha256(json.dumps([(str(p),file_digest(p)) for p in paths]).encode()).hexdigest()
                folder=(input_directory or self.store.project_dir(project_id)/'.input-cache'/key)/'similarity_path'
                folder.mkdir(parents=True,exist_ok=True)
                from .input_validation import merge_csv
                merged=merge_csv(paths,folder/'selected_similarity.csv','similarity')
                params['graph_inputs']=[dict(entries[0],file=str(merged),label='Similaridades mescladas',
                    compound_files=list(dict.fromkeys(p for e in entries for p in e['compound_files'])))]
            self._validate_parameters(project_id,dict(stage,parameters=params),cache_only=cache_only)
            return params
        if stage['operation']=='similarity':
            from .fingerprint_selection import generated_kind
            if pipeline is None:
                with self.store.connect() as db:
                    row=db.execute('SELECT pipeline FROM projects WHERE id=?',(project_id,)).fetchone()
                pipeline=json.loads(row[0])
            kind,_=generated_kind(stage,pipeline,results)
            if kind:params['fingerprint']=kind
        if stage['operation']=='docking_dock6' and self.dock6_path is not None:
            params['dock6_app_path']=str(self.dock6_path)
        bindings = dict(stage.get('bindings',{}))
        if stage['operation'] == 'consensus' and not params.get('base_input_path') and 'base_input_path' not in bindings:
            if 'base_vina_path' in bindings:
                bindings['base_input_path'] = bindings['base_vina_path']
            elif params.get('base_vina_path'):
                params['base_input_path'] = params['base_vina_path']
        for field,group in bindings.items():
            refs=input_sources(group)
            if (stage.get('_process_all_inputs') or
                    len(refs)>1 or any('asset' in ref for ref in refs) or
                    any(ref.get('selector','auto')!='auto' for ref in refs)):
                path,filename,selected=self._materialize_inputs(project_id,stage,field,refs,results,input_directory)
                params[field]=str(path)
                if stage.get('_individual_input') and field=='base_input_path' and stage['operation'] in ('prepare_structures','redocking','docking_vina','docking_dock6'):
                    raw={p.stem for p in selected if p.suffix.lower()=='.pdb' and p.name.count('.')==1}
                    receptors={p.name.split('.',1)[0] for p in selected if '.dockprep.' in p.name}
                    def matches(record):return str(record[0]) in raw or f'{record[0]}_{record[3]}' in receptors
                    metadata=path/params['target'].replace(' ','')/'pdb_codes.csv'
                    records=[]
                    if metadata.exists():
                        with metadata.open(encoding='utf-8-sig',newline='') as stream:
                            records=[[r['PDB_CODE'],r['LIGAND'],int(r['RESNUM']),r['CHAIN'],
                                      float(r.get('RESOLUTION') or 0)] for r in csv.DictReader(stream)]
                    key='pdb_codes' if stage['operation'] in ('prepare_structures','redocking') else 'pdb_code'
                    configured=params.get(key)
                    if configured:
                        configured=configured if isinstance(configured[0],(list,tuple)) else [configured]
                        chosen=[r for r in configured if matches(r)]
                        if not chosen and stage['operation']!='redocking':chosen=[r for r in records if matches(r)]
                        if not chosen:raise ValueError('Os registros PDB configurados não correspondem à estrutura selecionada.')
                        params[key]=chosen if key=='pdb_codes' else chosen[0][:4]
                        if stage['operation']=='redocking':
                            from .redocking_config import pair_key
                            keys={pair_key(r) for r in chosen}
                            params['preparation_pairs']={k:v for k,v in params.get('preparation_pairs',{}).items() if k in keys}
                if filename:
                    if stage['operation'] in ('admet','expand_similar_compounds'):params['input_file']=filename
                    elif stage['operation']=='retrieve_zinc':params['filename']=filename
                    elif stage['operation']=='similarity':params['filename']=filename
                    elif stage['operation']=='fingerprints':params['files']=[filename]
                    elif field=='base_selected_mols':params['mol_filename']=Path(filename).stem
                if item is not None:item.setdefault('input_files',{})[field]=[str(p) for p in selected]
                continue
            binding=refs[0]
            if 'stage' in binding:
                upstream = results.get(binding['stage'])
                if not upstream:
                    raise ValueError('Uma entrada depende de uma etapa desativada ou sem resultados.')
                path = select_input(upstream,field,binding.get('selector','auto'))
                if field=='base_input_path' and stage['operation']=='expand_similar_compounds':
                    original=next((Path(p) for p in upstream if 'ChEMBL' in Path(p).parts),None)
                    if original is None:
                        raise ValueError('A expansão requer uma etapa de recuperação com os downloads ChEMBL originais.')
                    path=Path(*original.parts[:original.parts.index('ChEMBL')])
                if field=='base_input_path' and stage['operation']=='admet' and binding.get('selector','auto')!='auto':
                    params['input_file'] = Path(binding['selector']).name
                if field=='base_input_path' and stage['operation']=='similarity' and binding.get('selector','auto')!='auto':
                    params['filename']=Path(binding['selector']).name
                if field=='base_input_path' and stage['operation']=='fingerprints' and binding.get('selector','auto')!='auto':
                    params['files']=[Path(binding['selector']).name]
            else:
                with self.store.connect() as db:
                    asset = db.execute('SELECT path FROM assets WHERE project_id=? AND id=?', (project_id,binding['asset'])).fetchone()
                if asset is None:
                    raise AccessDenied('Arquivo não autorizado.')
                file = self.store.scoped_path(project_id,asset['path'])
                path = file.parent
                if field == 'base_input_path' and stage['operation'] in ('admet','retrieve_zinc','similarity'):
                    params['input_file' if stage['operation'] == 'admet' else 'filename'] = file.name
                if field=='base_input_path' and stage['operation']=='fingerprints':
                    params['files']=[file.name]
            params[field] = str(self.store.scoped_path(project_id,path))
            if item is not None:
                selected=[p for p in upstream if Path(p).parent==path and (binding.get('selector','auto')=='auto' or Path(p).name==Path(binding['selector']).name)]
                item.setdefault('input_files',{})[field]=selected
            if field=='base_selected_mols' and stage['operation'] in ('docking_vina','docking_dock6'):
                selector=binding.get('selector','auto')
                if selector!='auto':
                    params['mol_filename']=Path(selector).stem
                elif 'asset' in binding:
                    params['mol_filename']=file.stem
                elif params.get('mol_filename')=='molecules' and not (path/'molecules.csv').exists():
                    compounds=[p for p in path.glob('*.csv') if {'canonical_smiles','molecule_chembl_id'}.issubset(_csv_columns(p))]
                    compounds.sort(key=lambda p:(p.name!='compounds.csv',bool(re.search(r'_(BBB|HIA)',p.stem)),p.name))
                    if compounds:
                        params['mol_filename']=compounds[0].stem
            if field == 'base_input_path' and stage['operation'] == 'admet' and 'input_file' not in params:
                compounds = [p for p in path.glob('*.csv') if {'canonical_smiles','molecule_chembl_id'}.issubset(_csv_columns(p))]
                compounds.sort(key=lambda p: (p.name not in ('compounds.csv','molecules.csv'),p.name))
                if compounds:
                    params['input_file'] = compounds[0].name
        for field in PATH_FIELDS:
            if params.get(field):
                path = Path(params[field])
                params[field] = str(self.store.scoped_path(project_id,path if path.is_absolute() else self.store.project_dir(project_id) / path))
        self._validate_parameters(project_id,dict(stage,parameters=params),cache_only=cache_only)
        if stage['operation']=='redocking' and not cache_only:
            from .redocking_config import validate_structure_pairs
            validate_structure_pairs(Path(params['base_input_path'])/params['target'].replace(' ',''),
                params['pdb_codes'], params.get('preparation_pairs') or {}, prepared=not params.get('prepare_complex',True))
        if stage['operation']=='retrieve_zinc':
            from urllib.parse import urlsplit
            with (Path(params['base_input_path']) / params['filename']).open() as stream:
                for line in stream:
                    if not line.strip():
                        continue
                    url = urlsplit(line.strip())
                    if url.scheme!='https' or url.hostname not in ('zinc.docking.org','zinc15.docking.org','zinc20.docking.org','files.docking.org') or url.username or url.password or url.port not in (None,443):
                        raise ValueError('O arquivo ZINC deve conter endereços HTTPS dos servidores ZINC.')
        return params

    def _import(self, project_id, params, output):
        output.mkdir(parents=True)
        target = params.get('target','MeuAlvo')
        if not re.fullmatch('[A-Za-z0-9_-]{1,80}',target):
            raise ValueError('O nome da pasta do alvo deve conter apenas letras, números, _ e -.')
        destination = output / target if params['kind'] in ('structures','prepared_structures') else output
        if params['kind'] == 'prepared_structures':
            destination = destination / 'Prepared'
        destination.mkdir(parents=True,exist_ok=True)
        names = set()
        for asset_id in params['asset_ids']:
            with self.store.connect() as db:
                asset = db.execute('SELECT path,name FROM assets WHERE project_id=? AND id=?', (project_id,asset_id)).fetchone()
            if asset is None:
                raise AccessDenied('Arquivo não autorizado.')
            source = self.store.scoped_path(project_id,asset['path'])
            if source.name in names:
                raise ValueError('Há dois arquivos com o mesmo nome na importação.')
            names.add(source.name)
            file = destination / source.name
            if source.name == 'pdb_codes.csv' and params['kind'] == 'prepared_structures':
                file = destination.parent / source.name
            shutil.copyfile(source,file)
            if params['kind'] == 'compounds' and file.suffix.lower() == '.csv':
                self._normalize_compounds(file)
        return [str(p) for p in output.rglob('*') if p.is_file()]

    @staticmethod
    def _normalize_compounds(path):
        from tempfile import NamedTemporaryFile
        aliases = {'smiles':'canonical_smiles','Canonical_SMILES':'canonical_smiles','name':'molecule_chembl_id'}
        temporary = None
        try:
            with path.open(encoding='utf-8-sig',newline='') as stream, NamedTemporaryFile(mode='w',dir=path.parent,delete=False,newline='',encoding='utf-8') as out:
                temporary = Path(out.name)
                reader = csv.DictReader(stream)
                fields = [aliases.get(k,k) for k in reader.fieldnames or []]
                if 'canonical_smiles' not in fields:
                    raise ValueError('O CSV de compostos precisa de uma coluna smiles ou canonical_smiles.')
                if 'molecule_chembl_id' not in fields:
                    fields.insert(0,'molecule_chembl_id')
                if len(set(fields)) != len(fields):
                    raise ValueError('O CSV contém colunas duplicadas após a conversão.')
                writer = csv.DictWriter(out,fieldnames=fields)
                writer.writeheader()
                for index,row in enumerate(reader,1):
                    result = {aliases.get(k,k):v for k,v in row.items()}
                    if not result.get('molecule_chembl_id'):
                        from rdkit import Chem
                        import hashlib
                        smiles=Chem.MolToSmiles(Chem.MolFromSmiles(result['canonical_smiles']))
                        result['molecule_chembl_id']='USER_'+hashlib.sha256(smiles.encode()).hexdigest()[:12]
                    identifier=result['molecule_chembl_id']
                    if not identifier or not re.fullmatch('[A-Za-z0-9_.+-]{1,100}',identifier) or '..' in identifier:
                        raise ValueError('Use identificadores moleculares com letras, números, _, -, + e ponto, sem caminhos.')
                    writer.writerow(result)
            temporary.replace(path)
        finally:
            if temporary:
                temporary.unlink(missing_ok=True)

    def _variants(self, project_id, stage, results=None):
        """One job per selected file; distinct input ports form combinations."""
        if stage['operation']=='redocking' and results is not None:
            from .bindings import pack
            stage=copy.deepcopy(stage)
            refs=[]
            prepared=not stage['parameters'].get('prepare_complex',True)
            records=stage['parameters'].get('pdb_codes') or []
            identities={f'{r[0]}_{r[3]}.dockprep.pdbqt' if prepared else f'{r[0]}.pdb' for r in records}
            for ref in input_sources(stage.get('bindings',{}).get('base_input_path',{})):
                if 'stage' not in ref or ref.get('selector','auto')!='auto':
                    refs.append(ref);continue
                files=[Path(p) for p in results.get(ref['stage'],[])]
                for path in files:
                    if path.name not in identities:continue
                    selector=path.name
                    for length in range(1,len(path.parts)):
                        selector='/'.join(path.parts[-length:])
                        if sum(p.as_posix().endswith('/'+selector) for p in files)==1:break
                    refs.append(dict(ref,selector=selector))
            if not refs:raise ValueError('Selecione arquivos de estruturas correspondentes aos pares de redocking configurados.')
            stage['bindings']['base_input_path']=pack(refs)
        fields=[];groups=[];shared={};asset_names={}
        bindings=stage.get('bindings',{})
        for field,group in bindings.items():
            if stage['operation']=='consensus' and field=='base_input_path' and group==bindings.get('base_vina_path'):continue
            refs=input_sources(group)
            primary=[];metadata=[]
            for ref in refs:
                name=Path(ref.get('selector','')).name
                if 'asset' in ref:
                    with self.store.connect() as db:
                        asset=db.execute('SELECT path FROM assets WHERE project_id=? AND id=?',(project_id,ref['asset'])).fetchone()
                    if asset is None:raise AccessDenied('Arquivo não autorizado.')
                    name=Path(asset['path']).name
                    asset_names[ref['asset']]=name
                (metadata if name in ('pdb_codes.csv','centers.csv') else primary).append(ref)
            fields.append(field);groups.append(primary or metadata);shared[field]=metadata if primary else []
        from .bindings import pack
        emitted=False
        for selection in product(*groups):
            if stage['operation']=='redocking' and 'base_input_path' in fields:
                ref=dict(zip(fields,selection))['base_input_path']
                filename=asset_names.get(ref.get('asset'),Path(ref.get('selector','')).name)
                if filename not in ('','auto'):
                    records=stage['parameters'].get('pdb_codes') or []
                    identities={str(r[0]) if filename.endswith('.pdb') else f'{r[0]}_{r[3]}' for r in records}
                    if filename.split('.',1)[0] not in identities:continue
            if stage['operation']=='docking_dock6' and results is not None and {'base_input_path','base_selected_mols','base_vina_path'}<=set(fields):
                from .docking_inputs import dock6_variant_matches
                if not dock6_variant_matches(stage,dict(zip(fields,selection)),results,self.store,project_id):continue
            if stage['operation']=='consensus' and {'base_vina_path','base_dock6_path'}<=set(fields):
                pair={field:ref for field,ref in zip(fields,selection)}
                def identity(ref):
                    filename=asset_names.get(ref.get('asset'),Path(ref.get('selector','')).name)
                    if filename in ('','auto'):return None
                    return filename.removesuffix('_scored.mol2').removesuffix('.pdbqt').removesuffix('.lig')
                vina,dock6=identity(pair['base_vina_path']),identity(pair['base_dock6_path'])
                if vina and dock6 and vina!=dock6:continue
            variant=copy.deepcopy(stage);variant['input_processing']='merge';variant['_individual_input']=True
            labels=[]
            for field,ref in zip(fields,selection):
                variant['bindings'][field]=pack([ref]+shared[field])
                labels.append(Path(ref.get('selector','entrada') if 'stage' in ref else asset_names[ref['asset']]).stem)
            if stage['operation']=='consensus' and bindings.get('base_input_path')==bindings.get('base_vina_path') and 'base_vina_path' in bindings:
                variant['bindings']['base_input_path']=copy.deepcopy(variant['bindings']['base_vina_path'])
            emitted=True
            yield variant,' + '.join(labels) or stage['name']
        if not emitted and stage['operation']=='redocking':
            raise ValueError('Selecione arquivos de estruturas correspondentes aos pares de redocking configurados.')
        if not emitted and stage['operation']=='docking_dock6':
            raise ValueError('Selecione receptores, compostos e poses Vina correspondentes para executar DOCK6.')
        if not emitted and stage['operation']=='consensus':
            raise ValueError('Selecione resultados Vina e DOCK6 com os mesmos identificadores de receptor e composto.')

    def _execute_stage(self, run_id, project_id, user_id, item, stages, run_dir, results, state):
        stage=execution_stage(item['configuration'],item.get('input_mode'))
        if stage.get('input_processing')!='individual' or stage.get('provided_results') or stage['operation']=='import_results':
            return self._execute_single(run_id,project_id,user_id,item,stages,run_dir,results,state)
        artifacts=[];inputs={};excluded=0
        item['batches']=[]
        for index,(variant,label) in enumerate(self._variants(project_id,stage,results),1):
            if state['event'].is_set():raise InterruptedError('Execução cancelada.')
            slug=re.sub(r'[^A-Za-z0-9_.+-]','_',label)[:90].strip('.') or 'entrada'
            directory=run_dir/stage['id']/'individual'/f'{index:03d}_{slug}'
            item['input_files']={}
            item['current_batch']={'index':index,'label':label}
            item['progress']=None
            try:
                produced=self._execute_single(run_id,project_id,user_id,item,stages,run_dir,results,state,variant,directory)
            except (AccessDenied,InterruptedError):raise
            except Exception as exc:raise RuntimeError(f'Arquivo {index} ({label}): {exc}') from exc
            batch=dict(label=label,artifacts=produced,input_files=copy.deepcopy(item.get('input_files',{})),
                       log_path=item.get('log_path'),job_id=item.get('job_id'),configuration=variant)
            item['batches'].append(batch)
            for field,paths in batch['input_files'].items():inputs.setdefault(field,[]).extend(paths)
            excluded+=item.get('excluded_records',0)
            artifacts.extend(produced)
            item['artifacts']=list(artifacts)
        item['input_files']={f:list(dict.fromkeys(paths)) for f,paths in inputs.items()}
        item['excluded_records']=excluded
        item.pop('current_batch',None)
        return artifacts

    def _execute_single(self, run_id, project_id, user_id, item, stages, run_dir, results, state, stage=None, stage_dir=None):
        stage = stage or execution_stage(item['configuration'],item.get('input_mode'))
        stage_dir = stage_dir or run_dir / stage['id']
        stage_dir.mkdir(parents=True)
        if stage['operation'] == 'import_results' or stage.get('provided_results'):
            params=stage.get('provided_results') or stage['parameters']
            item['provided']=bool(stage.get('provided_results'))
            artifacts=self._import(project_id,params,stage_dir / 'artifacts')
            if stage['operation']=='admet':
                from .visualizations import egg_view,write_view
                import pandas as pd
                for path in list(artifacts):
                    if path.endswith('.csv'):
                        file=Path(path)
                        view=egg_view(pd.read_csv(file),title='ADMET · resultados fornecidos')
                        destination=file.with_suffix('.biomol-view.json')
                        write_view(destination,view);artifacts.append(str(destination))
            return artifacts
        params = self._resolve(project_id,user_id,stage,results,stage_dir/'inputs',item,pipeline=execution_context(stages))
        if stage['operation']=='redocking':
            import sys
            from .redocking_config import validate_redocking_tools
            tool_path=str(Path(self.worker_python or sys.executable).parent)+os.pathsep+os.environ.get('PATH','')
            validate_redocking_tools(params.get('prepare_complex',True),tool_path)
        resource = materialize_templates(stage_dir / 'resources',stage.get('templates',{}))
        manager = JobManager(AppConfig(workspace=stage_dir, cpu_workers=self.cpu_workers,
            worker_python=self.worker_python, job_timeout=self.job_timeout, resource_dir=resource))
        try:
            with self._lock:
                state['manager'] = manager
                if state['event'].is_set():
                    raise InterruptedError('Execução cancelada.')
                with log_context(project_id=project_id, run_id=run_id, stage_id=stage['id']):
                    job = manager.submit(stage['operation'],params)
                state['job_id'] = job['id']
                item['log_path'] = job['log_path']
                item['job_id'] = job['id']
                self._update(run_id,'running',stages)
            while True:
                self.store._require_user(user_id,project_id,'editor')
                job = manager.get(job['id'])
                item['log_path'] = job['log_path']
                if job.get('progress') != item.get('progress'):
                    item['progress'] = job.get('progress')
                    self._update(run_id,'running',stages)
                if state['event'].is_set():
                    manager.cancel(job['id'])
                if job['status'] in TERMINAL:
                    break
                state['event'].wait(0.3)
            if job['status'] == 'cancelled':
                raise InterruptedError('Execução cancelada.')
            if job['status'] != 'succeeded':
                failure = RuntimeError(job['error'] or 'O estágio falhou. Consulte o log.')
                if (job.get('result') or {}).get('diagnostic'):
                    failure.error_code = job['result']['diagnostic']['error_code']
                    failure.action = job['result']['diagnostic']['action']
                raise failure
            from .molecule_quality import REPORT_NAME
            reports=list((stage_dir/'inputs').rglob(REPORT_NAME))
            item['excluded_records']=job['result'].get('details',{}).get('excluded_records',0)+sum(
                json.loads(path.read_text()).get('excluded_records',0) for path in reports)
            return job['result']['artifacts']+[str(p) for p in reports]
        finally:
            manager.close()
            with self._lock:
                state['manager'] = None
                state['job_id'] = None

    def _cache_key(self,project_id,user_id,stage,results,software,pipeline=None):
        if stage.get('input_processing')=='individual' and not stage.get('provided_results') and stage['operation']!='import_results':
            keys=[self._cache_key(project_id,user_id,variant,results,software,pipeline) for variant,_ in self._variants(project_id,stage,results)]
            return stage_key(stage,{'individual_keys':keys},{},software)
        parameters=stage.get('provided_results') or (stage['parameters'] if stage['operation']=='import_results' else self._resolve(project_id,user_id,stage,results,cache_only=True,pipeline=pipeline))
        inputs=input_manifests(stage,results)
        asset_ids=asset_references(stage)
        for asset_id in sorted(asset_ids):
            with self.store.connect() as db:
                asset=db.execute('SELECT path FROM assets WHERE project_id=? AND id=?',(project_id,asset_id)).fetchone()
            if asset is None:
                raise AccessDenied('Arquivo não autorizado.')
            inputs['asset:'+asset_id]=path_manifest(self.store.scoped_path(project_id,asset['path']))
        return stage_key(stage,parameters,inputs,software)

    def _cached_stage(self,project_id,user_id,stage,key,software):
        with self.store.connect() as db:
            previous=db.execute('SELECT id,stages FROM runs WHERE project_id=? ORDER BY created DESC',(project_id,)).fetchall()
        for run_id,encoded in previous:
            snapshot=json.loads(encoded)
            old=next((s for s in snapshot if s['id']==stage['id'] and s['status']=='succeeded' and s.get('artifacts')),None)
            if old is None:
                continue
            if old.get('outputs_removed'):
                continue
            try:
                for path in old['artifacts']:
                    self.store.scoped_path(project_id,path)
                manifest=old.get('artifact_manifest') or artifact_manifest(old['artifacts'])
                if set(manifest)!={str(Path(p).resolve()) for p in old['artifacts']}:
                    continue
                if not manifest_matches(manifest):
                    continue
                candidate_key=old.get('cache_key')
                if candidate_key and candidate_key!=key:
                    # The user's reuse choice can retain outputs across app
                    # upgrades or a different launcher. Verify the old request
                    # against its original provenance before comparing it with
                    # the current scientific request.
                    prior_results={s['id']:s['artifacts'] for s in snapshot if s['status']=='succeeded'}
                    configuration=execution_stage(old['configuration'],old.get('input_mode'))
                    context=execution_context([s for s in snapshot if s.get('configuration')])
                    provenance=old.get('cache_software')
                    if not provenance and old.get('software_version'):
                        provenance=old['software_version']+':'+software.split(':',1)[1]
                    if not provenance or self._cache_key(project_id,user_id,configuration,prior_results,provenance,pipeline=context)!=candidate_key:
                        continue
                    candidate_key=self._cache_key(project_id,user_id,configuration,prior_results,software,pipeline=context)
                if not candidate_key:
                    # Upgrade successful results produced before cache metadata
                    # existed, retaining their configuration and scoped inputs.
                    prior_results={s['id']:s['artifacts'] for s in snapshot if s['status']=='succeeded'}
                    # Legacy runs have no saved input hashes. Do not adopt their
                    # results if input files were edited after they completed.
                    configuration=old['configuration']
                    sources=stage_dependencies(configuration)
                    input_paths=[self.store.scoped_path(project_id,p) for source in sources for p in prior_results[source]]
                    for field in PATH_FIELDS:
                        if configuration['parameters'].get(field) and field not in configuration.get('bindings',{}):
                            path=self.store.scoped_path(project_id,configuration['parameters'][field])
                            input_paths.extend([path] if path.is_file() else [p for p in path.rglob('*') if p.is_file()])
                    for asset_id in asset_references(configuration):
                        with self.store.connect() as db:
                            asset=db.execute('SELECT path FROM assets WHERE project_id=? AND id=?',(project_id,asset_id)).fetchone()
                        if asset is None:
                            raise ValueError('Arquivo indisponível.')
                        input_paths.append(self.store.scoped_path(project_id,asset['path']))
                    if input_paths and (not old.get('finished_at') or any(p.stat().st_mtime>old['finished_at'] for p in input_paths)):
                        continue
                    candidate_key=self._cache_key(project_id,user_id,execution_stage(old['configuration'],old.get('input_mode')),prior_results,software,
                        pipeline=execution_context([s for s in snapshot if s.get('configuration')]))
                if candidate_key==key:
                    return dict(old,artifact_manifest=manifest,reused_from_run=run_id)
            except (OSError,ValueError,KeyError):
                continue
        return None

    def existing_results(self, token, project_id):
        """Report intact results persisted in this project's execution history."""
        self.store.project(token,project_id,'editor')
        found=set()
        for run in self.store.list_runs(token,project_id):
            for item in run['stages']:
                if item['id'] in found or item['status']!='succeeded' or not item.get('artifacts') or item.get('outputs_removed'):
                    continue
                try:
                    for path in item['artifacts']:self.store.scoped_path(project_id,path)
                    if manifest_matches(item.get('artifact_manifest') or artifact_manifest(item['artifacts'])):
                        found.add(item['id'])
                except (OSError,ValueError):continue
        return len(found)

    @staticmethod
    def _reuse_configuration(stage, previous):
        """Restore confirmed file choices only while the scientific setup matches."""
        parameters=dict(stage.get('parameters',{}));old_parameters=dict(previous.get('parameters',{}))
        if stage['operation']=='docking_dock6':
            parameters.pop('dock6_app_path',None);old_parameters.pop('dock6_app_path',None)
        if parameters!=old_parameters or any(stage.get(k,{})!=previous.get(k,{}) for k in ('templates','provided_results')):
            return None
        if stage['operation']!=previous['operation'] or stage_dependencies(stage)!=stage_dependencies(previous):return None
        bindings=stage.get('bindings',{});old_bindings=previous.get('bindings',{})
        if set(bindings)!=set(old_bindings):return None
        if stage.get('input_processing') and stage['input_processing']!=previous.get('input_processing','merge') and bindings:return None
        configuration=copy.deepcopy(stage)
        for field,group in bindings.items():
            refs=input_sources(group);old_refs=input_sources(old_bindings[field])
            identities=lambda values:{(r.get('stage'),r.get('asset')) for r in values}
            if identities(refs)!=identities(old_refs):return None
            if all('stage' in r and r.get('selector','auto')=='auto' for r in refs):
                configuration['bindings'][field]=copy.deepcopy(old_bindings[field])
            elif refs!=old_refs:return None
        if 'input_processing' not in configuration and 'input_processing' in previous:
            configuration['input_processing']=previous['input_processing']
        return configuration

    def _reuse_stage(self, project_id, user_id, stage, results, software, stages):
        with self.store.connect() as db:
            histories=db.execute('SELECT stages FROM runs WHERE project_id=? ORDER BY created DESC',(project_id,)).fetchall()
        for row in histories:
            for old in json.loads(row['stages']):
                if old['id']!=stage['id'] or old['status']!='succeeded' or not old.get('artifacts'):continue
                configuration=self._reuse_configuration(stage,old['configuration'])
                if configuration is None:continue
                try:
                    key=self._cache_key(project_id,user_id,execution_stage(configuration,'curated'),results,software,pipeline=execution_context(stages))
                    cached=self._cached_stage(project_id,user_id,configuration,key,software)
                except (OSError,ValueError,KeyError):continue
                if cached:return configuration,key,cached
        return None

    def _run(self, run_id, project_id, user_id, stages):
        state = self._active[run_id]
        results = {s['id']:s['artifacts'] for s in stages if s['status']=='succeeded'}
        try:
            software=implementation_digest() + ':' + str(self.worker_python or os.sys.executable)
            run_dir = self.store.project_dir(project_id) / 'runs' / run_id
            run_dir.mkdir(parents=True,exist_ok=True)
            self._update(run_id,'running',stages)
            by_id = {item['id']:item for item in stages}
            for item in stages:
                if state['event'].is_set():
                    raise InterruptedError('Execução cancelada.')
                self.store._require_user(user_id,project_id,'editor')
                if item['status'] in ('skipped','succeeded','failed'):
                    continue
                stage = item['configuration']
                dependencies = stage_dependencies(stage)
                item['input_mode']='curated'
                item['requires_curation']=bool(dependencies) and not (stage['operation']=='redocking' and stage['parameters'].get('pdb_codes'))
                unavailable = [by_id[source]['name'] for source in dependencies if by_id[source]['status'] != 'succeeded']
                if unavailable:
                    item.update(status='skipped', error='Esta etapa depende de blocos que não concluíram: ' + ', '.join(sorted(unavailable)))
                    self._update(run_id,'running',stages)
                    continue
                reusable=self._reuse_stage(project_id,user_id,stage,results,software,stages) if state.get('reuse_results',True) else None
                if reusable:
                    stage,key,cached=reusable
                    item.update(configuration=stage,curation_confirmed=True)
                if item.get('requires_curation') and not item.get('curation_confirmed') and not stage.get('provided_results'):
                    item.update(status='awaiting_input',error=None)
                    message=f'Configure o bloco “{item["name"]}” e escolha os arquivos de entrada para continuar o pipeline.'
                    with self._lock:
                        if state['event'].is_set():raise InterruptedError('Execução cancelada.')
                        self._update(run_id,'awaiting_input',stages,message)
                        self._active.pop(run_id,None)
                    return
                stage = execution_stage(stage,item.get('input_mode'))
                item['status'] = 'running'
                item['started_at'] = time.time()
                self._update(run_id,'running',stages)
                try:
                    if not reusable:
                        key=self._cache_key(project_id,user_id,stage,results,software,pipeline=execution_context(stages))
                        cached=self._cached_stage(project_id,user_id,stage,key,software) if state.get('reuse_results',True) else None
                    if cached:
                        artifacts=cached['artifacts']
                        item.update(reused=True,reused_from_run=cached['reused_from_run'],
                                    artifact_manifest=cached['artifact_manifest'])
                        if cached.get('log_path'):
                            item['log_path']=cached['log_path']
                    else:
                        artifacts = self._execute_stage(run_id,project_id,user_id,item,stages,run_dir,results,state)
                        item['artifact_manifest']=artifact_manifest(artifacts)
                    item['cache_key']=key
                    item['cache_software']=software
                    item['software_version']=software.split(':',1)[0]
                    if cached and cached.get('input_files'):item['input_files']=cached['input_files']
                    if cached and cached.get('provided'):item['provided']=True
                    if cached and cached.get('excluded_records'):item['excluded_records']=cached['excluded_records']
                    if cached and cached.get('batches'):item['batches']=cached['batches']
                    if state['event'].is_set():
                        raise InterruptedError('Execução cancelada.')
                    self.store._require_user(user_id,project_id,'editor')
                except (AccessDenied,InterruptedError):
                    raise
                except Exception as exc:
                    get_logger('backend').exception('Pipeline stage failed: %s', exc, extra={'event': 'stage.failed',
                        'run_id': run_id, 'project_id': project_id, 'stage_id': stage['id'], 'operation': stage['operation'], **diagnose_exception(exc)})
                    item.update(status='failed',error=str(exc),finished_at=time.time())
                    self._update(run_id,'running',stages)
                    continue
                item.update(status='succeeded',artifacts=artifacts,finished_at=time.time())
                results[stage['id']] = artifacts
                self._update(run_id,'running',stages)
            failed = [item for item in stages if item['status'] == 'failed']
            if failed:
                error = '\n'.join(f'{item["name"]}: {item["error"]}' for item in failed)
                self._update(run_id,'failed',stages,error)
            else:
                self._update(run_id,'succeeded',stages)
        except Exception as exc:
            get_logger('backend').exception('Pipeline failed: %s', exc, extra={'event': 'pipeline.failed', 'run_id': run_id, 'project_id': project_id, **diagnose_exception(exc)})
            cancelled = state['event'].is_set() or isinstance(exc,InterruptedError)
            for item in stages:
                if item['status'] in ('running','awaiting_input'):
                    item['status'] = 'cancelled' if cancelled else 'failed'
                    item['error'] = str(exc)
                    item['finished_at'] = time.time()
                elif item['status'] == 'queued':
                    item.update(status='skipped',error='A execução foi cancelada antes desta etapa.' if cancelled else
                                'A execução foi interrompida antes desta etapa: ' + str(exc))
            self._update(run_id,'cancelled' if cancelled else 'failed',stages,str(exc))
        finally:
            with self._lock:
                if self._active.get(run_id) is state:self._active.pop(run_id,None)

    def resume(self, token, run_id, configuration):
        """Confirm one pending stage without rerunning completed retrievals."""
        configuration=execution_stage(configuration,'curated')
        if configuration.get('operation')=='redocking' and not configuration.get('provided_results'):
            from .redocking_config import validate_pairs
            validate_pairs(configuration['parameters'].get('pdb_codes'),configuration['parameters'].get('preparation_pairs') or {})
        run = self.store.get_run(token,run_id)
        self.store.project(token,run['project_id'],'editor')
        user = self.store.user(token)
        with self._lock:
            if self._closed:raise RuntimeError('A aplicação está sendo encerrada.')
            run = self.store.get_run(token,run_id)
            if run['status']!='awaiting_input':raise ValueError('Esta execução não está aguardando seleção de arquivos.')
            stages = run['stages']
            pending = next(s for s in stages if s['status']=='awaiting_input')
            if configuration.get('id')!=pending['id'] or configuration.get('operation')!=pending['operation']:
                raise ValueError('Configure o bloco indicado pela execução.')
            if not configuration.get('enabled',True) or configuration.get('provided_results'):
                raise ValueError('Mantenha este bloco ativo e selecione suas entradas para executá-lo.')
            completed={s['id'] for s in stages if s['status']=='succeeded'}
            if not stage_dependencies(configuration).issubset(completed):
                raise ValueError('Escolha arquivos de blocos já concluídos nesta execução.')
            candidate=[configuration if s is pending else s['configuration'] for s in stages]
            validate_pipeline(candidate)
            context=[dict(c,_execution_batches=s.get('batches',[])) for c,s in zip(candidate,stages)]
            for group in configuration.get('bindings',{}).values():
                if any('stage' in ref and ref.get('selector','auto')=='auto' for ref in input_sources(group)):
                    raise ValueError('Escolha explicitamente os arquivos de entrada deste bloco antes de continuar.')
            self._validate_parameters(run['project_id'],configuration,partial=True)
            results={s['id']:s['artifacts'] for s in stages if s['status']=='succeeded'}
            if configuration.get('input_processing')=='individual':
                for variant,_ in self._variants(run['project_id'],configuration,results):
                    self._resolve(run['project_id'],user['id'],variant,results,cache_only=True,pipeline=context)
            else:self._resolve(run['project_id'],user['id'],configuration,results,cache_only=True,pipeline=context)
            pending.update(configuration=copy.deepcopy(configuration),name=configuration['name'],input_mode='curated',
                           status='queued',curation_confirmed=True,error=None)
            state={'event':threading.Event(),'manager':None,'job_id':None,'reuse_results':pending.get('reuse_results',True)}
            self._active[run_id]=state
            try:
                self._update(run_id,'queued',stages)
                self._executor.submit(self._run,run_id,run['project_id'],user['id'],stages)
            except Exception:
                self._active.pop(run_id,None)
                pending.update(status='awaiting_input',curation_confirmed=False)
                self._update(run_id,'awaiting_input',stages,'Não foi possível continuar. Tente novamente.')
                raise
        return self.store.get_run(token,run_id)

    def cancel(self, token, run_id):
        run = self.store.get_run(token,run_id)
        self.store.project(token,run['project_id'],'editor')
        with self._lock:
            state = self._active.get(run_id)
            if state:
                state['event'].set()
                if state['manager'] and state['job_id']:
                    state['manager'].cancel(state['job_id'])
            else:
                run=self.store.get_run(token,run_id)
                if run['status']=='awaiting_input':
                    for item in run['stages']:
                        if item['status']=='awaiting_input':item.update(status='cancelled',error='Execução cancelada.',finished_at=time.time())
                        elif item['status']=='queued':item.update(status='skipped',error='A execução foi cancelada antes desta etapa.')
                    self._update(run_id,'cancelled',run['stages'],'Execução cancelada.')

    def close(self):
        with self._lock:
            self._closed = True
            for state in self._active.values():
                state['event'].set()
        self._executor.shutdown(wait=True)
        self._owner.close()
