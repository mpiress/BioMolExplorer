"""Authorized result-file management and dynamic ADMET datasets."""
import csv
import io
import json
import math
import time
from pathlib import Path
from uuid import uuid4

from .workspace import AccessDenied
from .visualizations import MAX_VIEW_BYTES, MAX_GRAPH_BYTES, SUFFIX, egg_view, load_view

FILTERS = {'all': 'Todos os compostos da tabela', 'BBB+': 'Apenas BBB+',
           'BBB-': 'Apenas BBB−', 'HIA+': 'Apenas HIA+', 'HIA-':'Apenas HIA−',
           'selected':'Atendem à hipótese', 'rejected':'Não atendem à hipótese'}


class ResultFiles:
    def __init__(self, store):
        self.store = store

    def origins(self, token, pid, rid):
        from .result_origins import ResultOrigins
        def load(run_id):
            try:run=self.store.get_run(token,run_id)
            except ValueError:return None
            if run['project_id']!=pid:raise AccessDenied('Resultado de outro projeto.')
            return run
        run=load(rid)
        return ResultOrigins(run['stages'],rid,load_run=load)

    def stage(self, token, pid, rid, sid):
        run = self.store.get_run(token, rid)
        if run['project_id'] != pid:
            raise AccessDenied('Resultado de outro projeto.')
        stage = next((s for s in run['stages'] if s['id'] == sid), None)
        if stage is None:
            raise ValueError('Etapa não encontrada nesta execução.')
        return stage

    def files(self, token, pid, rid, sid):
        stage = self.stage(token, pid, rid, sid)
        from .artifacts import records
        from .artifact_naming import publication_metadata
        catalog={r['path']:r for r in records(stage,rid)}
        origins=self.origins(token,pid,rid)
        def information(path):
            record=catalog[str(path)]
            source=origins.describe(path,stage['id'],rid)
            return dict(record,**source,size=path.stat().st_size,description='Origem: '+source['origin']+' · '+record['selector'],
                        provenance='Origem: '+source['source_trace']+'\nExecução '+str(record['run_id'] or rid))
        return [information(p)
                for name in stage.get('artifacts', [])
                if (p := self.store.scoped_path(pid, name)).is_file()
                and not (stage['operation']=='retrieve_structures' and p.name in ('pdb_codes.csv','retrieval_report.json'))
                and not (stage['operation'] == 'admet' and p.suffix.lower() == '.json' and p.name!='molecule_exclusions.json'
                    and not publication_metadata(p).get('legacy_name',p.name).endswith('_admet_report.json'))]

    def preparation_failures(self,token,pid,rid,sid,filename):
        from .preparation_failures import REPORT_NAME
        stage=self.stage(token,pid,rid,sid);path=self.store.scoped_path(pid,filename)
        if path.name!=REPORT_NAME or str(path) not in stage.get('artifacts',[]):
            raise AccessDenied('Relatório não autorizado para esta etapa.')
        data=self.store.read_file(token,pid,str(path),MAX_VIEW_BYTES)
        return list(csv.DictReader(io.StringIO(data.decode('utf-8-sig'))))

    @staticmethod
    def _pdb_review_signature(path,ligands):
        from hashlib import sha256
        stat=path.stat()
        rows=sorted(tuple(str(r.get(key,'')) for key in ('PDB_CODE','LIGAND','RESNUM','CHAIN','RESOLUTION')) for r in ligands)
        return sha256(json.dumps([stat.st_size,stat.st_mtime_ns,rows]).encode()).hexdigest()

    def pdb_review_statuses(self,token,pid,rid):
        """Read each metadata table once; never parse every selected structure."""
        run=self.store.get_run(token,rid)
        if run['project_id']!=pid:raise AccessDenied('Resultado de outro projeto.')
        result={};tables={}
        for stage in run['stages']:
            if stage['operation']!='retrieve_structures':continue
            artifacts=set(stage.get('artifacts',[]))
            for value in artifacts:
                path=self.store.scoped_path(pid,value);metadata=path.parent/'pdb_codes.csv'
                if path.suffix.lower()!='.pdb' or str(metadata) not in artifacts or not path.is_file():continue
                if metadata not in tables:
                    data=self.store.read_file(token,pid,str(metadata),MAX_VIEW_BYTES)
                    grouped={}
                    for row in csv.DictReader(io.StringIO(data.decode('utf-8-sig'))):
                        grouped.setdefault(row['PDB_CODE'].upper(),[]).append(row)
                    tables[metadata]=grouped
                ligands=tables[metadata].get(path.stem.upper(),[])
                review=stage.get('pdb_reviews',{}).get(str(path),{})
                result[str(path)]=dict(ligands=ligands,reviewed=bool(review and
                    review.get('signature')==self._pdb_review_signature(path,ligands)),reviewed_at=review.get('reviewed_at'))
        return result

    def pdb_ligands(self, token, pid, rid, sid, filename):
        stage=self.stage(token,pid,rid,sid)
        path=self.store.scoped_path(pid,filename)
        if stage['operation']!='retrieve_structures' or str(path) not in stage.get('artifacts',[]) or path.suffix.lower()!='.pdb':
            raise AccessDenied('Estrutura não autorizada para esta etapa.')
        metadata=path.parent/'pdb_codes.csv'
        if str(metadata) not in stage.get('artifacts',[]):raise AccessDenied('Metadados PDB não autorizados.')
        data=self.store.read_file(token,pid,str(metadata),MAX_VIEW_BYTES)
        reader=csv.DictReader(io.StringIO(data.decode('utf-8-sig')))
        from hashlib import sha256
        from .pdb_ligands import analyze_ligands
        text=self.store.read_file(token,pid,str(path),MAX_VIEW_BYTES).decode('utf-8-sig')
        analysis=analyze_ligands(text,path.stem)
        available=[dict(row,RESNUM=str(row['RESNUM'])) for row in analysis['candidates']]
        available.sort(key=lambda r:(r['LIGAND'],int(r['RESNUM']),r['CHAIN'],r['ICODE']))
        ligands=[row for row in reader if row['PDB_CODE'].upper()==path.stem.upper()]
        review=stage.get('pdb_reviews',{}).get(str(path),{})
        return {'pdb_id':path.stem.upper(),'revision':sha256(data).hexdigest(),'resolution':analysis['resolution'],
                'ligands':ligands,'reviewed':bool(review and review.get('signature')==self._pdb_review_signature(path,ligands)),
                'available_ligands':available,
                'recommended_ligands':[r for r in available if r['RECOMMENDED']],
                'triage_available':str(path.parent/'ligand_candidates.csv') in stage.get('artifacts',[])}

    def set_pdb_ligands(self,token,pid,rid,sid,filename,records,revision):
        """Replace one structure's curated ligands, preserving other PDB rows."""
        import re
        from .storage import write_text
        from .stage_cache import artifact_manifest
        path=self.store.scoped_path(pid,filename)
        metadata=path.parent/'pdb_codes.csv'
        backup=None
        changes={}
        try:
            with self.store.change(token,pid,'editor','curate_pdb_ligands','Ligantes revisados: '+path.stem) as db:
                if db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running')",(pid,)).fetchone():
                    raise ValueError('Aguarde a execução terminar ou pausar antes de editar os ligantes.')
                context=self.pdb_ligands(token,pid,rid,sid,filename)
                if context['revision']!=revision:raise ValueError('Os ligantes foram alterados. Reabra a lista antes de salvar.')
                if not isinstance(records,list) or len(records)>10000:raise ValueError('Lista de ligantes inválida.')
                residues=set()
                for line in self.store.read_file(token,pid,str(path),MAX_VIEW_BYTES).decode().splitlines():
                    if line.startswith(('ATOM  ','HETATM')):
                        residues.add((line[17:20].strip().upper(),int(line[22:26]),line[21:22].strip()))
                curated=[];seen=set()
                with metadata.open(encoding='utf-8-sig',newline='') as stream:
                    reader=csv.DictReader(stream);fields=reader.fieldnames;existing=list(reader)
                for record in records:
                    ligand=str(record.get('LIGAND','')).strip().upper()
                    chain=str(record.get('CHAIN','')).strip()
                    try:number=int(record.get('RESNUM',''))
                    except (TypeError,ValueError):raise ValueError('Informe um número de resíduo inteiro.')
                    if not re.fullmatch(r'[A-Z0-9]{1,5}',ligand) or not re.fullmatch(r'[A-Za-z0-9]',chain):
                        raise ValueError('Informe código do ligante e cadeia válidos.')
                    key=(ligand,number,chain)
                    if key not in residues:raise ValueError('Este ligante/resíduo/cadeia não foi encontrado na estrutura baixada.')
                    if any(r.get('ICODE') for r in context['available_ligands']
                           if (r['LIGAND'],int(r['RESNUM']),r['CHAIN'])==key):
                        raise ValueError('Este resíduo possui código de inserção, não suportado na preparação atual. Selecione outra referência ou normalize a estrutura externamente.')
                    if key in seen:continue
                    seen.add(key)
                    row={k:'' for k in fields}
                    row.update(PDB_CODE=context['pdb_id'],LIGAND=ligand,RESNUM=str(number),CHAIN=chain)
                    if 'RESOLUTION' in fields:row['RESOLUTION']='' if context['resolution'] is None else str(context['resolution'])
                    curated.append(row)
                from .redocking_selection import sync_curated_pairs
                pipeline=json.loads(db.execute('SELECT pipeline FROM projects WHERE id=?',(pid,)).fetchone()[0])
                changes=sync_curated_pairs(pipeline,sid,path,
                    self.stage(token,pid,rid,sid)['artifacts'],
                    [row for row in existing if row['PDB_CODE'].upper()==context['pdb_id']],curated)
                if changes['updated_redocking_stages']:
                    db.execute('UPDATE projects SET pipeline=?,revision=revision+1,updated=? WHERE id=?',
                               (json.dumps(pipeline),time.time(),pid))
                rows=[row for row in existing if row['PDB_CODE'].upper()!=context['pdb_id']]+curated
                stream=io.StringIO(newline='');writer=csv.DictWriter(stream,fieldnames=fields)
                writer.writeheader();writer.writerows(rows)
                backup=metadata.read_text()
                write_text(stream.getvalue(),metadata)
                for run in db.execute('SELECT id,stages FROM runs WHERE project_id=?',(pid,)).fetchall():
                    stages=json.loads(run[1]);changed=False
                    for item in stages:
                        if str(metadata) in item.get('artifacts',[]):
                            item['artifact_manifest']=artifact_manifest(item['artifacts'])
                            item['curated_at']=time.time();changed=True
                            if str(path) in item.get('artifacts',[]):
                                item.setdefault('pdb_reviews',{})[str(path)]={
                                    'signature':self._pdb_review_signature(path,curated),'reviewed_at':time.time()}
                    if changed:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),run[0]))
        except Exception:
            if backup is not None:write_text(backup,metadata)
            raise
        return dict(self.pdb_ligands(token,pid,rid,sid,filename),**changes)

    def remove(self, token, pid, rid, sid, filename):
        """Audit deletion, invalidate producer caches and restore bytes on failure."""
        backup = path = None
        try:
            with self.store.change(token, pid, 'editor', 'remove_file', 'Arquivo removido: '+Path(filename).name) as db:
                if db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')", (pid,)).fetchone():
                    raise ValueError('Aguarde a execução terminar antes de remover arquivos.')
                stage = self.stage(token, pid, rid, sid)
                path = self.store.scoped_path(pid, filename)
                if str(path) not in stage.get('artifacts', []) or not path.is_file():
                    raise AccessDenied('Arquivo não autorizado para esta etapa.')
                backup_dir = self.store.project_dir(pid)/'.history'/'.removed'
                backup_dir.mkdir(parents=True, exist_ok=True)
                backup = backup_dir/uuid4().hex
                path.replace(backup)
                for row in db.execute('SELECT id,stages FROM runs WHERE project_id=?', (pid,)).fetchall():
                    stages = json.loads(row[1]); changed = False
                    for item in stages:
                        if str(path) in item.get('artifacts', []):
                            item['artifacts'].remove(str(path))
                            item['outputs_removed'] = True
                            item['curated_at'] = time.time()
                            changed = True
                    if changed:
                        db.execute('UPDATE runs SET stages=? WHERE id=?', (json.dumps(stages), row[0]))
        except Exception:
            if backup and backup.exists():
                backup.replace(path)
            raise
        self._cleanup_backup(backup)

    def remove_all(self, token, pid, rid, sid):
        """Remove a stage's artifacts atomically, including shared cache references."""
        backups=[]
        try:
            with self.store.change(token,pid,'editor','remove_stage_results','Resultados da etapa removidos') as db:
                if db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')",(pid,)).fetchone():
                    raise ValueError('Aguarde a execução terminar antes de remover arquivos.')
                stage=self.stage(token,pid,rid,sid)
                paths={str(self.store.scoped_path(pid,p)) for p in stage.get('artifacts',[])}
                backup_dir=self.store.project_dir(pid)/'.history'/'.removed'
                backup_dir.mkdir(parents=True,exist_ok=True)
                for filename in sorted(paths):
                    path=Path(filename)
                    if path.is_file():
                        backup=backup_dir/uuid4().hex
                        path.replace(backup);backups.append((path,backup))
                for row in db.execute('SELECT id,stages FROM runs WHERE project_id=?',(pid,)).fetchall():
                    stages=json.loads(row[1]);changed=False
                    for item in stages:
                        if paths.intersection(item.get('artifacts',[])):
                            item['artifacts']=[p for p in item['artifacts'] if p not in paths]
                            item['outputs_removed']=True;item['curated_at']=time.time();changed=True
                    if changed:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),row[0]))
        except Exception:
            for path,backup in reversed(backups):
                if backup.exists():backup.replace(path)
            raise
        for _,backup in backups:self._cleanup_backup(backup)
        return len(backups)

    def remove_project_results(self,token,pid):
        """Atomically remove every computed result while preserving user inputs."""
        backups=[]
        try:
            with self.store.change(token,pid,'editor','remove_project_results','Todos os resultados do projeto removidos') as db:
                if db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')",(pid,)).fetchone():
                    raise ValueError('Aguarde a execução terminar antes de remover arquivos.')
                rows=db.execute('SELECT id,stages FROM runs WHERE project_id=?',(pid,)).fetchall()
                decoded=[(run_id,json.loads(encoded)) for run_id,encoded in rows]
                # Imported/provided files belong to the project input library,
                # even when they are exposed as stage artifacts.
                result_paths={str(self.store.scoped_path(pid,name))
                    for _,stages in decoded for item in stages
                    if item.get('operation')!='import_results' and not item.get('provided')
                    and not item.get('configuration',{}).get('provided_results')
                    for name in item.get('artifacts',[])}
                protected={str(self.store.scoped_path(pid,row[0])) for row in
                    db.execute('SELECT path FROM assets WHERE project_id=?',(pid,)).fetchall()}
                protected.update(str(self.store.scoped_path(pid,name))
                    for _,stages in decoded for item in stages
                    if item.get('operation')=='import_results' or item.get('provided')
                    or item.get('configuration',{}).get('provided_results')
                    for name in item.get('artifacts',[]))
                paths=result_paths-protected
                backup_dir=self.store.project_dir(pid)/'.history'/'.removed'
                backup_dir.mkdir(parents=True,exist_ok=True)
                for filename in sorted(paths):
                    path=Path(filename)
                    if path.is_file():
                        backup=backup_dir/uuid4().hex
                        path.replace(backup);backups.append((path,backup))
                for run_id,stages in decoded:
                    changed=False
                    for item in stages:
                        computed=(item.get('operation')!='import_results' and not item.get('provided')
                            and not item.get('configuration',{}).get('provided_results'))
                        if computed and result_paths.intersection(item.get('artifacts',[])):
                            item['artifacts']=[p for p in item['artifacts'] if p not in result_paths]
                            item['outputs_removed']=True;item['curated_at']=time.time();changed=True
                    if changed:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),run_id))
        except Exception:
            for path,backup in reversed(backups):
                if backup.exists():backup.replace(path)
            raise
        for _,backup in backups:self._cleanup_backup(backup)
        return len(backups)

    @staticmethod
    def _cleanup_backup(backup):
        if backup is not None:
            try:
                backup.unlink(missing_ok=True)
            except OSError:
                # The audited change has committed; a cleanup failure must not
                # undo only its filesystem half. History excludes these backups.
                from .diagnostics import get_logger
                get_logger('backend').warning('Não foi possível limpar uma cópia temporária de arquivo removido.',exc_info=True)

    def admet_datasets(self, token, pid, rid, sid):
        stage = self.stage(token, pid, rid, sid)
        artifacts = [self.store.scoped_path(pid, p) for p in stage.get('artifacts', [])]
        result = []
        for path in artifacts:
            if path.suffix.lower() != '.csv' or not path.is_file():
                continue
            with path.open(encoding='utf-8-sig', newline='') as source:
                columns = set(csv.DictReader(source).fieldnames or [])
            if {'TPSA', 'WLOGP', 'molecule_chembl_id', 'canonical_smiles'} <= columns:
                result.append({'path': str(path), 'name': path.name})
        origins=self.origins(token,pid,rid)
        for item in result:
            item.update(origins.describe(item['path'],sid,rid))
            item['label']=item['display_name']
        return sorted(result,key=lambda item:(item['name'].endswith('_evaluated.csv'),len(item['name']),item['name']))

    @staticmethod
    def _admet_source_stem(path):
        import re
        from .artifact_naming import publication_metadata
        stem=Path(publication_metadata(path).get('legacy_name',Path(path).name)).stem
        return re.sub(r'_(?:evaluated|selected|rejected|alerts|BBB[+-]|HIA[+-]|rules_(?:lipinski-veber|lipinski|veber|fragment))$','',stem)

    def admet_alert_page(self,token,pid,rid,sid,filename,offset=0,limit=25,query=''):
        """One authorized, bounded review table per ADMET input, including history."""
        from .admet_config import ALERT_CATALOGS
        from .artifact_naming import output_role
        if type(offset) is not int or offset<0 or type(limit) is not int or not 1<=limit<=100:
            raise ValueError('Página de alertas inválida.')
        datasets=self.admet_datasets(token,pid,rid,sid)
        if filename not in {d['path'] for d in datasets}:raise AccessDenied('Tabela ADMET não autorizada.')
        current=Path(filename);stem=self._admet_source_stem(current)
        related=[Path(d['path']) for d in datasets if Path(d['path']).parent==current.parent and self._admet_source_stem(d['path'])==stem]
        def priority(path):
            role=output_role('admet',path,{})
            return 0 if role.startswith('alerts') else 1 if role=='evaluated' else 2 if role=='compounds' else 3
        path=min(related,key=priority)
        source=io.StringIO(self.store.read_file(token,pid,str(path),MAX_VIEW_BYTES+1).decode('utf-8-sig'))
        reader=csv.DictReader(source);fields=reader.fieldnames or []
        catalogs=[name for name in ALERT_CATALOGS if 'Alert_'+name in fields]
        legacy=not catalogs
        if legacy:
            config=self.stage(token,pid,rid,sid).get('configuration',{}).get('parameters',{})
            catalogs=[name for name in config.get('alert_catalogs',['PAINS','BRENK']) if name in ALERT_CATALOGS]
        rows=[];total=0;query=query.casefold().strip()
        for record in reader:
            identifier=record.get('molecule_chembl_id','');details=record.get('StructuralAlerts','')
            if query and query not in (identifier+' '+details+' '+record.get('canonical_smiles','')).casefold():continue
            if offset<=total<offset+limit:
                flags={name:(name+': ' in details if legacy else str(record.get('Alert_'+name,'')).lower() in ('true','1')) for name in catalogs}
                rows.append(dict(id=identifier,smiles=record.get('canonical_smiles',''),flags=flags,details=details))
            total+=1
        origins=self.origins(token,pid,rid)
        return dict(rows=rows,catalogs=catalogs,total=total,offset=offset,limit=limit,path=str(path),
                    label=origins.describe(str(path),sid,rid)['display_name'])

    def graph_datasets(self,token,pid,rid,sid):
        stage=self.stage(token,pid,rid,sid)
        result=[]
        for filename in stage.get('artifacts',[]):
            path=self.store.scoped_path(pid,filename)
            if not path.name.endswith(SUFFIX) or not path.is_file():continue
            if path.stat().st_size>MAX_GRAPH_BYTES:
                result.append({'path':str(path),'name':path.name+f' · acima do limite de {MAX_GRAPH_BYTES//(1024*1024)} MB','unavailable':True})
                continue
            model=load_view(self.store.read_file(token,pid,str(path),MAX_GRAPH_BYTES+1),max_bytes=MAX_GRAPH_BYTES)
            if model['kind']=='graph':
                result.append({'path':str(path),'name':model['title']})
        from collections import Counter
        titles=Counter(choice['name'] for choice in result)
        for index,choice in enumerate(result,1):
            if titles[choice['name']]>1:choice['name']=f"Análise {index:02d} · "+choice['name']
        origins=self.origins(token,pid,rid)
        for item in result:
            item.update(origins.describe(item['path'],sid,rid))
            item['label']=item['origin']+' — '+item['name']
        return result

    def graph_model(self,token,pid,rid,sid,filename):
        stage=self.stage(token,pid,rid,sid)
        path=self.store.scoped_path(pid,filename)
        if str(path) not in stage.get('artifacts',[]) or not path.name.endswith(SUFFIX):
            raise AccessDenied('Grafo não autorizado para esta etapa.')
        content=self.store.read_file(token,pid,str(path),MAX_GRAPH_BYTES+1)
        model=load_view(content,max_bytes=MAX_GRAPH_BYTES)
        if model['kind']!='graph':raise ValueError('Selecione um resultado de grafos.')
        model['dataset_origin']=self.origins(token,pid,rid).describe(path,sid,rid)
        source=model['dataset_origin']['source_label']
        if not model['title'].startswith(source):model['title']=source+' — '+model['title']
        # Opening old topology must never start an unbounded MCS search.
        import hashlib
        model['dashboard_context']={'project_id':pid,'run_id':rid,'stage_id':sid,'path':str(path),'sha256':hashlib.sha256(content).hexdigest()}
        return model

    def graph_sources(self,token,pid,rid):
        from .artifact_sources import records
        run=self.store.get_run(token,rid)
        if run['project_id']!=pid:raise AccessDenied('Resultado de outro projeto.')
        result=[]
        origins=self.origins(token,pid,rid)
        for stage in run['stages']:
            if stage.get('status')!='succeeded':continue
            for record in records(stage,rid):
                if Path(record['path']).suffix.lower()=='.csv':
                    item={key:record[key] for key in ('artifact_id','path','name','stage_id','stage_name','run_id','context')}
                    description=origins.describe(record['path'],stage['id'],rid)
                    item.update({key:description[key] for key in ('source_label','origin','display_name','source_trace')})
                    result.append(item)
        return result

    def graph_enrichment(self,token,pid,rid,sid,filename,source_id,expected_graph_sha256=None):
        from .molecular_network import enrich
        model=self.graph_model(token,pid,rid,sid,filename)
        if expected_graph_sha256 and expected_graph_sha256!=model['dashboard_context']['sha256']:
            raise ValueError('O resultado do grafo mudou. Atualize o dashboard antes de associar uma fonte.')
        source=next((s for s in self.graph_sources(token,pid,rid) if s['artifact_id']==source_id),None)
        if source is None:raise AccessDenied('Fonte molecular não autorizada para esta execução.')
        path=self.store.scoped_path(pid,source['path'])
        # Fingerprint sources can contain hundreds of MB of bits. Stream CSV
        # rows and retain only matching molecular metadata, never the bit column.
        import hashlib
        digest=hashlib.sha256()
        def lines():
            with path.open('rb') as stream:
                for index,raw in enumerate(stream):
                    digest.update(raw)
                    yield raw.decode('utf-8-sig' if index==0 else 'utf-8')
        result=enrich(model,lines(),source)
        result['source']=dict(source,content_sha256=digest.hexdigest())
        return result

    def graph_png(self,token,pid,rid,sid,filename):
        from caad.graph_results import report_png
        return report_png(self.graph_model(token,pid,rid,sid,filename))

    def graph_mcc_csv(self,token,pid,rid,sid,filename):
        model=self.graph_model(token,pid,rid,sid,filename)
        allowed=set(model['mcc']);records=[dict(n['properties'],molecule_chembl_id=n['id']) for n in model['nodes'] if n['id'] in allowed]
        fields=list(dict.fromkeys(['molecule_chembl_id','canonical_smiles','degree']+[k for row in records for k in row]))
        stream=io.StringIO(newline='');writer=csv.DictWriter(stream,fieldnames=fields)
        writer.writeheader();writer.writerows(records)
        return stream.getvalue().encode('utf-8')

    def remove_asset(self, token, pid, asset_id):
        backup=path=None
        try:
            with self.store.change(token,pid,'editor','remove_file','Arquivo de entrada removido') as db:
                if db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')",(pid,)).fetchone():
                    raise ValueError('Aguarde a execução terminar antes de remover arquivos.')
                from .bindings import asset_references
                pipeline=json.loads(db.execute('SELECT pipeline FROM projects WHERE id=?',(pid,)).fetchone()[0])
                if any(asset_id in asset_references(stage) for stage in pipeline):
                    raise ValueError('Este arquivo está associado a um bloco. Remova a associação nas configurações do bloco antes de excluir o arquivo.')
                row=db.execute('SELECT path FROM assets WHERE project_id=? AND id=?',(pid,asset_id)).fetchone()
                if row is None:raise AccessDenied('Arquivo não autorizado.')
                path=self.store.scoped_path(pid,row[0])
                if path.is_file():
                    backup_dir=self.store.project_dir(pid)/'.history'/'.removed';backup_dir.mkdir(parents=True,exist_ok=True)
                    backup=backup_dir/uuid4().hex;path.replace(backup)
                db.execute('DELETE FROM assets WHERE project_id=? AND id=?',(pid,asset_id))
        except Exception:
            if backup and backup.exists():backup.replace(path)
            raise
        self._cleanup_backup(backup)

    def admet_report(self,token,pid,rid,sid,filename):
        stage=self.stage(token,pid,rid,sid)
        if filename not in {c['path'] for c in self.admet_datasets(token,pid,rid,sid)}:
            raise AccessDenied('Tabela ADMET não autorizada.')
        path=self.store.scoped_path(pid,filename)
        from .artifact_naming import publication_metadata
        original=publication_metadata(path).get('legacy_name',path.name)
        stems=[Path(original).stem]
        for suffix in ('_evaluated','_selected','_rejected','_alerts','_rules_lipinski','_rules_veber','_rules_lipinski-veber','_rules_fragment','_BBB+','_BBB-','_HIA+','_HIA-'):
            if Path(original).stem.endswith(suffix):stems.append(Path(original).stem[:-len(suffix)])
        for stem in stems:
            for candidate in stage.get('artifacts',[]):
                report=Path(candidate)
                if report.parent==path.parent and publication_metadata(report).get('legacy_name',report.name)==stem+'_admet_report.json' and report.is_file():
                    return json.loads(self.store.read_file(token,pid,str(report),MAX_VIEW_BYTES+1))
        return None

    def admet_model(self, token, pid, rid, sid, filename, subset='all'):
        if subset not in FILTERS:
            raise ValueError('Filtro ADMET inválido.')
        choices = self.admet_datasets(token, pid, rid, sid)
        if filename not in {c['path'] for c in choices}:
            raise AccessDenied('Tabela ADMET não autorizada.')
        path = self.store.scoped_path(pid, filename)
        stage = self.stage(token, pid, rid, sid)
        model_path = str(path.with_name(path.stem+'_egg'+SUFFIX))
        if model_path in stage.get('artifacts', []) and Path(model_path).is_file():
            model = load_view(self.store.read_file(token, pid, model_path, MAX_VIEW_BYTES+1))
        else:
            import pandas as pd
            model = egg_view(pd.read_csv(io.BytesIO(self.store.read_file(token, pid, filename))), path.stem+' · BOILED-Egg')
        nodes = model['nodes']
        # User-supplied results may omit classification columns. Compute missing
        # classifications with the same evaluator used by the ADMET worker.
        if any(not n['properties'].get('BBB') or not n['properties'].get('HIA') for n in nodes):
            from caad.admet import MoleculeEvaluator
            evaluator = MoleculeEvaluator()
            for node in nodes:
                props = node['properties']
                if not props.get('HIA'):
                    props['HIA'] = evaluator.classify_hia(node['x'], node['y'])
                if not props.get('BBB'):
                    extra = evaluator.calculate_properties(props.get('canonical_smiles', ''))
                    if extra:
                        def descriptor(name):
                            value=props.get(name)
                            return value if type(value) in (int,float) and math.isfinite(value) else extra[name]
                        props['BBB'] = evaluator.classify_bbb(node['x'], node['y'],
                            descriptor('MW'),descriptor('HBD'),descriptor('RB'))
        if subset in ('selected','rejected'):
            expected=subset=='selected'
            nodes=[n for n in nodes if str(n['properties'].get('ADMET_selected','')).lower()==str(expected).lower()]
        elif subset != 'all':
            nodes = [n for n in nodes if n['properties'].get(subset[:3]) == subset]
        return dict(model,nodes=nodes,title=path.stem+' · '+FILTERS[subset],
                    dataset_origin=next(c for c in choices if c['path']==filename),
                    admet_report=self.admet_report(token,pid,rid,sid,filename))

    def admet_csv(self, token, pid, rid, sid, filename, subset):
        model = self.admet_model(token, pid, rid, sid, filename, subset)
        if subset == 'all':
            return self.store.read_file(token, pid, filename)
        # Existing subset CSVs contain only identifiers and SMILES. Downloads
        # from the explorer keep all properties shown by the interactive graph.
        fields = list(dict.fromkeys(['molecule_chembl_id', 'canonical_smiles', 'TPSA', 'WLOGP', 'BBB', 'HIA']+
            [k for n in model['nodes'] for k in n['properties']]))
        stream = io.StringIO(newline=''); writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(n['properties'] for n in model['nodes'])
        return stream.getvalue().encode('utf-8')

    def admet_png(self, token, pid, rid, sid, filename, subset):
        model = self.admet_model(token, pid, rid, sid, filename, subset)
        from .egg_plot import egg_png
        return egg_png(model)
