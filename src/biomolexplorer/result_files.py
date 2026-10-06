"""Authorized result-file management and dynamic ADMET datasets."""
import csv
import io
import json
import math
import time
from pathlib import Path
from uuid import uuid4

from .workspace import AccessDenied
from .visualizations import MAX_VIEW_BYTES, SUFFIX, egg_view, load_view

FILTERS = {'all': 'Todos os compostos avaliados', 'BBB+': 'Apenas BBB+',
           'BBB-': 'Apenas BBB−', 'HIA+': 'Apenas HIA+'}


class ResultFiles:
    def __init__(self, store):
        self.store = store

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
        return [{'path': str(p), 'name': p.name, 'size': p.stat().st_size}
                for name in stage.get('artifacts', [])
                if (p := self.store.scoped_path(pid, name)).is_file()
                and not (stage['operation']=='retrieve_structures' and p.name in ('pdb_codes.csv','retrieval_report.json'))
                and not (stage['operation'] == 'admet' and p.suffix.lower() == '.json' and p.name!='molecule_exclusions.json')]

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
        from Bio.PDB import parse_pdb_header
        resolution=parse_pdb_header(str(path)).get('resolution')
        available=set()
        for line in self.store.read_file(token,pid,str(path),MAX_VIEW_BYTES).decode().splitlines():
            if line.startswith('HETATM') and line[17:20].strip() not in ('HOH','WAT'):
                available.add((line[17:20].strip(),int(line[22:26]),line[21:22].strip()))
        return {'pdb_id':path.stem.upper(),'revision':sha256(data).hexdigest(),'resolution':resolution,
                'ligands':[row for row in reader if row['PDB_CODE'].upper()==path.stem.upper()],
                'available_ligands':[{'LIGAND':name,'RESNUM':str(number),'CHAIN':chain} for name,number,chain in sorted(available)]}

    def set_pdb_ligands(self,token,pid,rid,sid,filename,records,revision):
        """Replace one structure's curated ligands, preserving other PDB rows."""
        import re
        from .storage import write_text
        from .stage_cache import artifact_manifest
        path=self.store.scoped_path(pid,filename)
        metadata=path.parent/'pdb_codes.csv'
        backup=None
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
                    if key in seen:continue
                    seen.add(key)
                    row={k:'' for k in fields}
                    row.update(PDB_CODE=context['pdb_id'],LIGAND=ligand,RESNUM=str(number),CHAIN=chain)
                    if 'RESOLUTION' in fields:row['RESOLUTION']='' if context['resolution'] is None else str(context['resolution'])
                    curated.append(row)
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
                    if changed:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),run[0]))
        except Exception:
            if backup is not None:write_text(backup,metadata)
            raise
        return self.pdb_ligands(token,pid,rid,sid,filename)

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
        return result

    def graph_datasets(self,token,pid,rid,sid):
        stage=self.stage(token,pid,rid,sid)
        result=[]
        for filename in stage.get('artifacts',[]):
            path=self.store.scoped_path(pid,filename)
            if not path.name.endswith(SUFFIX) or not path.is_file():continue
            model=load_view(self.store.read_file(token,pid,str(path),MAX_VIEW_BYTES+1))
            if model['kind']=='graph':
                result.append({'path':str(path),'name':model['title']})
        from collections import Counter
        titles=Counter(choice['name'] for choice in result)
        for index,choice in enumerate(result,1):
            if titles[choice['name']]>1:choice['name']=f"Análise {index:02d} · "+choice['name']
        return result

    def graph_model(self,token,pid,rid,sid,filename):
        stage=self.stage(token,pid,rid,sid)
        path=self.store.scoped_path(pid,filename)
        if str(path) not in stage.get('artifacts',[]) or not path.name.endswith(SUFFIX):
            raise AccessDenied('Grafo não autorizado para esta etapa.')
        model=load_view(self.store.read_file(token,pid,str(path),MAX_VIEW_BYTES+1))
        if model['kind']!='graph':raise ValueError('Selecione um resultado de grafos.')
        # Old runs remain explorable without changing their files or cache hashes.
        if 'fragment' not in model:
            from caad.graph_results import common_fragment
            allowed=set(model['mcc'])
            model['fragment']=common_fragment([n['properties'].get('canonical_smiles') for n in model['nodes'] if n['id'] in allowed])
        return model

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
        if subset != 'all':
            nodes = [n for n in nodes if n['properties'].get(subset[:3]) == subset]
        return dict(model, nodes=nodes, title=path.stem+' · '+FILTERS[subset])

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
        # Object-oriented Agg avoids pyplot's global state in concurrent sessions.
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.patches import Ellipse
        figure = Figure(figsize=(11, 8)); FigureCanvasAgg(figure)
        ax = figure.subplots(); ax.set_facecolor('#EEF2F7')
        for center, width, height, color in [((75, 2), 150, 6, 'white'), ((42, 2.3), 94, 4.4, '#FFE066')]:
            ax.add_patch(Ellipse(center, width, height, facecolor=color, edgecolor='#64748B'))
        for classification, color in [('BBB+', '#DC2626'), ('BBB-', '#2552E8')]:
            nodes = [n for n in model['nodes'] if n['properties'].get('BBB') == classification]
            ax.scatter([n['x'] for n in nodes], [n['y'] for n in nodes], c=color, label=classification, s=35)
        ax.set(xlim=(min([0]+[n['x'] for n in model['nodes']]), max([200]+[n['x'] for n in model['nodes']])),
               ylim=(min([-2]+[n['y'] for n in model['nodes']]), max([7]+[n['y'] for n in model['nodes']])),
               xlabel='TPSA (Å²)', ylabel='WLOGP', title=model['title'])
        ax.legend(); figure.tight_layout()
        stream = io.BytesIO(); figure.savefig(stream, format='png', dpi=200)
        return stream.getvalue()
