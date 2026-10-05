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
                and not (stage['operation'] == 'admet' and p.suffix.lower() == '.json' and p.name!='molecule_exclusions.json')]

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
