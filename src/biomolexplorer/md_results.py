"""Authorized per-pair MD results; exclusions keep historical cache bytes intact."""
import json
from pathlib import Path
import re
from .result_files import ResultFiles
from .workspace import AccessDenied


class MDResults:
    def __init__(self,store):self.store=store;self.results=ResultFiles(store)

    def pairs(self,token,pid,rid,sid):
        stage=self.results.stage(token,pid,rid,sid)
        if stage['operation']!='dynamics_dock6':raise AccessDenied('Resultados de dinâmica indisponíveis.')
        current=next((s for s in self.store.project(token,pid)['pipeline'] if s['id']==sid),None)
        if current is None:return []
        params=current['parameters'];excluded=set(params.get('excluded_pairs') or [])
        selected=params.get('selected_pairs');rows=[]
        for name in stage.get('artifacts',[]):
            if Path(name).name!='dynamics_results.json':continue
            path=self.store.scoped_path(pid,name)
            for record in json.loads(self.store.read_file(token,pid,str(path))).get('pairs',[]):
                key=record.get('id','')
                if not re.fullmatch('[a-f0-9]{64}',key) or key in excluded or selected is not None and key not in selected:continue
                rows.append(dict(record,root=str(path.parent/'pairs'/key)))
        return rows

    def files(self,token,pid,rid,sid,key):
        row=next((r for r in self.pairs(token,pid,rid,sid) if r['id']==key),None)
        if row is None:raise AccessDenied('Par de dinâmica indisponível.')
        root=self.store.scoped_path(pid,row['root'])
        return [f for f in self.results.files(token,pid,rid,sid) if Path(f['path']).is_relative_to(root)]

    def dashboard(self,token,pid,rid,sid,key):
        files=self.files(token,pid,rid,sid,key)
        data={}
        for name in ('analysis.json','protocol.json','energies.json'):
            match=next((f for f in files if Path(f['path']).name==name),None)
            if match:data[name[:-5]]=json.loads(self.store.read_file(token,pid,match['path']))
        return data

    def comparison(self,token,pid,rid,sid):
        """Group by receptor content AND protocol; never compute a ranking."""
        from .calculation_cache import identity
        rows=[]
        for pair in self.pairs(token,pid,rid,sid):
            if pair['status']!='succeeded':continue
            data=self.dashboard(token,pid,rid,sid,pair['id'])
            protocol=data.get('protocol',{});analysis=data.get('analysis',{}).get('summary',{})
            receptor=protocol.get('inputs',{}).get('receptor',{}).get('sha256')
            if not receptor:continue
            group=identity(dict(receptor=receptor,settings=protocol.get('settings'),software=protocol.get('software')))
            rows.append(dict(receptor=pair['receptor'],ligand=pair['ligand'],group=group,
                **{k:analysis.get(k) for k in ('mean_protein_backbone_rmsd','mean_ligand_rmsd','mean_hydrogen_bonds','mean_protein_rg','mean_sasa')},
                amber_score=data.get('energies',{}).get('amber_score')))
        return rows

    def remove(self,token,pid,rid,sid,key):
        with self.store.change(token,pid,'editor','exclude_md_pair','Par de dinâmica removido: '+key[:12]) as db:
            if db.execute("SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')",(pid,)).fetchone():
                raise ValueError('Aguarde a execução terminar antes de remover o par de dinâmica.')
            if not any(r['id']==key for r in self.pairs(token,pid,rid,sid)):raise AccessDenied('Par de dinâmica indisponível.')
            row=db.execute('SELECT pipeline FROM projects WHERE id=?',(pid,)).fetchone()
            pipeline=json.loads(row[0]);stage=next(s for s in pipeline if s['id']==sid)
            excluded=set(stage['parameters'].get('excluded_pairs') or []);excluded.add(key)
            stage['parameters']['excluded_pairs']=sorted(excluded)
            db.execute('UPDATE projects SET pipeline=?,revision=revision+1 WHERE id=?',(json.dumps(pipeline),pid))
        return self.store.project(token,pid)
