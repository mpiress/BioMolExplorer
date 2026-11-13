"""Resolve a selected consensus pose for MD inside the authorized project job."""
from pathlib import Path

from .workspace import AccessDenied


class _JobResults:
    """Internal adapter: job authorization precedes MDInputs.resolve.

    This is not an HTTP/token service. Every historical run and artifact is
    still required to belong to the same project and to remain published.
    """
    def __init__(self,store):
        self.store=store;self._stage_snapshot=None

    def _allowed(self,pid,stage):
        return {str(self.store.scoped_path(pid,path)) for path in stage.get('artifacts',[])}

    def _record(self,path,index,version,*,root):
        import csv
        from itertools import islice
        from .result_state import input_digest
        if type(index) is not int or index<0:raise ValueError('Seleção de linha inválida.')
        if input_digest(path)!=version:raise ValueError('Atualize a tabela antes de visualizar a conformação.')
        with Path(path).open(encoding='utf-8-sig',newline='') as stream:
            row=next(islice(csv.DictReader(stream),index,index+1),None)
        if row is None:raise ValueError('A linha não está mais nesta tabela.')
        return row

    def _path(self,token,pid,rid,sid,filename):
        path=self.store.scoped_path(pid,filename)
        _,stage=self._stage(token,pid,rid,sid)
        if str(path) not in stage.get('artifacts',[]) or path.suffix.lower()!='.csv' or not path.is_file():
            raise AccessDenied('Tabela não autorizada para esta etapa.')
        return path

    def _stage(self,token,pid,rid,sid):
        with self.store.connect() as db:
            row=db.execute('''SELECT r.project_id,r.updated,r.status,coalesce(v.revision,0) AS revision
                FROM runs r LEFT JOIN result_read_versions v ON v.run_id=r.id WHERE r.id=?''',(rid,)).fetchone()
        if row is None or row['project_id']!=pid:raise AccessDenied('Resultado de outro projeto.')
        key=(pid,rid,sid,row['updated'],row['status'],row['revision'],None)
        if self._stage_snapshot and self._stage_snapshot[0]==key:return self._stage_snapshot[1:]
        from .result_catalog import read_stage
        stage=read_stage(self.store,rid,sid)
        if not stage or stage.get('status')!='succeeded' or stage.get('outputs_removed'):
            raise ValueError('Resultado deste alvo indisponível na origem.')
        run={'id':rid,'project_id':pid,'stages':[stage]}
        self._stage_snapshot=(key,run,stage)
        return run,stage


class ConsensusPoses:
    def __init__(self,store,pid):
        self.store=store;self.pid=pid;self.results=_JobResults(store)

    def ligand(self,receptor,ligand):
        row=ligand.get('compound') or {}
        if not row.get('target_scores'):return ligand
        from .consensus_targets import target_details
        targets=target_details(row)
        stem=Path(receptor['path']).name
        matches=[target for target in targets if stem==target+'.pdb' or stem.startswith(target+'.')]
        if not matches:
            # Legacy retrieval PDBs may omit the chain in the filename. Only
            # accept that spelling when the PDB identifies a unique target.
            matches=[target for target in targets if stem.casefold()==target.split('_',1)[0].casefold()+'.pdb']
        if len(matches)!=1:raise ValueError('A pose do consenso não possui um alvo correspondente ao receptor selecionado.')
        target=matches[0]
        if target==row.get('receptor_id'):return ligand
        table=self.store.scoped_path(self.pid,ligand['path'])
        relative=table.relative_to(self.store.project_dir(self.pid)).parts
        if len(relative)<4 or relative[0]!='runs':
            raise ValueError('Resultado deste alvo indisponível na origem.')
        rid,sid=relative[1:3]
        from .result_state import input_digest
        version=input_digest(table)
        index=ligand.get('source_row')
        if type(index) is not int or index<0:raise ValueError('Seleção de linha inválida.')
        actual=self.results._record(table,index,version,root=self.store.project_dir(self.pid))
        if any(actual.get(key)!=row.get(key) for key in ('molecule_chembl_id','receptor_id','target_scores')):
            raise ValueError('Os resultados foram alterados. Revise a confirmação antes de salvar.')
        from .consensus_targets import target_route
        route=target_route(self.results,None,self.pid,rid,sid,str(table),index,version,target,'dock6_pose')
        if route is None:raise ValueError('Os resultados foram alterados. Revise a confirmação antes de salvar.')
        child,source_rid,source_sid,selection=route
        source=Path(selection['table'])
        selected=child._record(source,selection['index'],selection['version'],root=self.store.project_dir(self.pid))
        pose=Path(selected['conformer_file'])
        pose=self.store.scoped_path(self.pid,pose if pose.is_absolute() else source.parent/pose)
        _,producer=child._stage(None,self.pid,source_rid,source_sid)
        if str(pose) not in producer.get('artifacts',[]) or not pose.is_file():
            raise ValueError('Resultado deste alvo indisponível na origem.')
        pose=str(pose)
        # A prepared pose retains its original coordinates and all provenance.
        # pair_id stays stable; molecular_identity hashes the resolved pose.
        compound=dict(actual,receptor_id=target,conformer_file=pose,dock6_pose=pose,engine='dock6')
        for field in ('pose_file','prepared_mol2','prepared_pdbqt'):compound.pop(field,None)
        compound['preparation_metadata']=''
        return dict(ligand,compound=compound,name=str(compound['molecule_chembl_id'])+' · '+target+' · DOCK6')
