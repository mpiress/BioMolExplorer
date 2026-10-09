"""Authorized score tables with audited row removal and actual-pose previews."""
import csv
from pathlib import Path
from .compound_tables import CompoundTables
from .input_validation import columns
from .workspace import AccessDenied


class DockingResults(CompoundTables):
    def tables(self,token,pid,rid,sid):
        _,stage=self._stage(token,pid,rid,sid)
        result=[]
        for name in stage.get('artifacts',[]):
            path=self.store.scoped_path(pid,name)
            if path.suffix=='.csv' and path.is_file() and {'molecule_chembl_id','canonical_smiles','conformer_file'}<=columns(path):
                fields=columns(path)
                if 'score' in fields or {'vina','dock6'}<=fields:
                    result.append(dict(path=str(path),name=path.name,integrated=False))
        # Keep a summary for every batch, without duplicated receptor tables.
        from .docking_data import result_tables
        authoritative={str(p) for p in result_tables(t['path'] for t in result)}
        result=[t for t in result if t['path'] in authoritative]
        return result

    def pose(self,token,pid,rid,sid,table,index,version,column='conformer_file'):
        path=self._path(token,pid,rid,sid,table)
        with path.open('rb') as stream:
            if self._digest(stream)!=version:raise ValueError('Atualize a tabela antes de visualizar a conformação.')
        with path.open(encoding='utf-8-sig',newline='') as stream:rows=list(csv.DictReader(stream))
        if not 0<=index<len(rows) or column not in ('conformer_file','vina_pose','dock6_pose'):raise AccessDenied('Conformação não autorizada.')
        value=rows[index].get(column)
        if not value:raise ValueError('A conformação selecionada não está disponível.')
        pose=self.store.scoped_path(pid,Path(value) if Path(value).is_absolute() else path.parent/value)
        _,stage=self._stage(token,pid,rid,sid)
        if str(pose) not in stage.get('artifacts',[]) or not pose.is_file():raise AccessDenied('Conformação não autorizada para esta etapa.')
        return str(pose)
