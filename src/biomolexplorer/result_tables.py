"""Read-only inspection of any published pipeline CSV, using the same pager."""
from pathlib import Path
from .compound_tables import CompoundTables


class ResultTables(CompoundTables):
    requires_compounds=False
    available_statuses=('succeeded','failed','interrupted','cancelled','skipped')

    def tables(self,token,pid,rid,sid):
        _,stage=self._stage(token,pid,rid,sid)
        return [dict(path=str(path),name=path.name,integrated=False) for value in stage.get('artifacts',[])
                if (path:=self.store.scoped_path(pid,value)).suffix.lower()=='.csv' and path.is_file()]

    def page(self,*args,**kwargs):
        result=super().page(*args,**kwargs);result['can_edit']=False
        return result

    def remove(self,*args,**kwargs):
        raise ValueError('Esta tabela é somente leitura.')
