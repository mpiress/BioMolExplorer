"""One header interaction and direction indicator for all result tables."""
import inspect
from biomolexplorer.table_sort import ordered


class TableSort:
    def __init__(self,key=None,ascending=True):
        self.key,self.ascending=key,ascending

    def order(self,rows,getter=None,numeric=()):
        if self.key is None:return list(rows)
        return ordered(rows,lambda row:getter(row,self.key) if getter else row.get(self.key),self.ascending,self.key in numeric)

    def bind(self,table,fields,changed,numeric=()):
        """None denotes an action column, which has no data to sort."""
        for column,key in zip(table.columns,fields):
            if key is None:continue
            column.numeric=key in numeric
            def handler(key):
                async def sort(e):
                    self.ascending=not self.ascending if self.key==key else True
                    self.key=key;self.indicator(table,fields)
                    result=changed()
                    if inspect.isawaitable(result):await result
                return sort
            column.on_sort=handler(key)
        self.indicator(table,fields)

    def indicator(self,table,fields):
        table.sort_column_index=fields.index(self.key) if self.key in fields else None
        table.sort_ascending=self.ascending
