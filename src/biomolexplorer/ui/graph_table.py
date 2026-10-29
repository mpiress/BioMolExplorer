"""Paged node/edge laboratory with compound inspection and complete CSV export."""
import csv
import io
import json
import flet as ft
from .localization import verbatim
from .table_sort import TableSort
from biomolexplorer.table_sort import is_numeric,ordered


class GraphTable:
    PAGE_SIZE=50

    def __init__(self,explorer):
        self.explorer=explorer;self.offset=0;self.sorting=TableSort()
        self.kind=ft.Dropdown(label='Tabela',value='nodes',width=160,on_select=self.search,
            options=[ft.DropdownOption(key='nodes',text='Nós'),ft.DropdownOption(key='edges',text='Arestas')])
        self.query=ft.TextField(label='Buscar na tabela',expand=True,on_submit=self.search)
        self.order=ft.Dropdown(label='Ordenar',value='id',width=200,on_select=self.preset,
            options=[ft.DropdownOption(key='id',text='Identificador / extremos'),ft.DropdownOption(key='value',text='Métrica / escore decrescente')])
        self.caption=ft.Text(size=12);self.body=ft.Column(spacing=4)

    def preset(self,e=None):
        self.sorting.key=None;self.search(e)

    def rows(self):
        explorer=self.explorer;view=explorer.view;graph=explorer.graph()
        query=(self.query.value or '').casefold()
        if self.kind.value=='nodes':
            values=explorer.analysis.metrics(view.mode,explorer.metric)[0]
            groups=explorer.analysis.community_cache.get(view.mode,({},None,''))[0]
            rows=[dict(compound_id=node,metric=explorer.metric,metric_value=values.get(node),
                community=groups.get(node),properties=view.nodes[node]['properties']) for node in graph]
            if query:rows=[row for row in rows if query in json.dumps(row,ensure_ascii=False).casefold()]
            rows=ordered(rows,lambda r:r.get('metric_value') if self.order.value=='value' else r['compound_id'],self.order.value!='value',self.order.value=='value')
        else:
            rows=[dict(source=a,target=b,**attrs) for a,b,attrs in graph.edges(data=True)]
            if query:rows=[row for row in rows if query in json.dumps(row,ensure_ascii=False).casefold()]
            rows=ordered(rows,lambda r:r.get('value') if self.order.value=='value' else r['source']+' '+r['target'],self.order.value!='value',self.order.value=='value')
        fields={key for row in rows for key in row}
        if self.sorting.key not in fields:self.sorting.key=None
        numeric={key for key in fields if key not in ('compound_id','source','target') and is_numeric(r.get(key) for r in rows)}
        return self.sorting.order(rows,numeric=numeric)

    def refresh(self):
        rows=self.rows();self.offset=min(self.offset,max(0,(len(rows)-1)//self.PAGE_SIZE*self.PAGE_SIZE))
        self.caption.value=f'{self.explorer.view.mode} · subgrafo filtrado · {len(rows)} registros · página {self.offset//self.PAGE_SIZE+1} / {max(1,(len(rows)+49)//50)}'
        fields=['compound_id','metric','metric_value','community','properties'] if self.kind.value=='nodes' else ['source','target']+sorted({key for row in rows for key in row}-{'source','target'})
        table=ft.DataTable(columns=[ft.DataColumn(verbatim(ft.Text(key))) for key in fields],rows=[],column_spacing=24)
        numeric={key for key in fields if key not in ('compound_id','source','target') and is_numeric(r.get(key) for r in rows)}
        self.sorting.bind(table,fields,self.search,numeric)
        async def choose(e,node):
            view=self.explorer.view
            await view.ui.guard(lambda:view.select(node,show_popup=True))
        for row in rows[self.offset:self.offset+self.PAGE_SIZE]:
            cells=[]
            for key in fields:
                value=row.get(key)
                if key in ('compound_id','source','target'):
                    async def open_node(e,node=value):await choose(e,node)
                    control=ft.TextButton(str(value),on_click=open_node)
                else:
                    text=json.dumps(value,ensure_ascii=False,sort_keys=True) if isinstance(value,dict) else str(value) if value is not None else '—'
                    control=ft.Text(text,max_lines=2,overflow=ft.TextOverflow.ELLIPSIS,tooltip=text,selectable=True)
                cells.append(ft.DataCell(verbatim(control)))
            table.rows.append(ft.DataRow(cells=cells))
        self.body.controls=[ft.Row([table],scroll=ft.ScrollMode.AUTO)] if rows else [ft.Text('Nenhum registro corresponde aos filtros.')]

    def search(self,e=None):
        self.offset=0;self.refresh();self.explorer.view.ui.page.update()

    def page(self,direction):
        self.offset=max(0,self.offset+direction*self.PAGE_SIZE);self.refresh();self.explorer.view.ui.page.update()

    def csv(self):
        rows=self.rows();stream=io.StringIO()
        if self.kind.value=='nodes':
            rows=[{**{k:v for k,v in row.items() if k!='properties'},**{'property.'+k:v for k,v in row['properties'].items()}} for row in rows]
        fields=sorted({key for row in rows for key in row})
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        writer.writerows({key:json.dumps(value,ensure_ascii=False) if isinstance(value,(dict,list)) else value for key,value in row.items()} for row in rows)
        return stream.getvalue().encode('utf-8')

    async def export(self,e=None):
        async def action():
            view=self.explorer.view
            await view.ui.call(view.ui.store.project,view.token,view.project_id)
            if not view.active or view.ui.token!=view.token:return
            data=await view.ui.call(self.csv)
            if view.active and view.ui.token==view.token:
                await view.ui.picker.save_file(file_name=f'network_{view.mode}_{self.kind.value}.csv',src_bytes=data)
        await self.explorer.view.ui.guard(action)

    def build(self):
        self.refresh()
        return ft.Column([ft.Row([self.kind,self.order]),ft.Row([self.query,ft.TextButton('Buscar',on_click=self.search)]),self.caption,
            ft.Row([ft.TextButton('Anterior',on_click=lambda e:self.page(-1)),ft.TextButton('Próxima',on_click=lambda e:self.page(1)),
                ft.TextButton('Exportar tabela completa CSV',on_click=self.export)],wrap=True),self.body])
