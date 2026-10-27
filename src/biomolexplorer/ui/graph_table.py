"""Paged node/edge laboratory with compound inspection and complete CSV export."""
import csv
import io
import json
import flet as ft
from .processing import save_file
from .localization import verbatim


class GraphTable:
    PAGE_SIZE=50

    def __init__(self,explorer):
        self.explorer=explorer;self.offset=0
        self.kind=ft.Dropdown(label='Tabela',value='nodes',width=160,on_select=self.search,
            options=[ft.DropdownOption(key='nodes',text='Nós'),ft.DropdownOption(key='edges',text='Arestas')])
        self.query=ft.TextField(label='Buscar na tabela',expand=True,on_submit=self.search)
        self.order=ft.Dropdown(label='Ordenar',value='id',width=200,on_select=self.search,
            options=[ft.DropdownOption(key='id',text='Identificador / extremos'),ft.DropdownOption(key='value',text='Métrica / escore decrescente')])
        self.caption=ft.Text(size=12);self.body=ft.Column(spacing=4)

    def rows(self):
        explorer=self.explorer;view=explorer.view;graph=explorer.graph()
        query=(self.query.value or '').casefold()
        if self.kind.value=='nodes':
            values=explorer.analysis.metrics(view.mode,explorer.metric)[0]
            groups=explorer.analysis.community_cache.get(view.mode,({},None,''))[0]
            rows=[dict(compound_id=node,metric=explorer.metric,metric_value=values.get(node),
                community=groups.get(node),properties=view.nodes[node]['properties']) for node in graph]
            if query:rows=[row for row in rows if query in json.dumps(row,ensure_ascii=False).casefold()]
            rows.sort(key=(lambda row:(-(row['metric_value'] or 0),row['compound_id'])) if self.order.value=='value' else lambda row:row['compound_id'])
        else:
            rows=[dict(source=a,target=b,**attrs) for a,b,attrs in graph.edges(data=True)]
            if query:rows=[row for row in rows if query in json.dumps(row,ensure_ascii=False).casefold()]
            rows.sort(key=(lambda row:(-float(row.get('value') or 0),row['source'],row['target'])) if self.order.value=='value' else lambda row:(row['source'],row['target']))
        return rows

    def refresh(self):
        rows=self.rows();self.offset=min(self.offset,max(0,(len(rows)-1)//self.PAGE_SIZE*self.PAGE_SIZE))
        self.caption.value=f'{self.explorer.view.mode} · subgrafo filtrado · {len(rows)} registros · página {self.offset//self.PAGE_SIZE+1} / {max(1,(len(rows)+49)//50)}'
        self.body.controls=[]
        async def choose(e,node):
            view=self.explorer.view
            await view.ui.guard(lambda:view.select(node,show_popup=True))
        for row in rows[self.offset:self.offset+self.PAGE_SIZE]:
            if self.kind.value=='nodes':
                node=row['compound_id'];props=row['properties']
                async def open_node(e,node=node):await choose(e,node)
                self.body.controls.append(ft.Container(padding=8,border=ft.Border.all(1,'#E2E8F0'),content=ft.Column([
                    ft.Row([verbatim(ft.TextButton(node,on_click=open_node)),ft.Text(f'{row["metric"]}: {row["metric_value"] if row["metric_value"] is not None else "—"} · C{row["community"] or "—"}',size=12)],wrap=True),
                    verbatim(ft.Text(json.dumps(props,ensure_ascii=False,sort_keys=True),size=11,selectable=True))])))
            else:
                async def open_source(e,node=row['source']):await choose(e,node)
                async def open_target(e,node=row['target']):await choose(e,node)
                self.body.controls.append(ft.Row([verbatim(ft.TextButton(row['source'],on_click=open_source)),ft.Text('↔'),
                    verbatim(ft.TextButton(row['target'],on_click=open_target)),verbatim(ft.Text(' · '.join(f'{k}: {v}' for k,v in row.items() if k not in ('source','target')),size=12))],wrap=True))
        if not rows:self.body.controls=[ft.Text('Nenhum registro corresponde aos filtros.')]

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
                await save_file(view.ui,file_name=f'network_{view.mode}_{self.kind.value}.csv',src_bytes=data)
        await self.explorer.view.ui.guard(action)

    def build(self):
        self.refresh()
        return ft.Column([ft.Row([self.kind,self.order]),ft.Row([self.query,ft.TextButton('Buscar',on_click=self.search)]),self.caption,
            ft.Row([ft.TextButton('Anterior',on_click=lambda e:self.page(-1)),ft.TextButton('Próxima',on_click=lambda e:self.page(1)),
                ft.TextButton('Exportar tabela completa CSV',on_click=self.export)],wrap=True),self.body])
