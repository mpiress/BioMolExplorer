"""Lazy, inline stage results: compound tables, ADMET explorer and file tables."""
from pathlib import Path
import flet as ft
from biomolexplorer.compound_tables import CompoundTables
from biomolexplorer.result_files import ResultFiles,FILTERS
from .compound_table import CompoundTableViewer
from .file_table import FileTable
from .results_viewer import ResultsViewer


class StageResults:
    def __init__(self,ui,project_id,run_id,stage,writable):
        self.ui,self.project_id,self.run_id,self.stage=ui,project_id,run_id,stage
        self.token=ui.token;self.active=True;self.loaded=False;self.expanded=False;self.sequence=0
        self.loading_results=False
        self.writable=writable;self.child=None
        self.service=ResultFiles(ui.store)
        self.root=ft.Column([ft.Text('Expanda para consultar os resultados.',size=13,color='#64748B')],spacing=20)
        self.loading=ft.ProgressBar(visible=False)
        self.graph=ft.Column(spacing=16)

    def valid(self):
        return self.active and self.ui.token==self.token and self.ui.current and self.ui.current['id']==self.project_id

    def close(self):
        self.active=False;self.sequence+=1
        if hasattr(self,'graph_files'):self.graph_files.active=False
        if self.child:
            self.child.active=False
            if hasattr(self.child,'selection_version'):self.child.selection_version+=1
            if hasattr(self.child,'preview_version'):self.child.preview_version+=1

    async def expand(self,e):
        self.expanded=str(e.data).lower()=='true'
        if self.expanded and not self.loaded:await self.ui.guard(self.load)

    async def load(self):
        if not self.valid() or self.loaded or self.loading_results:return
        self.loading_results=True
        self.root.controls=[ft.ProgressBar(),ft.Text('Carregando resultados…',size=13)]
        self.ui.page.update()
        try:
            operation=self.stage['operation']
            if operation in ('retrieve_compounds','expand_similar_compounds') and self.stage['status']=='succeeded':
                tables=await self.ui.call(CompoundTables(self.ui.store).tables,self.token,self.project_id,self.run_id,self.stage['id'])
                if not self.valid():return
                if tables:
                    self.child=CompoundTableViewer(self.ui,self.project_id,self.run_id,self.stage['id'],tables,inline=True)
                    self.root.controls=[self.child.build()];self.ui.page.update()
                    await self.child.load();self.loaded=True;return
            if operation=='admet':
                datasets=await self.ui.call(self.service.admet_datasets,self.token,self.project_id,self.run_id,self.stage['id'])
                if not self.valid():return
                if datasets:
                    self.dataset=ft.Dropdown(label='Arquivo de resultados ADMET',value=datasets[0]['path'],expand=True,
                        options=[ft.DropdownOption(key=d['path'],text=d['name']) for d in datasets],on_select=self.change_admet)
                    self.subset=ft.Dropdown(label='Compostos exibidos',value='all',expand=True,
                        options=[ft.DropdownOption(key=k,text=v) for k,v in FILTERS.items()],on_select=self.change_admet)
                    self.root.controls=[ft.ResponsiveRow([ft.Container(ft.Row([self.dataset]),col={'xs':12,'md':6}),
                        ft.Container(ft.Row([self.subset]),col={'xs':12,'md':6})],spacing=24,run_spacing=20),self.loading,self.graph]
                    self.ui.page.update();await self.change_admet(None);self.loaded=True;return
            if operation=='graphs':
                datasets=await self.ui.call(self.service.graph_datasets,self.token,self.project_id,self.run_id,self.stage['id'])
                if not self.valid():return
                if datasets:
                    self.dataset=ft.Dropdown(label='Análise de grafos · entrada e origem',value=datasets[0]['path'],expand=True,
                        options=[ft.DropdownOption(key=d['path'],text=d['name']) for d in datasets],on_select=self.change_graph)
                    files=await self.ui.call(self.service.files,self.token,self.project_id,self.run_id,self.stage['id'])
                    active=any(r['status'] in ('queued','running','awaiting_input') for r in await self.ui.call(self.ui.store.list_runs,self.token,self.project_id))
                    if not self.valid():return
                    async def remove(item):
                        await self.ui.call(self.service.remove,self.token,self.project_id,self.run_id,self.stage['id'],item['path'])
                        self.graph_file_records=[f for f in self.graph_file_records if f['path']!=item['path']]
                        if item['path']==self.dataset.value:
                            if self.child:self.child.active=False;self.child.selection_version+=1
                            self.graph.controls=[ft.Text('Visualização removida. Selecione outra análise ou execute a etapa novamente.')]
                    async def preview(item):await self.ui.preview_artifact(self.project_id,item['path'])
                    descriptions={'Molecules':'Compostos do MCC','maxcomp':'Arestas do MCC','similarity':'Similaridade calculada',
                        'centroids':'Fragmento molecular comum','plots':'Apresentação do MCC e graus'}
                    self.graph_file_records=[dict(f,description=descriptions.get(Path(f['path']).parent.name,'Resultado'))
                        for f in files if not f['name'].endswith('.biomol-view.json')]
                    self.graph_files=FileTable(self.ui,self.project_id,[],
                        self.writable and not active,remove=remove,preview=preview)
                    self.root.controls=[ft.Row([self.dataset]),ft.Text('Cada opção corresponde a um arquivo de entrada. As relações e os compostos de análises diferentes são apresentados separadamente.',size=12,color='#64748B'),
                        self.loading,self.graph,ft.ExpansionTile(title=ft.Text('Arquivos desta análise'),controls=[ft.Container(self.graph_files.build(),padding=20)])]
                    self.ui.page.update();await self.change_graph(None);self.loaded=True;return
            files=await self.ui.call(self.service.files,self.token,self.project_id,self.run_id,self.stage['id'])
            if not self.valid():return
            active=any(r['status'] in ('queued','running','awaiting_input') for r in await self.ui.call(self.ui.store.list_runs,self.token,self.project_id))
            if not self.valid():return
            async def remove(item):
                await self.ui.call(self.service.remove,self.token,self.project_id,self.run_id,self.stage['id'],item['path'])
                self.stage['artifacts']=[p for p in self.stage.get('artifacts',[]) if p!=item['path']]
            async def preview(item):await self.ui.preview_artifact(self.project_id,item['path'])
            self.child=FileTable(self.ui,self.project_id,files,self.writable and not active,remove=remove,preview=preview)
            self.root.controls=[self.child.build()] if files else [ft.Text('Nenhum arquivo de resultado disponível.',size=13,color='#64748B')]
            self.loaded=True;self.ui.page.update()
        except Exception:
            if self.valid():
                self.root.controls=[ft.Text('Não foi possível abrir estes resultados. Atualize a etapa para tentar novamente.',color='#B91C1C')]
                self.ui.page.update()
            raise
        finally:
            self.loading_results=False

    async def change_admet(self,e):
        async def action():
            if not self.valid():return
            self.sequence+=1;sequence=self.sequence
            if self.child:self.child.active=False;self.child.selection_version+=1
            self.loading.visible=True;self.graph.controls=[];self.ui.page.update()
            filename,subset=self.dataset.value,self.subset.value
            try:
                model=await self.ui.call(self.service.admet_model,self.token,self.project_id,self.run_id,self.stage['id'],filename,subset)
            finally:
                if self.valid() and sequence==self.sequence:self.loading.visible=False;self.ui.page.update()
            if not self.valid() or sequence!=self.sequence:return
            async def download(kind):
                async def save():
                    if not self.valid() or sequence!=self.sequence:return
                    data=await self.ui.call(self.service.admet_png if kind=='png' else self.service.admet_csv,
                        self.token,self.project_id,self.run_id,self.stage['id'],filename,subset)
                    name=Path(filename).stem+('' if subset=='all' else '_'+subset)+('_egg.png' if kind=='png' else '.csv')
                    if self.valid() and sequence==self.sequence:await self.ui.picker.save_file(file_name=name,src_bytes=data)
                await self.ui.guard(save)
            async def png(e):await download('png')
            async def csv(e):await download('csv')
            buttons=[ft.Button('Baixar gráfico EGG',icon=ft.Icons.DOWNLOAD,on_click=png),
                     ft.Button('Baixar arquivo selecionado',icon=ft.Icons.DOWNLOAD,on_click=csv)]
            self.child=ResultsViewer(self.ui,self.project_id,model,popup=True,actions=buttons)
            self.graph.controls=[self.child.build()];self.ui.page.update()
        if e is None:await action()
        else:await self.ui.guard(action)

    async def change_graph(self,e):
        async def action():
            if not self.valid():return
            self.sequence+=1;sequence=self.sequence;filename=self.dataset.value
            if self.child:self.child.active=False;self.child.selection_version+=1
            self.loading.visible=True;self.graph.controls=[];self.ui.page.update()
            try:
                model=await self.ui.call(self.service.graph_model,self.token,self.project_id,self.run_id,self.stage['id'],filename)
            finally:
                if self.valid() and sequence==self.sequence:self.loading.visible=False;self.ui.page.update()
            if not self.valid() or sequence!=self.sequence:return
            identifier=Path(filename).name.removesuffix('.biomol-view.json')
            self.graph_files.files=[f for f in self.graph_file_records if Path(f['name']).stem==identifier]
            self.graph_files.offset=0;self.graph_files.redraw()
            async def download(kind):
                async def save():
                    if not self.valid() or sequence!=self.sequence:return
                    data=await self.ui.call(self.service.graph_png if kind=='png' else self.service.graph_mcc_csv,
                        self.token,self.project_id,self.run_id,self.stage['id'],filename)
                    name=Path(filename).name.removesuffix('.biomol-view.json')+('_mcc_apresentacao.png' if kind=='png' else '_mcc.csv')
                    if self.valid() and sequence==self.sequence:await self.ui.picker.save_file(file_name=name,src_bytes=data)
                await self.ui.guard(save)
            async def png(e):await download('png')
            async def csv(e):await download('csv')
            buttons=[ft.Button('Baixar apresentação MCC',icon=ft.Icons.DOWNLOAD,on_click=png),
                     ft.Button('Baixar compostos do MCC',icon=ft.Icons.DOWNLOAD,on_click=csv)]
            self.child=ResultsViewer(self.ui,self.project_id,model,popup=True,actions=buttons)
            self.graph.controls=[self.child.build()];self.ui.page.update()
        if e is None:await action()
        else:await self.ui.guard(action)
