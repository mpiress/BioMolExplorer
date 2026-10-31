"""Lazy, inline stage results: compound tables, ADMET explorer and file tables."""
from pathlib import Path
import flet as ft
from .processing import save_file
from biomolexplorer.compound_tables import CompoundTables
from biomolexplorer.result_files import ResultFiles,FILTERS
from .compound_table import CompoundTableViewer
from .file_table import FileTable
from .results_viewer import ResultsViewer
from .localization import verbatim


class StageResults:
    def __init__(self,ui,project_id,run_id,stage,writable):
        self.ui,self.project_id,self.run_id,self.stage=ui,project_id,run_id,stage
        self.token=ui.token;self.active=True;self.loaded=False;self.expanded=False;self.sequence=0
        self.loading_results=False
        self.writable=writable;self.child=None;self.checkpoint_children=[]
        self.service=ResultFiles(ui.store)
        self.root=ft.Column([ft.Text('Expanda para consultar os resultados.',size=13,color='#64748B')],spacing=20)
        self.loading=ft.Container(visible=False)
        self.graph=ft.Column(spacing=16)

    def valid(self):
        return self.active and self.ui.token==self.token and self.ui.current and self.ui.current['id']==self.project_id

    def close(self):
        self.active=False;self.sequence+=1
        for child in self.checkpoint_children:child.close()
        if hasattr(self,'graph_files'):self.graph_files.active=False
        if self.child:
            if hasattr(self.child,'close'):self.child.close()
            self.child.active=False
            if hasattr(self.child,'selection_version'):self.child.selection_version+=1
            if hasattr(self.child,'preview_version'):self.child.preview_version+=1

    async def expand(self,e):
        self.expanded=str(e.data).lower()=='true'
        if self.expanded and not self.loaded:await self.ui.guard(self.load)

    async def load(self):
        if not self.valid() or self.loaded or self.loading_results:return
        self.loading_results=True
        self.root.controls=[ft.Text('Carregando resultados…',size=13)]
        self.ui.page.update()
        try:
            parts=self.stage.get('_result_sources',[])
            if parts:
                for child in self.checkpoint_children:child.close()
                self.checkpoint_children=[]
                active=any(r['status'] in ('queued','running','awaiting_input') for r in await self.ui.call(self.ui.store.list_runs,self.token,self.project_id))
                if not self.valid():return
                controls=[ft.Text('Resultados preservados em checkpoints distintos. Cada arquivo mantém sua origem e versão de processamento.',size=13,color='#64748B')]
                for part in parts:
                    child=StageResults(self.ui,self.project_id,part['publication_run_id'],dict(part,_checkpoint_subset=True),self.writable)
                    self.checkpoint_children.append(child)
                    await child.load()
                    if not self.valid():return
                    controls.append(child.root)
                    if self.writable:
                        async def remove_part(e,rid=part['publication_run_id']):
                            if self.valid():await self.ui.guard(lambda:self.ui.remove_stage_results(self.project_id,rid,self.stage['id']))
                        controls.append(ft.TextButton('Remover resultados deste conjunto',icon=ft.Icons.DELETE_OUTLINE,icon_color='#B91C1C',disabled=active,on_click=remove_part))
                self.root.controls=controls;self.loaded=True;self.ui.page.update();return
            operation='checkpoint_subset' if self.stage.get('_checkpoint_subset') else self.stage['operation']
            if self.stage.get('outputs_removed') and not self.stage.get('artifacts'):
                self.root.controls=[ft.Text('Resultados excluídos. Execute o bloco para gerar novamente. O histórico e os logs foram preservados.')]
                self.loaded=True;self.ui.page.update();return
            if self.stage.get('results_stale'):
                self.root.controls=[ft.Text('As entradas destes resultados foram removidas ou alteradas. Execute o bloco para atualizar a seleção.')]
                self.loaded=True;self.ui.page.update();return
            if operation=='dynamics_dock6':
                from biomolexplorer.md_results import MDResults
                from .md_results import MDResultsTable
                rows=await self.ui.call(MDResults(self.ui.store).pairs,self.token,self.project_id,self.run_id,self.stage['id'])
                if not self.valid():return
                self.child=MDResultsTable(self,rows)
                self.root.controls=[self.child.build()];self.loaded=True;self.ui.page.update();return
            if operation in ('docking_vina','docking_dock6','consensus'):
                if self.stage.get('skip_reason'):
                    self.root.controls=[ft.Text(self.stage['skip_reason'])]
                    self.loaded=True;self.ui.page.update();return
                if self.stage['status']=='succeeded':
                    from biomolexplorer.docking_results import DockingResults
                    from .docking_results import DockingResultsTable
                    tables=await self.ui.call(DockingResults(self.ui.store).tables,self.token,self.project_id,self.run_id,self.stage['id'])
                    if not self.valid():return
                    if tables:
                        self.child=DockingResultsTable(self.ui,self.project_id,self.run_id,self.stage['id'],tables,inline=True)
                        self.root.controls=[self.child.build()];await self.child.load()
                        self.loaded=True;return
            if operation=='redocking' and self.stage['status']=='succeeded':
                from biomolexplorer.redocking_results import RedockingResults
                from .redocking_results import RedockingResultsTable
                simulations=await self.ui.call(RedockingResults(self.ui.store).simulations,
                    self.token,self.project_id,self.run_id,self.stage['id'])
                if not self.valid():return
                self.child=RedockingResultsTable(self,simulations)
                self.root.controls=[self.child.build()] if simulations else [ft.Text('Nenhum valor de RMSD disponível para esta etapa.')]
                self.loaded=True;self.ui.page.update();return
            if operation in ('retrieve_compounds','retrieve_pubchem','expand_similar_compounds') and self.stage['status']=='succeeded':
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
                    self.dataset=ft.Dropdown(label='Arquivo de resultados ADMET',value=datasets[0]['path'],expand=True,enable_filter=True,enable_search=True,
                        options=[verbatim(ft.DropdownOption(key=d['path'],text=d.get('label',d['name'])),'text') for d in datasets],on_select=self.change_admet)
                    self.subset=ft.Dropdown(label='Compostos exibidos',value='all',expand=True,
                        options=[ft.DropdownOption(key=k,text=v) for k,v in FILTERS.items()],on_select=self.change_admet)
                    self.root.controls=[ft.ResponsiveRow([ft.Container(ft.Row([self.dataset]),col={'xs':12,'md':6}),
                        ft.Container(ft.Row([self.subset]),col={'xs':12,'md':6})],spacing=24,run_spacing=20),self.loading,self.graph]
                    self.ui.page.update();await self.change_admet(None);self.loaded=True;return
            if operation=='graphs':
                datasets=await self.ui.call(self.service.graph_datasets,self.token,self.project_id,self.run_id,self.stage['id'])
                if not self.valid():return
                if datasets:
                    self.dataset=ft.Dropdown(label='Análise de grafos · entrada e origem',value=next((d['path'] for d in datasets if not d.get('unavailable')),None),expand=True,enable_filter=True,enable_search=True,
                        options=[verbatim(ft.DropdownOption(key=d['path'],text=d.get('label',d['name']),disabled=bool(d.get('unavailable'))),'text') for d in datasets],on_select=self.change_graph)
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
                    from biomolexplorer.graph_exports import subset_label
                    descriptions={'Molecules':'Compostos do grafo','maxcomp':'Arestas do MCC','similarity':'Similaridade calculada',
                        'centroids':'Fragmento molecular comum','plots':'Apresentação do MCC e graus'}
                    self.graph_file_records=[dict(f,description=(subset_label(f['path'],getattr(self.ui,'tr',None)) or descriptions.get(Path(f['path']).parent.name,'Resultado'))+' · '+f.get('origin',''))
                        for f in files if not f['name'].endswith('.biomol-view.json')]
                    self.graph_files=FileTable(self.ui,self.project_id,[],
                        self.writable and not active,remove=remove,preview=preview,run_id=self.run_id,stage_id=self.stage['id'])
                    self.root.controls=[ft.Row([self.dataset]),ft.Text('Cada opção corresponde a um arquivo de entrada. As relações e os compostos de análises diferentes são apresentados separadamente.',size=12,color='#64748B'),
                        self.loading,self.graph,ft.ExpansionTile(title=ft.Text('Arquivos desta análise'),controls=[ft.Container(self.graph_files.build(),padding=20)])]
                    if any(d.get('unavailable') for d in datasets):self.root.controls.insert(1,ft.Text('Arquivos acima de 32 MB permanecem disponíveis para download, mas não são carregados no dashboard.',size=12,color='#B45309'))
                    self.ui.page.update()
                    if self.dataset.value:await self.change_graph(None)
                    self.loaded=True;return
            files=await self.ui.call(self.service.files,self.token,self.project_id,self.run_id,self.stage['id'])
            if not self.valid():return
            active=any(r['status'] in ('queued','running','awaiting_input') for r in await self.ui.call(self.ui.store.list_runs,self.token,self.project_id))
            if not self.valid():return
            if self.stage.get('_checkpoint_subset'):
                allowed=set(self.stage.get('artifacts',[]));files=[item for item in files if item['path'] in allowed]
            async def remove(item):
                await self.ui.call(self.service.remove,self.token,self.project_id,self.run_id,self.stage['id'],item['path'])
                self.stage['artifacts']=[p for p in self.stage.get('artifacts',[]) if p!=item['path']]
            async def preview(item):await self.ui.preview_artifact(self.project_id,item['path'])
            extra=None
            if self.stage['operation']=='retrieve_structures':
                from .pdb_results import PDBActions
                running=any(r['status'] in ('queued','running') for r in await self.ui.call(self.ui.store.list_runs,self.token,self.project_id))
                if not self.valid():return
                self.pdb_actions=PDBActions(self,self.writable and not running)
                extra=self.pdb_actions.actions
            self.child=FileTable(self.ui,self.project_id,files,self.writable and not active,remove=remove,preview=preview,extra_actions=extra,run_id=self.run_id,stage_id=self.stage['id'])
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
            if self.child:
                if hasattr(self.child,'close'):self.child.close()
                else:self.child.active=False;self.child.selection_version+=1
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
                    source_name=Path((model.get('admet_report') or {}).get('source_input','')).stem
                    prefix=source_name+'__' if source_name else ''
                    name=prefix+Path(filename).stem+('' if subset=='all' else '_'+subset)+('_egg.png' if kind=='png' else '.csv')
                    if self.valid() and sequence==self.sequence:await save_file(self.ui,file_name=name,src_bytes=data)
                await self.ui.guard(save)
            async def png(e):await download('png')
            async def csv(e):await download('csv')
            async def alerts(e):
                from .admet_alerts import ADMETAlertTable
                await self.ui.guard(ADMETAlertTable(self,filename).open)
            buttons=[ft.Button('Baixar gráfico EGG',icon=ft.Icons.DOWNLOAD,on_click=png),
                     ft.Button('Baixar arquivo selecionado',icon=ft.Icons.DOWNLOAD,on_click=csv)]
            async def table(e):
                from .result_table import show_result_table
                if self.valid() and sequence==self.sequence:
                    await self.ui.guard(lambda:show_result_table(self.ui,self.project_id,self.run_id,self.stage['id'],filename))
            buttons.append(ft.Button('Ver tabela',icon=ft.Icons.TABLE_CHART_OUTLINED,on_click=table))
            if model.get('admet_report'):
                async def report(e):
                    import json
                    data=json.dumps(model['admet_report'],ensure_ascii=False,indent=2).encode()
                    await self.ui.guard(lambda:save_file(self.ui,file_name=Path(filename).stem+'_admet_report.json',src_bytes=data))
                buttons.append(ft.Button('Baixar critérios e métodos ADMET',icon=ft.Icons.DOWNLOAD,on_click=report))
            self.child=ResultsViewer(self.ui,self.project_id,model,popup=True,actions=buttons)
            self.graph.controls=[self.child.build()];self.ui.page.update()
        if e is None:await action()
        else:await self.ui.guard(action)

    async def change_graph(self,e):
        async def action():
            if not self.valid():return
            self.sequence+=1;sequence=self.sequence;filename=self.dataset.value
            if self.child:
                if hasattr(self.child,'close'):self.child.close()
                else:self.child.active=False;self.child.selection_version+=1
            self.loading.visible=True;self.graph.controls=[];self.ui.page.update()
            try:
                model=await self.ui.call(self.service.graph_model,self.token,self.project_id,self.run_id,self.stage['id'],filename)
            finally:
                if self.valid() and sequence==self.sequence:self.loading.visible=False;self.ui.page.update()
            if not self.valid() or sequence!=self.sequence:return
            identifier=Path(filename).name.removesuffix('.biomol-view.json')
            parent=Path(filename).parent.parent
            self.graph_files.files=[f for f in self.graph_file_records if (Path(f['name']).stem==identifier or Path(f['name']).stem.startswith(identifier+'__')) and Path(f['path']).is_relative_to(parent)]
            self.graph_files.offset=0;self.graph_files.redraw()
            async def download(kind):
                async def save():
                    if not self.valid() or sequence!=self.sequence:return
                    data=await self.ui.call(self.service.graph_png if kind=='png' else self.service.graph_mcc_csv,
                        self.token,self.project_id,self.run_id,self.stage['id'],filename)
                    name=Path(filename).name.removesuffix('.biomol-view.json')+('_mcc_apresentacao.png' if kind=='png' else '_mcc.csv')
                    if self.valid() and sequence==self.sequence:await save_file(self.ui,file_name=name,src_bytes=data)
                await self.ui.guard(save)
            async def png(e):await download('png')
            async def csv(e):await download('csv')
            buttons=[ft.Button('Baixar apresentação MCC',icon=ft.Icons.DOWNLOAD,on_click=png),
                     ft.Button('Baixar compostos do MCC',icon=ft.Icons.DOWNLOAD,on_click=csv)]
            child=await self.ui.call(ResultsViewer,self.ui,self.project_id,model,popup=True,actions=buttons)
            if not self.valid() or sequence!=self.sequence:child.close();return
            self.child=child
            self.graph.controls=[self.child.build()];self.ui.page.update()
            await child.initialize()
        if e is None:await action()
        else:await self.ui.guard(action)
