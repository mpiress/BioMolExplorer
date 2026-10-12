"""On-demand molecular distributions, source joins and linked scatter selection."""
from pathlib import Path
import flet as ft
import flet.canvas as canvas
from biomolexplorer.molecular_network import numeric_properties,histogram,association
from biomolexplorer.visualizations import PointIndex
from .network_charts import card,distribution


class MolecularDashboard:
    def __init__(self,explorer):
        self.explorer=explorer;self.view=explorer.view;self.sources={};self.linked_ids=None;self.initialized=False
        self.source=ft.Dropdown(label='Fonte molecular · execução atual',value='graph',options=[ft.DropdownOption(key='graph',text='Dados registrados no grafo')],on_select=self.load_source,enable_filter=True,enable_search=True,menu_height=280,expand=True)
        self.source_query=ft.TextField(label='Pesquisar nome ou prefixo do arquivo',on_change=self.search_sources)
        self.source_page=0;self.source_count=ft.Text(size=12)
        self.scope=ft.Dropdown(label='Conjunto de compostos',value='mcc',expand=True,options=[ft.DropdownOption(key='mcc',text='MCC · maior componente conectado'),ft.DropdownOption(key='full',text='Full Graph')],on_select=self.analyze)
        self.property=ft.Dropdown(label='Descritor molecular',options=[],width=230)
        self.axis=ft.Dropdown(label='Eixo X',value='topology',width=230,options=[ft.DropdownOption(key='topology',text='Centralidade selecionada')])
        self.lower=ft.TextField(label='Descritor mínimo',width=140);self.upper=ft.TextField(label='Descritor máximo',width=140)
        self.notice=ft.Text(size=12,color='#64748B');self.bars=ft.Column(spacing=6)
        self.comparison=ft.Column(spacing=6);self.summary=ft.Text(size=12)
        self.drawing=canvas.Canvas(height=300,expand=True,on_resize=self.resize_scatter,resize_interval=100)
        self.scatter_pairs=[];self.scatter_width=660
        self.index=PointIndex({});self.points={};self.request=0
        self.neighbor=ft.Dropdown(label='Comparar composto conectado',disabled=True,hint_text='Selecione um composto na rede',options=[],width=280,enable_search=True,enable_filter=True,menu_height=300)
        self.structures=ft.Row(wrap=True)

    async def initialize(self):
        if self.initialized:return
        view=self.view;context=view.data.get('dashboard_context');initial_request=self.request
        await self.catalog()
        if not view.active or view.ui.token!=view.token:return
        # Native graph exports contain a molecular CSV and an edge CSV with
        # the same stem. Only the exact molecular artifact is a safe default.
        stem=Path(context['path']).name.removesuffix('.biomol-view.json') if context else ''
        correlated=[source for source in self.sources.values() if source['stage_id']==context['stage_id']
            and source['name']==stem+'.csv' and 'Molecules' in Path(source['path']).parts] if context else []
        if self.request==initial_request and self.source.value=='graph':
            self.source.value=correlated[0]['artifact_id'] if len(correlated)==1 else 'graph'
        automatic_source=self.source.value
        self.search_sources();await self.load_source()
        if view.active and view.ui.token==view.token and self.source.value==automatic_source and self.source.value!='graph' and not (view.molecular_enrichment or {}).get('matched'):
            self.source.value='graph';await self.load_source()
            self.notice.value='O arquivo correlato não forneceu correspondências válidas; usando os dados registrados no grafo.'
        if view.active and view.ui.token==view.token:self.initialized=True

    async def catalog(self,e=None):
        async def action():
            context=self.view.data.get('dashboard_context')
            if not context:self.notice.value='Procedência da execução indisponível; apenas atributos registrados no grafo.';self.view.ui.page.update();return
            from biomolexplorer.result_files import ResultFiles
            sources=await self.view.ui.call(ResultFiles(self.view.ui.store).graph_sources,self.view.token,self.view.project_id,context['run_id'])
            if not self.view.active or self.view.ui.token!=self.view.token:return
            self.sources={s['artifact_id']:s for s in sources}
            self.source_page=0;self.search_sources()
            self.notice.value='Escolha uma única fonte. Arquivos de outros experimentos não são combinados automaticamente.'
            self.view.ui.page.update()
        await self.view.ui.guard(action)

    def search_sources(self,e=None):
        if e is not None:self.source_page=0
        query=(self.source_query.value or '').casefold()
        found=[source for source in self.sources.values() if query in source['name'].casefold()]
        self.source_page=min(self.source_page,max(0,(len(found)-1)//50))
        shown=found[self.source_page*50:(self.source_page+1)*50]
        selected=self.sources.get(self.source.value)
        if selected and selected not in shown:shown=[selected,*shown]
        options=[ft.DropdownOption(key='graph',text='Dados registrados no grafo')]
        for source in shown:
            label=source['name']+' · '+source['stage_name']+' · '+source['stage_id'][:8]
            context=source['context']
            target=context.get('target') or context.get('search_term') or context.get('configured_target')
            if target:label+=' · '+str(target)
            options.append(ft.DropdownOption(key=source['artifact_id'],text=label,content=ft.Text(label,tooltip=label)))
        self.source.options=options
        self.source_count.value=f'{len(found)} correspondências · página {self.source_page+1} / {max(1,(len(found)+49)//50)} · até 50 arquivos por página'+(' · fonte selecionada preservada' if selected and selected not in found else '')
        self.source.tooltip=selected['name'] if selected else 'Nome completo e origem nas opções'
        if e is not None:self.view.ui.page.update()

    def page_sources(self,direction):
        self.source_page=max(0,self.source_page+direction);self.search_sources();self.view.ui.page.update()

    def resize_scatter(self,e):
        if not self.view.active or e.width<120:return
        self.scatter_width=e.width;self.draw_scatter();self.drawing.update()

    def draw_scatter(self):
        pairs=self.scatter_pairs;self.drawing.shapes=[];self.points={}
        if len(pairs)>2000:
            self.drawing.shapes=[canvas.Text(15,15,'Acima de 2.000 pares; restrinja ao LCC.')]
        elif pairs:
            xmin=min(p[1] for p in pairs);xmax=max(p[1] for p in pairs);ymin=min(p[2] for p in pairs);ymax=max(p[2] for p in pairs)
            right=max(100,self.scatter_width-25)
            self.points={n:(45+(right-45)*(a-xmin)/max(1e-12,xmax-xmin),245-200*(b-ymin)/max(1e-12,ymax-ymin)) for n,a,b in pairs}
            self.drawing.shapes=[canvas.Line(45,40,45,245),canvas.Line(45,245,right,245),
                canvas.Text(45,265,f'X: {xmin:.3g}–{xmax:.3g}',style=ft.TextStyle(size=11)),
                canvas.Text(45,10,f'Y: {ymin:.3g}–{ymax:.3g}',style=ft.TextStyle(size=11))]+[
                    canvas.Circle(px,py,4,paint=ft.Paint(color='#2563EB')) for px,py in self.points.values()]
        self.index=PointIndex(self.points)

    async def load_source(self,e=None):
        async def action():
            if not self.view.active or self.view.ui.token!=self.view.token:return
            context=self.view.data.get('dashboard_context');selected=self.source.value
            if selected=='graph':
                self.request+=1;self.view.molecular_enrichment=None
                self.notice.value='Fonte inicial: dados registrados no grafo. Escolha outro arquivo para complementar as propriedades.'
                await self.analyze()
                if self.view.active and self.view.selected in self.view.displayed():await self.view.select(self.view.selected,show_popup=False)
                return
            if not context or selected not in self.sources:self.view.ui.notify('Selecione uma fonte molecular autorizada.');return
            self.request+=1;request=self.request
            from biomolexplorer.result_files import ResultFiles
            result=await self.view.ui.call(ResultFiles(self.view.ui.store).graph_enrichment,self.view.token,self.view.project_id,
                context['run_id'],context['stage_id'],context['path'],selected,context.get('sha256'))
            if request!=self.request or not self.view.active or self.view.ui.token!=self.view.token:return
            self.view.molecular_enrichment=result
            self.view.selection_version+=1
            self.notice.value=f'Associação ao Full Graph: {result["matched"]} correspondências · {len(result["conflicts"])} IDs conflitantes excluídos · {result["unmatched"]} ausentes. Valores originais do grafo têm prioridade.'
            self.refresh_options();await self.analyze()
            if self.view.selected in self.view.displayed():await self.view.select(self.view.selected,show_popup=False)
        await self.view.ui.guard(action)

    def refresh_options(self,numeric=None):
        if numeric is None:_,numeric=numeric_properties(self.view.data['nodes'],self.view.molecular_enrichment)
        fields=sorted(numeric)
        self.property.disabled=self.axis.disabled=not bool(fields)
        self.property.options=[ft.DropdownOption(key=f,text=f) for f in fields]
        if self.property.value not in fields:
            preferred=('MW','molecule_properties.full_mwt','molecule_properties.alogp','WLOGP','TPSA')
            self.property.value=next((key for key in preferred if key in fields),fields[0] if fields else None)
        self.axis.options=[ft.DropdownOption(key='topology',text='Centralidade selecionada')]+[ft.DropdownOption(key=f,text=f) for f in fields]
        if self.axis.value not in ['topology',*fields]:self.axis.value='topology'
        if not fields:self.notice.value='Nenhum descritor numérico disponível; topologia e estruturas permanecem acessíveis.'

    async def analyze(self,e=None):
        async def action():
            view=self.view
            if not view.active or view.ui.token!=view.token:return
            await view.ui.call(view.ui.store.project,view.token,view.project_id)
            self.request+=1;request=self.request;mode=self.scope.value
            _,numeric=await view.ui.call(numeric_properties,view.data['nodes'],view.molecular_enrichment)
            if request!=self.request or not view.active or view.ui.token!=view.token or self.scope.value!=mode:return
            self.refresh_options(numeric);scope=self.explorer.analysis.universe(mode)
            y={n:v for n,v in numeric.get(self.property.value,{}).items() if n in scope}
            x=(await view.ui.call(self.explorer.analysis.metrics,mode,self.explorer.metric))[0] if self.axis.value=='topology' else numeric.get(self.axis.value,{})
            if request!=self.request or not view.active or view.ui.token!=view.token or self.scope.value!=mode:return
            x={n:value for n,value in x.items() if n in scope}
            pairs,r=association(x,y)
            self.summary.value=f'{"MCC / LCC" if mode=="mcc" else "Full Graph"} · {self.property.value or "sem descritor"}: {len(y)} / {len(scope)} valores · Pearson r: '+(f'{r:.4g}' if r is not None else 'indisponível (menos de 3 pares ou variância zero)')+' · associação exploratória, sem inferência causal.'
            bins=histogram(y.values())
            chart=distribution('Distribuição molecular',[(f'{left:.3g}–{right:.3g}',count) for left,right,count in bins],color='#0D9488',x_label=self.property.value or 'Descritor')
            chart.col=12;self.bars.controls=[chart]
            groups={};mapping=(await view.ui.call(self.explorer.analysis.communities,mode))[0]
            if request!=self.request or not view.active or view.ui.token!=view.token or self.scope.value!=mode:return
            for node,value in y.items():
                group=('LCC' if node in self.explorer.analysis.lcc else 'Outros componentes')
                groups.setdefault(group,[]).append(value)
                if node in mapping:groups.setdefault(f'C{mapping[node]}',[]).append(value)
            self.comparison.controls=[ft.Text(f'{group} · n={len(values)} · média {sum(values)/len(values):.4g} · intervalo {min(values):.4g}–{max(values):.4g}') for group,values in list(groups.items())[:25]]
            if len(groups)>25:self.comparison.controls.append(ft.Text('Primeiros 25 grupos; isole uma comunidade para explorar seu perfil.',size=12))
            self.scatter_pairs=pairs;self.draw_scatter();view.ui.page.update()
        await self.view.ui.guard(action)

    async def filter_range(self,e=None):
        async def action():
            import math
            _,numeric=numeric_properties(self.view.data['nodes'],self.view.molecular_enrichment)
            lower=float(self.lower.value) if self.lower.value else -math.inf
            upper=float(self.upper.value) if self.upper.value else math.inf
            if math.isnan(lower) or math.isnan(upper) or lower>upper:raise ValueError('Intervalo molecular inválido.')
            scope=self.explorer.analysis.universe(self.scope.value)
            self.linked_ids={n for n,value in numeric.get(self.property.value,{}).items() if n in scope and lower<=value<=upper}
            await self.explorer.apply_filters()
        await self.view.ui.guard(action)

    async def click(self,e):
        if e.local_position is None:return
        node=self.index.nearest(e.local_position.x,e.local_position.y)
        if node:
            async def action():
                await self.explorer.reset_filters()
                if node not in self.view.displayed():
                    self.view.component_id=next(i for i,g in self.view.components.items() if node in g)
                    self.view.redraw()
                if node in self.view.displayed():await self.view.follow_relation(node)
            await self.view.ui.guard(action)

    async def compare(self,e=None):
        async def action():
            view=self.view;selected=view.selected
            if selected not in view.graph:view.ui.notify('Selecione um composto para explorar sua vizinhança.');return
            await view.ui.call(view.ui.store.project,view.token,view.project_id)
            neighbors=sorted(view.graph.neighbors(selected))
            self.neighbor.options=[ft.DropdownOption(key=n,text=n) for n in neighbors]
            if self.neighbor.value not in neighbors:
                self.neighbor.value=neighbors[0] if neighbors else None
            properties,_=numeric_properties([view.nodes[n] for n in [selected,*([self.neighbor.value] if self.neighbor.value else [])]],view.molecular_enrichment)
            from biomolexplorer.visualizations import molecule_image
            controls=[]
            for node in [selected,*([self.neighbor.value] if self.neighbor.value else [])]:
                image=await view.ui.call(molecule_image,properties[node].get('canonical_smiles'))
                controls.append(ft.Column([ft.Text(node),ft.Image(src=image,width=220,height=165) if image else ft.Text('Estrutura indisponível: SMILES ausente ou inválido.')]))
            if view.active and view.ui.token==view.token:self.structures.controls=controls;view.ui.page.update()
        await self.view.ui.guard(action)

    def build(self):
        self.property.width=None;self.property.expand=True;self.axis.width=None;self.axis.expand=True
        self.neighbor.width=None;self.neighbor.expand=True
        parameters=card('Dados e parâmetros',ft.Column([
            ft.Row([self.scope]),ft.Text('A análise inicia no MCC, o maior grupo de compostos conectados. A fonte fornece as propriedades; o conjunto define quais compostos são analisados.',size=12),
            ft.TextButton('Atualizar fontes desta execução',on_click=self.catalog),self.source_query,ft.Row([self.source]),self.source_count,
            ft.Row([ft.TextButton('Anterior',on_click=lambda e:self.page_sources(-1)),ft.TextButton('Próxima',on_click=lambda e:self.page_sources(1))],wrap=True),
            ft.TextButton('Carregar fonte selecionada',on_click=self.load_source),self.notice,
            ft.Row([self.property]),ft.Row([self.axis]),ft.TextButton('Analisar descritores',on_click=self.analyze),
            ft.Row([self.lower,self.upper],wrap=True),ft.TextButton('Selecionar intervalo na rede',on_click=self.filter_range)],
            spacing=10,horizontal_alignment=ft.CrossAxisAlignment.STRETCH),{'xs':12,'lg':4})
        scatter=card('Estrutura × topologia',ft.Column([ft.Container(ft.GestureDetector(content=self.drawing,on_tap_up=self.click),height=300),
            ft.Text('Eixos: X = parâmetro escolhido · Y = descritor molecular. Clique em um ponto para localizar o composto.',size=12)]),{'xs':12,'xl':6})
        results=ft.Container(col={'xs':12,'lg':8},content=ft.Column([self.summary,
            ft.ResponsiveRow([ft.Container(self.bars,col={'xs':12,'xl':6}),scatter],spacing=12,run_spacing=12),
            ft.ExpansionTile(title=ft.Text('Comparação entre componentes e comunidades'),controls=[self.comparison]),
            card('Comparação estrutural · compostos conectados',ft.Column([ft.Row([self.neighbor]),ft.TextButton('Comparar estruturas 2D',on_click=self.compare),self.structures]),12),
            ft.Text('Dados ausentes não são imputados. Intervalos filtram somente a vista.',size=12)],horizontal_alignment=ft.CrossAxisAlignment.STRETCH))
        return ft.ResponsiveRow([parameters,results],spacing=12,run_spacing=12)
