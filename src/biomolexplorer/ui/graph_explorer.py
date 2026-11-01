"""Analytical controls shared by inline and standalone molecular graph views."""
import json
import math
from collections import Counter
import flet as ft
from .processing import save_file
from biomolexplorer.artifact_naming import download_name
import flet.canvas as canvas
from biomolexplorer.network_analysis import METRICS,PALETTE
from .localization import verbatim
from .network_charts import card,distribution


class GraphExplorer:
    def __init__(self,viewer):
        self.view=viewer;self.analysis=viewer.analysis
        self.metric='degree';self.community=None;self.ego=None;self.hops=1
        self.aggregate_mode=False;self.positions=None;self.request=0;self.community_enabled=False
        self.render_notice='';self.filtered_graph=self.analysis.graph
        self.geometry_notice=''
        self.query=ft.TextField(label='Filtrar identificadores',on_change=self.apply_filters)
        self.bounds=[ft.TextField(label=label,width=125,on_submit=self.apply_filters) for label in
            ('Grau mínimo','Grau máximo','Centralidade mínima','Centralidade máxima','Escore mínimo','Escore máximo')]
        self.metric_choice=ft.Dropdown(label='Centralidade / tamanho',value='degree',expand=True,
            options=[ft.DropdownOption(key=key,text=METRICS[key].split(' · ')[0]) for key in METRICS],on_select=self.change)
        self.community_choice=ft.Dropdown(label='Comunidade topológica',value='all',expand=True,
            options=[ft.DropdownOption(key='all',text='Todas as comunidades')],on_select=self.apply_filters,
            enable_filter=True,enable_search=True,menu_height=300)
        self.hop_choice=ft.Dropdown(label='Raio topológico',value='1',width=140,
            options=[ft.DropdownOption(key=str(i),text=f'{i}-hop') for i in (1,2,3)],on_select=self.apply_filters)
        self.focus_only=ft.Switch(label='Somente composto selecionado',on_change=self.apply_filters)
        self.color_communities=ft.Switch(label='Colorir comunidades',on_change=self.change)
        self.aggregate_choice=ft.Switch(label='Rede agregada por comunidade',on_change=self.change)
        self.notice=ft.Text(size=12,color='#B45309');self.method=ft.Text(size=12,color='#64748B')
        self.stats=ft.Column(spacing=10);self.ranking=ft.Column(spacing=4)
        self.distributions=ft.Column(spacing=8);self.groups=ft.Column(spacing=6)
        self.component_page=0
        self.community_page=0
        self.view_key=None;self.visible_graph=None
        from .molecular_network_dashboard import MolecularDashboard
        self.molecular=MolecularDashboard(self)
        from .molecular_inspector import MolecularInspector
        self.inspector=MolecularInspector(self)
        self.rank_query=ft.TextField(label='Pesquisar composto no ranking',on_change=self.refresh_ranking)
        self.rank_order=ft.Dropdown(label='Ordenação',value='desc',width=160,options=[ft.DropdownOption(key='desc',text='Decrescente'),ft.DropdownOption(key='asc',text='Crescente')],on_select=self.refresh_ranking)
        self.rank_limit=ft.Dropdown(label='Resultados',value='20',width=120,options=[ft.DropdownOption(key=str(n),text=str(n)) for n in (10,20,50,100)],on_select=self.refresh_ranking)
        self.rank_caption=ft.Text(size=12)
        self.rank_scope=ft.Dropdown(label='Escopo do ranking',value='visible',expand=True,on_select=self.calculate_ranking,
            options=[ft.DropdownOption(key='visible',text='Grafo exibido · aplicar filtros'),ft.DropdownOption(key='universe',text='Universo Full / MCC · antes dos filtros')])
        self.rank_notice=ft.Text(size=12,color='#64748B')
        self.rank_context=None;self.rank_context_key=None;self.rank_request=0
        self.focus_mode=False;self.sidebar_hidden=False;self.tab_index=0
        self.primary=ft.ResponsiveRow(spacing=10,run_spacing=10)
        from .graph_appearance import GraphAppearance
        self.appearance=GraphAppearance(self)

    def graph(self):
        key=(self.view.mode,self.view.component_id,self.community_choice.value,self.query.value,
             tuple(c.value for c in self.bounds),self.metric,self.ego,self.hop_choice.value,
             None if self.molecular.linked_ids is None else frozenset(self.molecular.linked_ids))
        key=(*key,bool(self.focus_only.value))
        if key==self.view_key and self.visible_graph is not None:return self.visible_graph
        def number(control):
            text=(control.value or '').strip()
            if not text:return None
            value=float(text)
            if not math.isfinite(value):raise ValueError('Informe limites numéricos finitos.')
            return value
        bounds=[number(c) for c in self.bounds]
        for start in (0,2,4):
            if bounds[start] is not None and bounds[start+1] is not None and bounds[start]>bounds[start+1]:
                raise ValueError('O limite mínimo deve ser menor ou igual ao máximo.')
        self.community=None if self.community_choice.value=='all' else int(self.community_choice.value)
        self.hops=int(self.hop_choice.value)
        self.filtered_graph=self.analysis.filtered(self.view.mode,self.view.component_id,self.community,
            self.query.value or '',tuple(bounds[:2]),self.metric,tuple(bounds[2:4]),tuple(bounds[4:]),self.ego,self.hops)
        if self.molecular.linked_ids is not None:self.filtered_graph=self.filtered_graph.subgraph(self.molecular.linked_ids)
        if self.focus_only.value and self.ego is not None:self.filtered_graph=self.filtered_graph.subgraph([self.ego])
        n=len(self.filtered_graph);m=self.filtered_graph.number_of_edges()
        self.render_notice=''
        if n>self.analysis.DRAW_NODES or m>self.analysis.DRAW_EDGES:
            self.render_notice=f'Rede com {n} nós / {m} arestas: desenho otimizado; todos os compostos continuam disponíveis para clique, busca e investigação. Nenhum nó foi amostrado.'
        self.visible_graph=self.filtered_graph
        self.view_key=key
        return self.visible_graph

    def close(self):
        self.rank_request+=1;self.rank_context=None;self.rank_context_key=None
        self.inspector.close()
        page=self.view.ui.page
        if getattr(page,'on_keyboard_event',None) is getattr(self,'keyboard_handler',None):page.on_keyboard_event=getattr(self,'previous_keyboard',None)
        if getattr(page,'on_resize',None) is getattr(self,'resize_handler',None):page.on_resize=getattr(self,'previous_resize',None)
        self.analysis.cache.clear();self.analysis.layouts.clear();self.analysis.community_cache.clear();self.analysis.statistics_cache.clear()
        self.positions=None;self.visible_graph=None;self.filtered_graph=None
        self.ranking.controls=[];self.groups.controls=[];self.distributions.controls=[]
        self.molecular.request+=1;self.molecular.drawing.shapes=[];self.molecular.sources={}

    async def initialize(self):
        view=self.view
        if not view.active or view.ui.token!=view.token:return
        await view.ui.call(view.ui.store.project,view.token,view.project_id)
        await self.prepare()
        if not view.active or view.ui.token!=view.token:return
        view.redraw()
        await self.molecular.initialize()
        if view.active:view.ui.page.update()

    async def prepare(self):
        view=self.view;mode=view.mode
        await view.ui.call(self.analysis.metrics,mode,self.metric)
        await view.ui.call(self.analysis.communities,mode)
        await view.ui.call(self.analysis.statistics,mode,True)

    async def change(self,e=None):
        async def action():
            view=self.view
            if not view.active or view.ui.token!=view.token:return
            await view.ui.call(view.ui.store.project,view.token,view.project_id)
            self.request+=1;request=self.request;mode=view.mode
            previous_layout=view.layout_mode;previous_aggregate=self.aggregate_mode;previous_metric=self.metric
            metric=self.metric_choice.value;layout=view.layout_choice.value
            if metric=='degree' and view.color_choice.value=='density':metric='density'
            community_enabled=bool(self.color_communities.value or self.aggregate_choice.value or layout=='community')
            self.notice.value='Calculando análise topológica…';view.ui.page.update()
            await view.ui.call(self.analysis.metrics,mode,metric)
            await view.ui.call(self.analysis.communities,mode)
            positions=None;layout_notice=''
            if layout not in ('network','circular'):
                positions,layout_notice=await view.ui.call(self.analysis.layout,mode,layout)
            statistics=await view.ui.call(self.analysis.statistics,mode,True)
            if request!=self.request or not view.active or view.ui.token!=view.token or view.mode!=mode:return
            self.metric=metric;self.community_enabled=community_enabled
            if previous_metric!=metric and self.molecular.axis.value=='topology':
                self.molecular.scatter_pairs=[];self.molecular.draw_scatter()
                self.molecular.summary.value='Métrica alterada; clique em Analisar descritores para atualizar a associação.'
            self.metric_choice.value=metric
            self.aggregate_mode=bool(self.aggregate_choice.value)
            if layout not in ('network','circular') and positions is None:
                view.layout_choice.value=view.layout_mode
            else:self.positions=positions;view.layout_mode=layout
            self.notice.value=layout_notice
            self.layout_notice=layout_notice
            self.refresh_choices()
            view.color_mode=view.color_choice.value;view.labels=bool(view.label_choice.value)
            view.color_choice.disabled=metric not in ('degree','density')
            view.redraw();self.update_panels(statistics)
            if view.selected in view.points and self.inspector.properties:
                await self.inspector.update(view.selected,self.inspector.properties,view.selection_version)
            if previous_layout!=view.layout_mode or previous_aggregate!=self.aggregate_mode:await view.fit_view()
            if previous_metric!=metric and self.molecular.initialized and self.molecular.axis.value=='topology':await self.molecular.analyze()
            view.ui.page.update()
        await self.view.ui.guard(action)

    async def apply_filters(self,e=None):
        async def action():
            view=self.view
            if not view.active or view.ui.token!=view.token:return
            await view.ui.call(view.ui.store.project,view.token,view.project_id)
            self.graph() # Validate first; retain the last drawing on invalid input.
            before=set(view.points);view.selection_version+=1;view.redraw()
            if set(view.points)!=before:await view.fit_view()
            if view.selected is not None and view.selected not in view.points:view.clear_selection()
            await self.inspector.sync_filters()
            view.ui.page.update()
        await self.view.ui.guard(action)

    async def ego_selected(self,e=None):
        if self.view.selected is None:self.view.ui.notify('Selecione um composto para explorar sua vizinhança.');return
        self.aggregate_mode=False;self.aggregate_choice.value=False
        self.ego=self.view.selected;self.focus_only.value=False;await self.apply_filters()

    def clear_filters(self):
        """Clear every restricting state, including component and aggregation."""
        self.request+=1;self.view_key=None;self.visible_graph=None
        self.aggregate_mode=False;self.aggregate_choice.value=False
        self.query.value='';self.community_choice.value='all';self.ego=None
        self.community=None;self.view.component_id=None
        if hasattr(self.view,'component_choice'):self.view.component_choice.value='all'
        self.molecular.linked_ids=None
        self.focus_only.value=False
        for control in self.bounds:control.value=''

    async def reset_filters(self,e=None):
        self.clear_filters()
        await self.apply_filters()
        await self.view.fit_view()

    async def restore_full(self,e=None):
        # Explicit recovery also works when "full" is already selected.
        from types import SimpleNamespace
        await self.view.change_mode(SimpleNamespace(control=SimpleNamespace(value='full')))

    def update_panels(self,statistics=None):
        view=self.view;graph=self.analysis.universe(view.mode)
        stats=statistics or self.analysis.statistics(view.mode)
        values,notice=self.analysis.cache.get((view.mode,self.metric),({},''))
        mapping,modularity,community_notice=self.analysis.community_cache.get(view.mode,({},None,''))
        scope='Full Graph' if view.mode=='full' else 'LCC / MCC'
        self.method.value=scope+' · '+METRICS[self.metric]+' Todas as métricas ignoram escores como distância. Centralidade não indica atividade biológica.'
        origin=view.data.get('origin',{})
        relation=('Similaridade de fingerprints · '+str(origin.get('metric',''))+' · limiar '+str(origin.get('threshold','não registrado'))
            if origin.get('kind')=='fingerprints' else 'Relações recebidas: value é o escore registrado; método e significado dependem da origem.')
        self.method.value+=' '+relation
        self.notice.value=' '.join(filter(None,[view.opening_notice,self.render_notice,self.geometry_notice,notice,community_notice,getattr(self,'layout_notice','')]))
        indicators=[('Nodes',stats['nodes'],'Número de compostos no universo'),('Edges',stats['edges'],'Relações no universo'),
            ('Connected Components',stats['components'],'Componentes conexos'),('Average Degree',f'{stats["mean_degree"]:.4g}','Número médio de conexões'),
            ('Density',f'{stats["density"]:.4g}','Arestas / relações possíveis; sem unidade'),
            ('Average Clustering',f'{stats["clustering"]:.4g}' if stats['clustering'] is not None else '—','Não calculado ou custo indisponível quando —'),
            ('Isolates',stats['isolates'],'Compostos sem conexões'),('Communities',len(set(mapping.values())) if mapping else 0 if not len(graph) else '—','Louvain · grupos topológicos do universo')]
        indicators.extend((label,'—' if stats.get(key) is None else f'{stats[key]:.5g}',tip) for label,key,tip in [
            ('Diameter · LCC','diameter_lcc','Maior caminho mínimo do LCC; até 1.500 nós / 20.000 arestas'),
            ('Average Path Length · LCC','mean_path_lcc','Média dos caminhos mínimos no LCC; sem pesos; até 1.500 nós / 20.000 arestas'),
            ('Assortativity','assortativity','Correlação entre os graus nas extremidades das arestas')])
        self.primary.controls=[ft.Container(col=4,padding=12,bgcolor='#FFFFFF',border_radius=10,border=ft.Border.all(1,'#E2E8F0'),content=ft.Column([ft.Text(label,size=11),ft.Text(str(value),size=23,weight=ft.FontWeight.W_600)],spacing=4)) for label,value in [('Nodes',stats['nodes']),('Edges',stats['edges']),('Communities',len(set(mapping.values())) if mapping else 0 if not len(graph) else '—')]]
        self.stats.controls=[ft.ResponsiveRow([ft.Container(col={'xs':6,'md':3},padding=12,bgcolor='#FFFFFF',border_radius=10,
            border=ft.Border.all(1,'#E2E8F0'),tooltip=tip,content=ft.Column([ft.Text(label,size=11,color='#64748B'),ft.Text(str(value),size=23,weight=ft.FontWeight.W_600)],spacing=4))
            for label,value,tip in indicators],spacing=10,run_spacing=10),
            ft.Text(f'{scope} · LCC {len(self.analysis.lcc)} / {len(self.analysis.graph)} · Modularidade: '+(f'{modularity:.4g}' if modularity is not None else 'não calculada')+
                f' · Exibição individual: {len(view.points)} compostos. Indicadores e rankings usam o universo antes dos filtros.',size=12)]
        self.refresh_ranking()
        from biomolexplorer.molecular_network import histogram
        bins=histogram(values.values())
        counts=Counter(mapping.values())
        plots=[distribution('Connectivity · Degree Distribution',stats['degrees'].items(),x_label='Grau'),
            distribution('Connectivity · Component Size',stats['component_sizes'].items(),x_label='Compostos / componente',y_label='Componentes'),
            distribution('Centrality · '+self.metric,[(f'{left:.3g}–{right:.3g}',count) for left,right,count in bins],color='#7C3AED',x_label=self.metric)]
        if counts:plots.append(distribution('Community Structure · Size',sorted(Counter(counts.values()).items()),color='#0D9488',x_label='Compostos / comunidade',y_label='Comunidades'))
        self.distributions.controls=[ft.ResponsiveRow(plots,spacing=12,run_spacing=12)]
        self.groups.controls=[] # Community details live in the paged selector and analytical chart.
        self.refresh_choices()

    def ranking_context(self):
        """Keep ranking caches separate from the graph's appearance universe."""
        from biomolexplorer.network_analysis import NetworkAnalysis
        view=self.view
        graph=self.graph() if self.rank_scope.value=='visible' else self.analysis.universe(view.mode)
        key=(self.rank_scope.value,view.mode,self.view_key if self.rank_scope.value=='visible' else None)
        if key!=self.rank_context_key or self.rank_context is None:
            universe=self.analysis.universe(view.mode)
            if len(graph)==len(universe) and graph.number_of_edges()==universe.number_of_edges():
                self.rank_context=(self.analysis,view.mode)
            else:
                data=dict(nodes=[{'id':node} for node in graph],edges=[dict(source=a,target=b) for a,b in graph.edges()])
                self.rank_context=(NetworkAnalysis(data),'full')
            self.rank_context_key=key
        return self.rank_context,key

    def refresh_ranking(self,e=None):
        view=self.view;((analysis,mode),_)=self.ranking_context()
        metric=self.rank_metric.value if hasattr(self,'rank_metric') else self.metric
        if metric in ('degree','density') and (mode,metric) not in analysis.cache:analysis.metrics(mode,metric)
        values,notice=analysis.cache.get((mode,metric),({},'Clique em Calcular ranking para analisar este escopo.'))
        query=(self.rank_query.value or '').casefold();descending=self.rank_order.value=='desc'
        ordered=sorted(values.items(),key=lambda pair:(-pair[1] if descending else pair[1],pair[0]))
        positions={node:i+1 for i,(node,_) in enumerate(ordered)}
        shown=[(node,value) for node,value in ordered if query in node.casefold()][:int(self.rank_limit.value)]
        maximum=max((abs(v) for v in values.values()),default=1) or 1
        mapping=self.analysis.community_cache.get(view.mode,({},None,''))[0]
        graph=analysis.universe(mode)
        scope='Grafo exibido' if self.rank_scope.value=='visible' else 'Universo antes dos filtros'
        self.rank_caption.value=f'{scope} · {view.mode} · {metric} · {len(graph)} nós / {graph.number_of_edges()} arestas · {len(shown)} / {len(values)} compostos'
        self.rank_notice.value=notice
        self.ranking.controls=[]
        for node,value in shown:
            async def choose(e,node=node):
                if node in view.displayed():
                    await view.ui.guard(lambda:view.select(node,show_popup=False))
                else:
                    view.search.value=node;await view.find(None)
                if hasattr(self,'tabs'):
                    self.tabs.selected_index=0;self.switch_tab()
            component=next((i for i,g in view.components.items() if node in g),None)
            label=f'C{mapping[node]}' if node in mapping else f'Componente {component}'
            self.ranking.controls.append(ft.Container(padding=10,border_radius=8,bgcolor='#EFF6FF' if positions[node]<=3 else '#FFFFFF',
                content=ft.ResponsiveRow([ft.Container(col={'xs':12,'md':5},content=ft.Row([ft.Text(str(positions[node]),width=30),verbatim(ft.TextButton(node,on_click=choose,tooltip=node,expand=True))])),
                    ft.Column([ft.Text(f'{value:.6g} · {label}',size=12),ft.ProgressBar(value=min(1,abs(value)/maximum),color='#0D9488',bgcolor='#E2E8F0')],col={'xs':12,'md':7},spacing=6)],run_spacing=6)))
        if not shown:self.ranking.controls=[ft.Text('Nenhum composto corresponde à pesquisa ou a métrica está indisponível.')]
        if e is not None:view.ui.page.update()

    def switch_tab(self,e=None):
        self.tab_index=self.tabs.selected_index
        if self.tab_index==4:self.table.refresh()
        for i,body in enumerate(self.sections):body.expanded=i==self.tab_index
        self.view.ui.page.update()

    async def toggle_focus(self,e=None):
        self.focus_mode=not self.focus_mode
        self.inspector_card.visible=not self.focus_mode and not self.sidebar_hidden
        self.graph_card.col=12 if self.focus_mode or self.sidebar_hidden else {'xs':12,'lg':9}
        self.tabs_bar.visible=not self.focus_mode;self.overview.visible=not self.focus_mode
        self.configuration.visible=not self.focus_mode;self.footer.visible=not self.focus_mode;self.primary.visible=not self.focus_mode
        self.graph_box.height=max(360,min(1000,(getattr(self.view.ui.page,'height',None) or 1000)-220)) if self.focus_mode else max(320,min(660,(getattr(self.view.ui.page,'height',None) or 1000)*.56))
        self.expand_button.tooltip='Sair do modo expandido' if self.focus_mode else 'Expandir grafo'
        self.expand_button.icon=ft.Icons.CLOSE_FULLSCREEN if self.focus_mode else ft.Icons.OPEN_IN_FULL
        self.view.ui.page.update()

    async def toggle_sidebar(self,e=None):
        self.sidebar_hidden=not self.sidebar_hidden
        self.inspector_card.visible=not self.sidebar_hidden and not self.focus_mode
        self.graph_card.col=12 if self.sidebar_hidden or self.focus_mode else {'xs':12,'lg':9}
        self.view.ui.page.update()

    def refresh_choices(self):
        if not hasattr(self.view,'component_choice'):return
        components=list(self.view.components.items())
        self.component_page=min(self.component_page,max(0,(len(components)-1)//50))
        subset=components[self.component_page*50:(self.component_page+1)*50]
        self.view.component_choice.options=[ft.DropdownOption(key='all',text=f'Todos · {len(components)} componentes')]+[
            ft.DropdownOption(key=str(i),text=f'Componente {i} · {len(group)} compostos · {self.analysis.component_edges[i]} arestas') for i,group in subset]
        counts=Counter(self.analysis.community_cache.get(self.view.mode,({},None,''))[0].values())
        groups=sorted(counts)
        self.community_page=min(self.community_page,max(0,(len(groups)-1)//50))
        self.community_choice.options=[ft.DropdownOption(key='all',text=f'Todas · {len(groups)} comunidades')]+[
            ft.DropdownOption(key=str(c),text=f'Comunidade {c} · {counts[c]} nós') for c in groups[self.community_page*50:(self.community_page+1)*50]]

    def page_choices(self,kind,direction):
        attr=kind+'_page';setattr(self,attr,max(0,getattr(self,attr)+direction))
        self.refresh_choices();self.view.ui.page.update()

    def draw_aggregate(self):
        graph,notice=self.analysis.aggregate(self.view.mode)
        self.render_notice='Representação agregada de todo o universo, antes dos filtros: cada nó é uma comunidade; arestas contam relações entre comunidades.'
        from biomolexplorer.visualizations import separated_graph_layout
        positions=separated_graph_layout(graph,circular=True)
        view=self.view;shapes=[];view.raster.visible=False
        view.raster.src=''
        self.appearance.map(graph,{node:attrs['size'] for node,attrs in graph.nodes(data=True)},'degree',aggregate=True)
        view.scale_label.value='Agregado: tamanho = número de compostos; cor = comunidade; espessura = quantidade de relações entre comunidades.'
        xmin=min((p[0] for p in positions.values()),default=0);ymin=min((p[1] for p in positions.values()),default=0)
        dx=max((p[0] for p in positions.values()),default=0)-xmin
        dy=max((p[1] for p in positions.values()),default=0)-ymin
        scale=min(30,3986/max(1,dx),4006/max(1,dy))
        view.WIDTH=max(720,dx*scale+110);view.HEIGHT=max(460,dy*scale+90)
        coordinates={n:(55+(p[0]-xmin)*scale,45+(p[1]-ymin)*scale) for n,p in positions.items()}
        for a,b,attrs in graph.edges(data=True):
            shapes.append(canvas.Line(*coordinates[a],*coordinates[b],paint=ft.Paint(color='#94A3B8',stroke_width=min(6,1+math.log1p(attrs['count'])))))
        for node,attrs in graph.nodes(data=True):
            x,y=coordinates[node];radius=min(24,6+math.sqrt(attrs['size']))
            shapes.append(canvas.Circle(x,y,radius,paint=ft.Paint(color=PALETTE[(node-1)%len(PALETTE)])))
            if len(graph)<=200 or view.labels:shapes.append(canvas.Text(x+radius+4,y-8,f'C{node} · {attrs["size"]} nós',style=ft.TextStyle(size=12)))
        view.points={};from biomolexplorer.visualizations import PointIndex
        view.index=PointIndex(coordinates);view.base_shapes=shapes;view.drawing.shapes=shapes
        view.summary.value=self.render_notice
        view.scene.width=view.drawing.width=view.raster_layer.width=view.WIDTH
        view.scene.height=view.drawing.height=view.raster_layer.height=view.HEIGHT
        view.raster.visible=len(graph)>self.analysis.DRAW_NODES or graph.number_of_edges()>self.analysis.DRAW_EDGES
        if view.raster.visible:
            from .graph_raster import paint_graph
            view.raster.width,view.raster.height=view.WIDTH,view.HEIGHT
            view.raster.src=paint_graph(shapes,view.WIDTH,view.HEIGHT)
            view.drawing.shapes=[]
        self.update_panels()

    async def open_community(self,identifier):
        """An aggregate node opens the same compound view as every other source."""
        self.aggregate_mode=False;self.aggregate_choice.value=False
        self.query.value='';self.ego=None;self.focus_only.value=False
        self.view.component_id=None;self.view.component_choice.value='all'
        self.molecular.linked_ids=None
        for control in self.bounds:control.value=''
        groups=sorted(set(self.analysis.community_cache.get(self.view.mode,({},None,''))[0].values()))
        if identifier not in groups:return
        self.community_page=groups.index(identifier)//50
        self.refresh_choices();self.community_choice.value=str(identifier)
        await self.apply_filters()
        await self.view.fit_view()

    async def export(self,kind):
        async def action():
            view=self.view
            await view.ui.call(view.ui.store.project,view.token,view.project_id)
            if not view.active or view.ui.token!=view.token:return
            if kind=='csv':data=await view.ui.call(self.analysis.export_csv,view.mode,self.metric)
            elif kind=='png':data=await view.ui.call(self.image_png)
            elif kind in ('gexf','graphml'):
                if self.aggregate_mode:
                    view.ui.notify('Abra a rede individual para exportar os compostos e suas propriedades.');return
                data=await view.ui.call(self.graph_exchange,kind)
            elif kind=='json':
                stats=await view.ui.call(self.analysis.statistics,view.mode,True)
                aggregate,_=await view.ui.call(self.analysis.aggregate,view.mode) if self.community_enabled else (None,'')
                data=json.dumps(dict(universe=view.mode,metric=self.metric,method='undirected unweighted topology',statistics=stats,
                    dataset_origin=view.data.get('dataset_origin'),
                    appearance=dict(layout=view.layout_mode,color_attribute=self.appearance.attribute.value,palette=self.appearance.palette.value,
                        size=self.appearance.size.value,community_colors=bool(self.color_communities.value),aggregate=self.aggregate_mode,
                        show_edges=bool(self.appearance.edges.value),labels=view.labels,legend_title=self.appearance.title,
                        color_maximum=self.appearance.maximum),
                    communities=self.analysis.community_cache.get(view.mode,({},None,'')),
                    aggregate=None if aggregate is None else dict(nodes=list(aggregate.nodes(data=True)),edges=list(aggregate.edges(data=True))),
                    displayed_nodes=list(view.points),filters=[c.value for c in self.bounds],query=self.query.value,ego=self.ego,hops=self.hops,focus_only=bool(self.focus_only.value),
                    molecular_source=None if not view.molecular_enrichment else {k:view.molecular_enrichment[k] for k in ('source','matched','unmatched','conflicts','rule')},
                    scaffold=None if self.inspector.result is None else {k:v for k,v in self.inspector.result.items() if k!='image'},
                    molecular_view={'property':self.molecular.property.value,'axis':self.molecular.axis.value,
                        'lower':self.molecular.lower.value,'upper':self.molecular.upper.value,
                        'selected_ids':None if self.molecular.linked_ids is None else sorted(self.molecular.linked_ids)}),ensure_ascii=False).encode()
            else:data=await view.ui.call(self.svg)
            if view.active and view.ui.token==view.token:await save_file(view.ui,file_name=download_name(view.data.get('dataset_origin',{}),result_type='Graph-'+str(view.mode),extension=kind,qualifier=view.mode,operation='graphs'),src_bytes=data)
        await self.view.ui.guard(action)

    def svg(self):
        from xml.sax.saxutils import escape
        view=self.view;parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{view.WIDTH}" height="{view.HEIGHT}" viewBox="0 0 {view.WIDTH} {view.HEIGHT}">',
            '<rect width="100%" height="100%" fill="white"/>',f'<title>{escape(view.data["title"])} · {view.mode} · {self.metric}</title>']
        shapes,height=self.image_shapes()
        parts[0]=f'<svg xmlns="http://www.w3.org/2000/svg" width="{view.WIDTH}" height="{height}" viewBox="0 0 {view.WIDTH} {height}">'
        for shape in shapes:
            if isinstance(shape,canvas.Line):parts.append(f'<line x1="{shape.x1}" y1="{shape.y1}" x2="{shape.x2}" y2="{shape.y2}" stroke="{shape.paint.color}" stroke-width="{shape.paint.stroke_width}"/>')
            elif isinstance(shape,canvas.Circle):
                outline=shape.paint.style==ft.PaintingStyle.STROKE
                parts.append(f'<circle cx="{shape.x}" cy="{shape.y}" r="{shape.radius}" fill="{"none" if outline else shape.paint.color}" stroke="{shape.paint.color if outline else "none"}" stroke-width="{shape.paint.stroke_width}"/>')
            elif isinstance(shape,canvas.Text):parts.append(f'<text x="{shape.x}" y="{shape.y}" font-size="10">{escape(shape.value or "")}</text>')
            elif isinstance(shape,canvas.Rect):parts.append(f'<rect x="{shape.x}" y="{shape.y}" width="{shape.width}" height="{shape.height}" fill="none" stroke="{shape.paint.color}"/>')
        return ('\n'.join(parts)+'\n</svg>').encode()

    def image_shapes(self):
        # The raster canvas contains only the live selection overlay.
        # Always export the complete base and the exact same legend mapping.
        legend,height=self.appearance.export_shapes(self.view.HEIGHT+16,self.view.WIDTH)
        return [*self.view.base_shapes,*(self.view.drawing.shapes if self.view.raster.visible else []),*legend],height

    def image_png(self):
        from .graph_raster import paint_graph
        shapes,height=self.image_shapes()
        return paint_graph(shapes,self.view.WIDTH,height)

    def graph_exchange(self,kind):
        """Export the visible molecular subgraph with attributes and appearance."""
        import io
        import networkx as nx
        view=self.view;graph=nx.Graph(name=view.data['title'])
        graph.graph.update(universe=view.mode,layout=view.layout_mode,metric=self.metric,
            topology='undirected, unweighted; similarity is an attribute, not a distance')
        if view.data.get('dataset_origin'):
            graph.graph.update(dataset_origin=view.data['dataset_origin']['origin'],
                               source_trace=view.data['dataset_origin'].get('source_trace',''))
        mapping=self.analysis.community_cache.get(view.mode,({},None,''))[0]
        for node in view.points:
            attributes={}
            for key,value in view.nodes[node]['properties'].items():
                if value is None:continue
                attributes[str(key)]=value if type(value) in (str,bool,int,float) else json.dumps(value,ensure_ascii=False)
            attributes.update({key:values[node] for (mode,key),(values,_) in self.analysis.cache.items() if mode==view.mode and node in values})
            if node in mapping:attributes['modularity_class']=mapping[node]
            color=self.appearance.colors.get(node,'#94A3B8');x,y=view.points[node]
            attributes.update(x=float(x),y=float(y),color=color,size=float(view.radii[node]))
            if kind=='gexf':attributes['viz']=dict(position=dict(x=float(x),y=-float(y),z=0.),
                color=dict(r=int(color[1:3],16),g=int(color[3:5],16),b=int(color[5:7],16),a=1.),size=float(view.radii[node]))
            graph.add_node(node,**attributes)
        for a,b,attrs in self.graph().edges(data=True):
            if a in graph and b in graph:
                attributes={str(k):v if type(v) in (str,bool,int,float) else json.dumps(v,ensure_ascii=False) for k,v in attrs.items() if v is not None}
                # Do not let a provider's "weight" change topological semantics.
                if 'weight' in attributes:attributes['relationship_weight']=attributes.pop('weight')
                graph.add_edge(a,b,**attributes)
        stream=io.BytesIO()
        if kind=='gexf':nx.write_gexf(graph,stream)
        else:nx.write_graphml(graph,stream)
        return stream.getvalue()

    def build(self,tools,style):
        view=self.view
        self.update_panels()
        view.color_choice.label='Normalização do grau'
        view.color_choice.tooltip='Outras centralidades usam a métrica selecionada para cor e tamanho.'
        for choice in (view.component_choice,view.layout_choice):choice.enable_filter=True;choice.enable_search=True;choice.menu_height=300
        view.layout_choice.options += [ft.DropdownOption(key=key,text=label) for key,label in
            [('fa2','ForceAtlas2 · componentes'),('fa2_linlog','ForceAtlas2 · LinLog · componentes'),
             ('force','Fruchterman–Reingold (até 1.500 nós)'),('fr_components','Fruchterman–Reingold · componentes'),
             ('random','Random · semente 42'),('radial','Concentric · grau por componente'),
             ('kamada','Kamada–Kawai (até 400 nós)'),('hierarchical','Hierarquia BFS por grau'),('community','Comunidades topológicas')]]
        view.layout_choice.on_select=self.change
        view.color_choice.on_select=self.change;view.label_choice.on_change=self.change
        # The toolbar wraps on small screens; children need bounded widths.
        view.mode_choice.width=None;view.mode_choice.expand=True
        async def detect(e):self.color_communities.value=True;await self.change(e)
        async def csv(e):await self.export('csv')
        async def report(e):await self.export('json')
        async def svg(e):await self.export('svg')
        filters=ft.ExpansionTile(title=ft.Text('Filtros visuais e vizinhança'),controls=[
            ft.Container(padding=12,content=ft.Column([self.query,ft.Row(self.bounds,wrap=True),
                ft.Row([ft.TextButton('Aplicar filtros',on_click=self.apply_filters),ft.TextButton('Limpar filtros / voltar à rede',on_click=self.reset_filters)],wrap=True),
                ft.Row([self.hop_choice,ft.TextButton('Ego do composto selecionado',on_click=self.ego_selected)],wrap=True),
                self.focus_only,
                ft.Text('Filtros não alteram arquivos nem recalculam centralidades sobre o subgrafo.',size=12)]))])
        community=ft.ExpansionTile(title=ft.Text('Comunidades e representação agregada'),controls=[ft.Container(padding=12,content=ft.Column([
            ft.TextButton('Detectar comunidades (Louvain)',on_click=detect),ft.Row([self.community_choice]),
            ft.Row([ft.TextButton('Comunidades anteriores',on_click=lambda e:self.page_choices('community',-1)),ft.TextButton('Mais comunidades',on_click=lambda e:self.page_choices('community',1))],wrap=True),
            self.color_communities,self.aggregate_choice,ft.Text('Cores se repetem; os IDs C1, C2… distinguem os grupos. Agregação usa o universo inteiro, antes dos filtros.',size=12)]))])
        # Keep inputs inside bounded Rows; never place flex inputs in wrapping Rows.
        style=[item for item in style if not (isinstance(item,ft.Container) and item.gradient)]
        settings=ft.ExpansionTile(title=ft.Text('Layout, centralidade e componentes'),controls=[ft.Container(padding=12,content=ft.Column([
            *style[:-1],ft.Row([ft.TextButton('Componentes anteriores',on_click=lambda e:self.page_choices('component',-1)),ft.TextButton('Mais componentes',on_click=lambda e:self.page_choices('component',1))],wrap=True),
            ft.Row([self.metric_choice]),ft.Row([self.appearance.attribute,self.appearance.size]),ft.Row([self.appearance.palette]),
            ft.Row([self.appearance.edges,self.appearance.components],wrap=True),self.method],horizontal_alignment=ft.CrossAxisAlignment.STRETCH))])
        from .flow_canvas import icon_button
        self.expand_button=icon_button(icon=ft.Icons.OPEN_IN_FULL,tooltip='Expandir grafo',on_click=self.toggle_focus)
        controls=ft.Row([ft.TextButton('Restaurar grafo completo',on_click=self.restore_full),*tools[1:],self.expand_button,icon_button(icon=ft.Icons.VIEW_SIDEBAR,tooltip='Recolher / mostrar composto',on_click=self.toggle_sidebar),
            icon_button(icon=ft.Icons.DESELECT,tooltip='Limpar seleção',on_click=view.clear_selection)],wrap=True,spacing=4)
        self.graph_box=ft.Container(view.viewer,height=560,bgcolor='#F8FAFD',border_radius=12,clip_behavior=ft.ClipBehavior.HARD_EDGE)
        self.configuration=ft.Column([settings,filters,community],spacing=4)
        self.graph_card=card('Graph Explorer',ft.Column([ft.Row([view.mode_choice]),controls,view.scale_label,self.appearance.legend,self.graph_box,view.status,self.configuration],
            spacing=10,horizontal_alignment=ft.CrossAxisAlignment.STRETCH),{'xs':12,'lg':9})
        self.inspector_card=self.inspector.build()
        graph_section=ft.Column([ft.Row([view.search]),ft.ResponsiveRow([self.graph_card,self.inspector_card],spacing=12,run_spacing=12),self.primary],
            horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
        indicators=ft.Column([self.stats,ft.TextButton('Calcular indicadores',on_click=self.change),self.distributions],horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
        self.rank_metric=ft.Dropdown(label='Métrica do ranking',value=self.metric,expand=True,options=self.metric_choice.options,on_select=self.ranking_metric)
        self.rank_info=icon_button(icon=ft.Icons.INFO_OUTLINE,tooltip='Sobre a métrica do ranking',on_click=self.ranking_help)
        ranking=ft.Column([ft.Row([self.rank_scope,ft.TextButton('Calcular ranking',on_click=self.calculate_ranking),ft.TextButton('Exportar ranking CSV',on_click=self.export_ranking)]),self.rank_caption,self.rank_notice,ft.ResponsiveRow([ft.Container(ft.Row([self.rank_metric,self.rank_info]),col={'xs':12,'md':5}),
            ft.Container(self.rank_query,col={'xs':12,'md':7})]),ft.Row([self.rank_order,self.rank_limit],wrap=True),self.ranking],horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
        # Retain mounted viewers when navigating tabs: Flutter camera state must
        # survive selection through the ranking and returns from analytical views.
        from .graph_table import GraphTable
        self.table=GraphTable(self)
        self.sections=[ft.ExpansionTile(title=ft.Container(height=0),controls=[body],expanded=i==0,maintain_state=True,
            min_tile_height=0,tile_padding=0,controls_padding=0,dense=True,show_trailing_icon=False,
            animation_style=ft.AnimationStyle(duration=0,reverse_duration=0),expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH)
            for i,body in enumerate([graph_section,indicators,self.molecular.build(),ranking,self.table.build()])]
        self.tabs_bar=ft.TabBar(tabs=[ft.Tab(label=title) for title in ('Rede','Indicadores','Análise molecular','Ranking','Laboratório de dados')],scrollable=True)
        self.tabs=ft.Tabs(length=5,on_change=self.switch_tab,content=ft.Column([self.tabs_bar,*self.sections],spacing=0,horizontal_alignment=ft.CrossAxisAlignment.STRETCH))
        self.overview=self.stats
        self.footer=ft.Row([ft.TextButton('Exportar métricas CSV',on_click=csv),ft.TextButton('Relatório JSON',on_click=report),ft.TextButton('Vista SVG',on_click=svg),*view.actions],wrap=True)
        for kind,label in [('png','Vista PNG'),('gexf','Gephi · GEXF'),('graphml','GraphML')]:
            async def save(e,kind=kind):await self.export(kind)
            self.footer.controls.append(ft.TextButton(label,on_click=save))
        self.previous_keyboard=getattr(view.ui.page,'on_keyboard_event',None)
        async def keyboard(e):
            if str(e.key).lower()=='escape' and view.active and self.tab_index==0:
                if self.inspector.dialog is not None:return
                if self.focus_mode:await self.toggle_focus()
                else:view.clear_selection()
            elif self.previous_keyboard:
                import inspect
                result=self.previous_keyboard(e)
                if inspect.isawaitable(result):await result
        self.keyboard_handler=keyboard;view.ui.page.on_keyboard_event=keyboard
        self.previous_resize=getattr(view.ui.page,'on_resize',None)
        async def resized(e):
            if self.previous_resize:
                import inspect
                result=self.previous_resize(e)
                if inspect.isawaitable(result):await result
            if not view.active:return
            height=getattr(view.ui.page,'height',None) or 1000
            self.graph_box.height=max(360,min(1000,height-220)) if self.focus_mode else max(320,min(660,height*.56))
            self.inspector.resize_details()
            view.ui.page.update()
        self.resize_handler=resized;view.ui.page.on_resize=resized
        return ft.Column([view.summary,self.notice,self.tabs,self.footer],spacing=14,horizontal_alignment=ft.CrossAxisAlignment.STRETCH)

    async def ranking_metric(self,e):
        self.rank_metric.value=e.control.value
        await self.calculate_ranking()

    async def calculate_ranking(self,e=None):
        async def action():
            view=self.view
            if not view.active or view.ui.token!=view.token:return
            await view.ui.call(view.ui.store.project,view.token,view.project_id)
            self.rank_request+=1;request=self.rank_request
            (analysis,mode),key=self.ranking_context();metric=self.rank_metric.value
            self.rank_notice.value='Calculando ranking…';view.ui.page.update()
            await view.ui.call(analysis.metrics,mode,metric)
            if request!=self.rank_request or not view.active or view.ui.token!=view.token:return
            _,current_key=self.ranking_context()
            if key!=current_key or metric!=self.rank_metric.value:return
            self.refresh_ranking();view.ui.page.update()
        await self.view.ui.guard(action)

    async def export_ranking(self,e=None):
        async def action():
            view=self.view
            if not view.active or view.ui.token!=view.token:return
            await view.ui.call(view.ui.store.project,view.token,view.project_id)
            (analysis,mode),key=self.ranking_context();metric=self.rank_metric.value
            data=await view.ui.call(analysis.export_csv,mode,metric,f'{self.rank_scope.value}:{view.mode}')
            if not view.active or view.ui.token!=view.token:return
            _,current_key=self.ranking_context()
            if current_key!=key or metric!=self.rank_metric.value:return
            await save_file(view.ui,file_name=download_name(view.data.get('dataset_origin',{}),result_type='Ranking-'+str(metric),extension='csv',qualifier=[self.rank_scope.value,view.mode,metric],operation='graphs'),src_bytes=data)
        await self.view.ui.guard(action)

    async def ranking_help(self,e=None):
        from .ranking_help import HELP,toy_values
        from .feedback import close_dialog
        async def action():
            view=self.view
            if not view.active or view.ui.token!=view.token:return
            await view.ui.call(view.ui.store.project,view.token,view.project_id)
            metric=self.rank_metric.value or self.metric
            meaning,purpose,example,limits=HELP[metric]
            values=await view.ui.call(toy_values,metric)
            if not view.active or view.ui.token!=view.token:return
            controls=[]
            for heading,body in [('O que mede',meaning),('Como usar na análise',purpose),('Exemplo simples',example),('Critérios e limites',limits)]:
                controls.extend([ft.Text(heading,size=15,weight=ft.FontWeight.W_600),ft.Text(body,size=13)])
                if heading=='Exemplo simples':
                    controls.extend([verbatim(ft.Text('A — B — C',size=24,weight=ft.FontWeight.W_600)),
                        verbatim(ft.Text(' · '.join(f'{node}: {values[node]:.3f}' for node in ('A','B','C')),size=13))])
            controls.append(ft.Text('Exemplo ilustrativo, sem dados do projeto. O ranking usa o escopo escolhido: grafo exibido com filtros ou universo Full / MCC. As cores do grafo mantêm o universo indicado na legenda. Centralidade não indica potência, segurança ou atividade biológica.',size=12,color='#64748B'))
            page=view.ui.page
            dialog=ft.AlertDialog(title=ft.Text(METRICS[metric].split(' · ')[0],size=22),
                content=ft.Container(width=max(260,min(620,(getattr(page,'width',None) or 1000)-80)),
                    height=max(240,min(560,(getattr(page,'height',None) or 900)-180)),
                    content=ft.Column(controls,spacing=12,scroll=ft.ScrollMode.AUTO)),
                actions=[ft.TextButton('Fechar',on_click=lambda e:close_dialog(page,dialog))])
            page.show_dialog(dialog)
        await self.view.ui.guard(action)
