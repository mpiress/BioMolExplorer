"""Native graph/EGG exploration with the shared molecular 3D viewer."""
import json
import asyncio
import base64

import flet as ft
import flet.canvas as canvas

from biomolexplorer.diagnostics import get_logger
from biomolexplorer.visualizations import PointIndex, molecule_image, separated_graph_layout
from .feedback import close_dialog
from .zoom import zoomable_view
from .localization import verbatim
from biomolexplorer.egg_plot import bounds as egg_bounds,REGIONS,COLORS,tick_values

BLUE = '#2552E8'
INK = '#172B4D'
VIRIDIS = ['#440154','#3B528B','#21918C','#5EC962','#FDE725']
GRAPH_NODE_RADIUS = 4


def degree_color(value,maximum):
    from .graph_appearance import continuous_color
    return continuous_color(value,maximum)


class ResultsViewer:
    WIDTH, HEIGHT = 720, 460
    VIEWPORT_HEIGHT = 460

    def __init__(self, ui, project_id, data, popup=False, actions=None):
        self.ui, self.project_id, self.data = ui, project_id, data
        self.token=ui.token;self.active=True;self.popup=popup;self.actions=actions or []
        self.nodes = {node['id']: node for node in data['nodes']}
        self.molecular_enrichment=None
        self.component_id = None
        self.components = {}
        self.explorer=None
        if data['kind'] == 'graph':
            import networkx as nx
            from biomolexplorer.network_analysis import NetworkAnalysis
            self.analysis=NetworkAnalysis(data)
            self.graph=self.analysis.graph
            groups = sorted(nx.connected_components(self.graph), key=lambda group: (-len(group), min(group)))
            self.components = {index+1: group for index, group in enumerate(groups)}
            # Old saved results benefit from the new component layout as well.
            saved={key:(node['x'],node['y']) for key,node in self.nodes.items()}
            self.network_positions=(saved if data.get('layout')=='components' and len(set(saved.values()))==len(saved) else separated_graph_layout(self.graph,circular=True))
            self.circular_positions = None
        self.mode = 'full'
        self.opening_notice=''
        self.layout_mode='network';self.color_mode='degree';self.labels=len(self.nodes)<=20
        self.viewport_width=720;self.viewport_height=560;self.viewport_ready=False;self.camera_request=0
        self.camera_lock=asyncio.Lock()
        self.selected = None
        self.selection_version = 0
        self.hovered = None
        self.status = ft.Text('Passe o mouse sobre um nó. Clique para abrir o composto.', size=12, color='#64748B')
        self.details = ft.Column([ft.Text('Informações do composto', size=19, weight=ft.FontWeight.W_600),
                                  ft.Text('Selecione um nó ou busque seu código.', size=13, color='#64748B')], spacing=14, scroll=ft.ScrollMode.AUTO)
        self.tip_text = ft.Text('', color='#FFFFFF', size=12)
        self.tip = ft.Container(self.tip_text, visible=False, bgcolor=INK, padding=10, border_radius=8,
                                ignore_interactions=True, left=0, top=0)
        egg=data['kind']=='egg'
        self.drawing = canvas.Canvas(width=None if egg else self.WIDTH,height=None if egg else self.HEIGHT,
            expand=egg,on_resize=self.resize_egg if egg else None,resize_interval=100)
        self.raster=ft.Image(src='',visible=False,fit=ft.BoxFit.FILL)
        self.raster_layer=ft.Container(self.raster,ignore_interactions=True,width=self.WIDTH,height=self.HEIGHT,left=0,top=0)
        self.scene = ft.Stack([self.raster_layer,ft.GestureDetector(content=self.drawing, hover_interval=60,
                             on_hover=self.hover, on_exit=self.leave, on_tap_up=self.click,expand=egg), self.tip],
                             width=None if egg else self.WIDTH,height=None if egg else self.HEIGHT,expand=egg)
        self.viewer = zoomable_view(content=self.scene, constrained=egg,
                                          alignment=ft.Alignment.CENTER if egg else ft.Alignment.TOP_LEFT, expand=True)
        if not egg:self.viewer.on_size_change=self.resize_graph
        self.search = ft.TextField(label='Código do composto', hint_text='CHEMBL… ou seu identificador',
                                   expand=True, on_submit=self.find)
        self.summary = ft.Text('', size=13, color='#64748B')
        self.scale_label=ft.Text(size=12,color='#64748B')
        self.fragment_panel=ft.Column(spacing=14,visible=False)
        self.populate_fragment()
        if data['kind']=='graph':
            from .graph_explorer import GraphExplorer
            self.explorer=GraphExplorer(self)
        self.redraw()

    async def initialize(self):
        if self.explorer:await self.explorer.initialize()

    async def resize_graph(self,event):
        if not self.active or event.width<120 or event.height<120:return
        if self.explorer and hasattr(self.explorer,'graph_box'):
            page_height=getattr(self.ui.page,'height',None) or 1000
            target=max(360,min(1000,page_height-220)) if self.explorer.focus_mode else max(320,min(660,page_height*.56))
            if abs(self.explorer.graph_box.height-target)>1:
                self.explorer.graph_box.height=target;self.ui.page.update();return
        changed=not self.viewport_ready or abs(self.viewport_width-event.width)>1 or abs(self.viewport_height-event.height)>1
        self.viewport_width,self.viewport_height=event.width,event.height;self.viewport_ready=True
        if changed:await self.fit_view()

    async def fit_view(self,e=None):
        if not self.active:return
        self.camera_request+=1;request=self.camera_request
        async with self.camera_lock:
            if request!=self.camera_request or not self.active:return
            # Flush the new scene dimensions and raster before imperative camera
            # calls; otherwise a community-sized client scene is still mounted.
            self.ui.page.update()
            await self.viewer.reset()
            if self.data['kind']=='egg' or not self.viewport_ready:return
            points=list(self.points.values())
            if not points and self.explorer and self.explorer.aggregate_mode:
                points=[(s.x,s.y) for s in self.base_shapes if isinstance(s,canvas.Circle)]
            if not points:return
            xmin=min(x for x,y in points);xmax=max(x for x,y in points)
            ymin=min(y for x,y in points);ymax=max(y for x,y in points)
            scale=min((self.viewport_width-64)/max(48,xmax-xmin+24),(self.viewport_height-64)/max(48,ymax-ymin+24),2)
            scale=max(.001,scale)
            await self.viewer.zoom(scale)
            await self.viewer.pan(self.viewport_width/(2*scale)-(xmin+xmax)/2,self.viewport_height/(2*scale)-(ymin+ymax)/2)

    def clear_selection(self,e=None):
        self.selection_version+=1;self.selected=self.hovered=None;self.tip.visible=False
        self.details.controls=[ft.Text('Selecione um composto para ver suas informações.')]
        if self.explorer:self.explorer.inspector.clear()
        self.highlight();self.ui.page.update()

    def resize_egg(self,event):
        if not self.active or event.width<120 or event.height<120:return
        if abs(self.WIDTH-event.width)<1 and abs(self.HEIGHT-event.height)<1:return
        self.WIDTH,self.HEIGHT=event.width,event.height
        self.hovered=None;self.tip.visible=False;self.redraw();self.ui.page.update()

    def populate_fragment(self):
        if self.data['kind']!='graph':return
        fragment=self.data.get('fragment',{})
        if not fragment:
            self.fragment_panel.controls=[ft.Text('Fragmento comum não registrado neste resultado. A abertura não executa uma nova busca MCS.',size=12)]
            return
        status=fragment.get('status','missing_structures')
        notices={'complete':'Busca concluída · comum a todos os compostos do MCC.',
            'partial':'Tempo limite atingido. Este é o melhor fragmento encontrado; seu tamanho máximo não está confirmado.',
            'missing_structures':'SMILES ausentes: o fragmento não está disponível.',
            'empty':'Não há compostos no MCC.','no_common_fragment':'Não foi encontrado um fragmento sob os critérios escolhidos.'}
        controls=[ft.Text('Fragmento comum do MCC',size=17,weight=ft.FontWeight.W_600),
            ft.Text(notices[status],size=12,color='#B45309' if status=='partial' else '#64748B')]
        if fragment.get('image'):
            smiles = fragment.get('smiles', '')
            if 'smiles' not in fragment and fragment.get('smarts'):
                # Historical artifacts only stored the SMARTS query. Extract
                # their displayed fragment from an MCC reference structure.
                from caad.graph_results import fragment_molecule
                from rdkit import Chem
                mcc = set(self.data['mcc'])
                for identifier in self.nodes:
                    if identifier not in mcc:
                        continue
                    reference = self.nodes[identifier]['properties'].get('canonical_smiles')
                    molecule = fragment_molecule(fragment['smarts'], reference)
                    if molecule is not None:
                        smiles = Chem.MolToSmiles(molecule)
                        break
            controls += [ft.Image(src=base64.b64decode(fragment['image']),width=320,height=225,fit=ft.BoxFit.CONTAIN),
                ft.Text(f"{fragment['atoms']} átomos · {fragment['bonds']} ligações · {fragment['compounds']} compostos",size=12),
                ft.TextField(label='SMILES do fragmento',value=smiles,read_only=True,multiline=True,min_lines=2,max_lines=4),
                ft.Text('Estrutura extraída da molécula de referência do MCC.',size=12,color='#64748B')]
        self.fragment_panel.controls=controls

    def displayed(self):
        if self.explorer:
            graph=self.explorer.graph()
            return {key:self.nodes[key] for key in graph}
        if self.mode == 'mcc':
            allowed = set(self.data['mcc'])
            return {key: value for key, value in self.nodes.items() if key in allowed}
        if self.component_id is not None:
            allowed = self.components.get(self.component_id, set())
            return {key: value for key, value in self.nodes.items() if key in allowed}
        return self.nodes

    def coordinates(self, nodes):
        if self.data['kind'] == 'egg':
            xmin,xmax,ymin,ymax=egg_bounds(list(nodes.values()))
            # All points, including outliers, remain inside the plot area.
            dx, dy = max(1, xmax-xmin), max(1, ymax-ymin)
            self.project = lambda x, y: (65+(x-xmin)/dx*(self.WIDTH-105),
                                         self.HEIGHT-55-(y-ymin)/dy*(self.HEIGHT-90))
            self.bounds = (xmin, xmax, ymin, ymax)
        else:
            positions = self.network_positions
            if self.layout_mode=='circular':
                if self.circular_positions is None:
                    self.circular_positions = separated_graph_layout(self.graph, circular=True)
                positions = self.circular_positions
            elif self.layout_mode!='network' and self.explorer and self.explorer.positions is not None:
                positions=self.explorer.positions
            # A pending/cost-limited layout may only cover the previous scope.
            # Complete it from the persisted full layout instead of aborting paint.
            if not set(nodes)<=positions.keys():positions={**self.network_positions,**positions}
            nodes={key:dict(node,x=positions[key][0],y=positions[key][1]) for key,node in nodes.items()}
            xmin = min((node['x'] for node in nodes.values()), default=-1)
            xmax = max((node['x'] for node in nodes.values()), default=1)
            ymin = min((node['y'] for node in nodes.values()), default=-1)
            ymax = max((node['y'] for node in nodes.values()), default=1)
            dx, dy = xmax-xmin, ymax-ymin
            # Preserve aspect ratio and a minimum 30 px distance between nodes.
            # Dense networks grow beyond the viewport and remain pannable.
            scale = max(30, min((720-110)/max(1, dx), (460-90)/max(1, dy)))
            if self.explorer:
                self.explorer.geometry_notice=''
                if dx*scale+110>4096 or dy*scale+90>4096:
                    scale=min(scale,(4096-110)/max(1,dx),(4096-90)/max(1,dy))
                    self.explorer.geometry_notice='Layout limitado a 4.096 px por dimensão. Use zoom ou ego para inspecionar regiões densas.'
            self.WIDTH, self.HEIGHT = max(720, dx*scale+110), max(460, dy*scale+90)
            self.scene.width = self.drawing.width = self.WIDTH
            self.scene.height = self.drawing.height = self.HEIGHT
            self.raster_layer.width,self.raster_layer.height=self.WIDTH,self.HEIGHT
            self.project = lambda x, y: ((self.WIDTH-dx*scale)/2+(x-xmin)*scale,
                                         (self.HEIGHT-dy*scale)/2+(y-ymin)*scale)
        return {key: self.project(node['x'], node['y']) for key, node in nodes.items()}

    def egg_background(self):
        shapes=[];xmin,xmax,ymin,ymax=self.bounds
        for value in tick_values(xmin,xmax):
            x,y=self.project(value,ymin)
            shapes.append(canvas.Line(x,35,x,self.HEIGHT-55,paint=ft.Paint(color='#DCE5EF',stroke_width=.6)))
            shapes.append(verbatim(canvas.Text(x-12,y+12,f'{value:g}',style=ft.TextStyle(size=12,color=INK))))
        for value in tick_values(ymin,ymax):
            x,y=self.project(xmin,value)
            shapes.append(canvas.Line(65,y,self.WIDTH-40,y,paint=ft.Paint(color='#DCE5EF',stroke_width=.6)))
            shapes.append(verbatim(canvas.Text(14,y-7,f'{value:g}',style=ft.TextStyle(size=12,color=INK))))
        for cx,cy,rx,ry,color in REGIONS:
            x1,y1=self.project(cx-rx,cy+ry);x2,y2=self.project(cx+rx,cy-ry)
            shapes.append(canvas.Oval(x1,y1,x2-x1,y2-y1,paint=ft.Paint(color=color)))
            shapes.append(canvas.Oval(x1,y1,x2-x1,y2-y1,paint=ft.Paint(color='#94A3B8',stroke_width=1,style=ft.PaintingStyle.STROKE)))
        shapes += [canvas.Line(65,35,65,self.HEIGHT-55,paint=ft.Paint(color='#64748B',stroke_width=1)),
            canvas.Line(65,self.HEIGHT-55,self.WIDTH-40,self.HEIGHT-55,paint=ft.Paint(color='#64748B',stroke_width=1)),
            verbatim(canvas.Text(self.WIDTH/2-35,self.HEIGHT-24,'TPSA (Å²)',style=ft.TextStyle(size=14,color=INK))),
            verbatim(canvas.Text(12,8,'WLOGP',style=ft.TextStyle(size=14,color=INK)))]
        return shapes

    def redraw(self):
        if self.explorer and self.explorer.aggregate_mode:
            self.explorer.draw_aggregate();return
        nodes = self.displayed()
        self.points = self.coordinates(nodes)
        self.index = PointIndex(self.points)
        shapes = self.egg_background() if self.data['kind'] == 'egg' else []
        if self.data['kind'] == 'graph' and self.mode == 'full' and self.component_id is None and self.layout_mode in ('network','radial','fr_components','fa2','fa2_linlog') and len(self.components)<=200 and self.explorer.appearance.components.value:
            for component, identifiers in self.components.items():
                points = [self.points[identifier] for identifier in identifiers if identifier in self.points]
                if not points:continue
                x, y = min(point[0] for point in points)-18, min(point[1] for point in points)-26
                width = max(point[0] for point in points)-x+18
                height = max(point[1] for point in points)-y+18
                shapes.append(canvas.Rect(x, y, width, height, border_radius=10,
                    paint=ft.Paint(color='#E2E8F0', stroke_width=1, style=ft.PaintingStyle.STROKE)))
                shapes.append(canvas.Text(x+8, y+5, f'Componente {component}',
                    style=ft.TextStyle(size=10, color='#64748B')))
        edges = [edge for edge in self.data['edges'] if edge['source'] in nodes and edge['target'] in nodes]
        if self.explorer:
            edges=[dict(source=a,target=b,**attrs) for a,b,attrs in self.explorer.filtered_graph.edges(data=True) if a in nodes and b in nodes]
        for edge in edges:
            if self.explorer and not self.explorer.appearance.edges.value:continue
            a, b = self.points[edge['source']], self.points[edge['target']]
            shapes.append(canvas.Line(*a, *b, paint=ft.Paint(color='#CBD5E1', stroke_width=.6+min(1,abs(float(edge.get('value',0))))*.7)))
        degrees={key:int(node['properties'].get('degree',0)) for key,node in nodes.items()}
        values={key:degree if self.color_mode=='degree' else degree/max(1,len(nodes)-1) for key,degree in degrees.items()}
        if self.explorer:
            metric=self.explorer.metric if self.explorer.metric!='degree' else ('density' if self.color_mode=='density' else 'degree')
            values=self.analysis.metrics(self.mode,metric)[0]
        maximum=max(values.values(),default=0)
        mcc=set(self.data['mcc'])
        self.radii={}
        if self.explorer:
            colors,radii=self.explorer.appearance.map(nodes,values,metric)
        for identifier, point in self.points.items():
            props = nodes[identifier]['properties']
            color = COLORS.get(props.get('BBB'),COLORS['unknown']) if self.data['kind'] == 'egg' else degree_color(values.get(identifier,0),maximum)
            if self.explorer:
                color=colors[identifier]
            radius=5 if self.data['kind']=='egg' else GRAPH_NODE_RADIUS
            if self.explorer:radius=radii[identifier]
            self.radii[identifier]=radius
            if self.data['kind']=='graph' and self.mode=='full' and identifier in mcc:
                shapes.append(canvas.Circle(*point,radius=radius+1.5,paint=ft.Paint(color=BLUE,stroke_width=1,style=ft.PaintingStyle.STROKE)))
            if self.data['kind']=='egg' and props.get('BBB')=='BBB-':
                shapes.append(canvas.Rect(point[0]-radius,point[1]-radius,2*radius,2*radius,paint=ft.Paint(color=color)))
            else:shapes.append(canvas.Circle(*point, radius=radius, paint=ft.Paint(color=color)))
            if self.data['kind']=='egg':
                shapes.append(canvas.Circle(*point,radius=radius+1,paint=ft.Paint(color='#FFFFFF',stroke_width=.8,style=ft.PaintingStyle.STROKE)))
            if self.data['kind']=='graph' and self.labels:
                shapes.append(verbatim(canvas.Text(point[0]+radius+4,point[1]-8,identifier,style=ft.TextStyle(size=10,color=INK))))
        self.base_shapes = shapes
        self.raster.visible=bool(self.explorer and (len(nodes)>self.analysis.DRAW_NODES or len(edges)>self.analysis.DRAW_EDGES))
        if self.raster.visible:
            from .graph_raster import paint_graph
            self.raster.width,self.raster.height=self.WIDTH,self.HEIGHT
            self.raster.src=paint_graph(shapes,self.WIDTH,self.HEIGHT)
        else:self.raster.src=''
        self.highlight()
        self.summary.value = (f'{len(nodes)} compostos · {len(edges)} relações · {len(self.components)} componentes · MCC: {len(mcc)} de {len(self.nodes)} compostos' if self.data['kind'] == 'graph' else f'{len(nodes)} compostos · BBB+: {sum(n["properties"].get("BBB")=="BBB+" for n in nodes.values())} · BBB−: {sum(n["properties"].get("BBB")=="BBB-" for n in nodes.values())}')
        self.scale_label.value=('Grau: 0 → '+str(int(maximum))+' conexões' if self.color_mode=='degree' else f'Conectividade: 0 → {maximum:.3f} · grau / (n − 1)')
        if self.explorer:
            statistics=self.analysis.statistics(self.mode)
            self.summary.value=f'{"Full Graph" if self.mode=="full" else "LCC / MCC"} · {statistics["nodes"]} compostos · {statistics["edges"]} relações · {statistics["components"]} componentes · Exibição: {len(nodes)} compostos / {len(edges)} relações'
            self.scale_label.value=self.explorer.appearance.size_label+' · Cor: '+self.explorer.appearance.title+' · Espessura: magnitude do escore disponível.'
            if not values and nodes:self.scale_label.value+=f' {metric}: indisponível; tamanho fixo.'
            self.explorer.update_panels()
        self.fragment_panel.visible=self.explorer is None and self.mode=='mcc' and self.data['kind']=='graph'
        if not nodes:
            self.status.value = 'Nenhum composto disponível nesta visualização.'

    def highlight(self):
        shapes = [] if self.raster.visible else list(self.base_shapes)
        if self.selected in self.points:
            if self.data['kind'] == 'graph':
                network=self.explorer.filtered_graph if self.explorer else self.graph
                for neighbor in network.neighbors(self.selected):
                    if neighbor in self.points:
                        shapes.append(canvas.Line(*self.points[self.selected], *self.points[neighbor],
                            paint=ft.Paint(color=BLUE, stroke_width=2)))
                        shapes.append(canvas.Circle(*self.points[neighbor], radius=self.radii[neighbor]+3,
                            paint=ft.Paint(color=BLUE, stroke_width=1.5, style=ft.PaintingStyle.STROKE)))
            shapes.append(canvas.Circle(*self.points[self.selected], radius=self.radii[self.selected]+5,
                          paint=ft.Paint(color=INK, stroke_width=2, style=ft.PaintingStyle.STROKE)))
        self.drawing.shapes = shapes

    def hover(self, event):
        if not self.active or self.ui.token!=self.token:return
        position = event.local_position
        if position is None:
            return
        identifier = self.index.nearest(position.x, position.y)
        if self.explorer and self.explorer.aggregate_mode:
            self.tip.visible=identifier is not None
            if identifier is not None:
                self.tip_text.value=f'Comunidade C{identifier} · clique para investigar seus compostos'
                self.tip.left=min(position.x+12,self.WIDTH-190);self.tip.top=max(0,position.y-42)
            self.ui.page.update();return
        if identifier == self.hovered:
            return
        self.hovered = identifier
        self.tip.visible = identifier is not None
        if identifier:
            overlapping = self.index.nearby(*self.points[identifier])
            self.tip_text.value = identifier + (f' (+{len(overlapping)-1})' if len(overlapping)>1 else '')
            if self.data['kind']=='egg':
                props=self.nodes[identifier]['properties']
                self.tip_text.value+=f'\nTPSA: {self.nodes[identifier]["x"]:.2f} Å² · WLOGP: {self.nodes[identifier]["y"]:.2f}\n'+str(props.get('BBB','?'))+' · '+str(props.get('HIA','?'))
            self.tip.left = min(position.x+12, self.WIDTH-190)
            self.tip.top = max(0, position.y-42)
            self.status.value = 'Composto: '+identifier+' · clique para ver a estrutura e as propriedades.'
            if self.data['kind']=='graph':self.status.value+=f" · grau: {self.graph.degree(identifier)}"
            if self.explorer:
                values,_=self.analysis.metrics(self.mode,self.explorer.metric)
                community=self.analysis.community_cache.get(self.mode,({},None,''))[0].get(identifier)
                self.tip_text.value+='\n'+self.mode+' · '+self.explorer.metric+': '+str(values.get(identifier,'indisponível'))
                if community is not None:self.tip_text.value+=f' · C{community}'
        else:
            self.status.value = 'Passe o mouse sobre um nó. Clique para abrir o composto.'
        self.ui.page.update()

    def leave(self, event):
        self.hovered = None
        self.tip.visible = False
        self.ui.page.update()

    async def click(self, event):
        position = event.local_position
        if position is not None:
            identifier = self.index.nearest(position.x, position.y)
            if identifier:
                if self.explorer and self.explorer.aggregate_mode:
                    await self.ui.guard(lambda:self.explorer.open_community(identifier))
                else:await self.ui.guard(lambda: self.pick(identifier))
            elif self.active and self.ui.token==self.token:self.clear_selection()

    async def pick(self, identifier):
        await self.ui.call(self.ui.store.project, self.ui.token, self.project_id)
        overlapping = self.index.nearby(*self.points[identifier])
        if len(overlapping) == 1:
            return await self.select(identifier)
        choices=ft.Column(spacing=6,scroll=ft.ScrollMode.AUTO)
        query=ft.TextField(label='Pesquisar composto neste ponto')
        offset=0
        def populate(reset=False):
            nonlocal offset
            if reset:offset=0;choices.controls=[]
            found=[code for code in overlapping if (query.value or '').casefold() in code.casefold()]
            for compound_id in found[offset:offset+50]:
                async def select_compound(event, compound_id=compound_id):
                    close_dialog(self.ui.page,dialog)
                    await self.ui.guard(lambda:self.select(compound_id))
                choices.controls.append(verbatim(ft.TextButton(compound_id,on_click=select_compound)))
            offset+=50;more.visible=offset<len(found)
        def search(e):populate(True);self.ui.page.update()
        def load_more(e):populate();self.ui.page.update()
        query.on_change=search
        more=ft.TextButton('Carregar mais compostos',on_click=load_more)
        populate(True)
        dialog=ft.AlertDialog(title=ft.Text('Compostos neste ponto'),
            content=ft.Container(width=360,height=min(460,len(overlapping)*50+90),
                                 content=ft.Column([query,choices,more],scroll=ft.ScrollMode.AUTO)),
            actions=[ft.TextButton('Fechar',on_click=lambda e:close_dialog(self.ui.page,dialog))])
        self.ui.page.show_dialog(dialog)

    async def find(self, event):
        async def action():
            identifier = (self.search.value or '').strip()
            matches = [key for key in self.nodes if key.casefold() == identifier.casefold()]
            if not matches:
                raise ValueError('Código não encontrado nesta visualização. No grafo, confira se você está no MCC ou no grafo completo.')
            chosen = matches[0]
            component_focus=self.component_id
            if self.explorer:
                self.explorer.aggregate_mode=False;self.explorer.aggregate_choice.value=False
                await self.explorer.reset_filters()
            if chosen not in self.displayed() or component_focus is not None:
                self.mode = 'full'
                self.mode_choice.value = 'full'
                self.component_choice.disabled = False
                self.component_id = next(index for index, group in self.components.items() if chosen in group)
                self.component_choice.value = str(self.component_id)
                if self.explorer:
                    await self.explorer.prepare()
                    if self.layout_mode not in ('network','circular'):
                        self.explorer.positions,_=await self.ui.call(self.analysis.layout,self.mode,self.layout_mode)
                self.redraw()
                await self.viewer.reset()
            if self.explorer and chosen not in self.displayed():
                self.explorer.ego=chosen;self.explorer.hop_choice.value='1';self.redraw()
                if chosen not in self.displayed():
                    self.explorer.focus_only.value=True;self.redraw()
                    self.ui.notify('A vizinhança excede o orçamento. Exibindo somente o composto selecionado; ajuste os filtros para expandir.')
                else:self.ui.notify('Exibindo ego de 1 salto para manter a exploração dentro do orçamento.')
            await self.select(chosen)
            await self.center_on(chosen)
        await self.ui.guard(action)

    async def show_3d(self, identifier, version):
        def current():
            project = getattr(self.ui, 'current', None)
            return (self.active and self.ui.token == self.token
                    and version == self.selection_version and self.selected == identifier
                    and project is not None and project['id'] == self.project_id
                    and identifier in self.displayed())

        if not current():
            return
        await self.ui.call(self.ui.store.project, self.token, self.project_id)
        if not current():
            return
        node=self.displayed()[identifier]
        properties=node['properties']
        if self.molecular_enrichment:
            from biomolexplorer.molecular_network import numeric_properties
            properties=numeric_properties([node],self.molecular_enrichment)[0][identifier]
        smiles = properties.get('canonical_smiles')
        url = await self.ui.call(self.ui.compound_view_url, self.project_id, smiles, identifier, self.token)
        if not current():
            return
        await self.ui.call(self.ui.store.project, self.token, self.project_id)
        if not current():
            return
        await ft.UrlLauncher().launch_url(url, mode=ft.LaunchMode.EXTERNAL_APPLICATION,
                                         web_only_window_name='_blank')

    async def select(self, identifier, show_popup=None):
        popup=(self.popup and not self.explorer) if show_popup is None else show_popup
        if show_popup is None and self.explorer and (self.explorer.focus_mode or self.explorer.sidebar_hidden):popup=True
        # Membership may have changed since opening this dialog.
        self.selection_version += 1
        version, token = self.selection_version, self.ui.token
        await self.ui.call(self.ui.store.project, token, self.project_id)
        node = self.displayed().get(identifier)
        if node is None:
            raise ValueError('O composto não pertence à visualização atual.')
        if self.explorer:
            from biomolexplorer.molecular_network import numeric_properties
            node=dict(node,properties=numeric_properties([node],self.molecular_enrichment)[0][identifier])
            from biomolexplorer.compound_properties import compound_properties
            properties=await self.ui.call(compound_properties,node['properties'],True)
            node=dict(node,properties=properties)
        image = None
        structure_notice = 'Estrutura indisponível: SMILES ausente ou inválido.'
        try:
            image = await self.ui.call(molecule_image, node['properties'].get('canonical_smiles'))
        except ImportError:
            get_logger('frontend').warning('RDKit indisponível para renderizar a estrutura molecular.')
            structure_notice = 'A estrutura 2D está indisponível neste ambiente. Consulte o SMILES e as propriedades abaixo.'
        if version != self.selection_version or token != self.ui.token:
            return
        if not self.active or token!=self.token or (getattr(self.ui,'current',None) is not None and self.ui.current['id']!=self.project_id):return
        # Structure rendering runs off the UI thread; recheck after that await.
        await self.ui.call(self.ui.store.project, token, self.project_id)
        if version != self.selection_version or token != self.ui.token:
            return
        controls = [verbatim(ft.Text(identifier, size=22, weight=ft.FontWeight.W_700, selectable=True))]
        controls.append(ft.Image(src=image, width=280, height=205, fit=ft.BoxFit.CONTAIN) if image else ft.Text(structure_notice, size=12))
        async def open_3d(event):
            await self.ui.guard(lambda: self.show_3d(identifier, version))
        controls.append(ft.TextButton('Visualizar estrutura 3D', on_click=open_3d, disabled=not bool(image)))
        properties=node['properties']
        if self.explorer:
            from biomolexplorer.molecular_network import numeric_properties
            properties=numeric_properties([node],self.molecular_enrichment)[0][identifier]
            if self.molecular_enrichment:controls.append(ft.Text('Fonte molecular adicional: '+str(self.molecular_enrichment['source']['stage_name'])+' · '+str(self.molecular_enrichment['source']['run_id']),size=11))
        for key, value in properties.items():
            if key in ('molecule_properties','properties') and any(k.startswith(key+'.') for k in properties):continue
            from biomolexplorer.admet_config import DESCRIPTORS
            label = DESCRIPTORS.get(key,str(key)) if self.data['kind']=='egg' else str(key)
            if self.data['kind']=='graph' and key in ('degree','degree_density','component'):label+=' · registrado no Full Graph'
            if str(key).startswith('ML_'):
                definition=(self.data.get('admet_report') or {}).get('provider',{}).get('endpoint_definitions',{}).get(str(key)[3:],{})
                label=str(key)+' · '+definition.get('units','unidade do endpoint')
            if self.data['kind']=='graph' and self.molecular_enrichment:
                label+=' · '+('grafo' if self.nodes[identifier]['properties'].get(key.split('.')[0]) is not None else 'fonte selecionada')
            rendered = '—' if value is None else json.dumps(value, ensure_ascii=False) if isinstance(value, (dict, list)) else f'{value:.5g}' if isinstance(value,float) else str(value)
            controls.append(ft.Column([ft.Text(label, size=11, color='#64748B'),
                                       verbatim(ft.Text(rendered, size=13, selectable=True))], spacing=4))
        if self.data['kind'] == 'graph':
            if self.explorer:
                metric=self.explorer.metric;values,notice=self.analysis.metrics(self.mode,metric)
                mapping=self.analysis.community_cache.get(self.mode,({},None,''))[0]
                controls.append(ft.Text(f'{self.mode} · {metric}: '+str(values.get(identifier,'indisponível'))+' · Comunidade: '+str(mapping.get(identifier,'não calculada'))))
            neighbors = sorted(self.graph[identifier].items(), key=lambda item: (-item[1].get('value',0), item[0]))
            controls.append(ft.Text(f'Relações deste composto ({len(neighbors)})', size=15, weight=ft.FontWeight.W_600))
            def relation_button(neighbor,attributes):
                async def follow(event, target=neighbor):
                    if self.explorer:self.explorer.inspector.close_details()
                    elif popup:
                        close_dialog(self.ui.page, dialog)
                    await self.ui.guard(lambda: self.follow_relation(target))
                return ft.TextButton(f"{neighbor} · escore "+(f"{attributes['value']:.3f}" if 'value' in attributes else 'indisponível'), on_click=follow,tooltip=json.dumps(attributes,ensure_ascii=False))
            controls.extend(relation_button(n,a) for n,a in neighbors[:50])
            if len(neighbors)>50:
                offset=50
                async def more_relations(event):
                    nonlocal offset
                    await self.ui.call(self.ui.store.project,token,self.project_id)
                    if not self.active or version!=self.selection_version or token!=self.ui.token:return
                    index=controls.index(more_button)
                    controls[index:index]=[relation_button(n,a) for n,a in neighbors[offset:offset+50]]
                    offset=min(len(neighbors),offset+50)
                    more_button.content=f'Carregar mais relações ({offset} / {len(neighbors)})'
                    more_button.disabled=offset>=len(neighbors)
                    if self.explorer:self.explorer.inspector.refresh_relations(controls)
                    self.ui.page.update()
                more_button=ft.TextButton(f'Carregar mais relações (50 / {len(neighbors)})',on_click=more_relations)
                controls.append(more_button)
        self.details.controls = controls
        self.selected = identifier
        if self.explorer:await self.explorer.inspector.update(identifier,properties,version)
        if version!=self.selection_version or not self.active:return
        self.highlight()
        self.ui.page.update()
        if popup and self.explorer:
            await self.explorer.inspector.details()
        elif popup:
            dialog=ft.AlertDialog(title=ft.Text(identifier+' · molécula 2D',size=22),
                content=ft.Container(width=460,height=520,content=ft.Column(controls,spacing=14,scroll=ft.ScrollMode.AUTO)),
                actions=[ft.TextButton('Fechar molécula',on_click=lambda e:close_dialog(self.ui.page,dialog))])
            self.ui.page.show_dialog(dialog)

    async def change_mode(self, event):
        async def action():
            if not self.active or self.ui.token!=self.token:return
            await self.ui.call(self.ui.store.project, self.ui.token, self.project_id)
            self.mode = event.control.value
            self.mode_choice.value=self.mode
            if self.explorer:
                self.explorer.inspector.invalidate()
                self.explorer.inspector.scope.value=self.mode
                self.explorer.inspector.describe_scope()
                self.explorer.clear_filters()
                request=self.explorer.request;mode=self.mode
                self.explorer.molecular.drawing.shapes=[];self.explorer.molecular.bars.controls=[];self.explorer.molecular.comparison.controls=[]
                self.explorer.molecular.points={};self.explorer.molecular.scatter_pairs=[];self.explorer.molecular.index=PointIndex({})
                self.explorer.molecular.summary.value='Universo alterado; clique em Analisar descritores para atualizar os gráficos moleculares.'
                await self.explorer.prepare()
                if request!=self.explorer.request or mode!=self.mode or not self.active or self.ui.token!=self.token:return
                layout_notice=''
                if self.layout_mode not in ('network','circular'):
                    positions,layout_notice=await self.ui.call(self.analysis.layout,mode,self.layout_mode)
                else:positions=None
                if request!=self.explorer.request or mode!=self.mode or not self.active or self.ui.token!=self.token:return
                self.explorer.positions=positions
                self.explorer.layout_notice=layout_notice
                if self.layout_mode not in ('network','circular') and positions is None:
                    self.layout_mode=self.layout_choice.value='network'
            self.component_choice.disabled = self.mode == 'mcc'
            self.selection_version += 1
            self.selected = self.hovered = None
            self.tip.visible = False
            self.details.controls = [ft.Text('Selecione um composto para ver suas informações.')]
            self.status.value = 'Passe o mouse sobre um nó. Clique para abrir o composto.'
            self.redraw()
            if self.explorer:
                self.explorer.inspector.clear()
                await self.explorer.inspector.follow_mode()
            await self.fit_view()
            if self.explorer and self.explorer.molecular.initialized:
                self.explorer.molecular.scope.value=self.mode
                await self.explorer.molecular.analyze()
            self.ui.page.update()
        await self.ui.guard(action)

    async def center_on(self, identifier):
        if identifier in self.points:
            await self.viewer.reset()
            x, y = self.points[identifier]
            await self.viewer.pan(self.viewport_width/2-x, self.viewport_height/2-y)

    def close(self):
        self.active=False;self.selection_version+=1
        self.drawing.shapes=[];self.base_shapes=[];self.points={}
        self.raster.visible=False;self.raster.src=''
        self.details.controls=[];self.fragment_panel.controls=[]
        if self.explorer:
            self.explorer.request+=1;self.explorer.close()

    async def follow_relation(self, identifier):
        if self.explorer and identifier not in self.displayed():await self.explorer.reset_filters()
        if self.explorer and identifier not in self.displayed():
            self.search.value=identifier
            return await self.find(None)
        await self.select(identifier)
        await self.center_on(identifier)

    async def change_component(self, event):
        async def action():
            if not self.active or self.ui.token != self.token:
                return
            await self.ui.call(self.ui.store.project, self.token, self.project_id)
            self.component_id = None if event.control.value == 'all' else int(event.control.value)
            self.selection_version += 1
            self.selected = self.hovered = None
            self.tip.visible = False
            self.details.controls = [ft.Text('Selecione um composto para ver suas informações.')]
            self.redraw()
            if self.explorer:
                self.explorer.inspector.clear()
                await self.explorer.inspector.follow_mode()
            await self.fit_view()
            self.ui.page.update()
        await self.ui.guard(action)

    async def change_style(self,event):
        if self.explorer:return await self.explorer.change(event)
        async def action():
            if not self.active or self.ui.token!=self.token:return
            await self.ui.call(self.ui.store.project,self.token,self.project_id)
            self.layout_mode=self.layout_choice.value;self.color_mode=self.color_choice.value;self.labels=bool(self.label_choice.value)
            self.hovered=None;self.tip.visible=False;self.redraw()
            await self.viewer.reset();self.ui.page.update()
        await self.ui.guard(action)

    def build(self):
        from .flow_canvas import icon_button
        from .dataset_origin import origin_header
        heading=origin_header(self.data)
        async def zoom_in(e): await self.viewer.zoom(1.25)
        async def zoom_out(e): await self.viewer.zoom(.8)
        async def reset(e): await self.fit_view()
        async def fit(e): await self.fit_view()
        tools = [icon_button(icon=ft.Icons.ZOOM_IN, tooltip='Ampliar', on_click=zoom_in),
                 icon_button(icon=ft.Icons.ZOOM_OUT, tooltip='Reduzir', on_click=zoom_out),
                 ft.TextButton('Recentrar', on_click=reset), ft.TextButton('Ajustar à tela', on_click=fit)]
        if self.data['kind'] == 'graph':
            self.mode_choice = ft.Dropdown(label='Visualização', value=self.mode, width=300, on_select=self.change_mode,
                         options=[ft.DropdownOption(key='full', text='Grafo completo'),
                                  ft.DropdownOption(key='mcc', text='Máximo componente conectado (MCC)')])
            tools.insert(0, self.mode_choice)
            self.component_choice = ft.Dropdown(label='Navegar entre componentes', value='all', expand=True,
                disabled=self.mode=='mcc',
                on_select=self.change_component, options=[ft.DropdownOption(key='all', text='Todos os componentes')]+[
                    ft.DropdownOption(key=str(index), text=f'Componente {index} · {len(group)} compostos · {self.graph.subgraph(group).number_of_edges()} relações')
                    for index, group in self.components.items()])
            self.color_choice=ft.Dropdown(label='Cor dos nós',value='degree',expand=True,on_select=self.change_style,
                options=[ft.DropdownOption(key='degree',text='Grau · número de conexões'),ft.DropdownOption(key='density',text='Conectividade normalizada')])
            self.layout_choice=ft.Dropdown(label='Organização do grafo',value='network',expand=True,on_select=self.change_style,
                options=[ft.DropdownOption(key='network',text='Rede · componentes separados'),ft.DropdownOption(key='circular',text='Anéis · comparar os nós')])
            self.label_choice=ft.Switch(label='Exibir códigos no grafo',value=self.labels,on_change=self.change_style)
            style=[ft.Row([self.component_choice]),ft.ResponsiveRow([ft.Container(ft.Row([self.color_choice]),col={'xs':12,'md':6}),
                ft.Container(ft.Row([self.layout_choice]),col={'xs':12,'md':6})],spacing=20,run_spacing=16),
                ft.Row([self.label_choice],wrap=True),
                ft.Container(height=12,border_radius=6,gradient=ft.LinearGradient(colors=VIRIDIS,begin=ft.Alignment.CENTER_LEFT,end=ft.Alignment.CENTER_RIGHT)),self.scale_label]
            body=self.explorer.build(tools,style)
            return ft.Column([*heading,body],spacing=16) if heading else body
        else:
            style=[ft.Row([ft.Text('● BBB+ · laranja',color=COLORS['BBB+'],size=12),ft.Text('■ BBB− · azul',color=COLORS['BBB-'],size=12),
                ft.Text('Branco: região HIA · amarelo: região BBB',size=12,color=INK)],wrap=True,spacing=18)]
        description = 'Os componentes ficam separados. Escolha um componente para explorar seus vértices. Clique em um composto para destacar seus vizinhos e navegar pelas relações. Arraste o fundo e use o zoom.' if self.data['kind'] == 'graph' else 'EGG local: aproximação heurística de permeação passiva; não é uma predição SwissADME validada. A cor indica BBB, não P-gp. Pontos fora das regiões continuam visíveis. Clique para inspecionar propriedades e motivos da seleção.'
        report=self.data.get('admet_report')
        if report:
            from biomolexplorer.admet_config import PROFILES
            config=report.get('configuration',{});profile=PROFILES.get(config.get('admet_profile'),{})
            style.insert(0,ft.Text(profile.get('label','ADMET')+' · '+str(report.get('selected',0))+' / '+str(report.get('evaluated',0))+' compostos atendem à hipótese na avaliação original.',size=12,color=INK))
        left = ft.Container(col={'xs':12, 'xl':9} if self.data['kind']=='egg' else {'xs':12, 'md':8}, content=ft.Column([
            ft.Row(tools, wrap=True), self.summary,*style,
            ft.Semantics(label='Grafo molecular interativo' if self.data['kind']=='graph' else 'Gráfico EGG interativo',
                content=ft.Container(self.viewer, height=560 if self.data['kind']=='egg' else self.VIEWPORT_HEIGHT, bgcolor='#EEF2F7' if self.data['kind']=='egg' else '#F8FAFD',
                         border=ft.Border.all(1, '#DCE5F4'), border_radius=12, clip_behavior=ft.ClipBehavior.HARD_EDGE)),
            self.status, ft.Text(description, size=11, color='#64748B')], spacing=12))
        right = ft.Container(col={'xs':12, 'xl':3} if self.data['kind']=='egg' else {'xs':12, 'md':4}, padding=18, height=None if self.popup else 570, bgcolor='#FFFFFF',
                            border=ft.Border.all(1, '#DCE5F4'), border_radius=12,
                            content=ft.Column([ft.Row([self.search, icon_button(icon=ft.Icons.SEARCH, tooltip='Buscar composto', on_click=self.find)]),
                                               *([*self.actions,ft.Text('Passe o mouse para identificar um composto. Clique no ponto para abrir sua estrutura 2D.',size=13,color='#64748B'),self.fragment_panel] if self.popup else [self.fragment_panel,ft.Container(self.details, expand=True)])], spacing=18))
        body=ft.ResponsiveRow([left, right], spacing=20, run_spacing=20)
        return ft.Column([*heading,body],spacing=16) if heading else body
