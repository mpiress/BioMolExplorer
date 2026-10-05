"""Native Flet graph/EGG exploration without external pages or JavaScript."""
import json
import base64

import flet as ft
import flet.canvas as canvas

from biomolexplorer.diagnostics import get_logger
from biomolexplorer.visualizations import PointIndex, molecule_image, separated_graph_layout
from .feedback import close_dialog
from .zoom import zoomable_view
from .localization import verbatim

BLUE = '#2552E8'
INK = '#172B4D'
VIRIDIS = ['#440154','#3B528B','#21918C','#5EC962','#FDE725']
GRAPH_NODE_RADIUS = 4


def degree_color(value,maximum):
    fraction=min(1,max(0,value/max(1e-12,maximum)))
    position=fraction*(len(VIRIDIS)-1);left=min(len(VIRIDIS)-2,int(position));weight=position-left
    a,b=VIRIDIS[left],VIRIDIS[left+1]
    return '#'+''.join(f'{round(int(a[i:i+2],16)*(1-weight)+int(b[i:i+2],16)*weight):02X}' for i in (1,3,5))


class ResultsViewer:
    WIDTH, HEIGHT = 720, 460
    VIEWPORT_HEIGHT = 460

    def __init__(self, ui, project_id, data, popup=False, actions=None):
        self.ui, self.project_id, self.data = ui, project_id, data
        self.token=ui.token;self.active=True;self.popup=popup;self.actions=actions or []
        self.nodes = {node['id']: node for node in data['nodes']}
        self.component_id = None
        self.components = {}
        if data['kind'] == 'graph':
            import networkx as nx
            self.graph = nx.Graph()
            self.graph.add_nodes_from(self.nodes)
            self.graph.add_edges_from((edge['source'], edge['target'], {'value': edge.get('value', 0)})
                                     for edge in data['edges'])
            groups = sorted(nx.connected_components(self.graph), key=lambda group: (-len(group), min(group)))
            self.components = {index+1: group for index, group in enumerate(groups)}
            # Old saved results benefit from the new component layout as well.
            self.network_positions = ({key: (node['x'], node['y']) for key, node in self.nodes.items()}
                                      if data.get('layout') == 'components' else separated_graph_layout(self.graph))
            self.circular_positions = None
        self.mode = 'full'
        self.layout_mode='network';self.color_mode='degree';self.labels=len(self.nodes)<=20
        self.selected = None
        self.selection_version = 0
        self.hovered = None
        self.status = ft.Text('Passe o mouse sobre um nó. Clique para abrir o composto.', size=12, color='#64748B')
        self.details = ft.Column([ft.Text('Informações do composto', size=19, weight=ft.FontWeight.W_600),
                                  ft.Text('Selecione um nó ou busque seu código.', size=13, color='#64748B')], spacing=14, scroll=ft.ScrollMode.AUTO)
        self.tip_text = ft.Text('', color='#FFFFFF', size=12)
        self.tip = ft.Container(self.tip_text, visible=False, bgcolor=INK, padding=10, border_radius=8,
                                ignore_interactions=True, left=0, top=0)
        self.drawing = canvas.Canvas(width=self.WIDTH, height=self.HEIGHT)
        self.scene = ft.Stack([ft.GestureDetector(content=self.drawing, hover_interval=60,
                             on_hover=self.hover, on_exit=self.leave, on_tap_up=self.click), self.tip],
                             width=self.WIDTH, height=self.HEIGHT)
        self.viewer = zoomable_view(content=self.scene, constrained=False,
                                          alignment=ft.Alignment.TOP_LEFT, expand=True)
        self.search = ft.TextField(label='Código do composto', hint_text='CHEMBL… ou seu identificador',
                                   expand=True, on_submit=self.find)
        self.summary = ft.Text('', size=13, color='#64748B')
        self.scale_label=ft.Text(size=12,color='#64748B')
        self.fragment_panel=ft.Column(spacing=14,visible=False)
        self.populate_fragment()
        self.redraw()

    def populate_fragment(self):
        if self.data['kind']!='graph':return
        fragment=self.data.get('fragment',{})
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
        if self.mode == 'mcc':
            allowed = set(self.data['mcc'])
            return {key: value for key, value in self.nodes.items() if key in allowed}
        if self.component_id is not None:
            allowed = self.components.get(self.component_id, set())
            return {key: value for key, value in self.nodes.items() if key in allowed}
        return self.nodes

    def coordinates(self, nodes):
        if self.data['kind'] == 'egg':
            xmin = min(0, min((node['x'] for node in nodes.values()), default=0))
            xmax = max(200, max((node['x'] for node in nodes.values()), default=200))
            ymin = min(-2, min((node['y'] for node in nodes.values()), default=-2))
            ymax = max(7, max((node['y'] for node in nodes.values()), default=7))
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
            nodes={key:dict(node,x=positions[key][0],y=positions[key][1]) for key,node in nodes.items()}
            xmin = min((node['x'] for node in nodes.values()), default=-1)
            xmax = max((node['x'] for node in nodes.values()), default=1)
            ymin = min((node['y'] for node in nodes.values()), default=-1)
            ymax = max((node['y'] for node in nodes.values()), default=1)
            dx, dy = xmax-xmin, ymax-ymin
            # Preserve aspect ratio and a minimum 30 px distance between nodes.
            # Dense networks grow beyond the viewport and remain pannable.
            scale = max(30, min((720-110)/max(1, dx), (460-90)/max(1, dy)))
            self.WIDTH, self.HEIGHT = max(720, dx*scale+110), max(460, dy*scale+90)
            self.scene.width = self.drawing.width = self.WIDTH
            self.scene.height = self.drawing.height = self.HEIGHT
            self.project = lambda x, y: ((self.WIDTH-dx*scale)/2+(x-xmin)*scale,
                                         (self.HEIGHT-dy*scale)/2+(y-ymin)*scale)
        return {key: self.project(node['x'], node['y']) for key, node in nodes.items()}

    def egg_background(self):
        shapes = []
        xmin, xmax, ymin, ymax = self.bounds
        for cx, cy, rx, ry, color in [(75, 2, 75, 3, '#FFFFFF'), (42, 2.3, 47, 2.2, '#FFE066')]:
            x1, y1 = self.project(cx-rx, cy+ry)
            x2, y2 = self.project(cx+rx, cy-ry)
            shapes.append(canvas.Oval(x1, y1, x2-x1, y2-y1, paint=ft.Paint(color=color)))
        for value in range(0, int(xmax)+1, max(50, int((xmax-xmin)/5))):
            x, y = self.project(value, ymin)
            shapes.append(canvas.Text(x, y+12, str(value), style=ft.TextStyle(size=11, color=INK)))
        for i in range(6):
            value = ymin+(ymax-ymin)*i/5
            x, y = self.project(xmin, value)
            shapes.append(canvas.Text(12, y-6, f'{value:.1f}', style=ft.TextStyle(size=11, color=INK)))
        shapes += [canvas.Line(65, 35, 65, self.HEIGHT-55, paint=ft.Paint(color='#64748B', stroke_width=1)),
                   canvas.Line(65, self.HEIGHT-55, self.WIDTH-40, self.HEIGHT-55, paint=ft.Paint(color='#64748B', stroke_width=1)),
                   canvas.Text(self.WIDTH/2-30, self.HEIGHT-20, 'TPSA (Å²)', style=ft.TextStyle(size=13, color=INK)),
                   canvas.Text(10, 8, 'WLOGP', style=ft.TextStyle(size=13, color=INK))]
        return shapes

    def redraw(self):
        nodes = self.displayed()
        self.points = self.coordinates(nodes)
        self.index = PointIndex(self.points)
        shapes = self.egg_background() if self.data['kind'] == 'egg' else []
        if self.data['kind'] == 'graph' and self.mode == 'full' and self.component_id is None:
            for component, identifiers in self.components.items():
                points = [self.points[identifier] for identifier in identifiers]
                x, y = min(point[0] for point in points)-18, min(point[1] for point in points)-26
                width = max(point[0] for point in points)-x+18
                height = max(point[1] for point in points)-y+18
                shapes.append(canvas.Rect(x, y, width, height, border_radius=10,
                    paint=ft.Paint(color='#E2E8F0', stroke_width=1, style=ft.PaintingStyle.STROKE)))
                shapes.append(canvas.Text(x+8, y+5, f'Componente {component}',
                    style=ft.TextStyle(size=10, color='#64748B')))
        edges = [edge for edge in self.data['edges'] if edge['source'] in nodes and edge['target'] in nodes]
        for edge in edges:
            a, b = self.points[edge['source']], self.points[edge['target']]
            shapes.append(canvas.Line(*a, *b, paint=ft.Paint(color='#CBD5E1', stroke_width=.6+float(edge.get('value',0))*.7)))
        degrees={key:int(node['properties'].get('degree',0)) for key,node in nodes.items()}
        values={key:degree if self.color_mode=='degree' else degree/max(1,len(nodes)-1) for key,degree in degrees.items()}
        maximum=max(values.values(),default=0)
        mcc=set(self.data['mcc'])
        self.radii={}
        for identifier, point in self.points.items():
            props = nodes[identifier]['properties']
            color = ('#DC2626' if props.get('BBB') == 'BBB+' else BLUE) if self.data['kind'] == 'egg' else degree_color(values[identifier],maximum)
            radius=5 if self.data['kind']=='egg' else GRAPH_NODE_RADIUS
            self.radii[identifier]=radius
            if self.data['kind']=='graph' and self.mode=='full' and identifier in mcc:
                shapes.append(canvas.Circle(*point,radius=radius+1.5,paint=ft.Paint(color=BLUE,stroke_width=1,style=ft.PaintingStyle.STROKE)))
            shapes.append(canvas.Circle(*point, radius=radius, paint=ft.Paint(color=color)))
            if self.data['kind']=='graph' and self.labels:
                shapes.append(verbatim(canvas.Text(point[0]+radius+4,point[1]-8,identifier,style=ft.TextStyle(size=10,color=INK))))
        self.base_shapes = shapes
        self.highlight()
        self.summary.value = (f'{len(nodes)} compostos · {len(edges)} relações · {len(self.components)} componentes · MCC: {len(mcc)} de {len(self.nodes)} compostos' if self.data['kind'] == 'graph' else f'{len(nodes)} compostos · vermelho: BBB+ · azul: BBB−')
        self.scale_label.value=('Grau: 0 → '+str(int(maximum))+' conexões' if self.color_mode=='degree' else f'Conectividade: 0 → {maximum:.3f} · grau / (n − 1)')
        self.fragment_panel.visible=self.mode=='mcc' and self.data['kind']=='graph'
        if not nodes:
            self.status.value = 'Nenhum composto disponível nesta visualização.'

    def highlight(self):
        shapes = list(self.base_shapes)
        if self.selected in self.points:
            if self.data['kind'] == 'graph':
                for neighbor in self.graph.neighbors(self.selected):
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
        if identifier == self.hovered:
            return
        self.hovered = identifier
        self.tip.visible = identifier is not None
        if identifier:
            overlapping = self.index.nearby(*self.points[identifier])
            self.tip_text.value = identifier + (f' (+{len(overlapping)-1})' if len(overlapping)>1 else '')
            self.tip.left = min(position.x+12, self.WIDTH-190)
            self.tip.top = max(0, position.y-42)
            self.status.value = 'Composto: '+identifier+' · clique para ver a estrutura e as propriedades.'
            if self.data['kind']=='graph':self.status.value+=f" · grau: {self.nodes[identifier]['properties'].get('degree',0)}"
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
                await self.ui.guard(lambda: self.pick(identifier))

    async def pick(self, identifier):
        await self.ui.call(self.ui.store.project, self.ui.token, self.project_id)
        overlapping = self.index.nearby(*self.points[identifier])
        if len(overlapping) == 1:
            return await self.select(identifier)
        choices = []
        for compound_id in overlapping:
            async def select_compound(event, compound_id=compound_id):
                close_dialog(self.ui.page,dialog)
                await self.ui.guard(lambda:self.select(compound_id))
            choices.append(verbatim(ft.TextButton(compound_id, on_click=select_compound)))
        dialog=ft.AlertDialog(title=ft.Text('Compostos neste ponto'),
            content=ft.Container(width=360,height=min(360,len(choices)*50),
                                 content=ft.Column(choices,scroll=ft.ScrollMode.AUTO)),
            actions=[ft.TextButton('Fechar',on_click=lambda e:close_dialog(self.ui.page,dialog))])
        self.ui.page.show_dialog(dialog)

    async def find(self, event):
        async def action():
            identifier = (self.search.value or '').strip()
            matches = [key for key in self.nodes if key.casefold() == identifier.casefold()]
            if not matches:
                raise ValueError('Código não encontrado nesta visualização. No grafo, confira se você está no MCC ou no grafo completo.')
            chosen = matches[0]
            if chosen not in self.displayed():
                self.mode = 'full'
                self.mode_choice.value = 'full'
                self.component_choice.disabled = False
                self.component_id = next(index for index, group in self.components.items() if chosen in group)
                self.component_choice.value = str(self.component_id)
                self.redraw()
                await self.viewer.reset()
            await self.select(chosen)
            await self.center_on(chosen)
        await self.ui.guard(action)

    async def select(self, identifier):
        # Membership may have changed since opening this dialog.
        self.selection_version += 1
        version, token = self.selection_version, self.ui.token
        await self.ui.call(self.ui.store.project, token, self.project_id)
        node = self.displayed().get(identifier)
        if node is None:
            raise ValueError('O composto não pertence à visualização atual.')
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
        for key, value in node['properties'].items():
            label = str(key)
            rendered = '—' if value is None else json.dumps(value, ensure_ascii=False) if isinstance(value, (dict, list)) else str(value)
            controls.append(ft.Column([ft.Text(label, size=11, color='#64748B'),
                                       verbatim(ft.Text(rendered, size=13, selectable=True))], spacing=4))
        if self.data['kind'] == 'graph':
            neighbors = sorted(self.graph[identifier].items(), key=lambda item: (-item[1]['value'], item[0]))
            controls.append(ft.Text(f'Relações deste composto ({len(neighbors)})', size=15, weight=ft.FontWeight.W_600))
            for neighbor, attributes in neighbors:
                async def follow(event, target=neighbor):
                    if self.popup:
                        close_dialog(self.ui.page, dialog)
                    await self.ui.guard(lambda: self.follow_relation(target))
                controls.append(ft.TextButton(f"{neighbor} · similaridade {attributes['value']:.3f}", on_click=follow))
        self.details.controls = controls
        self.selected = identifier
        self.highlight()
        self.ui.page.update()
        if self.popup:
            dialog=ft.AlertDialog(title=ft.Text(identifier+' · molécula 2D',size=22),
                content=ft.Container(width=460,height=520,content=ft.Column(controls,spacing=14,scroll=ft.ScrollMode.AUTO)),
                actions=[ft.TextButton('Fechar molécula',on_click=lambda e:close_dialog(self.ui.page,dialog))])
            self.ui.page.show_dialog(dialog)

    async def change_mode(self, event):
        async def action():
            if not self.active or self.ui.token!=self.token:return
            await self.ui.call(self.ui.store.project, self.ui.token, self.project_id)
            self.mode = event.control.value
            self.component_choice.disabled = self.mode == 'mcc'
            self.selection_version += 1
            self.selected = self.hovered = None
            self.tip.visible = False
            self.details.controls = [ft.Text('Selecione um composto para ver suas informações.')]
            self.status.value = 'Passe o mouse sobre um nó. Clique para abrir o composto.'
            self.redraw()
            await self.viewer.reset()
            self.ui.page.update()
        await self.ui.guard(action)

    async def center_on(self, identifier):
        if identifier in self.points:
            await self.viewer.reset()
            x, y = self.points[identifier]
            await self.viewer.pan(360-x, self.VIEWPORT_HEIGHT/2-y)

    async def follow_relation(self, identifier):
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
            await self.viewer.reset()
            self.ui.page.update()
        await self.ui.guard(action)

    async def change_style(self,event):
        async def action():
            if not self.active or self.ui.token!=self.token:return
            await self.ui.call(self.ui.store.project,self.token,self.project_id)
            self.layout_mode=self.layout_choice.value;self.color_mode=self.color_choice.value;self.labels=bool(self.label_choice.value)
            self.hovered=None;self.tip.visible=False;self.redraw()
            await self.viewer.reset();self.ui.page.update()
        await self.ui.guard(action)

    def build(self):
        from .flow_canvas import icon_button
        async def zoom_in(e): await self.viewer.zoom(1.25)
        async def zoom_out(e): await self.viewer.zoom(.8)
        async def reset(e): await self.viewer.reset()
        async def fit(e):
            await self.viewer.reset()
            await self.viewer.zoom(min(720/self.WIDTH, self.VIEWPORT_HEIGHT/self.HEIGHT))
        tools = [icon_button(icon=ft.Icons.ZOOM_IN, tooltip='Ampliar', on_click=zoom_in),
                 icon_button(icon=ft.Icons.ZOOM_OUT, tooltip='Reduzir', on_click=zoom_out),
                 ft.TextButton('Recentrar', on_click=reset), ft.TextButton('Ajustar à tela', on_click=fit)]
        if self.data['kind'] == 'graph':
            self.mode_choice = ft.Dropdown(label='Visualização', value='full', width=300, on_select=self.change_mode,
                         options=[ft.DropdownOption(key='full', text='Grafo completo'),
                                  ft.DropdownOption(key='mcc', text='Máximo componente conectado (MCC)')])
            tools.insert(0, self.mode_choice)
            self.component_choice = ft.Dropdown(label='Navegar entre componentes', value='all', expand=True,
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
        else:style=[]
        description = 'Os componentes ficam separados. Escolha um componente para explorar seus vértices. Clique em um composto para destacar seus vizinhos e navegar pelas relações. Arraste o fundo e use o zoom.' if self.data['kind'] == 'graph' else 'Regiões branca e amarela: HIA e BBB. Classificações heurísticas da plataforma.'
        left = ft.Container(col={'xs':12, 'md':8}, content=ft.Column([
            ft.Row(tools, wrap=True), self.summary,*style,
            ft.Semantics(label='Grafo molecular interativo' if self.data['kind']=='graph' else 'Gráfico EGG interativo',
                content=ft.Container(self.viewer, height=self.VIEWPORT_HEIGHT, bgcolor='#EEF2F7' if self.data['kind']=='egg' else '#F8FAFD',
                         border=ft.Border.all(1, '#DCE5F4'), border_radius=12, clip_behavior=ft.ClipBehavior.HARD_EDGE)),
            self.status, ft.Text(description, size=11, color='#64748B')], spacing=12))
        right = ft.Container(col={'xs':12, 'md':4}, padding=18, height=None if self.popup else 570, bgcolor='#FFFFFF',
                            border=ft.Border.all(1, '#DCE5F4'), border_radius=12,
                            content=ft.Column([ft.Row([self.search, icon_button(icon=ft.Icons.SEARCH, tooltip='Buscar composto', on_click=self.find)]),
                                               *([*self.actions,ft.Text('Passe o mouse para identificar um composto. Clique no ponto para abrir sua estrutura 2D.',size=13,color='#64748B'),self.fragment_panel] if self.popup else [self.fragment_panel,ft.Container(self.details, expand=True)])], spacing=18))
        return ft.ResponsiveRow([left, right], spacing=20, run_spacing=20)
