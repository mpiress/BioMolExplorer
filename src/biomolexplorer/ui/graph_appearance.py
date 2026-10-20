"""One appearance mapping shared by canvas, legends and Gephi exports."""
from collections import Counter
import flet as ft
import flet.canvas as canvas
from biomolexplorer.network_analysis import METRICS,PALETTE

PALETTES={
    'viridis':['#440154','#3B528B','#21918C','#5EC962','#FDE725'],
    'green_yellow':['#166534','#4D7C0F','#A3A635','#FDE725'],
    'blue_red':['#2166AC','#67A9CF','#F7F7F7','#EF8A62','#B2182B'],
}
MISSING='#94A3B8'


def continuous_color(value,maximum,palette='viridis'):
    colors=PALETTES[palette]
    fraction=min(1,max(0,value/max(1e-12,maximum)))
    position=fraction*(len(colors)-1);left=min(len(colors)-2,int(position));weight=position-left
    a,b=colors[left],colors[left+1]
    return '#'+''.join(f'{round(int(a[i:i+2],16)*(1-weight)+int(b[i:i+2],16)*weight):02X}' for i in (1,3,5))


class GraphAppearance:
    def __init__(self,explorer):
        self.explorer=explorer;self.lines=[];self.swatches=[];self.category_page=0
        self.attribute=ft.Dropdown(label='Cor dos nós',value='metric',expand=True,on_select=explorer.change,
            options=[ft.DropdownOption(key=k,text=v) for k,v in [('metric','Ranking · métrica selecionada'),('degree','Ranking · grau'),('component','Partição · componente'),('source','Partição · base de origem')]])
        self.palette=ft.Dropdown(label='Paleta',value='viridis',expand=True,on_select=explorer.change,
            options=[ft.DropdownOption(key=k,text=v) for k,v in [('viridis','Viridis · roxo → amarelo'),('green_yellow','Verde → amarelo'),('blue_red','Azul → vermelho')]])
        self.size=ft.Dropdown(label='Tamanho dos nós',value='metric',expand=True,on_select=explorer.change,
            options=[ft.DropdownOption(key='metric',text='Ranking · métrica selecionada'),ft.DropdownOption(key='uniform',text='Uniforme')])
        self.edges=ft.Switch(label='Exibir arestas',value=True,on_change=explorer.change)
        self.components=ft.Switch(label='Contornos dos componentes',value=True,on_change=explorer.change)
        self.legend=ft.Column(spacing=5,horizontal_alignment=ft.CrossAxisAlignment.STRETCH)

    def map(self,nodes,values,metric,aggregate=False):
        view=self.explorer.view;analysis=self.explorer.analysis
        self.colors={};self.lines=[];self.swatches=[];self.maximum=max(values.values(),default=0)
        size_metric=metric;size_maximum=self.maximum
        radii={node:4+6*(max(0,values.get(node,0))/max(self.maximum,1e-12))**.5
               if self.size.value=='metric' else 5 for node in nodes}
        attribute=self.attribute.value
        categories=None
        if aggregate or self.explorer.color_communities.value:
            mapping=analysis.community_cache.get(view.mode,({},None,''))[0]
            categories={node:f'C{node if aggregate else mapping[node]}' for node in nodes if aggregate or node in mapping}
            all_categories=[f'C{c}' for c in sorted(set(mapping.values()))]
            self.title='Comunidade topológica · cores categóricas'
        elif attribute=='component':
            mapping={node:i for i,group in view.components.items() for node in group}
            categories={node:f'Componente {mapping[node]}' for node in nodes}
            all_categories=[f'Componente {i}' for i in view.components]
            self.title='Componente conexo · cores categóricas'
        elif attribute=='source':
            mapping={node:str(record['properties'].get('source') or 'Origem não informada') for node,record in view.nodes.items()}
            categories={node:mapping[node] for node in nodes}
            all_categories=sorted(set(mapping.values()))
            self.title='Base de origem · cores categóricas'
        if categories is not None:
            colors={label:PALETTE[i%len(PALETTE)] for i,label in enumerate(all_categories)}
            self.colors={node:colors.get(categories.get(node),MISSING) for node in nodes}
            counts=Counter(categories.values())
            if aggregate:counts=Counter({f'C{node}':int(values[node]) for node in nodes})
            self.swatches=[(label,colors[label],counts[label]) for label in all_categories if counts[label]]
            self.category_page=min(self.category_page,max(0,(len(self.swatches)-1)//30))
            self.lines=[f'{label}: {count} nós' for label,_,count in self.swatches]
            if len(all_categories)>len(PALETTE):self.lines.append('Cores se repetem; use o identificador do grupo.')
            self.legend.controls=[ft.Text(self.title,size=12),ft.Row([
                ft.Container(padding=4,content=ft.Row([ft.Container(width=12,height=12,bgcolor=color),ft.Text(f'{label} · {count}',size=11)],tight=True))
                for label,color,count in self.swatches[self.category_page*30:(self.category_page+1)*30]],wrap=True)]
            if len(self.swatches)>30:
                self.legend.controls.append(ft.Row([ft.TextButton('Legenda anterior',on_click=lambda e:self.page(-1)),
                    ft.Text(f'{self.category_page+1} / {(len(self.swatches)+29)//30}',size=11),ft.TextButton('Mais cores',on_click=lambda e:self.page(1))],wrap=True))
            if len(categories)<len(nodes):self.legend.controls.append(ft.Text('Cinza: comunidade indisponível.',size=11))
        else:
            if attribute=='degree':values=analysis.metrics(view.mode,'degree')[0];metric='degree'
            self.maximum=max(values.values(),default=0)
            self.colors={node:continuous_color(values[node],self.maximum,self.palette.value) if node in values else MISSING for node in nodes}
            self.title=METRICS[metric].split(' · ')[0]+' · '+view.mode+' · escala do universo antes dos filtros'
            self.has_values=bool(values)
            palette=PALETTES[self.palette.value] if self.maximum else [continuous_color(0,0,self.palette.value)]*2
            if not values:palette=[MISSING]*2
            self.lines=[self.title,f'0 → {self.maximum:.5g}', 'Cinza: métrica indisponível. Contorno azul: pertence ao MCC.']
            self.legend.controls=[ft.Text(self.title,size=12),ft.Container(height=12,border_radius=4,
                gradient=ft.LinearGradient(colors=palette,begin=ft.Alignment.CENTER_LEFT,end=ft.Alignment.CENTER_RIGHT)),
                ft.Row([ft.Text('0 · menor valor',size=11),ft.Text(f'{self.maximum/2:.5g}',size=11),ft.Text(f'{self.maximum:.5g} · maior valor',size=11)],alignment=ft.MainAxisAlignment.SPACE_BETWEEN)]
            if not values:self.legend.controls.append(ft.Text('Cinza: métrica indisponível.',size=11))
        self.size_label=('Número de compostos por comunidade' if aggregate else 'Tamanho: '+('uniforme' if self.size.value=='uniform' else f'{size_metric} · 0 → {size_maximum:.5g}'))
        self.legend.controls.append(ft.Text(self.size_label+' · Contorno azul: MCC.',size=11))
        return self.colors,radii

    def page(self,direction):
        self.category_page=max(0,min(max(0,(len(self.swatches)-1)//30),self.category_page+direction))
        self.explorer.view.redraw();self.explorer.view.ui.page.update()

    def export_shapes(self,y,width):
        shapes=[canvas.Text(20,y,self.title,style=ft.TextStyle(size=12,color='#172B4D'))]
        if self.swatches:
            for i,(label,color,count) in enumerate(self.swatches[:90]):
                yy=y+22+i*18
                shapes.extend([canvas.Rect(20,yy,12,12,paint=ft.Paint(color=color)),canvas.Text(40,yy,f'{label} · {count}',style=ft.TextStyle(size=11,color='#172B4D'))])
            end=y+50+min(90,len(self.swatches))*18
            if len(self.swatches)>90:
                shapes.append(canvas.Text(20,end-20,'Legenda: primeiras 90 categorias; IDs completos no relatório e no GEXF.',style=ft.TextStyle(size=11)))
            return shapes,end
        for i in range(100):
            color=(continuous_color(i if self.maximum else 0,99,self.palette.value) if self.has_values else MISSING)
            shapes.append(canvas.Rect(20+i*(width-40)/100,y+24,(width-40)/100+1,12,paint=ft.Paint(color=color)))
        shapes.extend([canvas.Text(20,y+42,'0 · menor valor',style=ft.TextStyle(size=11)),canvas.Text(width-170,y+42,f'{self.maximum:.5g} · maior valor',style=ft.TextStyle(size=11))])
        shapes.append(canvas.Text(20,y+64,'Cinza: indisponível · Contorno azul: MCC · '+self.size_label,style=ft.TextStyle(size=11)))
        return shapes,y+90
