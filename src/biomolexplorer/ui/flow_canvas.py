"""Native, zoomable Flet graph editor; connections are validated before mutation."""
import copy
from uuid import uuid4
import flet as ft
from .localization import verbatim, stage_control
import flet.canvas as canvas
from biomolexplorer import flow
from biomolexplorer.catalog import TITLES, new_stage
from biomolexplorer.bindings import sources,input_labels
from .zoom import zoomable_view

def icon_button(**kwargs):
    label=kwargs.get('tooltip','Ação')
    return ft.Semantics(container=True,label=label,button=True,on_tap=None if kwargs.get('disabled') else kwargs.get('on_click'),disabled=kwargs.get('disabled',False),exclude_semantics=True,content=ft.IconButton(**kwargs))

COLORS={'Entradas':'#2563EB','Recuperação':'#0D9488','Análise':'#7C3AED','Docking':'#D97706'}

class FlowCanvas:
    def __init__(self,ui,writable):
        self.ui,self.writable=ui,writable
        self.history=[]; self.future=[]; self.pending=None
        self.stack=ft.Stack(width=flow.CANVAS_WIDTH,height=flow.CANVAS_HEIGHT)
        self.viewer=zoomable_view(content=ft.DragTarget(group='stage-library',content=self.stack,on_accept=self.drop),
            constrained=False,alignment=ft.Alignment.TOP_LEFT,expand=True)
        self.selection=ft.Dropdown(label='Bloco selecionado',width=260,options=[],on_select=self.select_from_menu)
        self.status=ft.Text('Arraste um bloco. Conecte a saída à entrada. Duplo clique para configurar.',size=12,color='#64748B')
        self.refresh(False)

    @property
    def stages(self): return self.ui.current['pipeline']
    def remember(self):
        self.history.append(copy.deepcopy(self.stages)); self.history=self.history[-50:]; self.future.clear()
    def changed(self):
        self.ui.dirty=True; self.refresh()
        if hasattr(self.ui,'schedule_save'):self.ui.schedule_save()
    def refresh(self,update=True,nodes=True):
        height=max(1600,max((flow.position(s,i)['y']+flow.node_height(s)+120 for i,s in enumerate(self.stages)),default=0))
        self.stack.height=height
        if not hasattr(self,'grid') or self.grid.height!=height:
            grid=[]
            for x in range(0,flow.CANVAS_WIDTH,40): grid.append(canvas.Line(x,0,x,height,paint=ft.Paint(color='#E8EDF5',stroke_width=1)))
            for y in range(0,int(height),40): grid.append(canvas.Line(0,y,flow.CANVAS_WIDTH,y,paint=ft.Paint(color='#E8EDF5',stroke_width=1)))
            self.grid=canvas.Canvas(shapes=grid,width=flow.CANVAS_WIDTH,height=height)
        shapes=[]
        by_id={s['id']:s for s in self.stages}
        for index,target in enumerate(self.stages):
            pos=flow.position(target,index)
            for j,port in enumerate(flow.input_ports(target)):
              group=target.get('bindings',{}).get(port['field'])
              for reference in sources(group) if group else []:
                source=by_id.get(reference.get('stage'))
                if not source: continue
                start=flow.position(source,self.stages.index(source)); x=start['x']+flow.NODE_WIDTH-30; y=start['y']+83
                endx=pos['x']+28; endy=pos['y']+117+j*34
                shapes.append(canvas.Path([canvas.Path.MoveTo(x,y),canvas.Path.CubicTo(x+100,y,endx-100,endy,endx,endy)],
                    paint=ft.Paint(color=COLORS[TITLES[source['operation']][2]],stroke_width=3,style=ft.PaintingStyle.STROKE)))
        self.stack.controls=[self.grid,canvas.Canvas(shapes=shapes,width=flow.CANVAS_WIDTH,height=height)]+([self.node(s,i) for i,s in enumerate(self.stages)] if nodes else self.stack.controls[2:])
        self.selection.options=[stage_control(ft.DropdownOption(key=s['id'],text=s['name']),s,'text') for s in self.stages]
        self.selection.value=self.ui.selected if self.ui.selected in by_id else None
        if update: self.ui.page.update()
    def node(self,stage,index):
        color=COLORS[TITLES[stage['operation']][2]]; pos=flow.position(stage,index)
        async def select(e): self.ui.selected=stage['id']; self.refresh()
        async def configure(e): await self.ui.guard(lambda:self.ui.open_stage_dialog(stage['id']))
        async def results(e): await self.ui.guard(lambda:self.ui.open_stage_results(stage['id']))
        async def information(e): self.show_information(stage)
        def begin(e):
            if self.writable: self.remember()
        def move(e):
            if not self.writable or e.local_delta is None: return
            stage['position']=flow.bounded_position(pos['x']+e.local_delta.x,pos['y']+e.local_delta.y,stage)
            pos.update(stage['position']); self.ui.dirty=True
            node.left=pos['x']; node.top=pos['y']; self.refresh(nodes=False)
        async def output(e):
            if self.writable:
                self.pending=stage['id']; self.status.value='Escolha uma entrada compatível para '+stage['name']; self.ui.page.update()
        out=ft.Draggable(group='flow-link',data=stage['id'],max_simultaneous_drags=1 if self.writable else 0,
            content=icon_button(icon=ft.Icons.RADIO_BUTTON_CHECKED,icon_color=color,tooltip='Conectar saída de '+stage['name'],on_click=output),
            content_feedback=ft.Icon(ft.Icons.RADIO_BUTTON_CHECKED,color=color,size=28))
        out=ft.Semantics(container=True,exclude_semantics=True,label='Conectar saída de '+stage['name'],button=True,on_tap=output,content=out)
        def moved(e):
            if self.writable and hasattr(self.ui,'schedule_save'):self.ui.schedule_save()
        header=ft.GestureDetector(on_tap=select,on_double_tap=configure,on_pan_start=begin,on_pan_update=move,on_pan_end=moved,
            mouse_cursor=ft.MouseCursor.MOVE,content=ft.Container(width=flow.NODE_WIDTH-4,height=64,padding=12,bgcolor=color,border_radius=ft.BorderRadius.only(top_left=14,top_right=14),
            content=ft.Column([ft.Text(TITLES[stage['operation']][2].upper(),size=10,color='#FFFFFF'),stage_control(ft.Text(stage['name'],size=14,color='#FFFFFF',weight=ft.FontWeight.W_600,max_lines=1,overflow=ft.TextOverflow.ELLIPSIS),stage)],spacing=4)))
        rows=[ft.Row([ft.Text('Saída',size=12,color=color,expand=True),out],spacing=0,height=34)]
        for port in flow.input_ports(stage):
            field=port['field']
            async def accept(e,field=field): await self.link(e.src.data,stage['id'],field)
            async def click(e,field=field):
                if self.pending: await self.link(self.pending,stage['id'],field)
                else: await configure(e)
            linked=field in stage.get('bindings',{})
            dot=ft.DragTarget(group='flow-link',on_accept=accept,content=icon_button(icon=ft.Icons.LINK if linked else ft.Icons.ADD_CIRCLE_OUTLINE,
                icon_color=color,tooltip='Entrada · '+port['label']+' · '+stage['name'],disabled=not self.writable,on_click=click))
            dot=ft.Semantics(container=True,exclude_semantics=True,label='Entrada · '+port['label']+' · '+stage['name'],button=True,on_tap=click,content=dot)
            async def unlink(e,field=field): self.remember(); flow.disconnect(stage,field); self.changed()
            rows.append(ft.Row([dot,ft.Text(port['label'],size=11,expand=True),icon_button(icon=ft.Icons.LINK_OFF,icon_size=16,tooltip='Desconectar '+port['label'],on_click=unlink,disabled=not self.writable) if linked else ft.Container(width=30)],height=34,spacing=0))
        issues=flow.stage_issues(stage)
        labels=input_labels(stage,getattr(self.ui,'project_assets',[]),getattr(self.ui,'automatic_labels',{}))
        rows.append(ft.Container(height=34,padding=ft.Padding.only(left=8,right=8),content=ft.Text(
            'Entradas: '+('; '.join(labels) if labels else 'Nenhuma / coleta própria'),size=10,max_lines=2,
            overflow=ft.TextOverflow.ELLIPSIS,tooltip='\n'.join(labels))))
        rows.append(ft.Row([ft.Icon(ft.Icons.ERROR_OUTLINE if issues else ft.Icons.CHECK_CIRCLE_OUTLINE,size=14,color='#D97706' if issues else color),
            ft.Text(issues[0] if issues else 'Resultados fornecidos · concluído' if stage.get('provided_results') else 'Pronto para executar',size=10,expand=True),
            icon_button(icon=ft.Icons.INFO_OUTLINE,icon_size=18,tooltip='Informações do bloco',on_click=information),
            icon_button(icon=ft.Icons.INSIGHTS,icon_size=18,tooltip='Resultados de '+stage['name'],on_click=results),
            icon_button(icon=ft.Icons.TUNE,icon_size=18,tooltip='Configurar '+stage['name'],on_click=configure)],height=32,spacing=4))
        node=ft.Container(left=pos['x'],top=pos['y'],width=flow.NODE_WIDTH,
            border=ft.Border.all(2,color if self.ui.selected==stage['id'] else '#DDE5F0'),border_radius=16,bgcolor='#FFFFFF',
            content=ft.Column([header,ft.Container(ft.Column(rows,spacing=0),padding=ft.Padding.symmetric(horizontal=8))],spacing=0))
        return node
    async def select_from_menu(self,e):
        self.ui.selected=self.selection.value; self.refresh()
    async def configure_selected(self,e):
        if self.ui.selected: await self.ui.guard(lambda:self.ui.open_stage_dialog(self.ui.selected))
        else: self.ui.notify('Selecione um bloco para configurar.')
    async def keyboard(self,e):
        if not self.writable or self.ui.current is None or self.ui.tab!='Pipeline' or getattr(self.ui,'editing_stage',False) or getattr(self,'typing',False): return
        if e.key=='Escape': self.pending=None; self.status.value='Conexão cancelada.'; self.ui.page.update()
        elif e.ctrl and e.key.lower()=='z': await (self.redo(e) if e.shift else self.undo(e))
        elif e.ctrl and e.key.lower()=='y': await self.redo(e)
        elif e.key=='Delete': await self.delete(e)
    async def link(self,source,target,field):
        if not self.writable: return
        before=copy.deepcopy(self.stages)
        try: flow.connect(self.stages,source,target,field)
        except ValueError as exc:
            self.pending=None
            self.status.value=str(exc)
            self.ui.notify(str(exc)); return
        self.history.append(before); self.history=self.history[-50:]; self.future.clear(); self.pending=None
        self.status.value='Conexão criada. A execução respeitará as dependências.'; self.changed()
    def show_information(self,stage):
        from .block_help import block_information
        from .feedback import close_dialog
        description,inputs,outputs,note=block_information(stage)
        body=[ft.Text(description),ft.Text('O que recebe',weight=ft.FontWeight.W_600),
              *[ft.Text(value) for value in inputs],ft.Text('O que fornece',weight=ft.FontWeight.W_600),
              *[ft.Text(value) for value in outputs]]
        if note:body.append(ft.Text(note))
        dialog=ft.AlertDialog(title=ft.Text(TITLES[stage['operation']][0]),
            content=ft.Container(width=520,content=ft.Column(body,tight=True,spacing=12,scroll=ft.ScrollMode.AUTO)),
            actions=[ft.TextButton('Fechar',on_click=lambda e:close_dialog(self.ui.page,dialog))])
        self.ui.page.show_dialog(dialog)
    async def drop(self,e):
        if self.writable: await self.add(e.src.data,e.local_position.x,e.local_position.y)
    async def add(self,operation,x=None,y=None):
        if not self.writable: return
        if len(self.stages)>=100:
            self.ui.notify('Um pipeline pode conter até 100 blocos.'); return
        self.remember(); stage=new_stage(operation)
        if hasattr(self.ui,'tr'):stage['name']=self.ui.tr(stage['name'])
        if operation=='docking_dock6' and self.ui.service.dock6_path: stage['parameters']['dock6_app_path']=str(self.ui.service.dock6_path)
        stage['position']=flow.bounded_position(x,y,stage) if x is not None and y is not None else flow.available_position(self.stages,stage)
        self.stages.append(stage); self.ui.selected=stage['id']; self.changed()
    async def undo(self,e):
        if self.writable and self.history:
            self.future.append(copy.deepcopy(self.stages)); self.ui.current['pipeline']=self.history.pop(); self.changed()
    async def redo(self,e):
        if self.writable and self.future:
            self.history.append(copy.deepcopy(self.stages)); self.ui.current['pipeline']=self.future.pop(); self.changed()
    async def delete(self,e):
        if self.writable and self.ui.selected: self.remember(); flow.remove(self.stages,self.ui.selected); self.ui.selected=None; self.changed()
    async def duplicate(self,e):
        stage=self.ui.selected_stage()
        if self.writable and stage:
            if len(self.stages)>=100:
                self.ui.notify('Um pipeline pode conter até 100 blocos.'); return
            self.remember(); clone=copy.deepcopy(stage); clone['id']=uuid4().hex; clone['name']+=' (cópia)'; pos=flow.position(stage); clone['position']=flow.bounded_position(pos['x']+40,pos['y']+60,clone)
            self.stages.append(clone); self.ui.selected=clone['id']; self.changed()
    async def arrange(self,e):
        if self.writable: self.remember(); flow.arrange(self.stages); self.changed()
    async def zoom_in(self,e): await self.viewer.zoom(1.2)
    async def zoom_out(self,e): await self.viewer.zoom(.8)
    async def reset(self,e): await self.viewer.reset()
    def build(self):
        palette=[]; items=[]; sections=[]
        for category in COLORS:
            blocks=[]
            for operation,(title,help_,group) in TITLES.items():
                if operation=='expand_similar_compounds':continue # Retain legacy pipelines, offer the independent PubChem block.
                if group!=category: continue
                async def add(e,op=operation): await self.add(op)
                async def information(e,op=operation):self.show_information(new_stage(op))
                blocks.append(ft.Container(padding=8,bgcolor='#FFFFFF',border_radius=10,border=ft.Border.all(1,'#E2E8F0'),tooltip=help_,content=ft.Row([ft.Draggable(group='stage-library',data=operation,max_simultaneous_drags=1 if self.writable else 0,
                    content_feedback=ft.Container(ft.Text(title,color='#FFFFFF'),bgcolor=COLORS[category],padding=16,border_radius=12),
                    content=ft.Container(width=126,padding=4,content=ft.Row([ft.Icon(ft.Icons.DRAG_INDICATOR,size=16,color=COLORS[category]),ft.Text(title,size=12,expand=True)],spacing=4))),
                    icon_button(icon=ft.Icons.INFO_OUTLINE,icon_size=18,tooltip='Informações do bloco',on_click=information),
                    icon_button(icon=ft.Icons.ADD,icon_size=18,tooltip='Adicionar '+title,on_click=add,disabled=not self.writable)],spacing=0)))
                items.append((title,blocks[-1]))
            section=ft.ExpansionTile(title=ft.Text(category,size=13,color=COLORS[category],weight=ft.FontWeight.W_600),
                controls=blocks,expanded=False,maintain_state=True,dense=True,min_tile_height=44,
                shape=ft.RoundedRectangleBorder(radius=10,side=ft.BorderSide(0,ft.Colors.TRANSPARENT)),
                collapsed_shape=ft.RoundedRectangleBorder(radius=10,side=ft.BorderSide(0,ft.Colors.TRANSPARENT)),
                tile_padding=ft.Padding.symmetric(horizontal=4,vertical=4),controls_padding=ft.Padding.only(bottom=8),
                expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
                icon_color=COLORS[category],collapsed_icon_color=COLORS[category])
            sections.append(section);palette.append(section)
        search=ft.TextField(label='Buscar etapa',text_size=12)
        expanded_before_search=None
        def filter_palette(e):
            nonlocal expanded_before_search
            query=(search.value or '').strip().casefold()
            if query and expanded_before_search is None:
                expanded_before_search=[section.expanded for section in sections]
            for title,item in items: item.visible=query in (self.ui.tr(title) if hasattr(self.ui,'tr') else title).casefold()
            for index,section in enumerate(sections):
                section.visible=any(item.visible for item in section.controls)
                if query:section.expanded=section.visible
                elif expanded_before_search is not None:section.expanded=expanded_before_search[index]
            if not query:expanded_before_search=None
            self.ui.page.update()
        search.on_change=filter_palette
        search.on_focus=lambda e:setattr(self,'typing',True)
        search.on_blur=lambda e:setattr(self,'typing',False)

        bar=ft.Row([icon_button(icon=icon,tooltip=label,on_click=handler,disabled=disabled) for icon,label,handler,disabled in [
            (ft.Icons.UNDO,'Desfazer',self.undo,not self.writable),(ft.Icons.REDO,'Refazer',self.redo,not self.writable),
            (ft.Icons.CONTENT_COPY,'Duplicar bloco',self.duplicate,not self.writable),(ft.Icons.DELETE_OUTLINE,'Excluir bloco',self.delete,not self.writable),
            (ft.Icons.AUTO_FIX_HIGH,'Organizar blocos',self.arrange,not self.writable),(ft.Icons.ZOOM_IN,'Ampliar',self.zoom_in,False),
            (ft.Icons.ZOOM_OUT,'Reduzir',self.zoom_out,False),(ft.Icons.CENTER_FOCUS_STRONG,'Redefinir visão',self.reset,False)]],spacing=0)
        return ft.Column([ft.Row([ft.Button('Salvar',on_click=self.ui.event(self.ui.save_draft),disabled=not self.writable),
            ft.Button('Modelos prontos',on_click=self.ui.event(self.ui.preset_dialog),disabled=not self.writable),
            ft.Button('Executar pipeline',icon=ft.Icons.PLAY_ARROW,on_click=self.ui.event(self.ui.run_pipeline),disabled=not self.writable),
            ft.Button('Executar até seleção',on_click=self.ui.event(self.ui.run_pipeline,True),disabled=not self.writable),bar],wrap=True),
            ft.Row([self.selection,ft.Button('Configurar seleção',icon=ft.Icons.TUNE,on_click=self.configure_selected)],wrap=True),
            self.status,ft.Row([ft.Container(width=235,padding=14,bgcolor='#F8FAFC',border_radius=14,
                content=ft.Column([ft.Text('Biblioteca de etapas',weight=ft.FontWeight.W_600),search,ft.Text('Arraste ou use + para adicionar.',size=11,color='#64748B'),*palette],scroll=ft.ScrollMode.AUTO)),
                ft.Container(self.viewer,expand=True,bgcolor='#F5F7FC',border_radius=14,border=ft.Border.all(1,'#DDE5F0'))],expand=True,spacing=14,vertical_alignment=ft.CrossAxisAlignment.STRETCH)],expand=True,spacing=8)
