"""Footprint chart popup with pan, zoom and the original PDF download."""
import inspect
import io
import flet as ft

from biomolexplorer.docking_results import DockingResults
from biomolexplorer.pdf_preview import pdf_preview
from .feedback import close_dialog
from .localization import verbatim
from .zoom import zoomable_view


async def footprint_dialog(ui,pid,rid,sid,selection,name,valid,origin=None,receptor=None):
    if not valid():return
    token=ui.token
    service=DockingResults(ui.store)
    async def read():
        return await ui.call(service.footprint_pdf,token,pid,rid,sid,**selection)
    data=await read()
    if not valid():return
    active=True
    image_ratio=1.1
    previous_resize=getattr(ui.page,'on_resize',None)
    chart=ft.Container(bgcolor='#FFFFFF',border_radius=12,alignment=ft.Alignment.CENTER,
        border=ft.Border.all(1,'#E2E8F0'),clip_behavior=ft.ClipBehavior.HARD_EDGE,
        content=ft.Column([ft.ProgressRing(width=32,height=32),ft.Text('Carregando gráfico…',color='#64748B')],
            alignment=ft.MainAxisAlignment.CENTER,horizontal_alignment=ft.CrossAxisAlignment.CENTER))
    viewer=zoomable_view(content=ft.Image(src=b'',fit=ft.BoxFit.CONTAIN,
        semantics_label='Gráfico footprint de '+name),expand=True)
    async def zoom_in(e):
        if active and valid():await viewer.zoom(1.25)
    async def zoom_out(e):
        if active and valid():await viewer.zoom(.8)
    async def reset(e):
        if active and valid():await viewer.reset()
    tools=ft.Row([ft.IconButton(ft.Icons.ZOOM_OUT,tooltip='Reduzir',on_click=zoom_out),
        ft.IconButton(ft.Icons.ZOOM_IN,tooltip='Ampliar',on_click=zoom_in),
        ft.TextButton('Recentrar',icon=ft.Icons.CENTER_FOCUS_STRONG,on_click=reset)],spacing=2,visible=False)
    caption='Use a roda do mouse para ampliar e arraste para explorar o gráfico.'
    content=ft.Container(content=ft.Column([
        ft.Row([ft.Container(content=ft.Text('DOCK6',size=12,weight=ft.FontWeight.W_600,color='#0F766E'),
            bgcolor='#CCFBF1',border_radius=8,padding=8),
            verbatim(ft.Text(name,size=18,weight=ft.FontWeight.W_600,selectable=True))],wrap=True),
        ft.Row([ft.Text('Energias de interação por resíduo',size=13,color='#64748B'),
            *([verbatim(ft.Text('· '+receptor,size=13,color='#64748B'))] if receptor else [])],wrap=True,spacing=4),
        tools,chart,ft.Text(caption,size=12,color='#64748B'),
        *([ft.Text('Footprint calculado sobre a entrada minimizada, antes do docking final.',size=12,color='#B45309')]
          if origin!='docked_pose' else [])],expand=True,spacing=6))
    def resize():
        width=getattr(ui.page,'width',None) or 1200
        height=getattr(ui.page,'height',None) or 900
        available_width=max(160,width-80)
        available_height=max(160,height-170)
        # Size the popup around the chart's aspect ratio instead of reserving
        # a wide landscape box for the nearly square footprint page.
        overhead=180 if available_width<540 else 150
        if origin!='docked_pose':overhead+=28
        chart_height=max(100,available_height-overhead)
        chart_width=min(available_width,chart_height*image_ratio)
        content.width=chart.width=max(160,chart_width)
        chart.height=min(chart_height,chart.width/image_ratio)
        content.height=chart.height+overhead
        viewer.content.width=chart.width
        viewer.content.height=chart.height
    async def resized(e):
        if previous_resize:
            result=previous_resize(e)
            if inspect.isawaitable(result):await result
        if active:
            resize();ui.page.update()
    def dismissed(e=None):
        nonlocal active
        active=False
        if getattr(ui.page,'on_resize',None) is resized:ui.page.on_resize=previous_resize
    def close(e):
        dismissed();close_dialog(ui.page,dialog)
    async def save():
        if not active or not valid():return
        pdf=await read()
        if active and valid():await ui.picker.save_file(file_name=name+'_footprint.pdf',src_bytes=pdf)
    async def download(e):await ui.guard(save)
    dialog=ft.AlertDialog(title=ft.Row([ft.Icon(ft.Icons.BAR_CHART,color='#0F766E'),ft.Text('Footprint',size=22)],spacing=10),
        content=content,inset_padding=16,content_padding=ft.Padding.symmetric(horizontal=24,vertical=12),
        clip_behavior=ft.ClipBehavior.HARD_EDGE,on_dismiss=dismissed,
        actions=[ft.Button('Baixar PDF',icon=ft.Icons.DOWNLOAD,on_click=download),ft.TextButton('Fechar',on_click=close)])
    resize();ui.page.on_resize=resized;ui.page.show_dialog(dialog)
    try:
        image=await ui.call(pdf_preview,data)
    except ValueError:
        if active and valid():
            chart.content=ft.Column([ft.Icon(ft.Icons.PICTURE_AS_PDF,size=48,color='#64748B'),
                ft.Text('A prévia não está disponível. Use “Baixar PDF” para consultar o gráfico.')],
                alignment=ft.MainAxisAlignment.CENTER,horizontal_alignment=ft.CrossAxisAlignment.CENTER)
    else:
        if active and valid():
            from PIL import Image
            with Image.open(io.BytesIO(image)) as preview:image_ratio=preview.width/preview.height
            viewer.content.src=image;chart.content=viewer;tools.visible=True
            resize()
    if active and valid():ui.page.update()
    elif active:close(None)
