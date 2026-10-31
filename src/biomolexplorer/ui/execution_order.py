"""Explicit, persistent scheduling order without editing scientific state."""
import flet as ft
from biomolexplorer.execution_order import move
from .feedback import close_dialog
from .localization import default_stage_names, verbatim


async def show_order(ui,continue_selected=None,selected_stage_id=None):
    token,pid=ui.token,ui.current['id']
    preference=await ui.call(ui.store.execution_order,token,pid)
    project=await ui.call(ui.store.project,token,pid)
    if ui.token!=token or not ui.current or ui.current['id']!=pid:return
    stages=project['pipeline'];by_id={s['id']:s for s in stages};order=preference['order'][:]
    rows=ft.Column(spacing=8,scroll=ft.ScrollMode.AUTO)
    def shift(identifier,direction):
        order[:]=move(stages,order,identifier,direction)
        redraw();ui.page.update()
    def redraw():
        rows.controls=[]
        for index,identifier in enumerate(order):
            stage=by_id[identifier]
            name=ui.tr(stage['name']) if stage['name'] in default_stage_names() else stage['name']
            rows.controls.append(ft.Container(padding=10,bgcolor='#F1F5F9',border_radius=10,
                content=ft.Row([ft.Text(str(index+1),width=32),ft.Container(verbatim(ft.Text(name)),expand=True),
                    ft.IconButton(ft.Icons.ARROW_UPWARD,tooltip=ui.tr('Mover acima'),disabled=move(stages,order,identifier,-1)==order,on_click=lambda e,key=identifier:shift(key,-1)),
                    ft.IconButton(ft.Icons.ARROW_DOWNWARD,tooltip=ui.tr('Mover abaixo'),disabled=move(stages,order,identifier,1)==order,on_click=lambda e,key=identifier:shift(key,1))])))
    async def save(e):
        if ui.token!=token or not ui.current or ui.current['id']!=pid:return
        await ui.call(ui.store.save_execution_order,token,pid,order[:])
        close_dialog(ui.page,dialog)
        if ui.token!=token or not ui.current or ui.current['id']!=pid:return
        if continue_selected is not None:
            ui.submitting=True
            try:await ui._submit_pipeline(continue_selected,selected_stage_id=selected_stage_id)
            finally:ui.submitting=False
    dialog=ft.AlertDialog(modal=True,title=ft.Text('Ordem de execução'),
        content=ft.Container(width=min(640,max(240,(ui.page.width or 1100)-80)),height=min(480,max(160,(ui.page.height or 900)-250)),
            content=ft.Column([ft.Text('Mover um bloco leva seus dependentes. As dependências sempre executam primeiro. A ordem se aplica às próximas execuções.'),ft.Container(rows,expand=True)],spacing=16)),
        actions=[ft.TextButton('Cancelar',on_click=lambda e:close_dialog(ui.page,dialog)),
            ft.Button('Salvar ordem',on_click=ui.event_action(save,configuration_dialog=lambda:dialog,processing_title='Salvando ordem de execução…'))])
    redraw();ui.page.show_dialog(dialog)
