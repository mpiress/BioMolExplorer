"""Paginated preparation exclusions beside docking results."""
import flet as ft
from .localization import verbatim
from .feedback import close_dialog


def show_preparation_failures(ui,records):
    token=ui.token;project_id=ui.current['id'];offset=0;limit=25
    body=ft.Column(spacing=12,scroll=ft.ScrollMode.AUTO,height=390)
    count=ft.Text(size=12,color='#64748B')
    previous=ft.TextButton('Página anterior');next_button=ft.TextButton('Próxima página')
    def render():
        body.controls=[]
        for row in records[offset:offset+limit]:
            body.controls.append(ft.Container(padding=12,border_radius=10,bgcolor='#F8FAFC',content=ft.Column([
                verbatim(ft.Text(row['compound_identifier'],size=15,weight=ft.FontWeight.W_600)),
                verbatim(ft.Text(row['reason'],size=12,selectable=True)),
                verbatim(ft.Text(row['preparation_step']+' · '+row['preparation_context'],size=11,color='#64748B',selectable=True)),
                verbatim(ft.Text(row['submitted_smiles'],size=11,selectable=True))],spacing=6)))
        count.value=f'{len(records)} compostos não utilizados · Página {offset//limit+1}'
        previous.disabled=offset==0;next_button.disabled=offset+limit>=len(records)
    def move(direction):
        def click(e):
            nonlocal offset
            if ui.token!=token or not ui.current or ui.current['id']!=project_id:return
            offset=max(0,offset+direction*limit);render();ui.page.update()
        return click
    previous.on_click=move(-1);next_button.on_click=move(1);render()
    dialog=ft.AlertDialog(title=ft.Text('Compostos não utilizados'),
        content=ft.Container(width=min(820,max(280,(ui.page.width or 1000)-100)),content=ft.Column([
            ft.Text('Estes compostos foram ignorados. Os demais seguiram para a preparação e o docking.',size=13),
            body,ft.Row([count,previous,next_button],wrap=True)],tight=True,spacing=12)),
        actions=[ft.TextButton('Fechar',on_click=lambda e:close_dialog(ui.page,dialog))])
    ui.page.show_dialog(dialog)
