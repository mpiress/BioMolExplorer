"""Paginated preparation exclusions beside docking results."""
import flet as ft
from .localization import verbatim
from .feedback import close_dialog
from .table_sort import TableSort


def show_preparation_failures(ui,records):
    token=ui.token;project_id=ui.current['id'];offset=0;limit=25
    body=ft.Column(spacing=12,scroll=ft.ScrollMode.AUTO,height=390)
    count=ft.Text(size=12,color='#64748B')
    previous=ft.TextButton('Página anterior');next_button=ft.TextButton('Próxima página')
    sorting=TableSort()
    def render():
        fields=list(dict.fromkeys(['compound_identifier','reason','preparation_step','preparation_context','submitted_smiles']+[key for row in records for key in row]))
        table=ft.DataTable(columns=[ft.DataColumn(verbatim(ft.Text(key))) for key in fields],rows=[],column_spacing=24)
        sorting.bind(table,fields,sort_changed,numeric={'input_row'})
        for row in sorting.order(records,numeric={'input_row'})[offset:offset+limit]:
            table.rows.append(ft.DataRow(cells=[ft.DataCell(verbatim(ft.Text(str(row.get(key) or '—'),
                max_lines=3,tooltip=str(row.get(key) or ''),selectable=True))) for key in fields]))
        body.controls=[ft.Row([table],scroll=ft.ScrollMode.AUTO)]
        count.value=f'{len(records)} compostos não utilizados · Página {offset//limit+1}'
        previous.disabled=offset==0;next_button.disabled=offset+limit>=len(records)
    def sort_changed():
        nonlocal offset
        if ui.token!=token or not ui.current or ui.current['id']!=project_id:return
        offset=0;render();ui.page.update()
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
