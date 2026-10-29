"""All CSV schemas share the score-table pagination and atomic rendering path."""
import flet as ft
from biomolexplorer.result_tables import ResultTables
from biomolexplorer.table_sort import is_numeric
from .compound_table import CompoundTableViewer
from .localization import verbatim
from .feedback import close_dialog


class ResultTableViewer(CompoundTableViewer):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs);self.service=ResultTables(self.ui.store)

    async def change_table(self,e):
        self.sorting.key=None
        await super().change_table(e)

    def build(self):return super().build(inspect_dataset=False)

    def render_rows(self,result):
        if not result['columns']:
            self.table.columns=[ft.DataColumn(ft.Text('—'))];self.table.rows=[];self.table.sort_column_index=None
            return
        fields=result['columns'];numeric={key for key in fields if key not in ('molecule_chembl_id','canonical_smiles','PDB_CODE') and is_numeric(r['properties'].get(key) for r in result['rows'])}
        molecular='canonical_smiles' in fields
        self.table.columns=[ft.DataColumn(verbatim(ft.Text(key))) for key in fields]+([ft.DataColumn(ft.Text('Visualização'))] if molecular else [])
        self.sorting.bind(self.table,fields+([None] if molecular else []),self.sort_changed,numeric)
        rows=[]
        for row in result['rows']:
            cells=[ft.DataCell(verbatim(ft.Text(str(row['properties'].get(key) or '—'),max_lines=2,
                    overflow=ft.TextOverflow.ELLIPSIS,tooltip=str(row['properties'].get(key) or '')[:2000],selectable=True))) for key in fields]
            if molecular:
                async def preview(e,row=row):await self.ui.guard(lambda:self.preview(row,False))
                cells.append(ft.DataCell(ft.TextButton('2D',on_click=preview,disabled=not bool(row['smiles']))))
            rows.append(ft.DataRow(cells=cells))
        self.table.rows=rows


async def show_result_table(ui,pid,rid,sid,path):
    token=ui.token
    tables=await ui.call(ResultTables(ui.store).tables,token,pid,rid,sid)
    if ui.token!=token or not ui.current or ui.current['id']!=pid:return
    selected=next((row for row in tables if row['path']==path),None)
    if selected is None:return
    viewer=ResultTableViewer(ui,pid,rid,sid,[selected],inline=True)
    def close(e=None):viewer.active=False;close_dialog(ui.page,dialog)
    dialog=ft.AlertDialog(title=verbatim(ft.Text(selected['name'])),on_dismiss=lambda e:setattr(viewer,'active',False),
        content=ft.Container(width=max(240,min(1150,(ui.page.width or 1300)-100)),
            height=max(250,min(650,(ui.page.height or 900)-230)),content=ft.Column([viewer.build()],scroll=ft.ScrollMode.AUTO)),
        actions=[ft.TextButton('Fechar',on_click=close)])
    ui.page.show_dialog(dialog)
    await viewer.load()
