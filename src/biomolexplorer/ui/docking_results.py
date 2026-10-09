"""Paginated scores and pose previews shared by Vina, DOCK6 and consensus."""
import flet as ft
from biomolexplorer.docking_results import DockingResults
from .compound_table import CompoundTableViewer
from .localization import verbatim


class DockingResultsTable(CompoundTableViewer):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.service=DockingResults(self.ui.store)

    async def load(self):
        await super().load()
        if not self.valid() or self.data is None:return
        properties=self.data['rows'][0]['properties'] if self.data['rows'] else {}
        scores=[k for k in ('score','vina','dock6','z-score','min-max') if k in properties]
        self.table.columns=[ft.DataColumn(ft.Text(n)) for n in ('Código do composto','Alvo',*scores,'SMILES','Ações')]
        rows=[]
        for row in self.data['rows']:
            fields=row['properties']
            async def remove(e,row=row):await self.ui.guard(lambda:self.confirm_remove(row))
            actions=[]
            for column,label in (('vina_pose','3D Vina'),('dock6_pose','3D DOCK6')) if 'vina_pose' in fields else (('conformer_file','3D'),):
                async def preview(e,row=row,column=column):await self.ui.guard(lambda:self.pose_preview(row,column))
                actions.append(ft.TextButton(label,on_click=preview,disabled=not bool(fields.get(column))))
            actions.append(ft.TextButton('Remover',icon=ft.Icons.DELETE_OUTLINE,on_click=remove,disabled=not self.data['can_edit']))
            values=[row['id'],fields.get('receptor_id','—'),*[fields.get(k,'—') for k in scores],row['smiles'] or '—']
            rows.append(ft.DataRow(cells=[ft.DataCell(verbatim(ft.Text(str(v),selectable=True))) for v in values]+
                                     [ft.DataCell(ft.Row(actions,spacing=2))]))
        self.table.rows=rows;self.ui.page.update()

    async def pose_preview(self,row,column):
        if not self.valid():return
        sequence=self.sequence
        pose=await self.ui.call(self.service.pose,self.token,self.project_id,self.run_id,self.stage_id,
            self.selector.value,row['index'],self.data['version'],column)
        if self.valid() and sequence==self.sequence:await self.ui.preview_artifact(self.project_id,pose)
