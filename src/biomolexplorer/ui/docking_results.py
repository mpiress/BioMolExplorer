"""Paginated scores and pose previews shared by Vina, DOCK6 and consensus."""
import flet as ft
from biomolexplorer.docking_results import DockingResults
from .compound_table import CompoundTableViewer
from .localization import verbatim


class DockingResultsTable(CompoundTableViewer):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.service=DockingResults(self.ui.store)
        self.sort_by=None;self.ascending=True

    def page_options(self):
        return {'sort_by':self.sort_by,'ascending':self.ascending} if self.sort_by else {}

    async def change_table(self,e):
        self.sort_by=None;self.table.sort_column_index=None
        await super().change_table(e)

    async def sort_rows(self,key,index):
        if not self.valid():return
        self.ascending=not self.ascending if self.sort_by==key else True
        self.sort_by=key;self.offset=0
        self.table.sort_column_index=index;self.table.sort_ascending=self.ascending
        await self.ui.guard(self.load)

    async def load(self):
        await super().load()
        if not self.valid() or self.data is None:return
        properties=self.data.get('columns',[]) or (self.data['rows'][0]['properties'] if self.data['rows'] else {})
        consensus='vina' in properties and 'dock6' in properties
        scores=['vina','dock6','normalized_score'] if consensus else [k for k in ('score','vina','dock6','z-score','min-max') if k in properties]
        columns=[('molecule_chembl_id','Código do composto')]+([] if consensus else [('receptor_id','Alvo')])+[
            (key,{'vina':'Score Vina','dock6':'Score DOCK6','normalized_score':'Score normalizado'}.get(key,key)) for key in scores]+([] if consensus else [('canonical_smiles','SMILES')])
        self.table.columns=[]
        for index,(key,label) in enumerate(columns):
            async def sort(e,key=key,index=index):await self.sort_rows(key,index)
            self.table.columns.append(ft.DataColumn(ft.Text(label),numeric=key in scores,on_sort=sort))
        self.table.columns.append(ft.DataColumn(ft.Text('Ações')))
        rows=[]
        for row in self.data['rows']:
            fields=row['properties']
            async def remove(e,row=row):await self.ui.guard(lambda:self.confirm_remove(row))
            actions=[]
            for column,label in (('vina_pose','3D Vina'),('dock6_pose','3D DOCK6')) if 'vina_pose' in fields else (('conformer_file','3D'),):
                async def preview(e,row=row,column=column):await self.ui.guard(lambda:self.pose_preview(row,column))
                actions.append(ft.TextButton(label,on_click=preview,disabled=not bool(fields.get(column))))
            async def contacts(e,row=row):await self.ui.guard(lambda:self.contacts(row))
            actions.append(ft.TextButton('Resíduos',icon=ft.Icons.SCIENCE_OUTLINED,on_click=contacts,
                disabled=not bool(fields.get('vina_pose') or fields.get('conformer_file'))))
            if fields.get('engine')=='dock6' or fields.get('dock6_pose') or fields.get('footprint_file'):
                async def footprint(e,row=row):await self.ui.guard(lambda:self.footprint(row))
                actions.append(ft.TextButton('Footprint',icon=ft.Icons.BAR_CHART,on_click=footprint,disabled=not row.get('footprint_available',False)))
            actions.append(ft.TextButton('Remover',icon=ft.Icons.DELETE_OUTLINE,on_click=remove,disabled=not self.data['can_edit']))
            values=[row['id']]+([] if consensus else [fields.get('receptor_id','—')])+[fields.get(k,'—') for k in scores]+([] if consensus else [row['smiles'] or '—'])
            rows.append(ft.DataRow(cells=[ft.DataCell(verbatim(ft.Text(str(v),selectable=True))) for v in values]+
                                     [ft.DataCell(ft.Row(actions,spacing=2))]))
        self.table.rows=rows;self.ui.page.update()

    async def pose_preview(self,row,column):
        if not self.valid():return
        sequence=self.sequence
        selection=dict(table=self.selector.value,index=row['index'],version=self.data['version'],column=column)
        if self.valid() and sequence==self.sequence:
            await self.ui.preview_docking(self.project_id,self.run_id,self.stage_id,'docking',selection,self.token)

    async def contacts(self,row):
        from .docking_scene import residue_dialog
        sequence=self.sequence
        selection=dict(table=self.selector.value,index=row['index'],version=self.data['version'],
            column='vina_pose' if row['properties'].get('vina_pose') else 'conformer_file')
        await residue_dialog(self.ui,self.project_id,
            lambda:self.service.scene(self.token,self.project_id,self.run_id,self.stage_id,**selection),
            lambda:self.valid() and sequence==self.sequence,lambda:self.pose_preview(row,selection['column']))

    async def footprint(self,row):
        if not self.valid():return
        from .footprint import footprint_dialog
        self.preview_version+=1;version=self.preview_version
        selection=dict(table=self.selector.value,index=row['index'],version=self.data['version'])
        await footprint_dialog(self.ui,self.project_id,self.run_id,self.stage_id,selection,row['id'],
            lambda:self.valid() and version==self.preview_version,
            origin=row['properties'].get('footprint_origin'),receptor=row['properties'].get('receptor_id'))
