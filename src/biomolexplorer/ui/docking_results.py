"""Paginated scores and pose previews shared by Vina, DOCK6 and consensus."""
import flet as ft
from biomolexplorer.docking_results import DockingResults
from .compound_table import CompoundTableViewer
from .localization import verbatim


class DockingResultsTable(CompoundTableViewer):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.service=DockingResults(self.ui.store)

    async def change_table(self,e):
        self.sorting.key=None;self.table.sort_column_index=None
        await super().change_table(e)

    def render_rows(self,result):
        properties=self.data.get('columns',[]) or (self.data['rows'][0]['properties'] if self.data['rows'] else {})
        consensus='vina' in properties and 'dock6' in properties
        scores=['vina','dock6','normalized_score'] if consensus else [k for k in ('score','vina','dock6','z-score','min-max') if k in properties]
        columns=[('molecule_chembl_id','Código do composto')]+([] if consensus else [('receptor_id','Alvo')])+[
            (key,{'vina':'Score Vina','dock6':'Score DOCK6','normalized_score':'Score normalizado'}.get(key,key)) for key in scores]+([] if consensus else [('canonical_smiles','SMILES')])
        self.table.columns=[]
        for index,(key,label) in enumerate(columns):
            self.table.columns.append(ft.DataColumn(ft.Text(label),numeric=key in scores))
        self.table.columns.append(ft.DataColumn(ft.Text('Ações')))
        self.sorting.bind(self.table,[key for key,label in columns]+[None],self.sort_changed,numeric=scores)
        rows=[]
        for row in self.data['rows']:
            fields=row['properties']
            async def remove(e,row=row):await self.ui.guard(lambda:self.confirm_remove(row))
            actions=[]
            for column,label in (('vina_pose','3D Vina'),('dock6_pose','3D DOCK6')) if 'vina_pose' in fields else (('conformer_file','3D'),):
                async def preview(e,row=row,column=column):await self.ui.guard(lambda:self.pose_preview(row,column))
                url=self.view_url(row,column) if fields.get(column) else None
                actions.append(ft.TextButton(label,on_click=None if url else preview,disabled=not bool(fields.get(column)),
                    url=ft.Url(url,target=ft.UrlTarget.BLANK) if url else None))
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
        self.table.rows=rows

    def view_url(self,row,column):
        if not getattr(self.ui.page,'web',False):return None
        selection=dict(table=self.selector.value,index=row['index'],version=self.data['version'],column=column)
        return self.ui.docking_view_url(self.project_id,self.run_id,self.stage_id,'docking',selection,self.token,name=str(row['id'])+' · Docking')

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
            lambda:self.valid() and sequence==self.sequence,lambda:self.pose_preview(row,selection['column']),
            open_3d_url=self.view_url(row,selection['column']))

    async def footprint(self,row):
        if not self.valid():return
        from .footprint import footprint_dialog
        self.preview_version+=1;version=self.preview_version
        selection=dict(table=self.selector.value,index=row['index'],version=self.data['version'])
        await footprint_dialog(self.ui,self.project_id,self.run_id,self.stage_id,selection,row['id'],
            lambda:self.valid() and version==self.preview_version,
            origin=row['properties'].get('footprint_origin'),receptor=row['properties'].get('receptor_id'))
