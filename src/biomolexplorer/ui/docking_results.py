"""Paginated scores and pose previews shared by Vina, DOCK6 and consensus."""
import flet as ft
import math
import json
from biomolexplorer.docking_results import DockingResults
from .compound_table import CompoundTableViewer
from .localization import verbatim
from .feedback import close_dialog
from .settings_layout import fit_dialog


def relevance_text(value):
    try:
        number=float(value)
    except (ValueError,TypeError):
        return '—' if value in (None,'') else str(value)
    return f'{number:.3f}' if math.isfinite(number) else '—'


def multiple_targets(fields):
    try:
        if int(fields.get('target_count') or 0)>1:return True
        targets=json.loads(fields.get('target_scores') or '[]')
        return len({item['receptor_id'] for item in targets})>1
    except (ValueError,TypeError,KeyError):return False


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
        relevance='consensus_model' in properties
        scores=['vina','dock6','normalized_score'] if consensus else [k for k in ('score','vina','dock6','z-score','min-max') if k in properties]
        columns=[('molecule_chembl_id','Código do composto')]+([] if consensus else [('receptor_id','Alvo')])+[
            (key,({'vina':'Relevância Vina','dock6':'Relevância DOCK6','normalized_score':'Consenso'} if relevance else
                  {'vina':'Score Vina','dock6':'Score DOCK6','normalized_score':'Score normalizado'}).get(key,key)) for key in scores]+([] if consensus else [('canonical_smiles','SMILES')])
        self.table.columns=[]
        for index,(key,label) in enumerate(columns):
            self.table.columns.append(ft.DataColumn(ft.Text(label),numeric=key in scores))
        self.table.columns.append(ft.DataColumn(ft.Text('Ações')))
        self.sorting.bind(self.table,[key for key,label in columns]+[None],self.sort_changed,numeric=scores)
        rows=[]
        for row in self.data['rows']:
            fields=row['properties']
            multi=consensus and multiple_targets(fields)
            async def remove(e,row=row):await self.ui.guard(lambda:self.confirm_remove(row))
            actions=[]
            for column,label in (('vina_pose','3D Vina'),('dock6_pose','3D DOCK6')) if 'vina_pose' in fields or multi else (('conformer_file','3D'),):
                async def preview(e,row=row,column=column):await self.ui.guard(lambda:self.pose_preview(row,column))
                url=self.view_url(row,column) if fields.get(column) and not multi else None
                actions.append(ft.TextButton(label,on_click=None if url else preview,disabled=not (multi or bool(fields.get(column))),
                    url=ft.Url(url,target=ft.UrlTarget.BLANK) if url else None))
            async def contacts(e,row=row):await self.ui.guard(lambda:self.contacts(row))
            actions.append(ft.TextButton('Resíduos',icon=ft.Icons.SCIENCE_OUTLINED,on_click=contacts,
                disabled=not (multi or bool(fields.get('vina_pose') or fields.get('conformer_file')))))
            if multi or fields.get('engine')=='dock6' or fields.get('dock6_pose') or fields.get('footprint_file'):
                async def footprint(e,row=row):await self.ui.guard(lambda:self.footprint(row))
                actions.append(ft.TextButton('Footprint',icon=ft.Icons.BAR_CHART,on_click=footprint,disabled=not (multi or row.get('footprint_available',False))))
            actions.append(ft.TextButton('Remover',icon=ft.Icons.DELETE_OUTLINE,on_click=remove,disabled=not self.data['can_edit']))
            values=[row['id']]+([] if consensus else [fields.get('receptor_id','—')])+[
                relevance_text(fields.get(k)) if relevance else fields.get(k,'—') for k in scores]+([] if consensus else [row['smiles'] or '—'])
            rows.append(ft.DataRow(cells=[ft.DataCell(verbatim(ft.Text(str(v),selectable=True))) for v in values]+
                                     [ft.DataCell(ft.Row(actions,spacing=2))]))
        self.table.rows=rows

    def view_url(self,row,column,target=None):
        if not getattr(self.ui.page,'web',False):return None
        selection=dict(table=self.selector.value,index=row['index'],version=self.data['version'],column=column)
        if target is not None:selection['target']=target
        return self.ui.docking_view_url(self.project_id,self.run_id,self.stage_id,'docking',selection,self.token,name=str(row['id'])+' · Docking')

    async def choose_target(self,row,action,column='vina_pose'):
        if not self.valid():return
        sequence=self.sequence;table=self.selector.value;version=self.data['version']
        choices=await self.ui.call(self.service.target_choices,self.token,self.project_id,self.run_id,self.stage_id,
            table,row['index'],version,column=column,action='footprint' if action=='footprint' else 'pose',validate=False)
        if not self.valid() or sequence!=self.sequence or table!=self.selector.value or version!=self.data['version']:return
        active=True
        options={choice['id']:choice for choice in choices}
        selector=ft.Dropdown(label='Alvo',hint_text='Pesquise e selecione um alvo',
            enable_filter=True,enable_search=True,menu_height=280,
            options=[ft.DropdownOption(key=choice['id'],text=choice['label']) for choice in choices])
        proceed=ft.TextButton('Continuar',disabled=True)
        dialog=ft.AlertDialog(title=ft.Text('Selecionar alvo'),
            content=ft.Container(ft.Column([verbatim(ft.Text(str(row['id']),selectable=True)),selector],tight=True),width=480),
            actions=[ft.TextButton('Cancelar',on_click=lambda e:cancel(e)),proceed])
        def current():
            return active and self.valid() and sequence==self.sequence and table==self.selector.value and version==self.data['version']
        def dismissed(e):
            nonlocal active
            active=False
        def cancel(e):
            dismissed(e)
            close_dialog(self.ui.page,dialog)
        def selected(e):
            if not current():return
            choice=options.get(selector.value)
            proceed.disabled=not choice or not choice['available']
            url=self.view_url(row,column,selector.value) if action=='pose' and not proceed.disabled else None
            proceed.url=ft.Url(url,target=ft.UrlTarget.BLANK) if url else None
            self.ui.page.update()
        async def execute():
            nonlocal active
            if not current() or proceed.disabled or selector.value not in options:return
            target=selector.value
            active=False;close_dialog(self.ui.page,dialog)
            if proceed.url:return
            if action=='footprint':await self.footprint(row,target)
            elif action=='contacts':await self.contacts(row,target)
            else:await self.pose_preview(row,column,target)
        async def confirm(e):
            await self.ui.guard(execute)
        selector.on_select=selected
        proceed.on_click=confirm
        dialog.on_dismiss=dismissed
        fit_dialog(dialog,self.ui.page)
        self.ui.page.show_dialog(dialog)

    async def pose_preview(self,row,column,target=None):
        if not self.valid():return
        if target is None and multiple_targets(row['properties']):
            await self.choose_target(row,'pose',column);return
        sequence=self.sequence
        selection=dict(table=self.selector.value,index=row['index'],version=self.data['version'],column=column)
        if target is not None:selection['target']=target
        if self.valid() and sequence==self.sequence:
            await self.ui.preview_docking(self.project_id,self.run_id,self.stage_id,'docking',selection,self.token)

    async def contacts(self,row,target=None):
        from .docking_scene import residue_dialog
        sequence=self.sequence
        if target is None and multiple_targets(row['properties']):
            await self.choose_target(row,'contacts','vina_pose' if row['properties'].get('vina_pose') else 'conformer_file');return
        selection=dict(table=self.selector.value,index=row['index'],version=self.data['version'],
            column='vina_pose' if row['properties'].get('vina_pose') else 'conformer_file')
        if target is not None:selection['target']=target
        await residue_dialog(self.ui,self.project_id,
            lambda:self.service.scene(self.token,self.project_id,self.run_id,self.stage_id,**selection),
            lambda:self.valid() and sequence==self.sequence,lambda:self.pose_preview(row,selection['column'],target),
            open_3d_url=self.view_url(row,selection['column'],target),
            naming_context=next((t for t in self.tables if t['path']==self.selector.value),None))

    async def footprint(self,row,target=None):
        if not self.valid():return
        if target is None and multiple_targets(row['properties']):
            await self.choose_target(row,'footprint','dock6_pose');return
        from .footprint import footprint_dialog
        self.preview_version+=1;version=self.preview_version
        selection=dict(table=self.selector.value,index=row['index'],version=self.data['version'])
        if target is not None:selection['target']=target
        await footprint_dialog(self.ui,self.project_id,self.run_id,self.stage_id,selection,row['id'],
            lambda:self.valid() and version==self.preview_version,
            origin=row['properties'].get('footprint_origin'),receptor=target or row['properties'].get('receptor_id'))
