"""Structured receptor and ligand selection shared by Vina and DOCK6."""
from pathlib import Path

import flet as ft

from .artifact_picker import ArtifactPicker
from .localization import verbatim


class DockingFilePicker(ArtifactPicker):
    def __init__(self,rows,origins,page=None,extras=None,writable=True,ligands=False):
        self.ligands=ligands
        self.table=ft.DataTable(columns=[],rows=[],column_spacing=20,
            data_row_min_height=72,data_row_max_height=112,heading_row_color='#EFF6FF',
            border=ft.Border.all(1,'#E2E8F0'),border_radius=12)
        self.table_host=ft.Container(content=ft.Row([self.table],scroll=ft.ScrollMode.ALWAYS))
        super().__init__(rows,origins,page,extras,writable)
        original=self.control.controls
        self.search.label='Pesquisar arquivo de ligantes' if ligands else 'Pesquisar receptor, PDB ou alvo'
        self.control.controls=[*original[:7],original[7],self.count,self.table_host,original[10]]

    def redraw(self):
        super().redraw()
        chosen=[(i,flag,ref) for i,(flag,ref) in enumerate(self.rows) if flag.value]
        self.offset=min(self.offset,max(0,((len(chosen)-1)//self.PAGE_SIZE)*self.PAGE_SIZE))
        headings=('Arquivo de ligantes','Origem','Triagem de compostos','Ações') if self.ligands else (
            'PDB / receptor','Alvo / origem','Estado','Ações')
        self.table.columns=[ft.DataColumn(ft.Text(label)) for label in headings]
        table_rows=[]
        for index,flag,ref in chosen[self.offset:self.offset+self.PAGE_SIZE]:
            def remove(e,index=index):self.set_selected([index],False)
            name=Path(getattr(flag,'physical_file',flag.label)).name
            origin=self.origins.get(self.origin(ref),self.origin(ref))
            available=getattr(flag,'available',True)
            extras=list(self.extra_controls(flag))
            choices=[control for control in extras if isinstance(control,ft.Column)]
            buttons=[control for control in extras if control not in choices]
            def configure(e,choices=choices,name=name):
                def close(event):
                    self.page.pop_dialog();self.update()
                dialog=ft.AlertDialog(title=verbatim(ft.Text(name),'value'),
                    content=ft.Column(choices,tight=True,spacing=12,width=480),
                    actions=[ft.TextButton('Concluir',on_click=close)],scrollable=True)
                self.page.show_dialog(dialog)
            detail=ft.Text('Todos os compostos' if self.ligands else 'Receptor selecionado',size=12,color='#475569')
            if choices:
                buttons.append(ft.IconButton(ft.Icons.FILTER_LIST,tooltip='Selecionar compostos',
                    on_click=configure,disabled=not self.writable or not self.page))
                choice=getattr(choices[0],'controls',[])
                selected=next((control.value for control in choice if isinstance(control,ft.Dropdown)),None)
                if selected:detail.value=selected
            status=ft.Row([ft.Icon(ft.Icons.CHECK_CIRCLE_OUTLINE if available else ft.Icons.ERROR_OUTLINE,
                color='#15803D' if available else '#B91C1C',size=18),
                ft.Text('Disponível' if available else 'Indisponível',size=12,
                    color='#15803D' if available else '#B91C1C')],spacing=8)
            actions=ft.Row([*buttons,ft.IconButton(ft.Icons.REMOVE_CIRCLE_OUTLINE,tooltip='Remover arquivo',
                on_click=remove,disabled=not self.writable)],spacing=4)
            cells=[verbatim(ft.Text(name,weight=ft.FontWeight.W_600,max_lines=2,
                    overflow=ft.TextOverflow.ELLIPSIS,tooltip=getattr(flag,'physical_file',flag.label)),'value'),
                verbatim(ft.Text(origin,max_lines=2,overflow=ft.TextOverflow.ELLIPSIS,tooltip=origin),'value'),
                detail if self.ligands else status,actions]
            table_rows.append(ft.DataRow(color='#FFFFFF' if available else '#FEF2F2',
                cells=[ft.DataCell(control) for control in cells]))
        self.table.rows=table_rows
        self.table_host.visible=bool(chosen)
        self.count.value=f'{len(chosen)} arquivo(s) selecionado(s)'
