"""Selected molecular files: bounded name table with always-visible actions."""
import flet as ft
from pathlib import Path

from .artifact_picker import ArtifactPicker
from .pair_table import PairTable
from .localization import verbatim


class MolecularFilePicker(ArtifactPicker):
    def __init__(self,rows,origins,page=None,extras=None,writable=True):
        self.pair_table=PairTable(('Identificador / arquivo','Ações'),page,min_width=240,pin_actions=True)
        super().__init__(rows,origins,page,extras,writable)

    def redraw(self):
        super().redraw()
        selected=[(index,flag) for index,(flag,_) in enumerate(self.rows) if flag.value]
        rendered=[]
        for index,flag in selected[self.offset:self.offset+self.PAGE_SIZE]:
            def remove(event,index=index):self.set_selected([index],False)
            actions=list(self.extra_controls(flag))+[ft.IconButton(ft.Icons.REMOVE_CIRCLE_OUTLINE,
                tooltip='Remover arquivo',on_click=remove,disabled=not self.writable)]
            identifier=Path(getattr(flag,'physical_file',flag.label)).name
            label=identifier+'\n'+self.label(flag) if identifier!=self.label(flag) else identifier
            name=ft.Text(label,width=180,max_lines=2,overflow=ft.TextOverflow.ELLIPSIS,
                tooltip=getattr(flag,'physical_file',flag.label))
            rendered.append(ft.DataRow(cells=[ft.DataCell(verbatim(name,'value','tooltip')),
                ft.DataCell(ft.Row(actions,spacing=4))]))
        self.pair_table.set_rows(rendered)
        self.selected.controls=[self.pair_table.control] if rendered else [ft.Text('Nenhum arquivo selecionado.',size=12)]
