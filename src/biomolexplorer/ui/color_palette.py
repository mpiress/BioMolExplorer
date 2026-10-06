"""Accessible project color swatches; hexadecimal codes stay in storage."""
import flet as ft
from biomolexplorer.workspace import COLORS

NAMES = ['Verde', 'Índigo', 'Âmbar', 'Rosa', 'Azul', 'Violeta']


class ColorPalette:
    def __init__(self, page, value=None):
        self.page=page
        self.value=value if value in COLORS else COLORS[0]
        self.swatches=[]
        for color,name in zip(COLORS,NAMES):
            chip=ft.IconButton(icon=ft.Icons.CHECK if color==self.value else ft.Icons.CIRCLE,
                icon_color='#FFFFFF',bgcolor=color,tooltip=name,
                on_click=lambda e,c=color:self.choose(c))
            self.swatches.append(chip)
        self.control=ft.Column([ft.Text('Cor'),ft.Row(self.swatches,wrap=True,spacing=12)],spacing=8)

    def choose(self,color):
        if color not in COLORS:raise ValueError('Cor inválida.')
        self.value=color
        for option,chip in zip(COLORS,self.swatches):
            chip.icon=ft.Icons.CHECK if option==color else ft.Icons.CIRCLE
        self.page.update()
