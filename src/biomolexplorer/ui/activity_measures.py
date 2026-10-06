"""Searchable, bounded activity vocabulary; selections survive filtering/paging."""
import json
from pathlib import Path
import flet as ft
from .localization import verbatim

VOCABULARY=Path(__file__).resolve().parents[1]/'resources/crawlers/activity_types.json'


class ActivityMeasures:
    PAGE_SIZE=24

    def __init__(self,page,selected=(),writable=True):
        self.page=page;self.writable=writable;self.selected=set(selected);self.initial=list(selected);self.offset=0
        self.types=json.loads(VOCABULARY.read_text(encoding='utf-8'))['types']
        self.search=ft.TextField(label='Buscar medida de atividade',hint_text='Ex.: IC50, Ki, Inhibition',on_change=self.filter)
        self.grid=ft.ResponsiveRow(spacing=12,run_spacing=4)
        self.summary=ft.Text();self.position=ft.Text()
        self.previous=ft.IconButton(ft.Icons.CHEVRON_LEFT,tooltip='Página anterior',on_click=lambda e:self.move(-1))
        self.next=ft.IconButton(ft.Icons.CHEVRON_RIGHT,tooltip='Próxima página',on_click=lambda e:self.move(1))
        self.custom=ft.TextField(label='Outra medida (nome exato na ChEMBL)',disabled=not writable,expand=True)
        self.control=ft.Column([ft.Text('Medidas de atividade'),self.summary,self.search,self.grid,
            ft.Row([self.previous,self.position,self.next,ft.TextButton('Limpar seleção',on_click=self.clear,disabled=not writable)],wrap=True),
            ft.Row([self.custom,ft.TextButton('Adicionar medida',on_click=self.add,disabled=not writable)]),
            ft.Text('Sem seleção: qualquer medida. Outras medidas podem exigir outra unidade e permitir valores sem pChEMBL.',size=12)],spacing=12)
        self.draw()

    def values(self):
        for cell in self.grid.controls:
            check=cell.content
            if check.value:self.selected.add(check.label)
            else:self.selected.discard(check.label)
        return [v for v in dict.fromkeys(self.initial+self.types+sorted(self.selected)) if v in self.selected]

    def filtered(self):
        query=(self.search.value or '').casefold().strip()
        # Selected values first makes review/removal convenient, including custom measures.
        values=list(dict.fromkeys(sorted(self.selected)+self.types))
        return [v for v in values if query in v.casefold()]

    def draw(self):
        values=self.filtered();self.offset=min(self.offset,max(0,(len(values)-1)//self.PAGE_SIZE*self.PAGE_SIZE))
        controls=[]
        for value in values[self.offset:self.offset+self.PAGE_SIZE]:
            def changed(e,value=value):
                if e.control.value:self.selected.add(value)
                else:self.selected.discard(value)
                self.summary.value='Selecionadas: '+str(len(self.selected));self.page.update()
            check=verbatim(ft.Checkbox(label=value,value=value in self.selected,on_change=changed,disabled=not self.writable),'label')
            controls.append(ft.Container(check,col={'xs':12,'sm':6,'md':4}))
        self.grid.controls=controls
        self.summary.value='Selecionadas: '+str(len(self.selected))
        self.position.value=f'{min(self.offset+1,len(values))}–{min(self.offset+self.PAGE_SIZE,len(values))} / {len(values)}'
        self.previous.disabled=self.offset==0;self.next.disabled=self.offset+self.PAGE_SIZE>=len(values)

    def filter(self,e):self.offset=0;self.draw();self.page.update()
    def move(self,direction):self.offset=max(0,self.offset+direction*self.PAGE_SIZE);self.draw();self.page.update()
    def clear(self,e):
        if self.writable:self.selected.clear();self.draw();self.page.update()

    def add(self,e):
        value=(self.custom.value or '').strip()
        if value and self.writable:
            self.selected.add(value);self.custom.value='';self.search.value='';self.offset=0;self.draw();self.page.update()
