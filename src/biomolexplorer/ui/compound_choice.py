"""Bounded presentation of compound IDs; all compounds remain selectable."""
import flet as ft
from .localization import verbatim


class CompoundChoice:
    PAGE_SIZE=100

    def __init__(self,choice,page):
        self.choice,self.page=choice,page
        self.codes=[];self.offset=0
        self.query=ft.TextField(label='Pesquisar código do composto',on_change=self.search)
        self.previous=ft.TextButton('Anterior',on_click=lambda e:self.move(-1))
        self.next=ft.TextButton('Próxima',on_click=lambda e:self.move(1))
        self.count=ft.Text(size=12)
        self.navigation=ft.Row([self.previous,self.count,self.next],wrap=True)
        self.control=ft.Column([self.query,choice,self.navigation],spacing=12,
            horizontal_alignment=ft.CrossAxisAlignment.STRETCH)

    def set_records(self,records):
        self.codes=sorted(records);self.offset=0;self.query.value='';self.refresh()

    def refresh(self):
        query=(self.query.value or '').casefold()
        found=[code for code in self.codes if query in code.casefold()]
        self.offset=min(self.offset,max(0,(len(found)-1)//self.PAGE_SIZE*self.PAGE_SIZE))
        keys=found[self.offset:self.offset+self.PAGE_SIZE]
        if self.choice.value and self.choice.value not in keys:keys=[*keys,self.choice.value]
        self.choice.options=[ft.DropdownOption(key='',text='Todos os compostos')]+[
            verbatim(ft.DropdownOption(key=code,text=code)) for code in keys]
        self.query.visible=self.navigation.visible=len(self.codes)>self.PAGE_SIZE
        self.previous.disabled=self.offset==0
        self.next.disabled=self.offset+self.PAGE_SIZE>=len(found)
        self.count.value=f'{self.offset+1 if found else 0}–{min(self.offset+self.PAGE_SIZE,len(found))} / {len(found)}'

    def search(self,e=None):
        self.offset=0;self.refresh()
        if self.page:self.page.update()

    def move(self,direction):
        self.offset=max(0,self.offset+direction*self.PAGE_SIZE);self.refresh()
        if self.page:self.page.update()
