"""Shared, client-animated processing splash with independent modal routes."""
from uuid import uuid4
import flet as ft
from .branding import BRAND_NAVY,processing_emblem

class ProcessingSplash:
    def __init__(self,page,saving=False,*,title=None,description=None):
        self.page=page;self.saving=saving;self.requests=set();self.messages={};self.dialog=None
        self.default_title=title or ('Salvando configurações…' if saving else 'Carregando configurações…')
        self.default_description=description or ('Validando e salvando as configurações do bloco. Aguarde a conclusão.' if saving else 'Lendo informações do bloco. Aguarde um instante.')
        self.control=ft.Container(left=0,top=0,right=0,bottom=0,padding=24,alignment=ft.Alignment.CENTER,bgcolor='#66091D48',on_click=lambda e:None,content=self._card(self.default_title,self.default_description))

    def _card(self,title,description):
        emblem=processing_emblem(label=title);self.ring=emblem._processing_ring
        self.title=ft.Text(title,size=16,weight=ft.FontWeight.W_600,color=BRAND_NAVY,text_align=ft.TextAlign.CENTER)
        self.description=ft.Text(description,size=13,color='#64748B',text_align=ft.TextAlign.CENTER)
        return ft.Container(width=320,padding=24,border_radius=20,bgcolor='#FFFFFF',content=ft.Column([emblem,self.title,self.description],tight=True,spacing=16,horizontal_alignment=ft.CrossAxisAlignment.CENTER))

    def show(self,request,title=None,description=None,*,modal=False,on_cancel=None):
        title=title or self.default_title;description=description or self.default_description
        self.requests.add(request);self.messages[request]=(title,description,on_cancel)
        if self.saving or modal or self.dialog is not None:
            if self.dialog is None or not self.dialog.open:
                self.dialog=ft.AlertDialog(key='processing-'+uuid4().hex,modal=True,content=self._card(title,description),content_padding=0,bgcolor=ft.Colors.TRANSPARENT,elevation=0,barrier_color='#66091D48')
            self.title.value=title;self.description.value=description;self.ring.semantics_label=title
            self.dialog.actions=[ft.TextButton('Interromper processamento',on_click=on_cancel)] if on_cancel else []
            if self.saving:self.control=self.dialog
            if not self.dialog.open:self.page.show_dialog(self.dialog)
            else:self.page.update()
            return
        if getattr(self.page,'overlay',None) is None:self.page.overlay=[]
        if not any(control is self.control for control in self.page.overlay):
            self.control.content=None;self.control.content=self._card(title,description)
            if hasattr(self.page,'localize'):self.page.localize(self.control)
            self.page.overlay.append(self.control)
        self.page.update()

    def hide(self,request):
        self.requests.discard(request);self.messages.pop(request,None)
        if not self.requests:self.close()
        elif self.messages:
            title,description,cancel=list(self.messages.values())[-1]
            self.title.value=title;self.description.value=description;self.ring.semantics_label=title
            if self.dialog is not None:self.dialog.actions=[ft.TextButton('Interromper processamento',on_click=cancel)] if cancel else []
            self.page.update()

    def close(self):
        self.requests.clear();self.messages.clear()
        if self.dialog is not None:self.dialog.open=False;self.page.update();return
        overlay=getattr(self.page,'overlay',[])
        retained=[control for control in overlay if control is not self.control]
        if len(retained)!=len(overlay):overlay[:]=retained;self.page.update()

ConfigurationLoading=ProcessingSplash
