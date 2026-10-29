"""Standard branded, client-animated splash for BioMolExplorer processing."""
import flet as ft

from .branding import BRAND_NAVY, processing_emblem


class ConfigurationLoading:
    def __init__(self, page, saving=False, *, title=None, description=None):
        self.page = page
        self.saving = saving
        self.requests = set()
        title=title or ('Salvando configurações…' if saving else 'Carregando configurações…')
        description=description or ('Validando e salvando as configurações do bloco. Aguarde a conclusão.' if saving else 'Lendo informações do bloco. Aguarde um instante.')
        emblem = ft.Stack([
            ft.ProgressRing(width=112, height=112, stroke_width=4, color=BRAND_BLUE,
                            semantics_label=title),
            logo(64),
        ], width=112, height=112, alignment=ft.Alignment.CENTER)
        card = ft.Container(width=320, padding=24, border_radius=20, bgcolor='#FFFFFF',
            content=ft.Column([
                emblem,
                ft.Text(title, size=16, weight=ft.FontWeight.W_600,
                        color=BRAND_NAVY, text_align=ft.TextAlign.CENTER),
                ft.Text(description, size=13,
                        color='#64748B', text_align=ft.TextAlign.CENTER),
            ], tight=True, spacing=16, horizontal_alignment=ft.CrossAxisAlignment.CENTER))
        self.control = ft.Container(left=0, top=0, right=0, bottom=0, padding=24,
            alignment=ft.Alignment.CENTER, bgcolor='#66091D48',
            on_click=lambda e: None, content=card)
        if saving:
            # Dialog routes sit above page overlays. Saving must cover the
            # configuration dialog, rather than animate invisibly behind it.
            self.control=ft.AlertDialog(modal=True,content=card,content_padding=0,
                bgcolor=ft.Colors.TRANSPARENT,elevation=0,barrier_color='#66091D48')

    def _card(self,title='Carregando configurações…',description='Lendo informações do bloco. Aguarde um instante.'):
        emblem=processing_emblem(label=title)
        self.ring=emblem._processing_ring
        self.title=ft.Text(title, size=16, weight=ft.FontWeight.W_600,
                        color=BRAND_NAVY, text_align=ft.TextAlign.CENTER)
        self.description=ft.Text(description, size=13,
                        color='#64748B', text_align=ft.TextAlign.CENTER)
        return ft.Container(width=320, padding=24, border_radius=20, bgcolor='#FFFFFF',
            content=ft.Column([
                emblem,
                self.title,self.description,
            ], tight=True, spacing=16, horizontal_alignment=ft.CrossAxisAlignment.CENTER))

    def show(self, request, title='Carregando configurações…',
             description='Lendo informações do bloco. Aguarde um instante.', *, modal=False,on_cancel=None):
        self.requests.add(request)
        if self.saving:
            self.page.show_dialog(self.control)
            return
        if getattr(self.page, 'overlay', None) is None:
            self.page.overlay = []
        if not any(control is self.control for control in self.page.overlay):
            # Flet skips structurally equal assignments; clear the old child
            # before installing the fresh card used by this overlay route.
            card=self._card(title,description)
            self.control.content=None
            self.control.content=card
            self.page.localize(self.control)
            self.page.overlay.append(self.control)
        self.page.update()

    def hide(self, request):
        self.requests.discard(request)
        self.messages.pop(request,None)
        if not self.requests:
            self.close()
        elif self.messages:
            title,description,cancel=list(self.messages.values())[-1]
            self.title.value=title;self.description.value=description
            self.ring.semantics_label=title
            if self.dialog is not None:
                self.dialog.actions=[ft.TextButton('Interromper processamento',on_click=cancel)] if cancel else []
            self.page.update()

    def close(self):
        self.requests.clear()
        if self.saving:
            self.control.open=False
            self.page.update()
            return
        overlay = getattr(self.page, 'overlay', [])
        retained = [control for control in overlay if control is not self.control]
        if len(retained) != len(overlay):
            overlay[:] = retained
            self.page.update()


# Existing configuration callers share the same standard component.
ConfigurationLoading=ProcessingSplash
