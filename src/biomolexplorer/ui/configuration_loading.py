"""Standard branded, client-animated splash for BioMolExplorer processing."""
import flet as ft

from .branding import BRAND_NAVY, processing_emblem


class ProcessingSplash:
    def __init__(self, page):
        self.page = page
        self.requests = set()
        self.messages = {}
        self.dialog = None
        self.control = ft.Container(left=0, top=0, right=0, bottom=0, padding=24,
            alignment=ft.Alignment.CENTER, bgcolor='#66091D48',
            on_click=lambda e: None, content=self._card())

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
        self.messages[request]=(title,description,on_cancel)
        self.title.value=title;self.description.value=description
        self.ring.semantics_label=title
        if modal or self.dialog is not None:
            # Page overlays are painted below native dialogs. Saving from the
            # settings dialog needs its own modal route above the editor.
            if self.dialog is None:
                overlay=getattr(self.page,'overlay',[])
                if any(control is self.control for control in overlay):
                    overlay[:]=[control for control in overlay if control is not self.control]
                    self.page.update()
                # A previous dialog may still be animating out on the client.
                # Give the new route its own controls instead of reparenting them.
                card=self._card(title,description)
                self.dialog = ft.AlertDialog(modal=True, content=card,
                    content_padding=0, bgcolor='#FFFFFF')
                if on_cancel:self.dialog.actions=[ft.TextButton('Interromper processamento',on_click=on_cancel)]
                self.page.show_dialog(self.dialog)
            else:
                self.dialog.actions=[ft.TextButton('Interromper processamento',on_click=on_cancel)] if on_cancel else []
                self.page.update()
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
        self.messages.clear()
        if self.dialog is not None:
            from .feedback import close_dialog
            close_dialog(self.page, self.dialog)
            self.dialog = None
        overlay = getattr(self.page, 'overlay', [])
        retained = [control for control in overlay if control is not self.control]
        if len(retained) != len(overlay):
            overlay[:] = retained
            self.page.update()


# Existing configuration callers share the same standard component.
ConfigurationLoading=ProcessingSplash
