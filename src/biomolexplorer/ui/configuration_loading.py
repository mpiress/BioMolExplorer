"""Centered, client-animated feedback while block configurations are prepared."""
import flet as ft

from .branding import BRAND_BLUE, BRAND_NAVY, logo


class ConfigurationLoading:
    def __init__(self, page, saving=False):
        self.page = page
        self.saving = saving
        self.requests = set()
        title='Salvando configurações…' if saving else 'Carregando configurações…'
        description='Validando e salvando as configurações do bloco. Aguarde a conclusão.' if saving else 'Lendo informações do bloco. Aguarde um instante.'
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

    def show(self, request):
        self.requests.add(request)
        if self.saving:
            self.page.show_dialog(self.control)
            return
        if getattr(self.page, 'overlay', None) is None:
            self.page.overlay = []
        if not any(control is self.control for control in self.page.overlay):
            self.page.localize(self.control)
            self.page.overlay.append(self.control)
            self.page.update()

    def hide(self, request):
        self.requests.discard(request)
        if not self.requests:
            self.close()

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
