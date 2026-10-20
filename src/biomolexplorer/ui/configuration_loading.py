"""Centered, client-animated feedback while block configurations are prepared."""
import flet as ft

from .branding import BRAND_BLUE, BRAND_NAVY, logo


class ConfigurationLoading:
    def __init__(self, page):
        self.page = page
        self.requests = set()
        emblem = ft.Stack([
            ft.ProgressRing(width=112, height=112, stroke_width=4, color=BRAND_BLUE,
                            semantics_label='Carregando configurações…'),
            logo(64),
        ], width=112, height=112, alignment=ft.Alignment.CENTER)
        card = ft.Container(width=320, padding=24, border_radius=20, bgcolor='#FFFFFF',
            content=ft.Column([
                emblem,
                ft.Text('Carregando configurações…', size=16, weight=ft.FontWeight.W_600,
                        color=BRAND_NAVY, text_align=ft.TextAlign.CENTER),
                ft.Text('Lendo informações do bloco. Aguarde um instante.', size=13,
                        color='#64748B', text_align=ft.TextAlign.CENTER),
            ], tight=True, spacing=16, horizontal_alignment=ft.CrossAxisAlignment.CENTER))
        self.control = ft.Container(left=0, top=0, right=0, bottom=0, padding=24,
            alignment=ft.Alignment.CENTER, bgcolor='#66091D48',
            on_click=lambda e: None, content=card)

    def show(self, request):
        self.requests.add(request)
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
        overlay = getattr(self.page, 'overlay', [])
        retained = [control for control in overlay if control is not self.control]
        if len(retained) != len(overlay):
            overlay[:] = retained
            self.page.update()
