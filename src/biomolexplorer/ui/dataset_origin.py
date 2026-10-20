"""Visible source identity shared by the ADMET and molecular graph dashboards."""
import flet as ft
from .localization import verbatim


def origin_header(model):
    origin=model.get('dataset_origin')
    if not origin:return []
    return [ft.Container(padding=16,bgcolor='#EFF6FF',border_radius=12,
        content=ft.Column([ft.Text('Origem dos dados',size=12,color='#475569'),
            verbatim(ft.Text(origin['source_label'],size=18,weight=ft.FontWeight.W_600,color='#172B4D')),
            verbatim(ft.Text(origin['origin'],size=12,color='#475569')),
            ft.ExpansionTile(title=ft.Text('Arquivo e caminho no pipeline',size=12),
                controls=[verbatim(ft.Text(origin.get('source_trace') or origin.get('display_name',''),
                    size=12,selectable=True))],dense=True)],spacing=6))]
