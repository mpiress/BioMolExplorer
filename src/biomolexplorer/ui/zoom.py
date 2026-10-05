"""Shared zoom limits without viewport-dependent restrictions on zooming out."""
import flet as ft

MIN_SCALE = .001
MAX_SCALE = 100.0


def zoomable_view(content, **kwargs):
    # Flutter's default zero boundary clamps the scale to the viewport size,
    # even when min_scale permits further reduction. All four margins must
    # be infinite to remove this independent restriction.
    return ft.InteractiveViewer(
        content=content, min_scale=MIN_SCALE, max_scale=MAX_SCALE,
        boundary_margin=ft.Margin.all(float('inf')), **kwargs)
