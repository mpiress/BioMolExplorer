"""Original project logos bundled for web, desktop and installed-package use."""
from functools import lru_cache
from pathlib import Path
import flet as ft

BRAND_BLUE='#2552E8'
BRAND_NAVY='#091D48'

@lru_cache(maxsize=4)
def image_bytes(name):
    return (Path(__file__).parents[1]/'resources'/'branding'/name).read_bytes()

def logo(size=48):
    return ft.Image(src=image_bytes('logo.png'),width=size,height=size,fit=ft.BoxFit.CONTAIN,semantics_label='Logo BioMolExplorer')

def wordmark(width=240):
    return ft.Image(src=image_bytes('wordmark.png'),width=width,height=width*209/662,fit=ft.BoxFit.CONTAIN,semantics_label='BioMolExplorer')

def flag(language):
    label='Português (Brasil)' if language=='pt' else 'English (United States)'
    return ft.Image(src=image_bytes('br.svg' if language=='pt' else 'us.svg'),
        width=28,height=20,fit=ft.BoxFit.CONTAIN,semantics_label=label)
