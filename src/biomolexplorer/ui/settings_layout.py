"""Shared, width-driven configuration layout. No fixed-height field rows."""
import flet as ft

FIELD_GAP=20
SECTION_GAP=24
CARD_PADDING=20
PARAMETER_HELP={
    'exhaustiveness':'Esforço de busca do Vina. Valores maiores ampliam a amostragem e o custo computacional.',
    'num_modes':'Número máximo de poses de docking solicitadas ao Vina.',
    'pH':'pH usado na preparação e protonação molecular.',
    'radius':'Raio do ambiente atômico (Morgan) ou da seleção espacial, conforme o método do bloco.',
    'morgan_n_bits':'Comprimento do vetor Morgan em bits. Afeta colisões e a representação usada na similaridade.',
    'threshold':'Limiar da seleção. A unidade e a faixa dependem do método escolhido.',
    'seed':'Semente pseudoaleatória; mantenha o valor para comparar execuções.',
    'density':'Densidade de pontos da superfície molecular para a preparação DOCK6.',
    'distance':'Distância de seleção do sítio em angstroms.',
    'chunk_size':'Quantidade de registros por lote de processamento.',
}


def field(control, *, full=False):
    """Give floating labels breathing room and wrap switches on narrow screens."""
    full=full or len(str(getattr(control,'label','') or ''))>40
    if isinstance(control,ft.Dropdown):
        full=full or any(len(str(option.text or option.key))>36 for option in control.options)
    if isinstance(control,(ft.TextField,ft.Dropdown)):
        control.width=None
        control.expand=True
        control.content_padding=ft.Padding.symmetric(horizontal=14,vertical=18)
        if isinstance(control,ft.Dropdown):
            control.menu_height=320
            control.enable_filter=True
            control.enable_search=True
    if isinstance(control,(ft.Switch,ft.Checkbox)):
        control.expand=None
        label=getattr(control,'_settings_label',None) or control.label
        control._settings_label=label
        if label:
            control.label=None
            content=ft.Row([ft.Semantics(label=label,content=control),ft.Text(label,expand=True)],spacing=12)
        else:content=control
        return ft.Container(content,col=12,
            padding=ft.Padding.symmetric(vertical=8),visible=control.visible)
    return ft.Container(control,col=12 if full else {'xs':12,'md':6},
        padding=ft.Padding.only(top=12,bottom=8),visible=control.visible)


def grid(controls, *, full=False):
    return ft.ResponsiveRow([field(c,full=full or getattr(c,'data',None) in ('wide','input')
        or not isinstance(c,(ft.TextField,ft.Dropdown,ft.Switch,ft.Checkbox))) for c in controls],
        spacing=FIELD_GAP,run_spacing=12)


def section(title,controls,description=None,*,expanded=None):
    body=ft.Column(([ft.Text(description,size=12,color='#64748B')] if description else [])+
        [grid(controls)],spacing=16,horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
    if expanded is not None:
        return ft.ExpansionTile(title=ft.Text(title,size=16,weight=ft.FontWeight.W_600),
            controls=[body],expanded=expanded,controls_padding=CARD_PADDING,tile_padding=16,
            expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH)
    return ft.Container(padding=CARD_PADDING,bgcolor='#FFFFFF',border_radius=16,
        border=ft.Border.all(1,'#E2E8F0'),content=ft.Column([
            ft.Text(title,size=17,weight=ft.FontWeight.W_600,color='#172B4D'),body],spacing=16))


def dialog_size(page):
    """Fit both small laptops and narrow windows; reserve header/actions space."""
    return dict(width=max(220,min(960,(page.width or 1160)-64)),
        height=max(140,min(660,(page.height or 900)-240)))


def add_field(tile,control,*,first=False):
    cells=tile.controls[0].controls[-1].controls
    cell=field(control,full=True)
    cells.insert(0,cell) if first else cells.append(cell)


def fit_dialog(dialog,page):
    """Shared content-sized dialogs with viewport-bounded overflow.

    AlertDialog keeps its natural height while it fits.  ``scrollable`` only
    activates the Material dialog viewport when that natural content exceeds
    the available screen, avoiding the large empty body formerly visible in
    short confirmation dialogs.
    """
    if not isinstance(dialog,ft.AlertDialog):return
    dialog.inset_padding=ft.Padding.symmetric(horizontal=24,vertical=20)
    dialog.scrollable=True
    dialog.actions_padding=ft.Padding.symmetric(horizontal=20,vertical=16)
    dialog.actions_alignment=ft.MainAxisAlignment.END
    dialog.actions_overflow_button_spacing=12
    content=dialog.content
    if isinstance(content,ft.Container):
        if getattr(content,'_settings_auto_size',False):
            limit=dialog_size(page)
            content.width,content.height=limit['width'],None
            return
        limit=dict(width=max(180,(page.width or 1160)-64),height=max(100,(page.height or 900)-240))
        original=getattr(content,'_settings_requested_size',None)
        if original is None:
            original=(content.width,content.height);content._settings_requested_size=original
        width,height=original
        if isinstance(width,(int,float)):content.width=min(width,limit['width'])
        if isinstance(height,(int,float)):content.height=min(height,limit['height'])
        if isinstance(content.content,ft.Column) and content.content.scroll is None and isinstance(height,(int,float)):
            content.content.scroll=ft.ScrollMode.AUTO
    if isinstance(dialog.title,ft.Column):
        dialog.title.horizontal_alignment=ft.CrossAxisAlignment.STRETCH
