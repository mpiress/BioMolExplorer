"""Native Flet login presentation; credentials and authentication stay in app.py.

Motion uses client-side implicit transforms, one finite reveal task and a slow
animation-end transition. No frame loop, provider request or chemistry engine.
"""
import asyncio

import flet as ft
import flet.canvas as canvas

from .branding import BRAND_BLUE, BRAND_NAVY, logo

MINT = '#88E4D0'
MUTED = '#64748B'
EASE = ft.AnimationCurve.EASE_OUT_CUBIC
SITE_URL = 'https://mpiress.github.io/BioMolExplorer/'
HERO_NAVY = '#07172D'
HERO_PANEL = '#0F2945'
CYAN = '#00B8E6'


class LoginExperience:
    """Responsive artwork/form composition with independently cancellable motion."""
    def __init__(self, ui, register, name, email, password, confirmation, error, submit, switch):
        self.ui = ui
        self.register = register
        self.active = True
        self.motion_enabled = True
        self.entry_task = None
        self.drift_direction = False
        self.groups = []
        self.submit = submit
        self.error = error
        self.error.expand = True
        self.error.size = 12
        self.error_box = ft.Container(
            ft.Semantics(live_region=True, content=ft.Row([
                ft.Icon(ft.Icons.ERROR_OUTLINE, size=18, color='#B4233A'), error,
            ], spacing=10, vertical_alignment=ft.CrossAxisAlignment.START)),
            bgcolor='#FFF1F3', border=ft.Border.all(1, '#FFD6DE'), border_radius=12,
            padding=12, visible=False, animate_opacity=160,
        )
        for field in (name, email, password, confirmation):
            field.width = None
            field.text_size = 14
            field.border = ft.OutlineInputBorder(border_radius=12, side=ft.BorderSide(1, '#DDE4F0'))
            field.border_color = '#DDE4F0'
            field.focused_border_color = BRAND_BLUE
            field.focused_border_width = 1.5
            field.filled = True
            field.fill_color = '#F8FAFE'
            field.content_padding = ft.Padding.symmetric(horizontal=16, vertical=17)
        name.prefix_icon = ft.Icons.PERSON_OUTLINE
        email.prefix_icon = ft.Icons.ALTERNATE_EMAIL
        email.hint_text = 'nome@instituicao.org'
        password.prefix_icon = confirmation.prefix_icon = ft.Icons.LOCK_OUTLINE
        submit.height = 54
        submit.style = ft.ButtonStyle(
            bgcolor={ft.ControlState.DEFAULT: BRAND_BLUE, ft.ControlState.HOVERED: '#1D44CC', ft.ControlState.DISABLED: '#A7B9E4'},
            color='#FFFFFF', icon_color='#FFFFFF', padding=16,
            shape=ft.RoundedRectangleBorder(radius=13), animation_duration=180,
            text_style=ft.TextStyle(size=14, weight=ft.FontWeight.W_600),
        )
        self.submit_shell = ft.Container(submit, border_radius=13, animate_scale=180,
            animate=ft.Animation(180, EASE), on_hover=self.hover_button)
        self.motion_button = ft.IconButton(ft.Icons.MOTION_PHOTOS_ON_OUTLINED,
            icon_size=17, icon_color=MUTED, tooltip='Reduzir animações',
            on_click=self.toggle_motion, width=34, height=34)
        header = ft.Row([
            ft.Container(logo(30), width=48, height=48, padding=8,
                bgcolor='#EDF2FF', border_radius=14, border=ft.Border.all(1, '#E2E9FC')),
            ft.Text('ESPAÇO DE PESQUISA', size=9, weight=ft.FontWeight.W_600,
                color=MUTED, expand=True, text_align=ft.TextAlign.RIGHT,
                style=ft.TextStyle(letter_spacing=1.3)),
        ], spacing=12)
        heading = ft.Column([
            ft.Text('Crie seu workspace.' if register else 'Bem-vindo de volta.',
                size=30, weight=ft.FontWeight.W_700, color=BRAND_NAVY),
            ft.Text('Organize seus estudos e colabore com sua equipe.' if register else
                'Entre para continuar sua pesquisa.', size=13, color=MUTED),
        ], spacing=8)
        fields = ([name] if register else []) + [email, password] + ([confirmation] if register else [])
        controls = [self.reveal_group(header), self.reveal_group(heading),
                    *[self.reveal_group(field) for field in fields]]
        if register:
            controls.append(self.reveal_group(ft.Text('Use uma senha com pelo menos 10 caracteres.', size=11, color=MUTED)))
        controls += [self.error_box, self.reveal_group(self.submit_shell), self.reveal_group(switch),
            self.reveal_group(ft.Row([
                ft.Icon(ft.Icons.LOCK_OUTLINE, size=15, color=MUTED),
                ft.Text('Seus projetos e dados permanecem sob seu controle.', size=11, color=MUTED, expand=True),
                self.motion_button,
            ], spacing=9))]
        self.card = ft.Container(
            ft.Column(controls, spacing=19, scroll=ft.ScrollMode.AUTO,
                horizontal_alignment=ft.CrossAxisAlignment.STRETCH),
            padding=36, bgcolor='#F5FFFFFF', blur=ft.Blur(12, 12), border_radius=26,
            animate=ft.Animation(200, EASE),
            border=ft.Border.all(1, '#FFFFFF'),
            shadow=[ft.BoxShadow(blur_radius=45, spread_radius=-12, color='#19234372', offset=ft.Offset(0, 16)),
                    ft.BoxShadow(blur_radius=1, color='#0D17345B', offset=ft.Offset(0, 1))],
        )
        self.form_frame = ft.Container(self.card, expand=True, alignment=ft.Alignment.CENTER)
        self.hero = self.build_hero()
        self.shell = ft.Container(ft.Row([self.hero, self.form_frame], spacing=32,
            vertical_alignment=ft.CrossAxisAlignment.STRETCH))
        self.control = ft.Container(expand=True, alignment=ft.Alignment.CENTER,
            gradient=ft.LinearGradient(colors=['#EFF3FC', '#FAFBFE', '#EEF5F8'],
                begin=ft.Alignment.TOP_LEFT, end=ft.Alignment.BOTTOM_RIGHT),
            content=ft.Column([self.shell], alignment=ft.MainAxisAlignment.CENTER,
                horizontal_alignment=ft.CrossAxisAlignment.CENTER, scroll=ft.ScrollMode.AUTO, expand=True))
        self.resize()

    def reveal_group(self, content):
        group = ft.Container(content, animate_opacity=ft.Animation(360, EASE),
            animate_offset=ft.Animation(420, EASE))
        self.groups.append(group)
        return group

    def build_hero(self):
        self.glow = ft.Container(width=480, height=480, left=-180, top=-130,
            gradient=ft.RadialGradient(colors=['#403974D9', '#0011274C'], radius=.65),
            animate_offset=ft.Animation(9000, ft.AnimationCurve.EASE_IN_OUT),
            on_animation_end=self.drift)
        self.spotlight = ft.Container(width=400, height=400, left=150, top=230,
            gradient=ft.RadialGradient(colors=['#253BD9C1', '#000A1738'], radius=.65),
            animate_offset=ft.Animation(350, EASE))
        # Top/bottom bounds create a real presentation region.  The artwork can
        # then be vertically centered instead of matching its own exact height.
        self.art = ft.Container(self.scientific_visualization(), left=22, right=22, top=330, bottom=62,
            alignment=ft.Alignment.CENTER,
            animate_offset=ft.Animation(450, EASE))
        self.hero_title = ft.Text('Exploração de dados moleculares para', size=36, color='#FFFFFF', weight=ft.FontWeight.W_700)
        self.hero_accent = ft.Text('descoberta inteligente de fármacos', size=36, color=MINT, weight=ft.FontWeight.W_700)
        title = ft.Container(ft.Column([
            ft.Row([ft.Icon(ft.Icons.BIOTECH_OUTLINED, size=16, color=MINT),
                ft.Text('PLATAFORMA PARA DESCOBERTA DE FÁRMACOS', size=10, color='#C6D7EE',
                    weight=ft.FontWeight.W_600, style=ft.TextStyle(letter_spacing=1.7))], spacing=10),
            ft.Container(height=10),
            ft.Column([self.hero_title, self.hero_accent], spacing=0),
            ft.Text('Uma plataforma computacional integrada para exploração molecular e recuperação de informações, '
                'projetada para apoiar a descoberta de fármacos assistida por computador por meio de análise de '
                'similaridade molecular, triagem virtual, validação por redocking e inspeção 3D interativa de poses '
                'de docking e interações receptor-ligante.', size=11.5, color='#BDCCE2',
                style=ft.TextStyle(height=1.35)),
            self.capability_grid(),
        ], spacing=10), left=36, right=32, top=28)
        async def official_site(event):
            await ft.UrlLauncher().launch_url(SITE_URL,mode=ft.LaunchMode.EXTERNAL_APPLICATION,
                web_only_window_name='_blank')
        footer = ft.Container(ft.TextButton('Conheça o site oficial do BioMolExplorer',icon=ft.Icons.OPEN_IN_NEW,
                on_click=official_site,style=ft.ButtonStyle(color=MINT,padding=0)),
            left=34, right=34, bottom=18, padding=ft.Padding.symmetric(horizontal=14, vertical=4), border_radius=13,
            bgcolor='#0FFFFFFF', border=ft.Border.all(1, '#20FFFFFF'))
        visual = ft.Stack([self.glow, self.spotlight, self.art, title, footer], expand=True)
        self.gesture = ft.GestureDetector(visual, on_hover=self.pointer, on_exit=self.leave,
            hover_interval=80, expand=True)
        return ft.Container(self.gesture, border_radius=28,
            gradient=ft.LinearGradient(colors=[HERO_NAVY, '#0F213A', '#0B3043'],
                begin=ft.Alignment.TOP_LEFT, end=ft.Alignment.BOTTOM_RIGHT),
            clip_behavior=ft.ClipBehavior.ANTI_ALIAS,
            border=ft.Border.all(1, '#122E55'))

    def scientific_visualization(self):
        """Detailed native counterpart of the official site's carbon data field."""
        stroke=lambda color,width=1:ft.Paint(color=color,stroke_width=width,style=ft.PaintingStyle.STROKE)
        fill=lambda color:ft.Paint(color=color)
        shapes=[]
        # Molecular similarity network.
        network=[(28,77,7,'#5EEAD4'),(64,43,5,CYAN),(82,92,9,'#2463EB'),
            (122,55,6,'#14B8A6'),(146,108,4,CYAN),(49,132,5,'#2463EB'),(111,145,4,'#5EEAD4')]
        links=[(28,77,64,43),(28,77,49,132),(64,43,122,55),(82,92,122,55),
            (82,92,49,132),(122,55,146,108),(49,132,111,145),(111,145,146,108)]
        shapes.extend(canvas.Line(*line,paint=ft.Paint(color='#557B9B',stroke_width=1.1)) for line in links)
        shapes.extend(canvas.Circle(x,y,r,paint=fill(color)) for x,y,r,color in network)
        shapes.append(canvas.Text(18,16,'MOLECULAR NETWORK',style=ft.TextStyle(size=8,color='#78A8BF',weight=ft.FontWeight.W_600)))
        # Carbon reference with the same orbital/data-field language as the site.
        cx,cy=238,101
        shapes.extend([canvas.Circle(cx,cy,78,paint=stroke('#2B6795',1)),
            canvas.Circle(cx,cy,57,paint=stroke('#26C8C3',1.2)),
            canvas.Circle(cx,cy,36,paint=stroke(MINT,1.1)),
            canvas.Line(cx-78,cy,cx-113,cy-22,paint=stroke('#287E98',1)),
            canvas.Line(cx+78,cy,cx+108,cy-23,paint=stroke('#287E98',1)),
            canvas.Circle(cx-61,cy-49,4,paint=fill(MINT)),canvas.Circle(cx+63,cy+45,4,paint=fill(CYAN)),
            canvas.Rect(cx-31,cy-31,62,62,border_radius=16,paint=fill('#0A3C54')),
            canvas.Rect(cx-31,cy-31,62,62,border_radius=16,paint=stroke('#62E8D5',1)),
            canvas.Text(cx-10,cy-23,'C',style=ft.TextStyle(size=27,color=MINT,weight=ft.FontWeight.W_700)),
            canvas.Text(cx-18,cy+12,'CARBON',style=ft.TextStyle(size=7,color='#8FC7D5',weight=ft.FontWeight.W_600))])
        # Fingerprint matrix card.
        shapes.extend([canvas.Rect(345,10,224,84,border_radius=12,paint=fill('#9A0A2945')),
            canvas.Rect(345,10,224,84,border_radius=12,paint=stroke('#28516D',1)),
            canvas.Text(358,20,'MOLECULAR FINGERPRINT',style=ft.TextStyle(size=8,color='#79E6DA',weight=ft.FontWeight.W_600)),
            canvas.Text(358,34,'binary descriptor field',style=ft.TextStyle(size=7,color='#7198B1'))])
        active={1,4,7,12,18,21,24,29,31,38,43,47,52,58,63,68}
        muted={3,9,15,20,27,35,41,49,56,65}
        for index in range(72):
            column=index%18;row=index//18
            color='#29D9BC' if index in active else '#147B9B' if index in muted else '#183E5B'
            shapes.append(canvas.Rect(358+column*10.6,52+row*8.2,6.4,4.5,border_radius=1,paint=fill(color)))
        # Similarity profile card.
        shapes.extend([canvas.Rect(345,105,224,86,border_radius=12,paint=fill('#9A0A2945')),
            canvas.Rect(345,105,224,86,border_radius=12,paint=stroke('#28516D',1)),
            canvas.Text(358,115,'SIMILARITY PROFILE',style=ft.TextStyle(size=8,color='#79E6DA',weight=ft.FontWeight.W_600)),
            canvas.Text(358,129,'Tanimoto · conceptual',style=ft.TextStyle(size=7,color='#7198B1')),
            canvas.Line(358,177,555,177,paint=stroke('#316079',1))])
        profile=(14,25,38,49,57,53,44,30,18,24,37,48,63,58,43,31,22,35,55,67,59,42)
        for index,height in enumerate(profile):
            shapes.append(canvas.Rect(359+index*8.7,176-height*.55,5.4,height*.55,border_radius=2,
                paint=fill('#28DDBB' if index in (4,12,19) else '#138CA5')))
        # Transformation chips visually connect the scientific concepts.
        for x,y,width,label in [(18,168,142,'STRUCTURE  →  DESCRIPTOR'),(162,181,165,'FINGERPRINT  →  GRAPH')]:
            shapes.extend([canvas.Rect(x,y,width,18,border_radius=9,paint=fill('#C20A2D49')),
                canvas.Rect(x,y,width,18,border_radius=9,paint=stroke('#2C6781',1)),
                canvas.Text(x+9,y+4,label,style=ft.TextStyle(size=7,color='#7BC4D2',weight=ft.FontWeight.W_600))])
        field=canvas.Canvas(shapes=shapes,width=590,height=200)
        visual=ft.Container(field, padding=12,bgcolor='#22102643',
            border=ft.Border.all(1,'#263E637D'),border_radius=20)
        return ft.Semantics(container=True,label='Carbono conectado a dados moleculares, similaridade, grafos, ADMET e docking',
            exclude_semantics=True,content=visual)

    def feature_card(self,icon,label):
        return ft.Container(ft.Row([ft.Container(width=3,height=28,bgcolor=MINT,border_radius=2),
            ft.Icon(icon,size=16,color=CYAN),ft.Text(label,size=10,color='#E5F2FA',expand=True)],spacing=8),
            width=170,padding=ft.Padding.symmetric(horizontal=9,vertical=7),bgcolor='#B30B2946',
            border=ft.Border.all(1,'#27416A78'),border_radius=10)

    def capability_grid(self):
        """Compact counterpart of the eight capability tiles in the official hero."""
        capabilities = [
            ('PDB', 'Estruturas proteicas'), ('ChEMBL', 'Bioatividade'),
            ('PubChem', 'Compostos'), ('ZINC', 'Triagem virtual'),
            ('Similaridade', 'Fingerprints'), ('Docking', 'Vina / DOCK6'),
            ('ADMET', 'Perfil in silico'), ('Visualização', 'Inspeção 3D'),
        ]
        tiles = ft.Row([
            ft.Container(ft.Column([
                ft.Text(title, size=9, color=MINT, weight=ft.FontWeight.W_600),
                ft.Text(detail, size=7.5, color='#AFC3DA', max_lines=1),
            ], spacing=1, tight=True), width=137, height=39,
                padding=ft.Padding.symmetric(horizontal=9, vertical=5),
                bgcolor='#8F0B2946', border=ft.Border.all(1, '#31516D'), border_radius=8)
            for title, detail in capabilities
        ], spacing=6, run_spacing=6, wrap=True,
            alignment=ft.MainAxisAlignment.CENTER, run_alignment=ft.MainAxisAlignment.CENTER)
        # A bounded full-width region makes centering observable and keeps the
        # two rows aligned independently of their intrinsic control widths.
        return ft.Container(tiles, width=590, alignment=ft.Alignment.CENTER)

    def pipeline_preview(self):
        stages=[(ft.Icons.DOWNLOAD_OUTLINED,'Recuperação'),(ft.Icons.HUB_OUTLINED,'Similaridade'),
            (ft.Icons.FAVORITE_OUTLINE,'ADMET'),(ft.Icons.GPS_FIXED,'Docking')]
        controls=[]
        for index,(icon,label) in enumerate(stages):
            if index:controls.append(ft.Icon(ft.Icons.ARROW_FORWARD,size=13,color='#6C91A9'))
            controls.append(ft.Column([ft.Icon(icon,size=16,color=MINT),ft.Text(label,size=9,color='#D8E4F4')],
                spacing=4,horizontal_alignment=ft.CrossAxisAlignment.CENTER,tight=True))
        return ft.Row(controls,alignment=ft.MainAxisAlignment.SPACE_AROUND,spacing=7)

    def valid(self):
        return self.active and not self.ui.token and getattr(self.ui, '_login_experience', None) is self

    def resize(self):
        width = self.ui.page.width or 1440
        height = self.ui.page.height or 900
        compact = width < 650
        desktop = width >= 1050
        padding = 16 if compact else 32
        self.control.padding = padding
        self.shell.width = min(1280, max(240, width - 2 * padding))
        self.shell.height = max(240, min(790, height - 76 - 2 * padding))
        self.hero.visible = desktop
        self.hero.width = self.shell.width * .53 if desktop else 0
        self.card.width = 490 if desktop else min(510, self.shell.width)
        self.card.height = min(self.shell.height, 704 if self.register else 634 if self.error.value else 570)
        self.card.padding = 24 if compact else 36
        self.card.blur = None if not desktop else ft.Blur(12, 12)
        self.card.content.spacing = 15 if compact else 19
        self.hero_title.size = self.hero_accent.size = 31 if width < 1240 else 36
        self.art.top = 330
        self.art.bottom = 62
        self.art.height = None

    def hover_button(self, event):
        if not self.valid() or self.submit.disabled or not self.motion_enabled:
            return
        hovering = str(event.data).lower() == 'true'
        self.submit_shell.scale = 1.012 if hovering else 1
        self.submit_shell.shadow = ft.BoxShadow(blur_radius=22, spread_radius=-5,
            color='#552552E8', offset=ft.Offset(0, 8)) if hovering else None
        self.ui.page.update()

    def pointer(self, event):
        if not self.valid() or not self.motion_enabled or not self.hero.visible:
            return
        position = event.local_position
        dx = max(-1, min(1, (position.x / max(self.hero.width, 1) - .5) * 2))
        dy = max(-1, min(1, (position.y / max(self.shell.height, 1) - .5) * 2))
        self.art.offset = ft.Offset(dx * .018, dy * .018)
        self.spotlight.offset = ft.Offset(dx * .10, dy * .10)
        self.ui.page.update()

    def leave(self, event):
        if not self.valid():
            return
        self.art.offset = self.spotlight.offset = ft.Offset(0, 0)
        self.ui.page.update()

    def drift(self, event=None):
        if not self.valid() or not self.motion_enabled or not self.hero.visible:
            return
        self.drift_direction = not self.drift_direction
        self.glow.offset = ft.Offset(.3, .15) if self.drift_direction else ft.Offset(0, 0)
        self.ui.page.update()

    def set_busy(self, busy):
        self.submit.disabled = busy
        self.submit.icon = ft.ProgressRing(width=17, height=17, stroke_width=2, color='#FFFFFF') if busy else ft.Icons.ARROW_FORWARD
        if busy:
            self.submit_shell.scale = 1
            self.submit_shell.shadow = None
        self.error_box.visible = bool(self.error.value)
        self.resize()

    def toggle_motion(self, event):
        self.motion_enabled = not self.motion_enabled
        if self.entry_task is not None:
            self.entry_task.cancel()
            self.entry_task = None
        self.apply_motion_preference()
        if self.motion_enabled:
            self.drift()
        self.ui.page.update()

    def apply_motion_preference(self):
        self.motion_button.icon = ft.Icons.MOTION_PHOTOS_ON_OUTLINED if self.motion_enabled else ft.Icons.MOTION_PHOTOS_PAUSED_OUTLINED
        self.motion_button.tooltip = 'Reduzir animações' if self.motion_enabled else 'Ativar animações'
        for group in self.groups:
            group.opacity = 1
            group.offset = ft.Offset(0, 0)
            group.animate_opacity = ft.Animation(360, EASE) if self.motion_enabled else None
            group.animate_offset = ft.Animation(420, EASE) if self.motion_enabled else None
        self.glow.animate_offset = ft.Animation(9000, ft.AnimationCurve.EASE_IN_OUT) if self.motion_enabled else None
        self.art.animate_offset = ft.Animation(450, EASE) if self.motion_enabled else None
        self.spotlight.animate_offset = ft.Animation(350, EASE) if self.motion_enabled else None
        self.submit_shell.animate_scale = 180 if self.motion_enabled else None
        self.submit_shell.animate = ft.Animation(180, EASE) if self.motion_enabled else None
        self.card.animate = ft.Animation(200, EASE) if self.motion_enabled else None
        if not self.motion_enabled:
            self.glow.offset = self.art.offset = self.spotlight.offset = ft.Offset(0, 0)
            self.submit_shell.scale = 1
            self.submit_shell.shadow = None

    def prepare_entry(self):
        raw = getattr(self.ui.page, '_page', self.ui.page)
        if not isinstance(raw, ft.Page):
            return
        for group in self.groups:
            group.opacity = 0
            group.offset = ft.Offset(0, .06)

    def start(self):
        raw = getattr(self.ui.page, '_page', self.ui.page)
        if isinstance(raw, ft.Page):
            self.entry_task = self.ui.page.run_task(self.reveal)

    async def reveal(self):
        try:
            try:
                features = await asyncio.wait_for(ft.SemanticsService().get_accessibility_features(), timeout=.35)
                if features.disable_animations or features.reduce_motion or features.accessible_navigation:
                    self.motion_enabled = False
            except (TimeoutError, RuntimeError, AttributeError):
                self.motion_enabled = False
            if not self.valid():
                return
            if not self.motion_enabled:
                self.apply_motion_preference()
                self.ui.page.update()
                return
            for group in self.groups:
                if not self.valid():
                    return
                group.opacity = 1
                group.offset = ft.Offset(0, 0)
                self.ui.page.update()
                await asyncio.sleep(.045)
            self.drift()
        finally:
            # Form controls remain usable if motion is interrupted or unsupported.
            if self.valid():
                self.apply_motion_preference()
                self.ui.page.update()

    def stop(self):
        self.active = False
        if self.entry_task is not None:
            self.entry_task.cancel()
            self.entry_task = None
