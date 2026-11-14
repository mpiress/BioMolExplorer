"""BioMolExplorer V5 authentication presentation using bundled local assets."""
import flet as ft

from .branding import BRAND_BLUE, BRAND_NAVY, image_bytes

BLUE = '#075ACB'
MUTED = '#63748E'
BORDER = '#C9D4E3'


class LoginExperience:
    """Responsive V5 card; authentication remains owned by ``WorkspaceUI``."""

    def __init__(self, ui, register, name, email, password, confirmation, error, submit, switch):
        self.ui = ui
        self.register = register
        self.active = True
        self.entry_task = None
        self.groups = []
        self.submit = submit
        self.error = error
        self.name = name
        self.email = email
        self.password = password
        self.confirmation = confirmation

        for field in (name, email, password, confirmation):
            field.label = None
            field.border_radius = 28
            field.border_color = BORDER
            field.focused_border_color = BLUE
            field.focused_border_width = 1.5
            field.height = 70
            field.text_size = 18
            field.content_padding = ft.Padding.symmetric(horizontal=23, vertical=15)
            field.filled = True
            field.fill_color = '#FFFFFF'
        name.hint_text = 'Seu nome'
        name.prefix_icon = ft.Icons.PERSON_OUTLINE
        email.hint_text = 'Usuário ou e-mail'
        email.prefix_icon = ft.Icons.MAIL_OUTLINE
        password.hint_text = 'Senha'
        password.prefix_icon = ft.Icons.LOCK_OUTLINE
        confirmation.hint_text = 'Confirme a senha'
        confirmation.prefix_icon = ft.Icons.LOCK_OUTLINE

        submit.height = 70
        submit.icon = ft.Icons.LOGIN if not register else ft.Icons.PERSON_ADD_OUTLINED
        submit.style = ft.ButtonStyle(
            bgcolor={ft.ControlState.DEFAULT: BLUE, ft.ControlState.HOVERED: '#064EAE',
                     ft.ControlState.DISABLED: '#9DB8DC'},
            color='#FFFFFF', icon_color='#FFFFFF',
            shape=ft.RoundedRectangleBorder(radius=34),
            text_style=ft.TextStyle(size=19, weight=ft.FontWeight.W_700),
        )

        self.error.expand = True
        self.error.size = 13
        self.error_box = ft.Container(
            ft.Semantics(live_region=True, content=ft.Row([
                ft.Icon(ft.Icons.ERROR_OUTLINE, size=18, color='#B4233A'), error,
            ], spacing=8, vertical_alignment=ft.CrossAxisAlignment.START)),
            bgcolor='#FFF1F3', border=ft.Border.all(1, '#FFD6DE'), border_radius=12,
            padding=10, visible=False,
        )

        self.remember = ft.Checkbox(label='Manter conectado', value=False, active_color=BLUE)
        self.forgot = ft.TextButton('Esqueceu a senha?', on_click=self.begin_recovery,
            style=ft.ButtonStyle(color=BRAND_BLUE, padding=0))
        options = ft.Row([self.remember, self.forgot],
            alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
            vertical_alignment=ft.CrossAxisAlignment.CENTER)

        self.title = ft.Text('Criar conta' if register else 'Entrar', color=BRAND_NAVY,
            size=50, weight=ft.FontWeight.W_700)
        self.caption = ft.Text(
            'Crie sua conta para acessar o BioMolExplorer.' if register else
            'Acesse seu workspace do BioMolExplorer.', color=MUTED, size=19)
        fields = ([name] if register else []) + [email, password] + ([confirmation] if register else [])
        controls = [self.title, self.caption, ft.Container(height=28)]
        for index, field in enumerate(fields):
            controls.append(field)
            if index < len(fields) - 1:
                controls.append(ft.Container(height=12))
        if not register:
            controls.extend([ft.Container(height=5), options])
        else:
            controls.extend([ft.Container(height=4),
                ft.Text('Use uma senha com pelo menos 10 caracteres.', size=12, color=MUTED)])
        controls.extend([ft.Container(height=12), self.error_box, submit,
            ft.Container(height=13), ft.Divider(color='#DCE4EF', height=1),
            ft.Container(height=5), ft.Row([switch], alignment=ft.MainAxisAlignment.CENTER)])
        self.login_controls=controls
        self.form = ft.Column(controls, spacing=0, scroll=ft.ScrollMode.AUTO,
            alignment=ft.MainAxisAlignment.CENTER,
            horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
        self.right = ft.Container(self.form, bgcolor='#FFFFFF')

        self.panel_image = ft.Image(src=image_bytes('login_v5_scientific_panel.png'),
            fit=ft.BoxFit.FILL, semantics_label='Painel científico BioMolExplorer')
        self.left = ft.Container(self.panel_image, bgcolor='#083458')

        self.card = ft.Container(
            ft.Row([self.left, self.right], spacing=0), bgcolor='#FFFFFF', border_radius=35,
            clip_behavior=ft.ClipBehavior.HARD_EDGE,
            shadow=ft.BoxShadow(blur_radius=40, color='#250C4270', offset=ft.Offset(0, 17)))
        self.background = ft.Image(src=image_bytes('login_v5_background.png'), fit=ft.BoxFit.FILL,
            semantics_label='Fundo decorativo da autenticação')
        self.center = ft.Container(self.card, alignment=ft.Alignment.CENTER, expand=True)
        self.control = ft.Stack([self.background, self.center], expand=True)
        self.resize()

    def style_recovery_field(self,field):
        field.border_radius=28;field.border_color=BORDER;field.focused_border_color=BLUE
        field.focused_border_width=1.5;field.height=70;field.text_size=18
        field.content_padding=ft.Padding.symmetric(horizontal=23,vertical=15)
        field.filled=True;field.fill_color='#FFFFFF'
        return field

    def recovery_button(self,label,handler,icon=ft.Icons.ARROW_FORWARD):
        return ft.Button(label,icon=icon,on_click=handler,height=70,style=self.submit.style)

    def recovery_shell(self,title,caption,controls):
        self.recovery_error=ft.Text('',size=13,color='#B4233A')
        self.form.controls=[ft.Text(title,color=BRAND_NAVY,size=self.title.size,weight=ft.FontWeight.W_700),
            ft.Text(caption,color=MUTED,size=self.caption.size),ft.Container(height=28),*controls,
            ft.Container(height=8),self.recovery_error]
        self.ui.page.update()

    def begin_recovery(self,event=None):
        self.recovery_email=(self.email.value or '').strip()
        self.recovery_email_field=self.style_recovery_field(ft.TextField(
            hint_text='Endereço de e-mail',value=self.recovery_email,prefix_icon=ft.Icons.MAIL_OUTLINE,
            keyboard_type=ft.KeyboardType.EMAIL,autofocus=True))
        send=self.recovery_button('Enviar código',self.send_recovery_code,ft.Icons.SEND_OUTLINED)
        back=ft.TextButton('Voltar ao login',icon=ft.Icons.ARROW_BACK,on_click=self.back_to_login)
        self.recovery_shell('Recuperar senha','Informe o e-mail associado à sua conta.',
            [self.recovery_email_field,ft.Container(height=20),send,ft.Row([back],alignment=ft.MainAxisAlignment.CENTER)])

    async def send_recovery_code(self,event=None):
        try:
            code=await self.ui.call(self.ui.recovery_service.request,self.recovery_email_field.value or '')
            self.recovery_email=(self.recovery_email_field.value or '').strip().lower()
        except Exception as exc:
            self.recovery_error.value=str(exc);self.ui.page.update();return
        self.ui.notify(self.ui.tr('Código enviado para seu e-mail.') if code is None else
            self.ui.tr('Código de teste: {0}').format(code))
        self.show_code_step()

    def show_code_step(self):
        self.code_field=self.style_recovery_field(ft.TextField(hint_text='Código de 8 dígitos',
            prefix_icon=ft.Icons.PIN_OUTLINED,keyboard_type=ft.KeyboardType.NUMBER,max_length=8,autofocus=True))
        self.code_field.on_submit=self.validate_recovery_code
        validate=self.recovery_button('Validar código',self.validate_recovery_code,ft.Icons.VERIFIED_OUTLINED)
        resend=ft.TextButton('Reenviar código',on_click=self.resend_recovery_code)
        back=ft.TextButton('Voltar ao login',icon=ft.Icons.ARROW_BACK,on_click=self.back_to_login)
        self.recovery_shell('Verifique seu e-mail',
            self.ui.tr('Informe o código enviado para {0}.').format(self.recovery_email),
            [self.code_field,ft.Container(height=20),validate,
             ft.Row([resend,back],alignment=ft.MainAxisAlignment.SPACE_BETWEEN)])

    async def resend_recovery_code(self,event=None):
        try:code=await self.ui.call(self.ui.recovery_service.request,self.recovery_email)
        except Exception as exc:
            self.recovery_error.value=str(exc);self.ui.page.update();return
        self.ui.notify(self.ui.tr('Novo código enviado para seu e-mail.') if code is None else
            self.ui.tr('Novo código de teste: {0}').format(code))
        self.code_field.value='';self.recovery_error.value='';self.ui.page.update()

    async def validate_recovery_code(self,event=None):
        try:self.reset_grant=await self.ui.call(
            self.ui.recovery_service.validate,self.recovery_email,self.code_field.value or '')
        except Exception as exc:
            self.recovery_error.value=str(exc);self.ui.page.update();return
        self.show_password_step()

    def show_password_step(self):
        self.new_password=self.style_recovery_field(ft.TextField(hint_text='Nova senha',
            prefix_icon=ft.Icons.LOCK_OUTLINE,password=True,can_reveal_password=True,autofocus=True))
        self.new_confirmation=self.style_recovery_field(ft.TextField(hint_text='Confirme a nova senha',
            prefix_icon=ft.Icons.LOCK_OUTLINE,password=True,can_reveal_password=True))
        self.new_confirmation.on_submit=self.save_new_password
        save=self.recovery_button('Salvar nova senha',self.save_new_password,ft.Icons.SAVE_OUTLINED)
        back=ft.TextButton('Voltar ao login',icon=ft.Icons.ARROW_BACK,on_click=self.back_to_login)
        self.recovery_shell('Crie uma nova senha','Use uma senha com pelo menos 10 caracteres.',
            [self.new_password,ft.Container(height=12),self.new_confirmation,
             ft.Container(height=20),save,ft.Row([back],alignment=ft.MainAxisAlignment.CENTER)])

    async def save_new_password(self,event=None):
        if self.new_password.value!=self.new_confirmation.value:
            self.recovery_error.value=self.ui.tr('As senhas não coincidem.');self.ui.page.update();return
        try:await self.ui.call(self.ui.recovery_service.reset,self.recovery_email,
            self.reset_grant,self.new_password.value or '')
        except Exception as exc:
            self.recovery_error.value=str(exc);self.ui.page.update();return
        self.ui.notify(self.ui.tr('Senha atualizada. Entre com sua nova senha.'))
        self.ui.show_login(False,{'email':self.recovery_email})

    def back_to_login(self,event=None):
        self.ui.show_login(False,{'email':getattr(self,'recovery_email',self.email.value or '')})

    def resize(self):
        width = self.ui.page.width or 1440
        height = max(420, (self.ui.page.height or 900) - 76)
        self.background.width = width
        self.background.height = height
        compact = width < 900
        if compact:
            card_width = max(292, min(width - 28, 570))
            card_height = min(max(560, height - 28), 760 if self.register else 660)
            self.left.visible = False
            self.right.width = card_width
            self.right.height = card_height
            inset = 23 if card_width < 420 else 32
            self.right.padding = ft.Padding.symmetric(horizontal=inset, vertical=24)
            self.card.width = card_width
            self.card.height = card_height
            self.title.size = 40
            self.caption.size = 15
            density = .83
        else:
            card_width = min(1432, width - 72)
            card_height = min(754, max(610, height - 56))
            left_width = round(card_width * (841 / 1432))
            right_width = card_width - left_width
            self.left.visible = True
            self.left.width = left_width
            self.left.height = card_height
            self.panel_image.width = left_width
            self.panel_image.height = card_height
            self.right.width = right_width
            self.right.height = card_height
            inset = max(25, min(45, right_width * .077))
            self.right.padding = ft.Padding.symmetric(horizontal=inset, vertical=24)
            self.card.width = card_width
            self.card.height = card_height
            self.title.size = 50 if card_width >= 1170 else 42
            self.caption.size = 19 if card_width >= 1170 else 16
            density = max(.75, min(1.0, (card_height - 40) / 714))
        for field in (self.name, self.email, self.password, self.confirmation):
            field.height = max(58, round(70 * density))
            field.text_size = 18 if density >= .9 else 16
        self.submit.height = max(60, round(70 * density))
        self.center.height = height

    def set_busy(self, busy):
        self.submit.disabled = busy
        self.submit.icon = (ft.ProgressRing(width=18, height=18, stroke_width=2, color='#FFFFFF')
                            if busy else ft.Icons.PERSON_ADD_OUTLINED if self.register else ft.Icons.LOGIN)
        self.error_box.visible = bool(self.error.value)
        self.resize()

    def prepare_entry(self):
        pass

    def start(self):
        pass

    def stop(self):
        self.active = False
        if self.entry_task is not None:
            self.entry_task.cancel()
            self.entry_task = None

    def valid(self):
        return self.active and not self.ui.token and getattr(self.ui, '_login_experience', None) is self
