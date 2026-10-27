"""Authenticated Flet workspace and pipeline editor."""
import argparse
import asyncio
import copy
import json
import os
import secrets
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import flet as ft

from biomolexplorer.catalog import TITLES, PRESETS, new_stage, template_names
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.templates import RESOURCE_ROOT
from biomolexplorer.diagnostics import configure_logging, get_logger, diagnose_exception
from biomolexplorer.visualizations import SUFFIX, MAX_VIEW_BYTES, MAX_GRAPH_BYTES, load_view
from .branding import flag, wordmark, BRAND_BLUE, BRAND_NAVY
from .removal import RemovalProgress
from .feedback import readable_error,close_dialog
from biomolexplorer.workspace import WorkspaceStore, AccessDenied, COLORS
from .project_tools import ProjectTools
from .zoom import zoomable_view
from .localization import LocalizedPage, verbatim, stage_control
from .processing import processing,background,with_processing,pick_files,save_file

INK = '#172B4D'
MUTED = '#64748B'
TEAL = BRAND_BLUE
BG = '#F3F6FC'
LINE = '#E2E8F0'
STATUS = {'queued':'Aguardando', 'running':'Executando', 'succeeded':'Concluído', 'failed':'Falhou',
          'cancelled':'Cancelado', 'interrupted':'Interrompido', 'skipped':'Não executado',
          'awaiting_input':'Aguardando seleção de arquivos'}
KINDS = {'compounds':'Compostos (CSV)', 'structures':'Complexos PDB', 'prepared_structures':'Receptores preparados',
         'molecular_poses':'Poses de ligantes (MOL2)','receptor_pdb':'Receptores (PDB)','dynamics':'Resultados de dinâmica',
         'fingerprints':'Fingerprints', 'similarity':'Relações de similaridade', 'vina':'Resultados Vina',
         'dock6':'Resultados DOCK6', 'scores':'Scores consolidados (CSV)', 'visualization':'Gráficos e visualizações',
         'zinc_urls':'Lista de downloads ZINC (TXT/URI/SH)', 'other':'Outros arquivos'}


def text(value, size=14, color=INK, weight=None, **kwargs):
    return ft.Text(value,size=size,color=color,weight=weight,**kwargs)


def user_text(value, *args, **kwargs):
    return verbatim(text(value, *args, **kwargs))


def button(label,handler,primary=False,icon=None,disabled=False):
    return ft.Button(label,on_click=handler,icon=icon,disabled=disabled,
        style=ft.ButtonStyle(bgcolor=TEAL if primary else '#FFFFFF',color='#FFFFFF' if primary else INK,
                             padding=16,shape=ft.RoundedRectangleBorder(radius=10)))


def panel(content, padding=24, **kwargs):
    kwargs.setdefault('border', ft.Border.all(1,LINE))
    return ft.Container(content=content,padding=padding,bgcolor='#FFFFFF',border_radius=16,
        **kwargs)


class WorkspaceUI(ProjectTools):
    def __init__(self,page,store,service,language='pt',structure_viewer=None):
        self.page,self.store,self.service = LocalizedPage(page,language),store,service
        from biomolexplorer.pdb_view import StructureViewers
        self.structure_viewer=structure_viewer or StructureViewers(store)
        self.language = language
        self.token = None
        self.current = None
        self.selected = None
        self.tab = 'Pipeline'
        self.dirty = False
        self.base_pipeline=[]
        self._save_lock=asyncio.Lock()
        self.current_run = None
        self.polling = None
        self.run_progress = None
        self.selection_dialog = None
        self.selection_dialog_key = None
        self.auto_selection_key = None
        self.reuse_dialog = None
        self.reuse_dialog_key = None
        self.submitting = False
        self.force_execution = False
        self.compound_viewer = None
        self.stage_results = {}
        self.results_focus = None
        self.uploads = {}
        self.picker = ft.FilePicker(on_upload=self.upload_event)
        page.title = 'BioMolExplorer · Workspace'
        page.bgcolor = BG
        page.padding = 0
        page.theme = ft.Theme(color_scheme_seed=TEAL,font_family='Roboto')
        page.theme_mode = ft.ThemeMode.LIGHT

    async def call(self,function,*args,**kwargs):
        async with processing(self,immediate=False) as feedback:
            if feedback:await feedback.start(function,args)
            return await asyncio.to_thread(function,*args,**kwargs)

    def tr(self,message):
        translator=getattr(self.page,'translator',None)
        return translator(message) if translator else message

    def notify(self,message):
        self.page.show_dialog(ft.SnackBar(content=ft.Text(str(message)),bgcolor=INK))

    async def guard(self,action):
        token = self.token
        try:
            async with processing(self,immediate=False):
                await action()
        except Exception as exc:
            get_logger('frontend').exception('Interface action failed in tab %s: %s', self.tab, exc,
                extra={'event': 'ui.action_failed', 'project_id': self.current['id'] if self.current else None, **diagnose_exception(exc)})
            if self.token != token:
                return
            if isinstance(exc,AccessDenied) and self.token:
                try:
                    await self.call(self.store.user,self.token)
                except AccessDenied:
                    if self.token != token:
                        return
                    self.close_project_views()
                    self.token = None
                    self.current = None
                    self.show_login()
                else:
                    if self.token != token:
                        return
                    if self.current:
                        project = self.current
                        project_id = project['id']
                        try:
                            fresh = await self.call(self.store.project,self.token,project_id)
                        except AccessDenied:
                            if self.token != token or self.current is not project:
                                return
                            await self.show_workspace()
                        else:
                            if self.token != token or self.current is not project:
                                return
                            if fresh['role'] != self.current['role']:
                                tab = self.tab
                                self.close_project_views()
                                await self.open_project(project_id,tab)
                    else:
                        await self.show_workspace()
            if self.token is None or self.token == token:
                self.notify(str(exc))

    def close_project_views(self):
        """Discard private dialogs/drafts when leaving a project or losing access."""
        loading=getattr(self,'_configuration_loading',None)
        if loading:loading.close()
        self.dismiss_reuse_prompt()
        viewer = getattr(self,'result_viewer',None)
        if viewer is not None:
            viewer.selection_version += 1
        self.result_viewer = None
        compound_viewer = getattr(self,'compound_viewer',None)
        if compound_viewer is not None:
            compound_viewer.preview_version += 1
        self.compound_viewer = None
        for result in getattr(self,'stage_results',{}).values():result.close()
        self.stage_results={};self.results_focus=None
        while self.page.pop_dialog() is not None:
            pass
        self.editing_stage = False
        self._dialog_request = getattr(self,'_dialog_request',0)+1
        self._compound_reference_cache = None
        self.docking_catalog = None
        self.docking_compounds = {}
        self.dirty = False
        self.selected = None
        self.current_run = None
        self.run_progress = None
        self.selection_dialog = None
        self.selection_dialog_key = None
        self.auto_selection_key = None
        self.page.on_keyboard_event = None
        self.pending_file_selections = {}
        self._configuration_dialogs = {}
        self.settings_dialog=None;self.settings_run_id=None;self.settings_body=None

    def render(self,content):
        self.page.controls.clear()
        self.page.add(content)

    def top_bar(self, on_language=None):
        self.language=getattr(self,'language','pt')
        if not isinstance(self.page,LocalizedPage):self.page=LocalizedPage(self.page,self.language)
        async def change(e):
            self.language=e.control.data
            self.page.set_language(self.language)
            language_choice.content=flag(self.language)
            language_choice.tooltip.message=self.tr('Selecionar idioma')
            if on_language:
                on_language()
            else:
                self.page.update()
        language_choice=ft.PopupMenuButton(content=flag(self.language),
            tooltip=ft.Tooltip(message=self.tr('Selecionar idioma'),exclude_from_semantics=True),
            width=48,height=48,menu_position=ft.PopupMenuPosition.UNDER,
            items=[ft.PopupMenuItem(content=flag(code),data=code,
                tooltip=label,on_click=change) for code,label in
                [('pt','Português (Brasil)'),('en','English (United States)')]])
        brand=wordmark(190)
        update_notice=getattr(self,'update_notice',None)
        bar=ft.Container(height=76,padding=ft.Padding.symmetric(horizontal=32),bgcolor='#FFFFFF',
            border=ft.Border.only(bottom=ft.BorderSide(1,LINE)),
            content=ft.Row([brand,ft.Container(expand=True),
                *([update_notice.button()] if update_notice else []),
                ft.Semantics(container=True,label='Selecionar idioma',button=True,content=language_choice)],spacing=8))
        def resize():
            compact=(getattr(self.page,'width',None) or 1440)<650
            bar.padding=ft.Padding.symmetric(horizontal=16 if compact else 32)
            brand.width=160 if compact else 190
            brand.height=brand.width*209/662
            if update_notice:update_notice.refresh()
        resize()
        return bar,language_choice,resize

    def show_login(self,register=False,values=None):
        self.language=getattr(self,'language','pt')
        if not isinstance(self.page,LocalizedPage):self.page=LocalizedPage(self.page,self.language)
        name = ft.TextField(label='Seu nome',autofocus=register)
        email = ft.TextField(label='E-mail',keyboard_type=ft.KeyboardType.EMAIL,autofocus=not register)
        password = ft.TextField(label='Senha',password=True,can_reveal_password=True)
        confirmation = ft.TextField(label='Confirme a senha',password=True,can_reveal_password=True)
        for key,field in (('name',name),('email',email),('password',password),('confirmation',confirmation)):
            field.value=(values or {}).get(key,'')
        def change_language():
            self.show_login(register,dict(name=name.value,email=email.value,password=password.value,confirmation=confirmation.value))
        top_bar,language_choice,resize_bar=self.top_bar(change_language)
        error = text('',color='#DC2626')
        submit = button('Criar conta' if register else 'Entrar no workspace',None,True,ft.Icons.ARROW_FORWARD)
        async def authenticate(e):
            submit.disabled = True
            language_choice.disabled = True
            error.value = ''
            self.page.update()
            try:
                if register and password.value != confirmation.value:
                    raise ValueError('As senhas não coincidem.')
                token = await self.call(self.store.register,name.value or '',email.value or '',password.value or '') if register else \
                        await self.call(self.store.login,email.value or '',password.value or '')
                self.token = token
                password.value = confirmation.value = ''
                await self.show_workspace()
            except Exception as exc:
                get_logger('frontend').exception('Falha na autenticação')
                error.value = str(exc)
                submit.disabled = False
                language_choice.disabled = False
                self.page.update()
        submit.on_click = authenticate
        password.on_submit = authenticate
        for field in (name,email,password,confirmation):
            field.width=None
            field.border=ft.OutlineInputBorder(border_radius=12,side=ft.BorderSide(1,'#DCE5F4'))
            field.filled=True
            field.fill_color='#F6F8FC'
            field.content_padding=ft.Padding.symmetric(horizontal=18,vertical=18)
        email.prefix_icon=ft.Icons.ALTERNATE_EMAIL
        password.prefix_icon=confirmation.prefix_icon=ft.Icons.LOCK_OUTLINE
        form=[text('CRIE SEU ESPAÇO DE PESQUISA' if register else 'BEM-VINDO DE VOLTA',11,TEAL,ft.FontWeight.W_600),
              text('Comece sua exploração.' if register else 'Sua próxima descoberta\ncomeça aqui.',30,weight=ft.FontWeight.W_700),
              text('Organize seus estudos e colabore com sua equipe.' if register else 'Entre para continuar seus projetos em um ambiente privado.',14,MUTED),ft.Container(height=8)]
        if register: form.append(name)
        form += [email,password]
        if register: form += [confirmation,text('Use uma senha com pelo menos 10 caracteres.',12,MUTED)]
        form += [error,submit,ft.TextButton('Já tenho uma conta' if register else 'Ainda não tenho uma conta',on_click=lambda e:self.show_login(not register)),
                 ft.Row([ft.Icon(ft.Icons.LOCK_OUTLINE,size=14,color=MUTED),
                    ft.Text('Você decide com quem compartilhar seus projetos.',size=11,color=MUTED,expand=True)],spacing=8)]
        login_card=ft.Container(width=560,height=730,padding=ft.Padding.symmetric(horizontal=40,vertical=32),
            bgcolor='#FFFFFF',border_radius=28,border=ft.Border.all(1,'#E1E8F5'),
            content=ft.Column(form,spacing=18,scroll=ft.ScrollMode.AUTO,horizontal_alignment=ft.CrossAxisAlignment.STRETCH))
        compact=(self.page.width or 1440)<650
        login_card.height=max(240,min(680,(self.page.height or 900)-76-(32 if compact else 64)))
        login_card.padding=24 if compact else ft.Padding.symmetric(horizontal=40,vertical=32)
        login_card.content.spacing=14 if compact else 18
        login_view=ft.Container(expand=True,padding=16 if compact else 32,
            gradient=ft.LinearGradient(colors=['#EDF2FF','#F8FAFD'],begin=ft.Alignment.TOP_LEFT,end=ft.Alignment.BOTTOM_RIGHT),
            alignment=ft.Alignment.CENTER,
            content=ft.Column([login_card],
                expand=True,scroll=ft.ScrollMode.AUTO,horizontal_alignment=ft.CrossAxisAlignment.CENTER,alignment=ft.MainAxisAlignment.CENTER))
        async def resize_login(e):
            if self.token: return
            compact=(self.page.width or 1440)<650
            resize_bar()
            login_card.height=max(240,min(680,(self.page.height or 900)-76-(32 if compact else 64)))
            login_card.padding=24 if compact else ft.Padding.symmetric(horizontal=40,vertical=32)
            login_card.content.spacing=14 if compact else 18
            login_view.padding=16 if compact else 32
            self.page.update()
        self.page.on_resize=resize_login
        self.render(ft.Column([top_bar,login_view],expand=True,spacing=0))

    def shell(self,title,subtitle,body,actions=None,groups=None,sidebar_content=None):
        from .action_help import action_row,resize_help
        user = self.store.user(self.token)
        project=getattr(self,'current',None)
        navigation=groups if groups is not None else [action_row(self,'projects',
            button('Meus projetos',self.event(self.show_workspace),icon=ft.Icons.GRID_VIEW)),*(actions or [])]
        sidebar = ft.Container(width=280,bgcolor='#EEF3FB',padding=16,
            border=ft.Border.only(right=ft.BorderSide(1,LINE)),content=ft.Column([
                ft.Column([text('Projeto' if project else 'Workspace',16,INK,ft.FontWeight.W_600),
                    *navigation,*([ft.Divider(color=LINE),sidebar_content] if sidebar_content else [])],
                    expand=True,spacing=10,scroll=ft.ScrollMode.AUTO),
                ft.Divider(color=LINE),ft.Row([ft.Icon(ft.Icons.ACCOUNT_CIRCLE_OUTLINED,color=MUTED),
                    ft.Column([user_text(user['name'],13,weight=ft.FontWeight.W_600),
                        user_text(user['email'],11,MUTED,max_lines=1,overflow=ft.TextOverflow.ELLIPSIS)],expand=True,spacing=2),
                    ft.TextButton('Sair',on_click=self.logout)],spacing=8)],expand=True,spacing=10))
        heading=user_text(title,24,weight=ft.FontWeight.W_700) if project and title==project['name'] else text(title,28,weight=ft.FontWeight.W_700)
        heading.max_lines=1;heading.overflow=ft.TextOverflow.ELLIPSIS
        summary=user_text(subtitle,size=12,color=MUTED) if project and subtitle==project['description'] else text(subtitle,12,MUTED)
        summary.max_lines=1;summary.overflow=ft.TextOverflow.ELLIPSIS
        top_bar,_,resize_bar=self.top_bar()
        def resize():
            compact=(getattr(self.page,'width',None) or 1440)<900
            resize_bar()
            sidebar.visible=getattr(self,'sidebar_visible',None) if getattr(self,'sidebar_visible',None) is not None else not compact
            sidebar.width=min(280,max(180,(getattr(self.page,'width',None) or 1440)-120))
            main.padding=12 if project else 16 if compact else 26
            menu.tooltip='Ocultar menu lateral' if sidebar.visible else 'Mostrar menu lateral'
            resize_help(self.page,getattr(self,'help_dialog',None))
            settings=getattr(self,'settings_body',None)
            if settings is not None:
                from .settings_layout import dialog_size
                settings.width=dialog_size(self.page)['width']
                settings.height=None
            selection=getattr(self,'selection_dialog',None)
            if selection is not None and selection.open:
                from .file_selection import selection_dialog_size
                for key,value in selection_dialog_size(self.page,getattr(selection,'_redocking_selection',False)).items():
                    setattr(selection.content,key,value)
        def toggle_menu(e):
            self.sidebar_visible=not sidebar.visible;resize();self.page.update()
        menu=ft.IconButton(icon=ft.Icons.MENU,tooltip='Mostrar menu lateral',on_click=toggle_menu)
        header=ft.Row([menu,ft.Column([heading,summary],expand=True,spacing=3)],spacing=8)
        main=ft.Container(expand=True,padding=12,content=ft.Column([header,body],expand=True,spacing=12))
        async def on_resize(e):
            resize();self.page.update()
        resize();self.page.on_resize=on_resize
        return ft.Column([top_bar,ft.Row([sidebar,main],expand=True,spacing=0,
            vertical_alignment=ft.CrossAxisAlignment.STRETCH)],expand=True,spacing=0)

    def event(self,coroutine,*args):
        async def handler(e):
            await self.guard(lambda:coroutine(*args))
        return handler

    async def logout(self,e):
        if self.token:
            await self.call(self.store.logout,self.token)
        self.token = None
        self.close_project_views()
        self.current = None
        self.current_run = None
        self.show_login()

    @with_processing('Carregando projetos…','Preparando a lista de projetos. Aguarde um instante.')
    async def show_workspace(self,query='',archived=False):
        token = self.token
        self.close_project_views()
        self.current = None
        self.current_run = None
        projects = await self.call(self.store.list_projects,token,archived)
        if self.token != token:
            return
        invites = await self.call(self.store.invitations,token)
        if self.token != token:
            return
        search = ft.TextField(hint_text='Buscar por nome ou descrição',value=query,prefix_icon=ft.Icons.SEARCH,expand=True)
        archive = ft.Checkbox(label='Mostrar arquivados',value=archived)
        async def filter_projects(e):
            await self.guard(lambda:self.show_workspace(search.value or '',bool(archive.value)))
        search.on_submit = filter_projects
        archive.on_change = filter_projects
        cards = []
        for project in projects:
            if query and query.lower() not in (project['name']+' '+project['description']).lower():
                continue
            card = panel(ft.Column([
                ft.Row([ft.Container(width=12,height=12,bgcolor=project['color'],border_radius=6),user_text(project['name'],20,weight=ft.FontWeight.W_600)],wrap=True),
                user_text(project['description'],color=MUTED) if project['description'] else text('Seu próximo caminho de exploração começa aqui.',color=MUTED),
                text(f"{len(project['pipeline'])} etapas · " + self.tr({'owner':'Proprietário','editor':'Editor','viewer':'Leitor'}[project['role']]),12,MUTED),
                ft.Row([ft.TextButton('Abrir projeto',icon=ft.Icons.ARROW_FORWARD,on_click=self.event(self.open_project,project['id'])),
                    ft.TextButton('Histórico',icon=ft.Icons.HISTORY,on_click=self.event(self.history_dialog,project['id'])),
                    ft.TextButton('Exportar',icon=ft.Icons.DOWNLOAD,on_click=self.event(self.export_project_dialog,project['id'])),
                    *([ft.TextButton('Excluir projeto',icon=ft.Icons.DELETE_OUTLINE,
                        style=ft.ButtonStyle(color='#B91C1C'),on_click=self.event(self.delete_project_dialog,project['id']))]
                      if project['role']=='owner' else [])],wrap=True,alignment=ft.MainAxisAlignment.SPACE_BETWEEN)
                ],spacing=18),col={'sm':12,'md':6,'xl':4})
            cards.append(card)
        invitation_controls=[]
        for invitation in invites:
            invitation_controls.append(panel(ft.Row([ft.Icon(ft.Icons.MAIL_OUTLINE,color=TEAL),
                ft.Column([text('Convite para '+invitation['name'],weight=ft.FontWeight.W_600),
                           text('Por '+invitation['owner_name']+' · '+self.tr('Editor' if invitation['role']=='editor' else 'Leitor'),12,MUTED)],expand=True),
                button('Aceitar',self.event(self.accept_invite,invitation['id'],True,invitation['role']),True),
                button('Recusar',self.event(self.accept_invite,invitation['id'],False,invitation['role']))])))
        body=ft.Column([*invitation_controls,ft.Row([search,button('Buscar',filter_projects),archive]),
            ft.ResponsiveRow(cards,spacing=18,run_spacing=18) if cards else panel(ft.Column([
                ft.Icon(ft.Icons.FOLDER_OPEN,size=48,color=TEAL),text('Um espaço para sua próxima descoberta.',22),
                text('Crie um projeto ou aceite um convite para começar.',color=MUTED)]))],expand=True,scroll=ft.ScrollMode.AUTO,spacing=20)
        from .action_help import action_row
        self.render(self.shell('Seu workspace','Explore, organize e compartilhe suas pesquisas.',body,
            [action_row(self,'import',button('Importar projeto',self.event(self.import_project_dialog),icon=ft.Icons.UPLOAD_FILE)),
             action_row(self,'new',button('Novo projeto',self.event(self.project_dialog),True,ft.Icons.ADD))]))

    async def accept_invite(self,project_id,accepted,expected_role=None):
        try:
            await self.call(self.store.accept_invitation,self.token,project_id,accepted,expected_role)
        except ValueError:
            await self.show_workspace()
            raise
        await self.show_workspace()

    async def project_dialog(self,project=None):
        name = ft.TextField(label='Nome do projeto',value=project['name'] if project else '')
        description = ft.TextField(label='Descrição',value=project['description'] if project else '',multiline=True,min_lines=2,max_lines=4)
        from .color_palette import ColorPalette
        color = ColorPalette(self.page,project['color'] if project else None)
        token=self.token
        selection={'plan':None}
        async def choose_parent(path):
            plan=await self.call(self.store.project_destination,token,name.value or '',path)
            if self.token!=token:return
            async def accept(e=None):
                if self.token!=token:return
                selection['plan']=plan
                directory.value=plan['directory'];directory.error=None;self.page.update()
                if plan['replace']:close_dialog(self.page,confirmation)
            if plan['replace']:
                confirmation=ft.AlertDialog(modal=True,title=text('Substituir projeto existente?',24),
                    content=ft.Container(width=520,content=ft.Column([
                        text('Ao salvar, o projeto existente e todo o conteúdo desta pasta serão removidos permanentemente:'),
                        user_text(plan['directory'],12),
                        text('Entradas, resultados e histórico serão apagados. Esta ação não pode ser desfeita.',12,'#B91C1C')],tight=True,spacing=16)),
                    actions=[ft.TextButton('Manter pasta',on_click=lambda e:close_dialog(self.page,confirmation)),
                             button('Utilizar e substituir ao salvar',accept,True)])
                self.page.show_dialog(confirmation)
            else:
                await accept()
        directory,folder=self.folder_controls(project['directory'] if project else None,locked=bool(project),
            project_name=None if project else lambda:name.value, on_select=None if project else choose_parent)
        def renamed(e):
            if project:return
            selection['plan']=None
            directory.value='';directory.error='Escolha novamente a pasta principal após alterar o nome.'
            self.page.update()
        name.on_change=renamed
        async def save(e):
            async def action():
                if self.token!=token:return
                values=(name.value or '',description.value or '',color.value,project.get('tags',[]) if project else [])
                if project:
                    await self.call(self.store.update_project,self.token,project['id'],*values,archived=bool(project['archived']))
                    close_dialog(self.page,dialog)
                    await self.open_project(project['id'])
                else:
                    if not directory.value or not directory.value.strip():
                        directory.error='Informe ou escolha a pasta que receberá os arquivos e resultados.'
                        self.page.update()
                        raise ValueError('Informe a pasta do projeto.')
                    plan=selection['plan']
                    if not plan or plan['name']!=(name.value or '').strip():
                        raise ValueError('Escolha a pasta principal usando o botão de pasta após preencher o nome do projeto.')
                    created=await self.call(self.store.create_project_in_parent,token,*values,
                        parent=plan['parent'],confirmation=plan if plan['replace'] else None)
                    if self.token!=token:return
                    close_dialog(self.page,dialog)
                    await self.open_project(created['id'])
            await self.guard(action)
        dialog=ft.AlertDialog(title=text('Editar projeto' if project else 'Novo projeto',24),
            content=ft.Container(width=560,height=550,content=ft.Column([name,description,color.control,folder],
                scroll=ft.ScrollMode.AUTO,spacing=20,horizontal_alignment=ft.CrossAxisAlignment.STRETCH)),
            actions=[ft.TextButton('Cancelar',on_click=lambda e:close_dialog(self.page,dialog)),button('Salvar',save,True)])
        self.page.show_dialog(dialog)

    @with_processing('Abrindo projeto…','Carregando o pipeline e os arquivos do projeto. Aguarde um instante.')
    async def open_project(self,project_id,tab='Pipeline'):
        token = self.token
        project = await self.call(self.store.project,token,project_id)
        if self.token != token:
            return
        if getattr(self,'reuse_dialog_key',None)!=(token,project_id):self.dismiss_reuse_prompt()
        self.current = project
        self.base_pipeline=copy.deepcopy(project['pipeline'])
        self.dirty = False
        self.tab = tab
        stages = self.current['pipeline']
        if self.selected not in [s['id'] for s in stages]:
            self.selected = stages[0]['id'] if stages else None
        await self.draw_project()
        previous=getattr(self,'_project_poll',None)
        if previous and previous is not asyncio.current_task():previous.cancel()
        self._project_poll=self.page.run_task(self.poll_project,project_id)

    async def delete_project_dialog(self,project_id):
        token = self.token
        project=await self.call(self.store.project,token,project_id,'owner')
        runs=await self.call(self.store.list_runs,token,project_id)
        if self.token != token:
            return
        active=any(run['status'] in ('queued','running','awaiting_input') for run in runs)
        async def remove(e):
            if self.token != token:
                return
            async def action():
                await self.call(self.store.delete_project,token,project_id,expected_directory=project['directory'])
                if self.token != token:
                    return
                close_dialog(self.page,dialog)
                await self.show_workspace()
                self.notify('Projeto e pasta removidos permanentemente.')
            await progress.run(action)
        dialog=ft.AlertDialog(modal=True,title=text('Excluir projeto?',24),
            content=ft.Container(width=460,content=ft.Column([
                user_text(project['name'],18,weight=ft.FontWeight.W_600),
                text('O projeto será excluído do workspace e sua pasta será removida permanentemente, incluindo entradas, resultados e histórico.'),
                user_text(project['directory'],12),
                text('Há uma execução ativa. Conclua ou cancele a execução antes de excluir o projeto.' if active else
                     'Esta ação é permanente e não pode ser desfeita.',12,'#B91C1C')],tight=True,spacing=18)),
            actions=[ft.TextButton('Manter projeto',on_click=lambda e:close_dialog(self.page,dialog)),
                     button('Excluir projeto',remove,icon=ft.Icons.DELETE_OUTLINE,disabled=active)])
        progress=RemovalProgress(self,dialog,'Excluindo projeto…')
        self.page.show_dialog(dialog)

    @with_processing('Atualizando projeto…','Preparando as informações e os resultados exibidos. Aguarde um instante.')
    async def draw_project(self):
        project = self.current
        if project is None:
            return
        token, tab = self.token, self.tab
        # Recheck membership even when rendering an unsaved local draft.
        fresh = await self.call(self.store.project,token,project['id'])
        if self.token != token or self.current is not project or self.tab != tab:
            return
        current_response=fresh['revision']>=project['revision'] and fresh['updated']>=project['updated']
        if current_response and ((not self.dirty and not getattr(self,'editing_stage',False)) or fresh['role'] == 'viewer'):
            project.update(fresh)
            self.base_pipeline=copy.deepcopy(fresh['pipeline'])
            self.dirty = False
        else:
            project['role'] = fresh['role']
        if project['role'] != 'owner' and self.tab == 'Compartilhar':
            self.tab = 'Pipeline'
            tab = self.tab
        writable = project['role'] != 'viewer'
        nav = ft.Row([button(label,self.event(self.change_tab,label),primary=self.tab==label)
                      for label in ('Pipeline','Arquivos','Execuções')],wrap=True)
        if self.tab=='Pipeline':
            content = await self.pipeline_view(writable)
        elif self.tab=='Arquivos':
            content = await self.files_view(writable)
        elif self.tab=='Execuções':
            content = await self.runs_view(writable)
        else:
            content = await self.sharing_view()
        if self.token != token or self.current is not project or self.tab != tab:
            return
        from .action_help import action_row,action_group
        project_actions=[action_row(self,'workspace',button('Workspace',self.event(self.show_workspace),icon=ft.Icons.ARROW_BACK))]
        if writable:
            project_actions.append(action_row(self,'project',button('Projeto',self.event(self.project_dialog,project),icon=ft.Icons.TUNE)))
        if project['role']=='owner':
            project_actions.append(action_row(self,'share',button('Compartilhar',self.event(self.change_tab,'Compartilhar'),icon=ft.Icons.PERSON_ADD_ALT_1)))
        project_actions.append(action_row(self,'history',button('Histórico',self.event(self.history_dialog,project['id']),icon=ft.Icons.HISTORY)))
        execution_actions=[action_row(self,key,button(label,self.event(handler,*args),icon=icon,disabled=not writable))
            for key,label,handler,args,icon in [
                ('save','Salvar',self.save_draft,(),ft.Icons.SAVE_OUTLINED),
                ('templates','Modelos prontos',self.preset_dialog,(),ft.Icons.DASHBOARD_OUTLINED),
                ('run','Executar pipeline',self.run_pipeline,(),ft.Icons.PLAY_ARROW),
                ('selected','Executar até seleção',self.run_pipeline,(True,),ft.Icons.PLAY_CIRCLE_OUTLINE)]]
        groups=[action_group(self,'project','Navegação do projeto',ft.Icons.FOLDER_OPEN,project_actions),
            action_group(self,'execution','Pipeline e execução',ft.Icons.PLAY_CIRCLE_OUTLINE,execution_actions)]
        progress = getattr(self,'run_progress',None)
        extras = [progress.badge] if progress and progress.valid() else []
        self.render(self.shell(project['name'],project['description'] or 'Construa o pipeline que sua pesquisa precisa.',
            ft.Column([nav,*extras,content],expand=True,spacing=10),groups=groups,
            sidebar_content=self.flow_editor.library if self.tab=='Pipeline' else None))
        if self.tab=='Execuções':
            for result in list(self.stage_results.values()):
                if result.expanded and not result.loaded and result.valid():await result.load()

    async def change_tab(self,tab):
        if self.tab=='Pipeline':
            if self.dirty:
                await self.save_draft()
        self.tab = tab
        await self.draw_project()

    def selected_stage(self):
        return next((s for s in self.current['pipeline'] if s['id']==self.selected),None)

    async def pipeline_view(self,writable):
        from .flow_canvas import FlowCanvas
        self.project_assets=await self.call(self.store.assets,self.token,self.current['id'])
        runs=await self.call(self.store.list_runs,self.token,self.current['id'])
        self.automatic_labels={}
        for run in runs:
            for item in run['stages']:
                if item['status']=='succeeded' and item['id'] not in self.automatic_labels:
                    from biomolexplorer.artifact_choices import choices
                    options=choices(item)
                    self.automatic_labels[item['id']]='compounds.csv' if any(p.endswith('/compounds.csv') or p=='compounds.csv' for p in options) else next(iter(sorted(options)),item['name']).rsplit('/',1)[-1]
        self.flow_editor = FlowCanvas(self,writable)
        self.page.on_keyboard_event=self.flow_editor.keyboard
        return self.flow_editor.build(include_library=False)

    async def open_stage_dialog(self,stage_id,waiting_run_id=None,configuration=None):
        request=getattr(self,'_dialog_request',0)+1
        self._dialog_request=request
        if not hasattr(self,'_dialog_lock'):self._dialog_lock=asyncio.Lock()
        token,project=self.token,self.current
        async with processing(self,'Carregando configurações…',
                'Lendo informações do bloco. Aguarde um instante.',modal=False):
            async with self._dialog_lock:
                if self._dialog_request!=request or self.token!=token or self.current is not project:return
                await self._open_stage_dialog(stage_id,waiting_run_id,configuration,request)

    def _track_configuration_dialog(self,dialog,run_id,stage_id):
        key=(run_id,stage_id)
        dialog._configuration_key=key
        if not hasattr(self,'_configuration_dialogs'):self._configuration_dialogs={}
        self._configuration_dialogs.setdefault(key,[]).append(dialog)

    def _finish_configuration_dialogs(self,run_id,stage_id):
        """Finish the whole input/curation/settings flow, including older routes."""
        key=(run_id,stage_id)
        dialogs=getattr(self,'_configuration_dialogs',{}).pop(key,[])
        if getattr(self,'selection_dialog_key',None)==key:
            selection=getattr(self,'selection_dialog',None)
            if selection is not None and not any(d is selection for d in dialogs):dialogs.append(selection)
            self.selection_dialog=None;self.selection_dialog_key=None
        settings=getattr(self,'settings_dialog',None)
        if settings is not None and getattr(settings,'_configuration_key',None)==key:
            if not any(d is settings for d in dialogs):dialogs.append(settings)
            self.settings_dialog=None;self.settings_run_id=None;self.settings_body=None
        for dialog in dialogs:
            invalidate=getattr(dialog,'_configuration_invalidate',None)
            if invalidate:invalidate()
            for action in dialog.actions:action.disabled=True
            close_dialog(self.page,dialog)
        getattr(self,'pending_file_selections',{}).pop(key,None)
        self.editing_stage=any(getattr(self,name,None) is not None and getattr(self,name).open
            for name in ('selection_dialog','settings_dialog'))
        if not self.editing_stage:self._dialog_request=getattr(self,'_dialog_request',0)+1
        progress=getattr(self,'run_progress',None)
        if run_id and progress and progress.run['id']==run_id:
            progress.is_open=False;progress.closing=True;progress.reopen_requested=False
            close_dialog(self.page,progress.dialog)

    async def _open_stage_dialog(self,stage_id,waiting_run_id=None,configuration=None,request=None):
        from .guided import GuidedForm
        from biomolexplorer.pipeline import validate_pipeline
        token, project = self.token, self.current
        fresh=await self.call(self.store.project,token,project['id'])
        if self.token != token or self.current is not project or self._dialog_request!=request:
            return
        project['role']=fresh['role']
        self.selected=stage_id
        original=next((s for s in project['pipeline'] if s['id']==stage_id),None)
        if not original:
            return
        draft=copy.deepcopy(original)
        writable=self.current['role']!='viewer'
        assets=await self.call(self.store.assets,token,project['id'])
        if self.token != token or self.current is not project or self._dialog_request!=request:
            return
        runs=await self.call(self.store.list_runs,token,project['id'])
        if self.token != token or self.current is not project or self._dialog_request!=request:
            return
        if waiting_run_id:
            waiting_run=next((r for r in runs if r['id']==waiting_run_id and r['status']=='awaiting_input'),None)
            pending=next((s for s in waiting_run['stages'] if s['id']==stage_id and s['status']=='awaiting_input'),None) if waiting_run else None
            if pending is None:raise ValueError('Este bloco não está aguardando configuração.')
            draft=copy.deepcopy(pending['configuration'])
            # Selectors must reflect this execution's freshly produced files.
            runs=[waiting_run]+[r for r in runs if r['id']!=waiting_run_id]
        if configuration is not None:
            if configuration.get('id')!=stage_id or configuration.get('operation')!=original['operation']:
                raise ValueError('A seleção de arquivos pertence a outro bloco.')
            draft=copy.deepcopy(configuration)
        draft['parameters'].pop('verbose',None)
        from biomolexplorer.templates import quiet_template
        draft['templates']={name:quiet_template(name,text) for name,text in draft.get('templates',{}).items()}
        self.artifact_choices={}
        self.artifact_items={}
        for run in runs:
            for item in run['stages']:
                if item['id'] in self.artifact_choices or item['status']!='succeeded':continue
                from biomolexplorer.artifact_choices import choices
                self.artifact_choices[item['id']]=choices(item)
                self.artifact_items[item['id']]=dict(item,run_id=run['id'])
        for source in self.current['pipeline']:
            if source['operation']=='import_results' or source.get('provided_results'):
                params=source.get('provided_results') or source['parameters']
                self.artifact_choices.setdefault(source['id'],set()).update(a['name'] for a in assets if a['id'] in params.get('asset_ids',[]))
            elif source['operation'] in ('retrieve_compounds','retrieve_pubchem','retrieve_zinc'):
                self.artifact_choices.setdefault(source['id'],set()).add('compounds.csv')
        self.redocking_records={}
        if draft['operation']=='retrieve_pubchem':
            from biomolexplorer.pubchem_retrieval import reference_choices
            references={};seen=set()
            self.pubchem_references=references
            for run in runs:
                for item in run['stages']:
                    if self._dialog_request!=request or self.token!=token or self.current is not project:return
                    if item['status']!='succeeded' or item['id'] in seen:continue
                    seen.add(item['id']);records={}
                    for filename in item.get('artifacts',[]):
                        path=self.store.scoped_path(project['id'],filename)
                        if path.suffix.lower()=='.csv':records[str(path)]=await self.call(reference_choices,[path])
                    references['stage:'+item['id']]=records
            for asset in assets:
                if asset['kind'] not in ('compounds','other'):continue
                path=await self.call(self.store.asset_path,token,project['id'],asset['id'])
                references['asset:'+asset['id']]={str(path):await self.call(reference_choices,[path])}
            for source in project['pipeline']:
                if source['operation']=='import_results' or source.get('provided_results'):
                    params=source.get('provided_results') or source['parameters']
                    references['stage:'+source['id']]={p:r for a in params.get('asset_ids',[])
                        for p,r in references.get('asset:'+a,{}).items()}
            direct=draft['parameters'].get('base_input_path')
            if direct:
                folder=self.store.scoped_path(project['id'],direct)
                references['configured-path']={str(p):await self.call(reference_choices,[p]) for p in folder.glob('*.csv')}
            if self.token!=token or self.current is not project or self._dialog_request!=request:return
        if draft['operation'] in ('docking_vina','docking_dock6'):
            from .compound_catalog import CompoundCatalog
            paths={};seen=set()
            for run in runs:
                for item in run['stages']:
                    if self._dialog_request!=request or self.token!=token or self.current is not project:return
                    if item['status']!='succeeded' or item['id'] in seen:continue
                    seen.add(item['id'])
                    paths['stage:'+item['id']]=[self.store.scoped_path(project['id'],p)
                        for p in item.get('artifacts',[]) if Path(p).suffix.lower()=='.csv']
            for asset in assets:
                if asset['kind'] not in ('compounds','other'):continue
                path=await self.call(self.store.asset_path,token,project['id'],asset['id'])
                paths['asset:'+asset['id']]=[path]
            for source in project['pipeline']:
                if source['operation']=='import_results' or source.get('provided_results'):
                    params=source.get('provided_results') or source['parameters']
                    paths['stage:'+source['id']]=[p for a in params.get('asset_ids',[])
                        for p in paths.get('asset:'+a,[])]
            direct=draft['parameters'].get('base_selected_mols') or (draft['parameters'].get('base_vina_path') if draft['operation']=='docking_dock6' else None)
            if direct:
                folder=self.store.scoped_path(project['id'],direct)
                paths['configured-path']=await self.call(lambda:list(folder.glob('*.csv')))
            catalog=CompoundCatalog(self,paths,lambda:self.token==token and self.current is project
                and self._dialog_request==request)
            self.docking_catalog=catalog
            self.docking_compounds=catalog.references
            from biomolexplorer.bindings import sources
            group=draft.get('bindings',{}).get('base_selected_mols') or (draft.get('bindings',{}).get('base_vina_path') if draft['operation']=='docking_dock6' else None)
            origins=[('stage:'+r['stage']) if 'stage' in r else ('asset:'+r['asset'])
                for r in sources(group) if 'stage' in r or 'asset' in r] if group else []
            if direct:origins.append('configured-path')
            for origin in dict.fromkeys(origins):await catalog.load(origin)
        if draft['operation']=='redocking':
            from biomolexplorer.redocking_config import metadata_records
            from biomolexplorer.result_files import ResultFiles
            self.pdb_review_states={}
            seen=set()
            for run in runs:
                if any(item.get('operation')=='retrieve_structures' and item['id'] not in seen and item['status']=='succeeded' for item in run['stages']):
                    self.pdb_review_states.update(await self.call(ResultFiles(self.store).pdb_review_statuses,token,project['id'],run['id']))
                for item in run['stages']:
                    if self._dialog_request!=request or self.token!=token or self.current is not project:return
                    if item['status']!='succeeded' or item['id'] in seen:continue
                    seen.add(item['id'])
                    paths=[self.store.scoped_path(project['id'],p) for p in item.get('artifacts',[])]
                    self.redocking_records['stage:'+item['id']]=await self.call(metadata_records,paths)
            for asset in assets:
                path=await self.call(self.store.asset_path,token,project['id'],asset['id'])
                self.redocking_records['asset:'+asset['id']]=await self.call(metadata_records,[path])
            for source in project['pipeline']:
                if source['operation']=='import_results' or source.get('provided_results'):
                    params=source.get('provided_results') or source['parameters']
                    self.redocking_records['stage:'+source['id']]=[r for a in params.get('asset_ids',[]) for r in self.redocking_records.get('asset:'+a,[])]
            direct=draft['parameters'].get('base_input_path')
            if direct:
                folder=self.store.scoped_path(project['id'],direct)
                self.redocking_records['configured-path']=await self.call(metadata_records,list(folder.rglob('pdb_codes.csv')))
            if self.token!=token or self.current is not project or self._dialog_request!=request:return
        if self._dialog_request!=request or self.token!=token or self.current is not project:return
        if waiting_run_id in getattr(self,'_cancelled_run_ids',set()):return
        if not any(s['id']==stage_id for s in project['pipeline']):return
        mode='visual'
        reader=None
        from .settings_layout import dialog_size
        body=ft.Container(width=dialog_size(self.page)['width'],padding=16,bgcolor=BG,border_radius=16)
        body._settings_auto_size=True
        self.settings_body=body
        async def render_mode(new_mode):
            nonlocal draft,mode,reader
            if reader:
                draft=reader()
            mode=new_mode
            if mode=='visual':
                form=GuidedForm(self,draft,assets,writable)
                reader=form.read
                body.content=form.layout()
            else:
                title=ft.TextField(label='Nome do bloco',value=draft['name'],disabled=not writable)
                parameters=ft.TextField(label='Parâmetros completos (JSON)',value=json.dumps(draft['parameters'],indent=2,ensure_ascii=False),multiline=True,min_lines=8,disabled=not writable)
                bindings=ft.TextField(label='Conexões avançadas (JSON)',value=json.dumps(draft['bindings'],indent=2),multiline=True,min_lines=3,disabled=not writable)
                templates={name:ft.TextField(value=draft.get('templates',{}).get(name,(RESOURCE_ROOT/name).read_text()),multiline=True,min_lines=8,text_style=ft.TextStyle(font_family='monospace'),disabled=not writable) for name in template_names(draft['operation'])}
                body.content=ft.Column([title,parameters,bindings,*[ft.ExpansionTile(title=text(name),controls=[field]) for name,field in templates.items()]],scroll=ft.ScrollMode.AUTO,spacing=24,horizontal_alignment=ft.CrossAxisAlignment.STRETCH)
                def read_advanced():
                    result=copy.deepcopy(draft); result['name']=title.value or draft['name']; result['parameters']=json.loads(parameters.value); result['parameters'].pop('verbose',None); result['bindings']=json.loads(bindings.value)
                    result['templates']={name:field.value for name,field in templates.items() if field.value!=(RESOURCE_ROOT/name).read_text()}
                    return result
                reader=read_advanced
        async def visual(e):
            await self.guard(lambda:render_mode('visual')); self.page.update()
        async def advanced(e):
            await self.guard(lambda:render_mode('advanced')); self.page.update()
        applied=False
        async def apply(e):
            nonlocal applied
            if applied or self.token!=token or self.current is not project:return
            async with processing(self,'Aplicando configuração…',
                    'Validando e salvando as opções do bloco. Aguarde um instante.'):
                await self.call(self.store.project,self.token,self.current['id'],'editor')
                updated=reader()
                if updated['operation'] in ('redocking','prepare_structures','docking_vina','docking_dock6'):
                    from biomolexplorer.pipeline import execution_stage
                    updated=execution_stage(updated,'curated')
                from biomolexplorer.bindings import sources as binding_sources
                from biomolexplorer.flow import reset_input_configuration
                if not updated.get('provided_results'):
                    for field,group in original.get('bindings',{}).items():
                        before={json.dumps(ref,sort_keys=True) for ref in binding_sources(group)}
                        after={json.dumps(ref,sort_keys=True) for ref in binding_sources(updated.get('bindings',{}).get(field))}
                        # Redocking's reader has already reconciled its selected
                        # files, curated pairs and preparation options. Refining
                        # auto bindings into artifact references must not erase
                        # that plan immediately before resume/save.
                        if before-after and updated['operation']!='redocking':
                            reset_input_configuration(updated,field)
                if updated.get('provided_results'):
                    from biomolexplorer.input_validation import validate_bundle
                    provided=updated['provided_results']
                    paths=[await self.call(self.store.asset_path,self.token,self.current['id'],a) for a in provided['asset_ids']]
                    await self.call(validate_bundle,paths,provided['kind'],updated['operation'])
                if updated['operation']!='import_results' and not updated.get('provided_results'):
                    from biomolexplorer.operations import OPERATIONS, validate_operation
                    spec=OPERATIONS[updated['operation']]
                    parameters={key:'pendente' for key in spec.required}
                    parameters.update(updated['parameters'])
                    parameters.update({key:'pending-input' for key in updated.get('bindings',{})})
                    validate_operation(updated['operation'],parameters,
                        defer_redocking_selection=updated['operation']=='redocking',
                        defer_pubchem_selection=bool(updated.get('bindings',{})))
                candidate=[updated if s['id']==stage_id else s for s in self.current['pipeline']]
                validate_pipeline(candidate)
                if waiting_run_id:
                    resumed=await self.call(self.service.resume,token,waiting_run_id,updated)
                    applied=True
                    self.editing_stage=False
                    close_dialog(self.page,dialog)
                    self.current_run=resumed['id']
                    if getattr(self,'run_progress',None):self.run_progress.update(resumed)
                if hasattr(self,'flow_editor'):self.flow_editor.remember()
                self.current['pipeline']=candidate; self.dirty=True
                await self.save_draft(silent=True)
                applied=True
                self._finish_configuration_dialogs(waiting_run_id,stage_id)
                if waiting_run_id:
                    self.tab='Execuções'
                    await self.draw_project()
                    if not self.polling or self.polling.done():self.polling=self.page.run_task(self.poll_runs)
                else:
                    self.editing_stage=False
                    close_dialog(self.page,dialog)
                    if hasattr(self,'flow_editor'):self.flow_editor.refresh()
        await render_mode(mode)
        async def cancel_settings(e):
            if waiting_run_id and writable:
                await self.cancel_run(waiting_run_id)
            else:
                close_dialog(self.page,dialog)
        dialog=ft.AlertDialog(title=ft.Column([text(TITLES[original['operation']][2].upper(),11,TEAL,ft.FontWeight.W_600),stage_control(user_text(original['name'],24,weight=ft.FontWeight.W_700),original),text(TITLES[original['operation']][1],13,MUTED),ft.Row([
            ft.TextButton('Visual · formulário guiado',on_click=visual),ft.TextButton('Avançado · parâmetros e scripts',on_click=advanced)],wrap=True)],spacing=12),content=body,
            actions=[ft.TextButton('Cancelar execução' if waiting_run_id and writable else 'Cancelar' if writable else 'Fechar',on_click=self.event_action(cancel_settings) if waiting_run_id and writable else lambda e:close_dialog(self.page,dialog)),button('Aplicar e continuar pipeline' if waiting_run_id else 'Aplicar configuração',self.event_action(apply),True,disabled=not writable)])
        self._track_configuration_dialog(dialog,waiting_run_id,stage_id)
        self.settings_dialog=dialog
        self.settings_run_id=waiting_run_id
        self.editing_stage=True
        def dismissed(e):
            if getattr(self,'settings_dialog',None) is dialog:
                self.settings_dialog=None;self.settings_run_id=None
            if self._dialog_request==request:
                self.settings_body=None
                self.editing_stage=False
                self._dialog_request+=1
        dialog.on_dismiss=dismissed
        self.page.show_dialog(dialog)

    def event_action(self,handler):
        running=False
        async def wrapped(e):
            nonlocal running
            if running:return
            running=True
            control=getattr(e,'control',None)
            disabled=getattr(control,'disabled',False)
            try:
                if control is not None:
                    control.disabled=True
                    self.page.update()
                await self.guard(lambda:handler(e))
            finally:
                running=False
                if control is not None:
                    control.disabled=disabled
                    self.page.update()
        return wrapped

    def auto_bind(self,stage,previous):
        operation=stage['operation']
        wanted={
            'retrieve_pubchem': {'base_input_path':['retrieve_compounds','admet','graphs','import_results']},
            'expand_similar_compounds': {'base_input_path':['retrieve_compounds']},
            'admet': {'base_input_path':['graphs','retrieve_compounds','retrieve_pubchem','import_results']},
            'fingerprints': {'base_input_path':['admet','retrieve_compounds','retrieve_pubchem','import_results']},
            'similarity': {'base_input_path':['fingerprints','import_results']},
            'graphs': {'base_input_path':['retrieve_compounds','import_results'], 'similarity_path':['similarity']},
            'prepare_structures': {'base_input_path':['import_results','retrieve_structures']},
            'redocking': {'base_input_path':['retrieve_structures','import_results']},
            'docking_vina': {'base_input_path':['prepare_structures','redocking','import_results'], 'base_selected_mols':['admet','graphs','retrieve_compounds','retrieve_pubchem']},
            'docking_dock6': {'base_input_path':['prepare_structures','redocking','import_results'], 'base_selected_mols':['docking_vina','docking_dock6','admet','graphs','retrieve_compounds','retrieve_pubchem','import_results']},
            'consensus': {'base_input_path':['docking_vina'], 'base_vina_path':['docking_vina','import_results'], 'base_dock6_path':['docking_dock6','import_results']},
        }.get(operation,{})
        for field,candidates in wanted.items():
            def compatible(s):
                if s['operation'] not in candidates:
                    return False
                if s['operation']!='import_results':
                    return True
                from biomolexplorer.flow import output_types
                kinds=output_types(s)
                if field=='base_selected_mols' or operation in ('retrieve_pubchem','admet','fingerprints','graphs'):
                    return bool(kinds & {'compounds','vina','dock6'}) if field=='base_selected_mols' and operation in ('docking_vina','docking_dock6') else 'compounds' in kinds
                if field=='base_vina_path':
                    return 'vina' in kinds
                if field=='base_dock6_path':
                    return 'dock6' in kinds
                if operation in ('prepare_structures','redocking','docking_vina','docking_dock6'):
                    return bool(kinds & {'structures','prepared_structures'})
                return True
            source=next((s for s in reversed(previous) if compatible(s)),None)
            if source:
                stage['bindings'][field]={'stage':source['id'],'selector':'auto'}
                if field=='base_input_path' and 'target' in source['parameters'] and 'target' in stage['parameters']:
                    stage['parameters']['target']=source['parameters']['target']
                if operation=='expand_similar_compounds':
                    stage['parameters']['search_term']=source['parameters']['search_term']
        if operation=='docking_dock6' and self.service.dock6_path:
            stage['parameters']['dock6_app_path']=str(self.service.dock6_path)

    async def preset_dialog(self):
        choice=ft.Dropdown(label='Modelo de pipeline',value=next(iter(PRESETS)),options=[ft.DropdownOption(key=k,text=k) for k in PRESETS],width=440)
        async def apply(e):
            stages=[]
            for operation in PRESETS[choice.value]:
                stage=new_stage(operation)
                stage['name']=self.tr(stage['name'])
                if operation=='import_results' and choice.value=='Meus PDBs → preparação':
                    stage['parameters']['kind']='structures'
                self.auto_bind(stage,stages)
                if choice.value=='ChEMBL + PubChem → ADMET' and operation=='admet':
                    stage['bindings']['base_input_path']={'sources':[{'stage':s['id'],'selector':'auto'}
                        for s in stages if s['operation'] in ('retrieve_compounds','retrieve_pubchem')]}
                    stage['input_processing']='merge'
                stages.append(stage)
            from biomolexplorer.flow import arrange
            self.current['pipeline']+=stages
            arrange(self.current['pipeline'])
            self.selected=stages[0]['id']
            self.dirty=True
            self.page.pop_dialog()
            await self.draw_project()
        self.page.show_dialog(ft.AlertDialog(title=text('Comece com um caminho pronto',24),content=ft.Column([
            text('As etapas serão adicionadas ao pipeline atual.\nVocê poderá mudar tudo depois.',color=MUTED),choice],tight=True),
            actions=[ft.TextButton('Cancelar',on_click=lambda e:self.page.pop_dialog()),button('Adicionar modelo',self.event_action(apply),True)]))

    @with_processing('Salvando pipeline…','Gravando as configurações e o histórico do projeto. Aguarde um instante.')
    async def save_draft(self,silent=False):
        async with self._save_lock:
            await self._save_draft(silent)

    async def _save_draft(self,silent=False):
        token, project = self.token, self.current
        if project is None or project['role']=='viewer':
            return
        if not self.dirty:return
        sent=copy.deepcopy(project['pipeline'])
        saved=await self.call(self.store.save_pipeline,token,project['id'],sent,project['revision'],base_pipeline=self.base_pipeline)
        if self.token != token or self.current is not project:
            return
        previous_pipeline=project['pipeline']
        if project['pipeline']==sent:
            # Dialogs and scheduled work hold this project identity. Preserve
            # unchanged stage objects too, so canvas gestures stay attached.
            if saved['pipeline']==sent:
                saved=dict(saved,pipeline=project['pipeline'])
            project.update(saved);self.dirty=False
        else:
            from biomolexplorer.project_merge import merge_pipeline
            pending=merge_pipeline(sent,project['pipeline'],saved['pipeline'])
            project.update({key:value for key,value in saved.items() if key!='pipeline'})
            if pending!=project['pipeline']:project['pipeline']=pending
            self.schedule_save()
        self.base_pipeline=copy.deepcopy(saved['pipeline'])
        if self.tab=='Pipeline' and hasattr(self,'flow_editor') and project['pipeline'] is not previous_pipeline:
            self.flow_editor.refresh()
        if not silent:self.notify('Pipeline salvo.')

    @with_processing('Preparando execução…','Validando o pipeline e organizando os dados de entrada. Aguarde um instante.')
    async def run_pipeline(self,selected=False):
        if getattr(self,'submitting',False):
            return
        self.submitting = True
        try:
            await self._submit_pipeline(selected)
        finally:
            self.submitting = False

    def dismiss_reuse_prompt(self):
        dialog=getattr(self,'reuse_dialog',None)
        self.reuse_dialog=None
        self.reuse_dialog_key=None
        if dialog is not None and dialog.open:close_dialog(self.page,dialog)

    async def _submit_pipeline(self,selected=False,reuse_results=None,selected_stage_id=None):
        token, project_id = self.token, self.current['id']
        active_dialog=getattr(self,'reuse_dialog',None)
        if reuse_results is None and active_dialog is not None:
            if getattr(self,'reuse_dialog_key',None)==(token,project_id) and active_dialog.open:return
            self.dismiss_reuse_prompt()
        selected_stage = selected_stage_id or self.selected
        await self.save_draft()
        if self.token != token or not self.current or self.current['id'] != project_id:
            return
        if reuse_results is None and await self.call(self.service.existing_results,token,project_id):
            if self.token!=token or not self.current or self.current['id']!=project_id:return
            async def choose(reuse):
                if getattr(self,'reuse_dialog',None) is not dialog:return
                if self.token!=token or not self.current or self.current['id']!=project_id:
                    self.dismiss_reuse_prompt()
                    return
                if self.submitting:return
                self.dismiss_reuse_prompt()
                self.submitting=True
                try:await self._submit_pipeline(selected,reuse,selected_stage)
                finally:self.submitting=False
            def cancel(e):
                if getattr(self,'reuse_dialog',None) is dialog:self.dismiss_reuse_prompt()
            def dismissed(e):
                if getattr(self,'reuse_dialog',None) is dialog:
                    self.reuse_dialog=None
                    self.reuse_dialog_key=None
            dialog=ft.AlertDialog(modal=True,title=text('Reaproveitar dados e resultados?',24),
                content=ft.Container(width=560,content=text(
                    'Este projeto já possui resultados salvos. Deseja reaproveitá-los? '
                    'Etapas compatíveis com arquivos íntegros serão consideradas concluídas. '
                    'No redocking, receptores e pares compatíveis de execuções concluídas serão reaproveitados; '
                    'somente os pares novos ou alterados serão calculados e incorporados à tabela.')),
                actions=[ft.TextButton('Cancelar',on_click=cancel),
                    button('Não, executar novamente',self.event_action(lambda e:choose(False))),
                    button('Sim, reaproveitar',self.event_action(lambda e:choose(True)),True)],on_dismiss=dismissed)
            self.reuse_dialog=dialog
            self.reuse_dialog_key=(token,project_id)
            self.page.show_dialog(dialog)
            return
        if self.token!=token or not self.current or self.current['id']!=project_id:return
        run=await self.call(self.service.submit,token,project_id,[selected_stage] if selected else None,
                            reuse_results=reuse_results if reuse_results is not None else not getattr(self,'force_execution',False))
        if self.token != token or not self.current or self.current['id'] != project_id:
            return
        self.current_run=run['id']
        from .run_progress import RunProgress
        self.run_progress=RunProgress(self,run)
        self.run_progress.show()
        self.tab='Execuções'
        await self.draw_project()
        if not self.polling or self.polling.done():
            self.polling=self.page.run_task(self.poll_runs)

    async def files_view(self,writable):
        token,project_id=self.token,self.current['id']
        assets=await self.call(self.store.assets,token,project_id)
        kind=ft.Dropdown(label='Tipo dos arquivos enviados',value='compounds',options=[ft.DropdownOption(key=k,text=v) for k,v in KINDS.items()],width=320,disabled=not writable)
        async def upload(e):
            async def action():
                files=await pick_files(self,allow_multiple=True)
                if not files:
                    return
                project_id=self.current['id']
                for file in files:
                    if file.path and not self.page.web:
                        await self.call(self.store.import_local_file,self.token,project_id,file.path,kind.value)
                    else:
                        if file.name in self.uploads:
                            raise ValueError('Aguarde o envio anterior do arquivo '+file.name+'.')
                        ticket=await self.call(self.store.prepare_upload,self.token,project_id,file.name,kind.value)
                        self.uploads[file.name]=(ticket,project_id)
                        await self.picker.upload([ft.FilePickerUploadFile(name=file.name,id=file.id,
                            upload_url=self.page.get_upload_url(ticket,600))])
                if not self.page.web:
                    await self.draw_project()
                else:
                    self.notify('Envio iniciado. Os arquivos aparecerão ao concluir.')
            await self.guard(action)
        from .file_table import FileTable
        from biomolexplorer.result_files import ResultFiles
        active=any(r['status'] in ('queued','running','awaiting_input') for r in await self.call(self.store.list_runs,token,project_id))
        async def download(item):
            path=await self.call(self.store.asset_path,token,project_id,item['id'])
            data=await self.call(self.store.read_file,token,project_id,path)
            if self.token==token and self.current and self.current['id']==project_id:
                await save_file(self,file_name=item['name'],src_bytes=data)
        async def remove(item):await self.call(ResultFiles(self.store).remove_asset,token,project_id,item['id'])
        files=FileTable(self,project_id,[dict(a,description=KINDS[a['kind']]) for a in assets],writable and not active,
            download=download,remove=remove)
        return ft.Column([panel(ft.Column([text('Sua biblioteca de entradas',22,weight=ft.FontWeight.W_600),
            text('Use resultados próprios em qualquer etapa. Os arquivos ficam restritos aos membros deste projeto.',color=MUTED),
            ft.Row([kind,button('Enviar arquivos',upload,True,ft.Icons.UPLOAD_FILE,not writable)],wrap=True),
            text('Compostos: CSV com canonical_smiles/molecule_chembl_id ou smiles/name.\nPDBs próprios: envie as estruturas; configure ligante, resíduo e cadeia na preparação.\nReceptores preparados: inclua PDBQT, centers.csv e pdb_codes.csv para o protocolo Vina existente.',12,MUTED)],spacing=14)),
            panel(files.build())],expand=True,scroll=ft.ScrollMode.AUTO,spacing=12)

    async def upload_event(self,e):
        pending=self.uploads.get(e.file_name)
        if isinstance(pending,dict):
            future=pending['future']
            if future.done():return
            try:
                if e.error:raise ValueError('Falha no envio: '+e.error)
                if e.progress==1:
                    result=self.store.staging/pending['ticket'] if pending['archive'] else await self.call(self.store.finish_upload,pending['token'],pending['ticket'])
                    future.set_result(result)
            except Exception as exc:
                get_logger('frontend').exception('Falha no envio de arquivo')
                future.set_exception(exc)
            return
        async def action():
            if e.error:
                raise ValueError('Falha no envio: '+e.error)
            if e.progress==1 and e.file_name in self.uploads:
                ticket,project_id=self.uploads.pop(e.file_name)
                await self.call(self.store.finish_upload,self.token,ticket)
                if self.current and self.current['id']==project_id and self.tab=='Arquivos':
                    await self.draw_project()
        await self.guard(action)

    async def runs_view(self,writable):
        if not hasattr(self,'stage_results'):self.stage_results={}
        if not hasattr(self,'results_focus'):self.results_focus=None
        token = self.token
        project_id = self.current['id']
        runs=await self.call(self.store.list_runs,token,project_id)
        if self.token != token or not self.current or self.current['id'] != project_id:
            return ft.Column([])
        active=next((run for run in runs if run['status'] in ('queued','running','awaiting_input')),None)
        visible={(r['id'],s['id']) for r in runs for s in r['stages']}
        for key in list(self.stage_results):
            if key not in visible:self.stage_results.pop(key).close()
        if active and getattr(self,'run_progress',None) is None:
            from .run_progress import RunProgress
            self.run_progress=RunProgress(self,active)
        if active and (not self.polling or self.polling.done()):
            self.current_run=active['id']
            self.polling=self.page.run_task(self.poll_runs)
        progress=getattr(self,'run_progress',None)
        if progress and progress.valid():
            refreshed=next((run for run in runs if run['id']==progress.run['id']),None)
            if refreshed:progress.update(refreshed)
        focus_id=active['id'] if active else progress.run['id'] if progress and progress.valid() else runs[0]['id'] if runs else None
        content=[]
        history=[]
        for run in runs:
            stage_rows=[]
            for stage in run['stages']:
                from .stage_results import StageResults
                key=(run['id'],stage['id'])
                signature=(stage['status'],tuple(stage.get('artifacts',[])),writable,bool(active),
                    json.dumps(stage.get('artifact_manifest'),sort_keys=True),stage.get('curated_at'))
                result=self.stage_results.get(key)
                if result is None or result.signature!=signature:
                    expanded=result.expanded if result else False
                    if result:result.close()
                    result=StageResults(self,project_id,run['id'],stage,writable)
                    result.expanded=expanded
                    result.signature=signature;self.stage_results[key]=result
                if self.results_focus==key:result.expanded=True
                async def log(e,path=stage.get('log_path'),project_id=project_id):
                    await self.guard(lambda:self.open_stage_log(project_id,path))
                from biomolexplorer.molecule_quality import REPORT_NAME
                reports=[p for p in stage.get('artifacts',[]) if Path(p).name==REPORT_NAME]
                report_buttons=[ft.TextButton('Baixar relatório de exclusões'+(f' ({index+1})' if len(reports)>1 else ''),
                    icon=ft.Icons.DOWNLOAD,on_click=self.event(self.download_stage_report,project_id,path))
                    for index,path in enumerate(reports)]
                from .run_presentation import result_state
                display_status,stage_label=result_state(stage)
                if stage.get('excluded_records'):stage_label+=f' · {stage["excluded_records"]} registros excluídos'
                if stage.get('redocking_failures'):stage_label=self.tr(stage_label)+' · '+self.tr(f'{stage["redocking_failures"]} falhas na preparação de complexos')
                from .run_presentation import STATES,state_icon
                stage_details=[]
                from biomolexplorer.preparation_failures import REPORT_NAME as PREPARATION_REPORT
                failures=[p for p in stage.get('artifacts',[]) if Path(p).name==PREPARATION_REPORT]
                if failures:
                    path=min(failures,key=lambda p:len(Path(p).parts))
                    stage_details.append(ft.Row([
                        ft.TextButton('Compostos não utilizados',icon=ft.Icons.WARNING_AMBER_OUTLINED,
                            on_click=self.event(self.open_preparation_failures,project_id,run['id'],stage['id'],path)),
                        ft.TextButton('Baixar lista de compostos não utilizados',icon=ft.Icons.DOWNLOAD,
                            on_click=self.event(self.download_stage_report,project_id,path))],wrap=True))
                if writable and stage.get('artifacts'):
                    stage_details.append(ft.TextButton('Remover todos os resultados',icon=ft.Icons.DELETE_SWEEP_OUTLINED,
                        disabled=bool(active),on_click=self.event(self.remove_stage_results,project_id,run['id'],stage['id'])))
                if stage.get('selection_notice'):stage_details.append(text(stage['selection_notice'],color='#B45309'))
                if stage.get('error'):stage_details.append(text(readable_error(stage['error']),color='#B91C1C'))
                if report_buttons:stage_details.append(ft.Row(report_buttons,wrap=True))
                stage_details.append(result.root)
                if stage.get('log_path'):stage_details.append(ft.TextButton('Ver log',icon=ft.Icons.ARTICLE_OUTLINED,on_click=log))
                stage_rows.append(ft.ExpansionTile(title=stage_control(user_text(stage['name'],15,weight=ft.FontWeight.W_600),stage),
                    leading=state_icon(display_status),subtitle=text(stage_label,12,STATES[display_status][0]),
                    expanded=result.expanded,on_change=result.expand,maintain_state=True,
                    tile_padding=ft.Padding.symmetric(horizontal=8,vertical=4),
                    shape=ft.RoundedRectangleBorder(radius=12,side=ft.BorderSide(1,LINE)),
                    collapsed_shape=ft.RoundedRectangleBorder(radius=12,side=ft.BorderSide(1,LINE)),
                    controls=[ft.Container(content=ft.Column(stage_details,spacing=16),padding=16,bgcolor='#F8FAFC',border_radius=12)]))
            actions=[]
            if run['status'] in ('queued','running','awaiting_input') and writable:
                actions=[button('Cancelar execução',self.event(self.cancel_run,run['id']),icon=ft.Icons.STOP)]
                if run['status']=='awaiting_input':
                    pending=next(s for s in run['stages'] if s['status']=='awaiting_input')
                    actions.insert(0,button('Configurar '+pending['name'],self.event(self.configure_pending,run['id']),True,ft.Icons.EDIT))
            timestamp=datetime.fromtimestamp(run['created'],timezone.utc).strftime('%d/%m/%Y %H:%M UTC')
            label='Execução atual' if run['status'] in ('queued','running','awaiting_input') else 'Última execução' if run['id']==focus_id else 'Execução anterior'
            if run['status']=='succeeded':
                first=next((s for s in run['stages'] if s['status']=='succeeded' and s.get('artifacts')),None)
                if first:
                    async def explore(e,key=(run['id'],first['id'])):
                        self.results_focus=key;await self.guard(self.draw_project)
                    actions.append(button('Explorar resultados',explore,True,ft.Icons.INSIGHTS))
            from .run_presentation import outcome_overview
            details=[ft.Row([text(label,12,MUTED),text('·',12,MUTED),user_text(timestamp,12,MUTED)],spacing=4,wrap=True),outcome_overview(run,actions)]
            if run.get('error'):details.append(text(readable_error(run['error']),color='#B45309' if run['status']=='awaiting_input' else '#B91C1C'))
            details.extend([text('Resultados por etapa',14,INK,ft.FontWeight.W_600),*stage_rows])
            card=panel(ft.Column(details,spacing=16),padding=18)
            if run['id']==focus_id:
                content.append(card)
            else:
                history.append(card)
        if history:
            content.append(ft.ExpansionTile(title=text(f'Execuções anteriores ({len(history)})',16),
                subtitle=text('Histórico de resultados e diagnósticos. Estas execuções não fazem parte da execução atual.',12,MUTED),
                expanded=bool(self.results_focus and self.results_focus[0]!=focus_id),controls=history))
        if not content:
            content=[panel(ft.Column([text('Pronto para explorar?',24),text('Configure suas etapas e execute o pipeline. Os resultados e logs aparecerão aqui.',color=MUTED)]))]
        has_results=any(stage.get('artifacts') for run in runs for stage in run['stages']
                        if stage.get('operation')!='import_results' and not stage.get('provided'))
        tools=[button('Atualizar resultados',self.event(self.draw_project),icon=ft.Icons.REFRESH)]
        if writable:
            tools.append(button('Excluir todos os resultados',self.event(self.remove_all_results,project_id),
                icon=ft.Icons.DELETE_SWEEP_OUTLINED,disabled=bool(active) or not has_results))
        return ft.Column([ft.Row(tools,wrap=True),*content],expand=True,scroll=ft.ScrollMode.AUTO,spacing=16)

    async def remove_all_results(self,project_id):
        from biomolexplorer.result_files import ResultFiles
        token=self.token
        async def action():
            if self.token!=token:return
            await self.call(ResultFiles(self.store).remove_project_results,token,project_id)
            close_dialog(self.page,dialog)
            if self.current and self.current['id']==project_id:await self.draw_project()
        async def remove(e):await progress.run(action)
        dialog=ft.AlertDialog(modal=True,title=text('Excluir todos os resultados?',24),
            content=ft.Column([
                text('Todos os arquivos calculados de todas as execuções serão removidos.'),
                text('Arquivos enviados pelo usuário, configurações, logs e histórico serão preservados. A restauração continua disponível pelo histórico.',12),
                text('Esta ação afeta todo o projeto.',12,'#B91C1C')],tight=True,spacing=12),
            actions=[ft.TextButton('Manter resultados',on_click=lambda e:close_dialog(self.page,dialog)),
                button('Excluir todos',remove,icon=ft.Icons.DELETE_SWEEP_OUTLINED)])
        progress=RemovalProgress(self,dialog,'Excluindo todos os resultados…')
        self.page.show_dialog(dialog)

    async def remove_stage_results(self,project_id,run_id,stage_id):
        from biomolexplorer.result_files import ResultFiles
        token=self.token
        async def action():
            if self.token!=token:return
            await self.call(ResultFiles(self.store).remove_all,token,project_id,run_id,stage_id)
            close_dialog(self.page,dialog)
            if self.current and self.current['id']==project_id:await self.draw_project()
        async def remove(e):
            await progress.run(action)
        dialog=ft.AlertDialog(modal=True,title=text('Remover todos os resultados?',24),
            content=text('Os arquivos desta etapa serão removidos, incluindo referências reaproveitadas em outras execuções. O registro da execução e o log serão preservados. Você poderá restaurar os arquivos pelo histórico.'),
            actions=[ft.TextButton('Cancelar',on_click=lambda e:close_dialog(self.page,dialog)),
                button('Remover resultados',remove,icon=ft.Icons.DELETE_OUTLINE)])
        progress=RemovalProgress(self,dialog,'Excluindo resultados da etapa…')
        self.page.show_dialog(dialog)

    @with_processing('Carregando relatório de preparação…','Consultando os compostos e complexos não utilizados. Aguarde um instante.')
    async def open_preparation_failures(self,project_id,run_id,stage_id,path):
        from biomolexplorer.result_files import ResultFiles
        token=self.token
        rows=await self.call(ResultFiles(self.store).preparation_failures,token,project_id,run_id,stage_id,path)
        from .preparation_failures import show_preparation_failures
        if self.token==token and self.current and self.current['id']==project_id:show_preparation_failures(self,rows)

    @with_processing('Carregando log da etapa…','Lendo o registro de execução do bloco. Aguarde um instante.')
    async def open_stage_log(self,project_id,path):
        if not path:
            raise ValueError('O log estará disponível quando a etapa iniciar.')
        token = self.token
        data=await self.call(self.store.read_file,token,project_id,path,200000)
        if self.token != token or not self.current or self.current['id'] != project_id:
            return
        self.page.show_dialog(ft.AlertDialog(title=text('Log da etapa',24),content=ft.Container(width=750,height=450,
            content=ft.TextField(label='Conteúdo do log',value=data.decode('utf-8',errors='replace'),
                read_only=True,multiline=True,min_lines=18,max_lines=18,text_size=12,
                text_style=ft.TextStyle(font_family='monospace'))),
            actions=[ft.TextButton('Fechar log',on_click=lambda e:self.page.pop_dialog())]))

    @with_processing('Preparando download…','Lendo o arquivo para download. Aguarde um instante.')
    async def download_stage_report(self,project_id,path):
        token=self.token
        data=await self.call(self.store.read_file,token,project_id,path)
        if self.token==token and self.current and self.current['id']==project_id:
            await save_file(self,file_name=Path(path).name,src_bytes=data)

    @with_processing('Carregando compostos…','Preparando as tabelas de compostos disponíveis. Aguarde um instante.')
    async def open_compound_tables(self,project_id,run_id,stage_id):
        from biomolexplorer.compound_tables import CompoundTables
        from .compound_table import CompoundTableViewer
        token=self.token
        tables=await self.call(CompoundTables(self.store).tables,token,project_id,run_id,stage_id)
        if self.token!=token or not self.current or self.current['id']!=project_id:
            return
        if not tables:
            raise ValueError('Esta etapa não produziu tabelas de compostos para consulta.')
        viewer=CompoundTableViewer(self,project_id,run_id,stage_id,tables)
        self.compound_viewer=viewer
        width=max(280,min(1160,(self.page.width or 1440)-100))
        height=max(300,min(740,(self.page.height or 1000)-200))
        self.editing_stage=True
        def dismissed(e):
            if self.compound_viewer is viewer:
                viewer.preview_version+=1
                self.compound_viewer=None
                self.editing_stage=False
        dialog=ft.AlertDialog(title=text('Compostos recuperados',22),
            content=ft.Container(width=width,height=height,content=viewer.build()),on_dismiss=dismissed,
            actions=[ft.TextButton('Fechar tabela',on_click=lambda e:close_dialog(self.page,dialog))])
        self.page.show_dialog(dialog)
        await viewer.load()

    def pdb_view_url(self,project_id,path):
        from biomolexplorer.pdb_view import StructureViewers
        if not getattr(self,'structure_viewer',None):self.structure_viewer=StructureViewers(self.store)
        key=self.structure_viewer.issue(self.token,project_id,path,getattr(self,'language','pt'))
        return self.structure_viewer.url(key,web=self.page.web,page_url=self.page.url if self.page.web else None)

    def compound_view_url(self,project_id,smiles,name,token):
        from biomolexplorer.pdb_view import StructureViewers
        if not getattr(self,'structure_viewer',None):self.structure_viewer=StructureViewers(self.store)
        key=self.structure_viewer.issue_compound(token,project_id,smiles,name,getattr(self,'language','pt'))
        return self.structure_viewer.url(key,web=self.page.web,page_url=self.page.url if self.page.web else None)

    def docking_view_url(self,project_id,run_id,stage_id,kind,selection,token=None,name=None):
        """Issue a navigation URL for a button's client-side user gesture."""
        token=token or self.token
        if token!=self.token or not self.current or self.current['id']!=project_id:return None
        from biomolexplorer.pdb_view import StructureViewers
        if not getattr(self,'structure_viewer',None):self.structure_viewer=StructureViewers(self.store)
        key=self.structure_viewer.reserve_result(token,project_id,run_id,stage_id,kind,selection,
            getattr(self,'language','pt'),name=name)
        return self.structure_viewer.url(key,web=self.page.web,page_url=self.page.url if self.page.web else None)

    async def preview_docking(self,project_id,run_id,stage_id,kind,selection,token=None,document=False):
        token=token or self.token
        if token!=self.token:return
        from biomolexplorer.pdb_view import StructureViewers
        if not getattr(self,'structure_viewer',None):self.structure_viewer=StructureViewers(self.store)
        if document:
            key=await self.call(self.structure_viewer.issue_document,token,project_id,run_id,stage_id,selection,getattr(self,'language','pt'))
        else:
            key=await self.call(self.structure_viewer.issue_result,token,project_id,run_id,stage_id,kind,selection,getattr(self,'language','pt'))
        if token!=self.token or not self.current or self.current['id']!=project_id:return
        url=self.structure_viewer.url(key,web=self.page.web,page_url=self.page.url if self.page.web else None)
        await ft.UrlLauncher().launch_url(url,mode=ft.LaunchMode.EXTERNAL_APPLICATION,web_only_window_name='_blank')

    @with_processing('Preparando visualização…','Lendo os dados e preparando sua visualização. Aguarde um instante.')
    async def preview_artifact(self,project_id,path):
        token = self.token
        if Path(path).suffix.lower() in ('.pdb','.pdbqt','.mol2','.sdf'):
            url=await self.call(self.pdb_view_url,project_id,path)
            if self.token!=token or not self.current or self.current['id']!=project_id:return
            await ft.UrlLauncher().launch_url(url,mode=ft.LaunchMode.EXTERNAL_APPLICATION,web_only_window_name='_blank')
            return
        limit=MAX_GRAPH_BYTES if str(path).endswith(SUFFIX) else MAX_VIEW_BYTES
        data=await self.call(self.store.read_file,token,project_id,path,limit+1)
        if self.token != token or not self.current or self.current['id'] != project_id:
            return
        if len(data)>limit:
            raise ValueError(f'Visualização maior que o limite de {limit//(1024*1024)} MB. Baixe o arquivo para consultá-lo localmente.')
        width=max(320,min(1160,(self.page.width or 1440)-100))
        height=max(300,min(690,(self.page.height or 1000)-220))
        result=None
        if str(path).endswith(SUFFIX):
            from .results_viewer import ResultsViewer
            model=await self.call(load_view,data,max_bytes=limit)
            if self.token != token or not self.current or self.current['id'] != project_id:
                return
            # Raw uploaded JSON cannot assert a trusted run/artifact context.
            model.pop('dashboard_context',None)
            viewer=await self.call(ResultsViewer,self,project_id,model)
            if self.token!=token or not self.current or self.current['id']!=project_id:viewer.close();return
            self.result_viewer=viewer
            result=viewer
            content=ft.Column([self.result_viewer.build()],scroll=ft.ScrollMode.AUTO)
            await viewer.initialize()
            title=model['title']
        elif Path(path).suffix.lower() in ('.png','.jpg','.jpeg'):
            title=Path(path).name
            from .flow_canvas import icon_button
            viewer=zoomable_view(content=ft.Image(src=data,fit=ft.BoxFit.CONTAIN),expand=True)
            async def zoom_in(e): await viewer.zoom(1.25)
            async def zoom_out(e): await viewer.zoom(.8)
            async def reset(e): await viewer.reset()
            content=ft.Column([text('Use o mouse para ampliar e arrastar o gráfico.',12,MUTED),
                ft.Row([icon_button(icon=ft.Icons.ZOOM_IN,tooltip='Ampliar',on_click=zoom_in),
                        icon_button(icon=ft.Icons.ZOOM_OUT,tooltip='Reduzir',on_click=zoom_out),
                        ft.TextButton('Recentrar',on_click=reset)]),viewer],expand=True)
        else:
            raise ValueError('Este arquivo não possui visualizador.')
        previous_editing=getattr(self,'editing_stage',False)
        self.editing_stage=True
        def dismissed(e):
            self.editing_stage=previous_editing
            if result:
                result.close()
                if getattr(self,'result_viewer',None) is result:self.result_viewer=None
        self.page.show_dialog(ft.AlertDialog(title=text(title,22),content=ft.Container(width=width,height=height,content=content),
                                            on_dismiss=dismissed,
                                            actions=[ft.TextButton('Fechar',on_click=lambda e:self.page.pop_dialog())]))

    @with_processing('Carregando resultados da etapa…','Organizando as tabelas, estruturas e gráficos do bloco. Aguarde um instante.')
    async def open_stage_results(self,stage_id):
        token = self.token
        project_id=self.current['id']
        runs=await self.call(self.store.list_runs,token,project_id)
        if self.token != token or not self.current or self.current['id'] != project_id:
            return
        for run in runs:
            stage=next((s for s in run['stages'] if s['id']==stage_id and s['status']=='succeeded'),None)
            if stage:
                self.results_focus=(run['id'],stage_id)
                await self.change_tab('Execuções')
                result=self.stage_results.get(self.results_focus)
                self.results_focus=None
                if result:await result.load()
                return
        raise ValueError('Execute este bloco para gerar seus gráficos. Resultados antigos continuam disponíveis em Execuções.')

    @with_processing('Cancelando execução…','Interrompendo os processos e atualizando o estado da execução. Aguarde um instante.')
    async def cancel_run(self,run_id):
        token,project=self.token,self.current
        await self.call(self.service.cancel,token,run_id)
        run=await self.call(self.store.get_run,token,run_id)
        if self.token!=token or self.current is not project or not project or project['id']!=run['project_id']:return
        if not hasattr(self,'_cancelled_run_ids'):self._cancelled_run_ids=set()
        self._cancelled_run_ids.add(run_id)
        key=getattr(self,'selection_dialog_key',None)
        if key and key[0]==run_id:
            dialog=getattr(self,'selection_dialog',None)
            if dialog:close_dialog(self.page,dialog)
            self.selection_dialog=None;self.selection_dialog_key=None
            self.editing_stage=False
        if getattr(self,'settings_run_id',None)==run_id:
            dialog=getattr(self,'settings_dialog',None)
            if dialog:close_dialog(self.page,dialog)
            self.settings_dialog=None;self.settings_run_id=None;self.settings_body=None
            self._dialog_request=getattr(self,'_dialog_request',0)+1
            self.editing_stage=False
        drafts=getattr(self,'pending_file_selections',{})
        for key in list(drafts):
            if key[0]==run_id:drafts.pop(key,None)
        progress=getattr(self,'run_progress',None)
        if progress and progress.valid() and progress.run['id']==run_id:progress.update(run)
        if run['status'] in ('queued','running','awaiting_input'):
            self.current_run=run_id
            if not self.polling or self.polling.done():self.polling=self.page.run_task(self.poll_runs)
        elif self.current_run==run_id:self.current_run=None
        await self.draw_project()
        self.notify('Execução cancelada.' if run['status']=='cancelled' else 'Cancelamento solicitado.')

    @with_processing('Preparando seleção de arquivos…','Organizando as entradas e os metadados do próximo bloco. Aguarde um instante.',foreground=True)
    async def configure_pending(self,run_id):
        from .file_selection import FileSelection,selection_dialog_size
        token,project=self.token,self.current
        run=await self.call(self.store.get_run,token,run_id)
        fresh=await self.call(self.store.project,token,run['project_id'])
        assets=await self.call(self.store.assets,token,run['project_id'])
        if self.token!=token or self.current is not project or not project or project['id']!=run['project_id']:return
        if run_id in getattr(self,'_cancelled_run_ids',set()):return
        if fresh['role']=='viewer':raise AccessDenied('Somente editores podem selecionar os arquivos da execução.')
        if run['status']!='awaiting_input':raise ValueError('A execução não está aguardando configuração.')
        pending=next(s for s in run['stages'] if s['status']=='awaiting_input')
        key=(run_id,pending['id'])
        if getattr(self,'selection_dialog_key',None)==key and getattr(self,'selection_dialog',None) and self.selection_dialog.open:return
        progress=getattr(self,'run_progress',None)
        if progress and progress.is_open:await progress.dismiss()
        if not hasattr(self,'pending_file_selections'):self.pending_file_selections={}
        drafts=self.pending_file_selections
        from .pdb_results import PDBActions
        from .stage_results import StageResults
        def file_actions(stage,path):
            if stage['operation']!='retrieve_structures':return []
            results=StageResults(self,run['project_id'],run_id,stage,True)
            return PDBActions(results,True,on_saved=lambda item,saved:
                form.apply_pdb_reviews({item['path']:saved})).actions({'path':str(path),'name':path.name})
        compound_references=None
        if pending['configuration']['operation'] in ('retrieve_pubchem','docking_vina','docking_dock6'):
            from biomolexplorer.pubchem_retrieval import reference_choices
            compound_references={}
            for item in run['stages']:
                if item['status']!='succeeded':continue
                for filename in item.get('artifacts',[]):
                    path=self.store.scoped_path(project['id'],filename)
                    if path.suffix.lower()=='.csv':compound_references[str(path)]=await self.call(reference_choices,[path])
            for asset in assets:
                if asset['kind'] not in ('compounds','other'):continue
                path=await self.call(self.store.asset_path,token,project['id'],asset['id'])
                compound_references['asset:'+asset['id']]=await self.call(reference_choices,[path])
            if self.token!=token or self.current is not project:return
        if run_id in getattr(self,'_cancelled_run_ids',set()):return
        reviews={}
        if pending['configuration']['operation']=='redocking':
            from biomolexplorer.result_files import ResultFiles
            reviews=await self.call(ResultFiles(self.store).pdb_review_statuses,token,project['id'],run_id)
            if self.token!=token or self.current is not project or run_id in getattr(self,'_cancelled_run_ids',set()):return
        form=FileSelection(run,pending,assets,state=drafts.get(key),file_actions=file_actions,
                           compound_references=compound_references,page=self.page)
        if reviews:form.apply_pdb_reviews(reviews)
        inactive=False
        def remember():drafts[key]=form.state()
        def dismiss(e=None):
            nonlocal inactive
            remember()
            inactive=True
            self.editing_stage=False
            close_dialog(self.page,dialog)
        async def apply(e):
            if inactive or self.token!=token or self.current is not project:return
            updated=form.read()
            if updated['operation']=='redocking':
                dismiss()
                await self.open_stage_dialog(pending['id'],waiting_run_id=run_id,configuration=updated)
                return
            resumed=await self.call(self.service.resume,token,run_id,updated)
            dismiss()
            self.current_run=resumed['id']
            if progress:progress.update(resumed)
            # Store the same inputs in the block so reopening its settings does
            # not require another selection. Keep other current block settings.
            if hasattr(self,'flow_editor'):self.flow_editor.remember()
            for stage in self.current['pipeline']:
                if stage['id']==pending['id']:
                    stage['bindings']=copy.deepcopy(updated['bindings'])
                    stage['input_processing']=updated.get('input_processing','individual')
                    if stage['operation']=='retrieve_pubchem':
                        for field in ('selection_mode','compound_id'):stage['parameters'][field]=updated['parameters'].get(field,'')
                    self.dirty=True
                    break
            await self.save_draft(silent=True)
            drafts.pop(key,None)
            await self.draw_project()
        async def configure(e):
            if inactive:return
            updated=form.read(require_selection=False)
            dismiss()
            await self.open_stage_dialog(pending['id'],waiting_run_id=run_id,configuration=updated)
        async def cancel_execution(e):
            if inactive or self.token!=token or self.current is not project:return
            await self.cancel_run(run_id)
            dismiss()
            drafts.pop(key,None)
        redocking=pending['configuration']['operation']=='redocking'
        dialog=ft.AlertDialog(modal=True,title=text('Selecionar arquivos · '+pending['name'],22),
            content=ft.Container(content=form.control,**selection_dialog_size(self.page,redocking)),
            actions=[ft.TextButton('Selecionar depois',on_click=dismiss),
                     ft.TextButton('Configurar etapa',on_click=self.event_action(configure)),
                     ft.TextButton('Cancelar execução',on_click=self.event_action(cancel_execution)),
                     button('Selecionar pares e configurar redocking' if pending['configuration']['operation']=='redocking' else 'Continuar com os arquivos selecionados',self.event_action(apply),True)])
        def invalidate():
            nonlocal inactive
            inactive=True
        dialog._configuration_invalidate=invalidate
        self._track_configuration_dialog(dialog,run_id,pending['id'])
        self.selection_dialog_key=key
        dialog._redocking_selection=redocking
        self.selection_dialog=dialog
        self.auto_selection_key=key
        self.editing_stage=True
        def dismissed(e):
            if not inactive and self.selection_dialog is dialog:
                remember()
                self.editing_stage=False
        dialog.on_dismiss=dismissed
        self.page.show_dialog(dialog)

    @background
    async def poll_runs(self):
        previous = None
        while self.token and self.current_run:
            token, run_id = self.token, self.current_run
            await asyncio.sleep(1)
            if self.token != token or self.current_run != run_id:
                continue
            try:
                run=await self.call(self.store.get_run,token,run_id)
                if self.token != token or self.current_run != run_id:
                    continue
                progress = getattr(self,'run_progress',None)
                if progress and progress.valid() and progress.run['id'] == run_id:
                    progress.update(run)
                    if run['status'] not in ('queued','running','awaiting_input'):
                        progress.show()
                    self.page.update()
                # Progress updates do not rebuild expanded result cards every second.
                snapshot = (run_id,run['status'],run.get('error'),json.dumps([
                    (s['id'],s['status'],s.get('log_path'),s.get('artifacts'),s.get('error')) for s in run['stages']]))
                if snapshot != previous and self.current and self.current['id']==run['project_id'] and self.tab=='Execuções':
                    await self.draw_project()
                if run['status']=='awaiting_input' and self.current and self.current['id']==run['project_id'] and self.current['role']!='viewer':
                    pending=next(s for s in run['stages'] if s['status']=='awaiting_input')
                    key=(run_id,pending['id'])
                    if getattr(self,'auto_selection_key',None)!=key and not getattr(self,'editing_stage',False):
                        await self.configure_pending(run_id)
                        self.auto_selection_key=key
                previous = snapshot
                if self.token == token and self.current_run == run_id and run['status'] not in ('queued','running','awaiting_input'):
                    failures=sum(s.get('redocking_failures',0) for s in run['stages'])
                    if failures:self.notify(f'Redocking: {failures} complexos ignorados por falha na preparação. Consulte os resultados destacados.')
                    self.current_run=None
            except Exception as exc:
                if self.token != token or self.current_run != run_id:
                    continue
                get_logger('frontend').exception('Falha ao acompanhar execução; run=%s',run_id)
                self.current_run=None
                self.notify(str(exc))
                if self.token:
                    await self.guard(self.show_workspace)

    async def sharing_view(self):
        members=await self.call(self.store.members,self.token,self.current['id'])
        email=ft.TextField(label='E-mail de quem você quer convidar',expand=True)
        role=ft.Dropdown(label='Permissão',value='viewer',options=[ft.DropdownOption(key='viewer',text='Leitor'),ft.DropdownOption(key='editor',text='Editor')],width=180)
        async def invite(e):
            async def action():
                await self.call(self.store.invite,self.token,self.current['id'],email.value or '',role.value)
                self.notify('Convite disponível no workspace do usuário.')
                await self.draw_project()
            await self.guard(action)
        rows=[]
        for member in members:
            rows.append(ft.Row([ft.Column([user_text(member['name'],16),user_text(member['email'],12,MUTED)],expand=True),
                text(self.tr('Editor' if member['role']=='editor' else 'Leitor')+(' · convite pendente' if not member['accepted'] else ''),12,MUTED),
                button('Revogar',self.event(self.revoke_member,member['id']))]))
        async def archive():
            p=self.current
            await self.call(self.store.update_project,self.token,p['id'],p['name'],p['description'],p['color'],p['tags'],not p['archived'])
            await self.show_workspace()
        async def confirm_delete():
            await self.delete_project_dialog(self.current['id'])
        return ft.Column([panel(ft.Column([text('Explore em equipe',24,weight=ft.FontWeight.W_600),
            text('Leitores consultam arquivos e resultados. Editores configuram etapas, enviam arquivos e executam análises.',color=MUTED),
            text('Convide uma conta já cadastrada. O convite é entregue dentro da aplicação.',12,MUTED),
            ft.Row([email,role,button('Convidar',invite,True,ft.Icons.PERSON_ADD_OUTLINED)]),*rows],spacing=18)),
            panel(ft.Column([text('Organização do projeto',20),ft.Row([button('Restaurar projeto' if self.current['archived'] else 'Arquivar projeto',self.event(archive)),
                button('Excluir projeto',self.event(confirm_delete))],wrap=True)]))],expand=True,scroll=ft.ScrollMode.AUTO,spacing=18)

    async def revoke_member(self,user_id):
        await self.call(self.store.revoke,self.token,self.current['id'],user_id)
        await self.draw_project()


def main():
    parser=argparse.ArgumentParser(description='BioMolExplorer workspace Flet')
    parser.add_argument('--data-dir',type=Path,default=Path.home()/'.local/share/biomolexplorer')
    parser.add_argument('--web',action='store_true')
    parser.add_argument('--no-browser',action='store_true',help='Serve a interface web sem abrir um navegador.')
    parser.add_argument('--host',default='127.0.0.1')
    parser.add_argument('--port',type=int,default=8550)
    parser.add_argument('--log-dir',type=Path,help='Diretório central dos logs da interface e dos workers.')
    parser.add_argument('--worker-python',type=Path)
    parser.add_argument('--dock6-path',type=Path)
    parser.add_argument('--language',choices=('en','pt'),default='en',help='Initial interface language; users can change it on the login screen.')
    args=parser.parse_args()
    if args.no_browser:
        os.environ['FLET_FORCE_WEB_SERVER']='true'
    if args.log_dir:
        os.environ['BIOMOL_LOG_DIR']=str(args.log_dir.expanduser().resolve())
    configure_logging('frontend')
    store=WorkspaceStore(args.data_dir)
    os.environ.setdefault('FLET_SECRET_KEY',secrets.token_urlsafe(48))
    os.environ.setdefault('FLET_MAX_UPLOAD_SIZE',str(store.max_upload_bytes))
    service=PipelineService(store,worker_python=args.worker_python,dock6_path=args.dock6_path)
    from biomolexplorer.updates import UpdateChecker
    updates=UpdateChecker(store.root)
    from biomolexplorer.pdb_view import StructureViewers
    viewers=StructureViewers(store)
    async def session(page):
        loop=asyncio.get_running_loop()
        def async_error(loop,context):
            exception=context.get('exception')
            get_logger('frontend').error('Falha assíncrona: %s',context.get('message'),exc_info=(type(exception),exception,exception.__traceback__) if exception else None)
            loop.default_exception_handler(context)
        loop.set_exception_handler(async_error)
        ui=WorkspaceUI(page,store,service,language=args.language,structure_viewer=viewers)
        from .updates import UpdateNotice
        ui.update_notice=UpdateNotice(ui,updates)
        ui.show_login()
        ui.update_notice.task=page.run_task(ui.update_notice.poll)
        async def disconnected(event):
            ui.update_notice.task.cancel()
        async def connected(event):
            if ui.update_notice.task.done():
                ui.update_notice.task=page.run_task(ui.update_notice.poll)
        page.on_disconnect=disconnected
        page.on_connect=connected
    try:
        if args.web or args.no_browser:
            from .web_host import run_web
            run_web(session,store,viewers,args.host,args.port,not args.no_browser,os.environ['FLET_SECRET_KEY'])
        else:
            ft.run(session,view=ft.AppView.FLET_APP,host=args.host,port=args.port,
                   upload_dir=str(store.staging),assets_dir=None)
    finally:
        viewers.close()
        service.close()


if __name__=='__main__':
    main()
