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
from biomolexplorer.diagnostics import configure_logging, get_logger
from biomolexplorer.visualizations import SUFFIX, MAX_VIEW_BYTES, load_view
from .branding import flag, wordmark, BRAND_BLUE, BRAND_NAVY
from .feedback import readable_error,close_dialog
from biomolexplorer.workspace import WorkspaceStore, AccessDenied, COLORS
from .project_tools import ProjectTools
from .zoom import zoomable_view
from .localization import LocalizedPage, verbatim

INK = '#172B4D'
MUTED = '#64748B'
TEAL = BRAND_BLUE
BG = '#F3F6FC'
LINE = '#E2E8F0'
STATUS = {'queued':'Aguardando', 'running':'Executando', 'succeeded':'Concluído', 'failed':'Falhou',
          'cancelled':'Cancelado', 'interrupted':'Interrompido', 'skipped':'Não executado',
          'awaiting_input':'Aguardando seleção de arquivos'}
KINDS = {'compounds':'Compostos (CSV)', 'structures':'Complexos PDB', 'prepared_structures':'Receptores preparados',
         'fingerprints':'Fingerprints', 'similarity':'Relações de similaridade', 'vina':'Resultados Vina',
         'dock6':'Resultados DOCK6', 'scores':'Scores consolidados (CSV)', 'visualization':'Gráficos e visualizações', 'other':'Outros arquivos'}


def text(value, size=14, color=INK, weight=None):
    return ft.Text(value,size=size,color=color,weight=weight)


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
        return await asyncio.to_thread(function,*args,**kwargs)

    def tr(self,message):
        translator=getattr(self.page,'translator',None)
        return translator(message) if translator else message

    def notify(self,message):
        self.page.show_dialog(ft.SnackBar(content=ft.Text(str(message)),bgcolor=INK))

    async def guard(self,action):
        token = self.token
        try:
            await action()
        except Exception as exc:
            get_logger('frontend').exception('Falha em ação da interface; project=%s tab=%s',self.current['id'] if self.current else None,self.tab)
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
        self.dirty = False
        self.selected = None
        self.current_run = None
        self.run_progress = None
        self.selection_dialog = None
        self.selection_dialog_key = None
        self.auto_selection_key = None
        self.page.on_keyboard_event = None
        self.pending_file_selections = {}

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
        bar=ft.Container(height=76,padding=ft.Padding.symmetric(horizontal=32),bgcolor='#FFFFFF',
            border=ft.Border.only(bottom=ft.BorderSide(1,LINE)),
            content=ft.Row([brand,ft.Container(expand=True),
                ft.Semantics(container=True,label='Selecionar idioma',button=True,content=language_choice)],spacing=8))
        def resize():
            compact=(getattr(self.page,'width',None) or 1440)<650
            bar.padding=ft.Padding.symmetric(horizontal=16 if compact else 32)
            brand.width=160 if compact else 190
            brand.height=brand.width*209/662
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

    def shell(self,title,subtitle,body,actions=None):
        user = self.store.user(self.token)
        sidebar = ft.Container(width=220,bgcolor=BRAND_NAVY,padding=22,content=ft.Column([
            text('Workspace',18,'#FFFFFF',ft.FontWeight.W_600),
            ft.Container(height=12),button('Meus projetos',self.event(self.show_workspace),icon=ft.Icons.GRID_VIEW),
            ft.Container(expand=True),user_text(user['name'],16,'#FFFFFF'),user_text(user['email'],12,'#94A3B8'),
            ft.TextButton('Sair',on_click=self.logout)],spacing=12))
        project=getattr(self,'current',None)
        heading=user_text(title,30,weight=ft.FontWeight.W_700) if project and title==project['name'] else text(title,30,weight=ft.FontWeight.W_700)
        summary=user_text(subtitle,color=MUTED) if project and subtitle==project['description'] else text(subtitle,color=MUTED)
        header = ft.Column([ft.Column([heading,summary]),ft.Row(actions or [],wrap=True)],spacing=16)
        main = ft.Container(expand=True,padding=30,content=ft.Column([header,ft.Container(height=8),body],expand=True,spacing=16))
        top_bar,_,resize_bar=self.top_bar()
        mobile_nav=ft.Container(padding=ft.Padding.symmetric(horizontal=16),bgcolor='#FFFFFF',
            content=ft.Row([ft.TextButton('Meus projetos',icon=ft.Icons.GRID_VIEW,on_click=self.event(self.show_workspace)),
                ft.Container(expand=True),ft.TextButton('Sair',on_click=self.logout)]))
        def resize():
            compact=(getattr(self.page,'width',None) or 1440)<900
            resize_bar()
            sidebar.visible=not compact
            mobile_nav.visible=compact
            main.padding=16 if compact else 30
        async def on_resize(e):
            resize();self.page.update()
        resize()
        self.page.on_resize=on_resize
        return ft.Column([top_bar,mobile_nav,ft.Row([sidebar,main],expand=True,spacing=0,
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
        self.render(self.shell('Seu workspace','Explore, organize e compartilhe suas pesquisas.',body,
            [button('Importar projeto',self.event(self.import_project_dialog),icon=ft.Icons.UPLOAD_FILE),
             button('Novo projeto',self.event(self.project_dialog),True,ft.Icons.ADD)]))

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
                self.page.pop_dialog()
                await self.show_workspace()
                self.notify('Projeto e pasta removidos permanentemente.')
            await self.guard(action)
        self.page.show_dialog(ft.AlertDialog(modal=True,title=text('Excluir projeto?',24),
            content=ft.Container(width=460,content=ft.Column([
                user_text(project['name'],18,weight=ft.FontWeight.W_600),
                text('O projeto será excluído do workspace e sua pasta será removida permanentemente, incluindo entradas, resultados e histórico.'),
                user_text(project['directory'],12),
                text('Há uma execução ativa. Conclua ou cancele a execução antes de excluir o projeto.' if active else
                     'Esta ação é permanente e não pode ser desfeita.',12,'#B91C1C')],tight=True,spacing=18)),
            actions=[ft.TextButton('Manter projeto',on_click=lambda e:self.page.pop_dialog()),
                     button('Excluir projeto',remove,icon=ft.Icons.DELETE_OUTLINE,disabled=active)]))

    async def draw_project(self):
        project = self.current
        if project is None:
            return
        token, tab = self.token, self.tab
        # Recheck membership even when rendering an unsaved local draft.
        fresh = await self.call(self.store.project,token,project['id'])
        if self.token != token or self.current is not project or self.tab != tab:
            return
        if (not self.dirty and not getattr(self,'editing_stage',False)) or fresh['role'] == 'viewer':
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
                      for label in ('Pipeline','Arquivos','Execuções','Compartilhar') if label!='Compartilhar' or project['role']=='owner'],wrap=True)
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
        actions = [button('Workspace',self.event(self.show_workspace),icon=ft.Icons.ARROW_BACK)]
        if writable:
            actions.append(button('Projeto',self.event(self.project_dialog,project),icon=ft.Icons.TUNE))
        if project['role']=='owner':
            actions.append(button('Compartilhar',self.event(self.change_tab,'Compartilhar'),icon=ft.Icons.PERSON_ADD_ALT_1))
        actions.append(button('Histórico',self.event(self.history_dialog,project['id']),icon=ft.Icons.HISTORY))
        progress = getattr(self,'run_progress',None)
        extras = [progress.badge] if progress and progress.valid() else []
        self.render(self.shell(project['name'],project['description'] or 'Construa o pipeline que sua pesquisa precisa.',
            ft.Column([nav,*extras,content],expand=True,spacing=18),actions))
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
        return self.flow_editor.build()

    async def open_stage_dialog(self,stage_id,waiting_run_id=None,configuration=None):
        from .guided import GuidedForm
        from biomolexplorer.pipeline import validate_pipeline
        token, project = self.token, self.current
        fresh=await self.call(self.store.project,token,project['id'])
        if self.token != token or self.current is not project:
            return
        project['role']=fresh['role']
        self.selected=stage_id
        original=self.selected_stage()
        if not original:
            return
        draft=copy.deepcopy(original)
        writable=self.current['role']!='viewer'
        assets=await self.call(self.store.assets,token,project['id'])
        if self.token != token or self.current is not project:
            return
        runs=await self.call(self.store.list_runs,token,project['id'])
        if self.token != token or self.current is not project:
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
        self.artifact_choices={}
        for run in runs:
            for item in run['stages']:
                if item['id'] in self.artifact_choices or item['status']!='succeeded':continue
                from biomolexplorer.artifact_choices import choices
                self.artifact_choices[item['id']]=choices(item)
        for source in self.current['pipeline']:
            if source['operation']=='import_results' or source.get('provided_results'):
                params=source.get('provided_results') or source['parameters']
                self.artifact_choices.setdefault(source['id'],set()).update(a['name'] for a in assets if a['id'] in params.get('asset_ids',[]))
            elif source['operation']=='retrieve_compounds':
                self.artifact_choices.setdefault(source['id'],set()).add('compounds.csv')
        mode='visual'
        reader=None
        body=ft.Container(width=min(960,(self.page.width or 1160)-100),height=max(300,min(620,(self.page.height or 900)-250)),padding=24,bgcolor=BG,border_radius=16)
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
                    result=copy.deepcopy(draft); result['name']=title.value or draft['name']; result['parameters']=json.loads(parameters.value); result['bindings']=json.loads(bindings.value)
                    result['templates']={name:field.value for name,field in templates.items() if field.value!=(RESOURCE_ROOT/name).read_text()}
                    return result
                reader=read_advanced
        async def visual(e):
            await self.guard(lambda:render_mode('visual')); self.page.update()
        async def advanced(e):
            await self.guard(lambda:render_mode('advanced')); self.page.update()
        async def apply(e):
            await self.call(self.store.project,self.token,self.current['id'],'editor')
            updated=reader()
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
                validate_operation(updated['operation'],parameters)
            candidate=[updated if s['id']==stage_id else s for s in self.current['pipeline']]
            validate_pipeline(candidate)
            if waiting_run_id:
                resumed=await self.call(self.service.resume,token,waiting_run_id,updated)
            if hasattr(self,'flow_editor'):self.flow_editor.remember()
            self.current['pipeline']=candidate; self.dirty=True
            await self.save_draft(silent=True)
            self.editing_stage=False
            close_dialog(self.page,dialog)
            if waiting_run_id:
                self.current_run=resumed['id']
                if getattr(self,'run_progress',None):self.run_progress.update(resumed)
                self.tab='Execuções'
                await self.draw_project()
                if not self.polling or self.polling.done():self.polling=self.page.run_task(self.poll_runs)
            elif hasattr(self,'flow_editor'):self.flow_editor.refresh()
        await render_mode(mode)
        dialog=ft.AlertDialog(title=ft.Column([text(TITLES[original['operation']][2].upper(),11,TEAL,ft.FontWeight.W_600),user_text(original['name'],24,weight=ft.FontWeight.W_700),text(TITLES[original['operation']][1],13,MUTED),ft.Row([
            ft.TextButton('Visual · formulário guiado',on_click=visual),ft.TextButton('Avançado · parâmetros e scripts',on_click=advanced)],wrap=True)],spacing=12),content=body,
            actions=[ft.TextButton('Cancelar' if writable else 'Fechar',on_click=lambda e:close_dialog(self.page,dialog)),button('Aplicar e continuar pipeline' if waiting_run_id else 'Aplicar configuração',self.event_action(apply),True,disabled=not writable)])
        self.editing_stage=True
        def dismissed(e): self.editing_stage=False
        dialog.on_dismiss=dismissed
        self.page.show_dialog(dialog)

    def event_action(self,handler):
        async def wrapped(e):
            await self.guard(lambda:handler(e))
        return wrapped

    def auto_bind(self,stage,previous):
        operation=stage['operation']
        wanted={
            'expand_similar_compounds': {'base_input_path':['retrieve_compounds']},
            'admet': {'base_input_path':['graphs','retrieve_compounds','import_results']},
            'fingerprints': {'base_input_path':['admet','retrieve_compounds','import_results']},
            'similarity': {'base_input_path':['fingerprints','import_results']},
            'graphs': {'base_input_path':['retrieve_compounds','import_results'], 'similarity_path':['similarity']},
            'prepare_structures': {'base_input_path':['import_results','retrieve_structures']},
            'redocking': {'base_input_path':['retrieve_structures','import_results']},
            'docking_vina': {'base_input_path':['prepare_structures','redocking','import_results'], 'base_selected_mols':['admet','graphs','retrieve_compounds']},
            'docking_dock6': {'base_input_path':['prepare_structures','redocking','import_results'], 'base_selected_mols':['admet','graphs','retrieve_compounds'], 'base_vina_path':['docking_vina']},
            'consensus': {'base_input_path':['docking_vina'], 'base_vina_path':['docking_vina','import_results'], 'base_dock6_path':['docking_dock6','import_results']},
        }.get(operation,{})
        for field,candidates in wanted.items():
            def compatible(s):
                if s['operation'] not in candidates:
                    return False
                if s['operation']!='import_results':
                    return True
                kind=s['parameters'].get('kind')
                if field=='base_selected_mols' or operation in ('admet','fingerprints','graphs'):
                    return kind=='compounds'
                if field=='base_vina_path':
                    return kind=='vina'
                if field=='base_dock6_path':
                    return kind=='dock6'
                if operation in ('prepare_structures','redocking','docking_vina','docking_dock6'):
                    return kind in ('structures','prepared_structures')
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
        self.base_pipeline=copy.deepcopy(saved['pipeline'])
        if project['pipeline']==sent:
            self.current=saved;self.dirty=False
        else:
            project['revision']=saved['revision'];project['updated']=saved['updated'];self.schedule_save()
        if self.tab=='Pipeline' and hasattr(self,'flow_editor'):
            self.flow_editor.refresh()
        if not silent:self.notify('Pipeline salvo.')

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
                    'Etapas novas ou alteradas serão processadas normalmente.')),
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
                files=await self.picker.pick_files(allow_multiple=True)
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
                await self.picker.save_file(file_name=item['name'],src_bytes=data)
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
                stage_label='Reaproveitado' if stage.get('reused') else 'Resultados fornecidos' if stage.get('provided') else STATUS[stage['status']]
                if stage.get('excluded_records'):stage_label+=f' · {stage["excluded_records"]} registros excluídos'
                stage_rows.append(ft.ExpansionTile(title=user_text(stage['name'],16),subtitle=text(stage_label,12,TEAL if stage['status']=='succeeded' else MUTED),
                    expanded=result.expanded,on_change=result.expand,
                    controls=[text(readable_error(stage.get('error','')),color='#DC2626'),ft.Row(report_buttons,wrap=True),result.root,ft.TextButton('Ver log',on_click=log)]))
            actions=[]
            if run['status'] in ('queued','running','awaiting_input') and writable:
                actions=[button('Cancelar execução',self.event(self.cancel_run,run['id']),icon=ft.Icons.STOP)]
                if run['status']=='awaiting_input':
                    pending=next(s for s in run['stages'] if s['status']=='awaiting_input')
                    actions.insert(0,button('Configurar '+pending['name'],self.event(self.configure_pending,run['id']),True,ft.Icons.EDIT))
            timestamp=datetime.fromtimestamp(run['created'],timezone.utc).strftime('%d/%m/%Y %H:%M UTC')
            label='Execução atual' if run['status'] in ('queued','running','awaiting_input') else 'Última execução' if run['id']==focus_id else 'Execução anterior'
            card=panel(ft.Column([text(label+' · '+timestamp,12,MUTED),ft.Row([text(STATUS[run['status']],20,weight=ft.FontWeight.W_600),
                text(run['id'][:8],12,MUTED),*actions],wrap=True),text(readable_error(run.get('error')),color=TEAL if run['status']=='awaiting_input' else '#DC2626'),*stage_rows],spacing=10))
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
        return ft.Column([button('Atualizar resultados',self.event(self.draw_project),icon=ft.Icons.REFRESH),*content],expand=True,scroll=ft.ScrollMode.AUTO,spacing=16)

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

    async def download_stage_report(self,project_id,path):
        token=self.token
        data=await self.call(self.store.read_file,token,project_id,path)
        if self.token==token and self.current and self.current['id']==project_id:
            await self.picker.save_file(file_name=Path(path).name,src_bytes=data)

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

    async def preview_artifact(self,project_id,path):
        token = self.token
        if Path(path).suffix.lower()=='.pdb':
            url=await self.call(self.pdb_view_url,project_id,path)
            if self.token!=token or not self.current or self.current['id']!=project_id:return
            await ft.UrlLauncher().launch_url(url,mode=ft.LaunchMode.EXTERNAL_APPLICATION,web_only_window_name='_blank')
            return
        data=await self.call(self.store.read_file,token,project_id,path,MAX_VIEW_BYTES+1)
        if self.token != token or not self.current or self.current['id'] != project_id:
            return
        if len(data)>MAX_VIEW_BYTES:
            raise ValueError('Visualização maior que o limite de 32 MB. Baixe o arquivo para consultá-lo localmente.')
        width=max(320,min(1160,(self.page.width or 1440)-100))
        height=max(300,min(690,(self.page.height or 1000)-220))
        if str(path).endswith(SUFFIX):
            from .results_viewer import ResultsViewer
            model=await self.call(load_view,data)
            if self.token != token or not self.current or self.current['id'] != project_id:
                return
            self.result_viewer=ResultsViewer(self,project_id,model)
            content=ft.Column([self.result_viewer.build()],scroll=ft.ScrollMode.AUTO)
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
        def dismissed(e): self.editing_stage=previous_editing
        self.page.show_dialog(ft.AlertDialog(title=text(title,22),content=ft.Container(width=width,height=height,content=content),
                                            on_dismiss=dismissed,
                                            actions=[ft.TextButton('Fechar',on_click=lambda e:self.page.pop_dialog())]))

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

    async def cancel_run(self,run_id):
        await self.call(self.service.cancel,self.token,run_id)
        self.notify('Cancelamento solicitado.')

    async def configure_pending(self,run_id):
        from .file_selection import FileSelection
        token,project=self.token,self.current
        run=await self.call(self.store.get_run,token,run_id)
        fresh=await self.call(self.store.project,token,run['project_id'])
        assets=await self.call(self.store.assets,token,run['project_id'])
        if self.token!=token or self.current is not project or not project or project['id']!=run['project_id']:return
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
            return PDBActions(results,True).actions({'path':str(path),'name':path.name})
        form=FileSelection(run,pending,assets,state=drafts.get(key),file_actions=file_actions)
        inactive=False
        def remember():drafts[key]=form.state()
        def dismiss(e=None):
            nonlocal inactive
            remember()
            inactive=True
            self.editing_stage=False
            close_dialog(self.page,dialog)
        async def apply(e):
            if self.token!=token or self.current is not project:return
            updated=form.read()
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
                    self.dirty=True
                    break
            await self.save_draft(silent=True)
            drafts.pop(key,None)
            await self.draw_project()
        async def configure(e):
            updated=form.read(require_selection=False)
            dismiss()
            await self.open_stage_dialog(pending['id'],waiting_run_id=run_id,configuration=updated)
        dialog=ft.AlertDialog(modal=True,title=text('Selecionar arquivos · '+pending['name'],22),
            content=ft.Container(content=form.control,width=min(860,(self.page.width or 1160)-100),
                                 height=max(300,min(580,(self.page.height or 900)-250))),
            actions=[ft.TextButton('Selecionar depois',on_click=dismiss),
                     ft.TextButton('Configurar etapa',on_click=self.event_action(configure)),
                     button('Continuar com os arquivos selecionados',self.event_action(apply),True)])
        self.selection_dialog_key=key
        self.selection_dialog=dialog
        self.auto_selection_key=key
        self.editing_stage=True
        def dismissed(e):
            if not inactive:
                remember()
                self.editing_stage=False
        dialog.on_dismiss=dismissed
        self.page.show_dialog(dialog)

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
        ui.show_login()
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
