"""Sidebar actions, contextual help, canvas placement and outcome presentation."""
import asyncio
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import AsyncMock

import flet as ft
from biomolexplorer.catalog import new_stage
from biomolexplorer import flow
from biomolexplorer.workspace import WorkspaceStore
from biomolexplorer.ui.app import WorkspaceUI
from biomolexplorer.ui.action_help import ACTION_HELP,help_dimensions,show_action_help,action_row
from biomolexplorer.ui.flow_canvas import FlowCanvas
from biomolexplorer.ui.localization import LocalizedPage,Translator
from biomolexplorer.ui.run_presentation import outcome_overview


def walk(root):
    if isinstance(root,(list,tuple)):
        for item in root:yield from walk(item)
    elif isinstance(root,ft.Control):
        yield root
        for attr in LocalizedPage._children:yield from walk(getattr(root,attr,None))


def caption(control):
    return control.content.value if isinstance(control.content,ft.Text) else control.content


class Page:
    width=1440; height=1000; web=False
    def __init__(self):self.controls=[];self.dialogs=[]
    def add(self,*controls):self.controls.extend(controls)
    def update(self,*args):pass
    def show_dialog(self,dialog):dialog.open=True;self.dialogs.append(dialog)
    def pop_dialog(self):return self.dialogs.pop() if self.dialogs else None
    def run_task(self,*args):return SimpleNamespace(done=lambda:False,cancel=lambda:None)


class FrontendNavigationTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.store=WorkspaceStore(self.tmp.name)
        self.token=self.store.register('Researcher','front@example.org','front-test-password')
        self.project=self.store.create_project(self.token,'Layout study')
        self.page=Page();self.ui=WorkspaceUI(self.page,self.store,SimpleNamespace(dock6_path=None))
        self.addCleanup(self.ui.structure_viewer.close)
        self.ui.token=self.token
        async def call(function,*args,**kwargs):return function(*args,**kwargs)
        self.ui.call=call

    def sidebar(self):return self.page.controls[0].controls[1].controls[0]

    def test_workspace_actions_live_in_sidebar_and_help_opens(self):
        asyncio.run(self.ui.show_workspace())
        buttons=[c for c in walk(self.sidebar()) if isinstance(c,ft.Button)]
        self.assertEqual([caption(c) for c in buttons],['Meus projetos','Importar projeto','Novo projeto'])
        icons=[c for c in walk(self.sidebar()) if isinstance(c,ft.IconButton) and c.icon==ft.Icons.INFO_OUTLINE]
        self.assertEqual(len(icons),3)
        for icon in icons:
            icon.on_click(None)
            self.assertGreater(self.page.dialogs[-1].content.height,0)
        main=self.page.controls[0].controls[1].controls[1]
        self.assertFalse(any(isinstance(c,ft.Button) and caption(c) in ('Importar projeto','Novo projeto') for c in walk(main)))

    def test_project_groups_collapse_persist_and_menu_reclaims_space(self):
        self.ui.current=self.project
        asyncio.run(self.ui.draw_project())
        sidebar=self.sidebar()
        groups=[c for c in walk(sidebar) if isinstance(c,ft.ExpansionTile)]
        self.assertEqual([g.title.value for g in groups[:2]],['Navegação do projeto','Pipeline e execução'])
        self.assertEqual(len(groups),6)
        groups[1].on_change(SimpleNamespace(data='false'))
        asyncio.run(self.ui.draw_project())
        sidebar=self.sidebar()
        groups=[c for c in walk(sidebar) if isinstance(c,ft.ExpansionTile)]
        self.assertFalse(groups[1].expanded)
        menu=next(c for c in walk(self.page.controls) if isinstance(c,ft.IconButton) and c.icon==ft.Icons.MENU)
        menu.on_click(None);self.assertFalse(sidebar.visible)
        menu.on_click(None);self.assertTrue(sidebar.visible)
        main=self.page.controls[0].controls[1].controls[1]
        self.assertFalse(any(isinstance(c,ft.Button) and caption(c)=='Executar pipeline' for c in walk(main)))
        self.assertGreaterEqual(self.ui.flow_editor.stack.width,4000)

    def test_viewer_cannot_run_or_edit_but_can_read_action_help(self):
        guest=self.store.register('Reader','reader-front@example.org','reader-test-password')
        self.store.invite(self.token,self.project['id'],'reader-front@example.org','viewer')
        self.store.accept_invitation(guest,self.project['id'])
        self.ui.token=guest;self.ui.current=self.store.project(guest,self.project['id'])
        asyncio.run(self.ui.draw_project())
        buttons=[c for c in walk(self.sidebar()) if isinstance(c,ft.Button)]
        self.assertFalse(any(caption(c) in ('Projeto','Compartilhar') for c in buttons))
        for name in ('Salvar','Modelos prontos','Executar pipeline','Executar até seleção'):
            self.assertTrue(next(c for c in buttons if caption(c)==name).disabled)
        info=next(c for c in walk(self.sidebar()) if isinstance(c,ft.IconButton) and 'Executar pipeline' in str(c.tooltip))
        info.on_click(None);self.assertEqual(self.page.dialogs[-1].title.value,'Executar pipeline')

    def test_help_scales_with_content_viewport_and_language(self):
        short=help_dimensions(self.page,['Uma ação simples.'])
        long=help_dimensions(self.page,['Explicação detalhada. '*100])
        self.assertLess(short[1],long[1]);self.assertLessEqual(long[1],self.page.height-210)
        self.page.width=360;self.page.height=560
        width,height=help_dimensions(self.page,['Explicação detalhada. '*100])
        self.assertLessEqual(width,360-64);self.assertLessEqual(height,350)
        translator=Translator('en')
        for key,values in ACTION_HELP.items():
            for value in values[1:]:self.assertNotEqual(translator(value),value,(key,value))
        self.ui.page.set_language('en')
        dialog=show_action_help(self.ui,'run')
        self.assertEqual(dialog.title.value,'Run pipeline')
        row=action_row(self.ui,'run',ft.Button('Executar pipeline'))
        self.ui.page.add(row)
        self.assertEqual(row.controls[1].tooltip,'Information · Run pipeline')
        self.assertNotIn('A configuração',dialog.content.content.controls[1].value)

    def test_duplicate_and_drop_find_free_space_without_moving_existing_blocks(self):
        original=new_stage('admet');original['position']={'x':60,'y':64}
        self.ui.current=dict(self.project,pipeline=[original]);self.ui.selected=original['id']
        self.ui.schedule_save=lambda:None
        editor=FlowCanvas(self.ui,True)
        asyncio.run(editor.duplicate(None))
        asyncio.run(editor.add('admet',60,64))
        self.assertEqual(original['position'],{'x':60,'y':64})
        boxes=[flow.position(s) for s in editor.stages]
        self.assertEqual(len({(p['x'],p['y']) for p in boxes}),3)
        for i,p in enumerate(boxes):
            for q in boxes[i+1:]:
                self.assertTrue(abs(p['x']-q['x'])>=flow.NODE_WIDTH+24 or abs(p['y']-q['y'])>=flow.node_height(original)+24)

    def test_complete_summary_counts_unique_files_and_reuse(self):
        run=dict(id='test-run',status='succeeded',created=100,updated=225,stages=[
            dict(status='succeeded',artifacts=['a.csv','b.png'],reused=True),
            dict(status='succeeded',artifacts=['b.png','c.csv'])])
        values=[c.value for c in walk(outcome_overview(run,[])) if isinstance(c,ft.Text)]
        self.assertIn('2 / 2',values);self.assertIn('3',values);self.assertIn('1',values)
        self.assertIn('00:02:05',values);self.assertIn('Explore os resultados de cada etapa abaixo.',values)
