"""Navigation, responsive login and collapsible stage-library interactions."""
import asyncio
import unittest
from types import SimpleNamespace

import flet as ft

from biomolexplorer.ui.app import WorkspaceUI
from biomolexplorer.ui.flow_canvas import FlowCanvas
from biomolexplorer.ui.localization import LocalizedPage,Translator


def walk(control):
    if isinstance(control,(list,tuple)):
        for child in control:yield from walk(child)
    elif isinstance(control,ft.Control):
        yield control
        for attribute in LocalizedPage._children:yield from walk(getattr(control,attribute,None))


class LayoutTests(unittest.TestCase):
    def page(self):
        page=SimpleNamespace(width=1200,height=900,controls=[],web=False,update=lambda:None)
        page.add=lambda *controls:page.controls.extend(controls)
        return page

    def test_login_bar_remains_visible_when_resized(self):
        page=self.page();ui=WorkspaceUI(page,None,None)
        ui.show_login(True)
        root=page.controls[0]
        bar,body=root.controls
        self.assertEqual(bar.height,76)
        self.assertTrue(any(isinstance(c,ft.PopupMenuButton) for c in walk(bar)))
        self.assertFalse(any(isinstance(c,ft.Dropdown) for c in walk(body)))
        page.width=360;page.height=640
        asyncio.run(page.on_resize(None))
        self.assertEqual(bar.content.controls[0].width,160)
        self.assertEqual(body.padding,16)
        self.assertLessEqual(body.content.controls[0].height,640-76-32)

    def test_workspace_language_change_preserves_unsaved_pipeline(self):
        page=self.page();ui=WorkspaceUI(page,None,None)
        ui.current={'pipeline':[{'name':'Meu estudo'}]};ui.dirty=True
        bar,menu,_=ui.top_bar()
        ui.page.add(bar,ft.Text('Biblioteca de etapas'))
        english=next(item for item in menu.items if item.data=='en')
        asyncio.run(english.on_click(SimpleNamespace(control=english)))
        self.assertEqual(page.controls[-1].value,'Stage library')
        self.assertEqual(ui.current['pipeline'],[{'name':'Meu estudo'}])
        self.assertTrue(ui.dirty)

    def test_library_search_expands_matches_and_restores_categories(self):
        ui=SimpleNamespace(current={'pipeline':[]},selected=None,page=self.page(),tr=Translator('en'),
            event=lambda *args:lambda e:None,save_draft=None,preset_dialog=None,run_pipeline=None)
        controls=list(walk(FlowCanvas(ui,True).build()))
        sections=[c for c in controls if isinstance(c,ft.ExpansionTile)]
        self.assertEqual(len(sections),4)
        self.assertTrue(all(not s.expanded for s in sections))
        search=next(c for c in controls if isinstance(c,ft.TextField))
        sections[0].expanded=True
        search.value='ADMET';search.on_change(None)
        matches=[s for s in sections if s.visible]
        self.assertEqual(len(matches),1)
        self.assertEqual(matches[0].title.value,'Análise')
        self.assertTrue(matches[0].expanded)
        self.assertEqual(sum(c.visible for c in matches[0].controls),1)
        search.value='does not exist';search.on_change(None)
        self.assertTrue(all(not s.visible for s in sections))
        search.value='';search.on_change(None)
        self.assertTrue(all(s.visible for s in sections))
        self.assertEqual([s.expanded for s in sections],[True,False,False,False])


if __name__=='__main__':unittest.main()
