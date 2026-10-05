"""The project import form retains selected packages across its event callbacks."""
import asyncio
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from biomolexplorer.workspace import WorkspaceStore
try:
    from biomolexplorer.ui.app import WorkspaceUI
except ImportError:WorkspaceUI=None


@unittest.skipIf(WorkspaceUI is None,'Install the ui extra')
class ProjectToolsUITests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)
        self.store=WorkspaceStore(self.root/'workspace')
        self.token=self.store.register('Owner','owner@example.org','owner-password-2026')
        self.project=self.store.create_project(self.token,'UI import')
        self.archive=self.store.export_project(self.token,self.project['id'])
        self.dialogs=[];self.messages=[]
        self.ui=object.__new__(WorkspaceUI)
        self.ui.page=SimpleNamespace(web=True,update=lambda:None,show_dialog=self.dialogs.append)
        self.ui.store=self.store;self.ui.token=self.token;self.ui.current=None
        async def call(function,*args,**kwargs):return function(*args,**kwargs)
        async def guard(action):await action()
        async def files(*args,**kwargs):return [self.archive]
        async def opened(pid):self.ui.current=self.store.project(self.token,pid)
        self.ui.call=call;self.ui.guard=guard;self.ui.pick_uploads=files;self.ui.open_project=opened;self.ui.notify=self.messages.append

    def tearDown(self):self.temp.cleanup()

    def test_new_project_requires_an_explicit_destination_folder(self):
        asyncio.run(self.ui.project_dialog())
        dialog=self.dialogs[-1]
        controls=dialog.content.content.controls
        controls[0].value='Novo estudo'
        directory=controls[-1].controls[0]
        self.assertEqual(directory.value,'')
        with self.assertRaisesRegex(ValueError,'Informe a pasta do projeto'):
            asyncio.run(dialog.actions[-1].on_click(None))
        self.assertIn('arquivos e resultados',directory.error)
        self.assertEqual(len(self.store.list_projects(self.token)),1)
        self.assertIsNone(self.ui.current)

    def test_new_project_uses_the_folder_selected_by_the_user(self):
        asyncio.run(self.ui.project_dialog())
        dialog=self.dialogs[-1]
        controls=dialog.content.content.controls
        controls[0].value='Novo estudo'
        destination=self.root/'experimentos'/'estudo com espaços'
        field=controls[-1].controls[0]
        self.assertTrue(field.read_only)
        self.ui.page.web=False
        async def selected_folder(**kwargs):return str(destination)
        self.ui.picker=SimpleNamespace(get_directory_path=selected_folder)
        asyncio.run(field.suffix_icon.on_tap(None))
        asyncio.run(dialog.actions[-1].on_click(None))
        self.assertEqual(self.ui.current['directory'],str(destination))
        self.assertTrue((destination/'project.json').is_file())
        self.assertFalse((self.store.projects_root/self.ui.current['id']).exists())
        self.assertFalse(dialog.open)

    def test_native_folder_picker_uses_the_selected_folder(self):
        destination=self.root/'pasta escolhida'
        self.ui.page.web=False
        async def selected_folder(**kwargs):return str(destination)
        self.ui.picker=SimpleNamespace(get_directory_path=selected_folder)
        field,folder=self.ui.folder_controls()
        field.error='Informe a pasta do projeto.'
        self.assertTrue(field.read_only)
        asyncio.run(field.suffix_icon.on_tap(None))
        self.assertEqual(field.value,str(destination))
        self.assertIsNone(field.error)

    def test_import_button_uses_package_selected_by_the_picker(self):
        asyncio.run(self.ui.import_project_dialog())
        dialog=self.dialogs[-1]
        controls=dialog.content.content.controls
        controls[-1].controls[0].value=str(self.root/'selected-destination')
        asyncio.run(controls[1].on_click(None))
        self.assertIn('Pacote selecionado',controls[2].value)
        asyncio.run(dialog.actions[-1].on_click(None))
        self.assertEqual(self.ui.current['directory'],str(self.root/'selected-destination'))
        self.assertTrue((self.root/'selected-destination'/'project.json').exists())
        self.assertFalse(dialog.open)

    def test_export_response_is_discarded_after_logout(self):
        async def logout_after_call(function,*args,**kwargs):
            value=function(*args,**kwargs)
            self.ui.token=None
            return value
        self.ui.call=logout_after_call
        asyncio.run(self.ui.export_project_dialog(self.project['id']))
        self.assertEqual(self.messages,[])

    def test_folder_picker_cancellation_preserves_selection(self):
        self.ui.page.web=False
        async def cancel(**kwargs):return None
        self.ui.picker=SimpleNamespace(get_directory_path=cancel)
        field,_=self.ui.folder_controls(str(self.root))
        asyncio.run(field.suffix_icon.on_tap(None))
        self.assertEqual(field.value,str(self.root))

    def test_existing_project_folder_cannot_be_changed(self):
        asyncio.run(self.ui.project_dialog(self.project))
        field=self.dialogs[-1].content.content.controls[-1].controls[0]
        self.assertTrue(field.read_only)
        self.assertTrue(field.suffix_icon.disabled)
        self.assertIsNone(field.suffix_icon.on_tap)
        self.assertEqual(field.value,self.project['directory'])

    def test_web_browser_navigates_creates_and_selects_folder(self):
        from biomolexplorer.ui.folder_browser import FolderBrowser
        selected=[]
        browser=FolderBrowser(self.ui,selected.append,str(self.root))
        asyncio.run(browser.open())
        self.assertEqual(browser.current,self.root)
        browser.name.value='estudo com espaços'
        asyncio.run(browser.create())
        self.assertEqual(browser.current,self.root/'estudo com espaços')
        asyncio.run(browser.select())
        self.assertEqual(selected,[str(self.root/'estudo com espaços')])
        self.assertFalse(browser.dialog.open)

    def test_web_browser_rejects_nonempty_folder_and_invalid_names(self):
        from biomolexplorer.ui.folder_browser import FolderBrowser
        selected=[]
        browser=FolderBrowser(self.ui,selected.append,str(self.root))
        asyncio.run(browser.open())
        asyncio.run(browser.select())
        self.assertEqual(selected,[])
        self.assertIn('vazia',browser.error.value)
        for name in ('../escape','a/b','a\\b','..',''):
            browser.name.value=name
            asyncio.run(browser.create())
            self.assertIn('válido',browser.error.value)

    def test_web_browser_discards_selection_after_logout(self):
        from biomolexplorer.ui.folder_browser import FolderBrowser
        selected=[]
        empty=self.root/'empty';empty.mkdir()
        browser=FolderBrowser(self.ui,selected.append,str(empty))
        asyncio.run(browser.open())
        self.ui.token=None
        asyncio.run(browser.select())
        self.assertEqual(selected,[])


if __name__=='__main__':unittest.main()
