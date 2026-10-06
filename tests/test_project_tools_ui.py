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
        destination.mkdir(parents=True)
        field.value=str(destination/'Novo estudo')
        asyncio.run(field.suffix_icon.on_tap(None))
        asyncio.run(self.dialogs[-1].actions[-1].on_click(None))
        asyncio.run(dialog.actions[-1].on_click(None))
        self.assertEqual(self.ui.current['directory'],str(destination/'Novo estudo'))
        self.assertTrue((destination/'Novo estudo'/'project.json').is_file())
        self.assertFalse((self.store.projects_root/self.ui.current['id']).exists())
        self.assertFalse(dialog.open)

    def test_name_is_required_before_opening_folder_browser(self):
        asyncio.run(self.ui.project_dialog());dialog=self.dialogs[-1]
        field=dialog.content.content.controls[-1].controls[0]
        with self.assertRaisesRegex(ValueError,'nome de projeto válido'):
            asyncio.run(field.suffix_icon.on_tap(None))
        self.assertEqual(self.dialogs,[dialog])

    def test_replacement_requires_confirmation_and_waits_for_save(self):
        base=self.root/'research';base.mkdir()
        old=self.store.create_project(self.token,'Existing',directory=base/'Study')
        (base/'Study'/'data.txt').write_text('old data')
        (base/'unrelated.txt').write_text('keep')
        asyncio.run(self.ui.project_dialog());dialog=self.dialogs[-1]
        controls=dialog.content.content.controls;controls[0].value='Study'
        field=controls[-1].controls[0];field.value=str(base/'Study')
        asyncio.run(field.suffix_icon.on_tap(None))
        asyncio.run(self.dialogs[-1].actions[-1].on_click(None))
        confirmation=self.dialogs[-1]
        self.assertIn('Substituir',confirmation.title.value)
        self.assertEqual(self.store.project(self.token,old['id'])['name'],'Existing')
        asyncio.run(confirmation.actions[-1].on_click(None))
        self.assertTrue((base/'Study'/'data.txt').exists())
        asyncio.run(dialog.actions[-1].on_click(None))
        self.assertFalse((base/'Study'/'data.txt').exists())
        self.assertEqual((base/'unrelated.txt').read_text(),'keep')
        self.assertEqual(self.ui.current['directory'],str(base/'Study'))

    def test_deletion_confirmation_removes_the_associated_folder(self):
        folder=Path(self.project['directory']);(folder/'results.txt').write_text('data')
        self.ui.page.pop_dialog=lambda:self.dialogs[-1].__setattr__('open',False)
        async def workspace():self.ui.current=None
        self.ui.show_workspace=workspace
        asyncio.run(self.ui.delete_project_dialog(self.project['id']));dialog=self.dialogs[-1]
        contents=dialog.content.content.controls
        self.assertIn('permanentemente',contents[1].value)
        self.assertEqual(contents[2].value,str(folder))
        self.assertTrue(folder.exists())
        asyncio.run(dialog.actions[-1].on_click(None))
        self.assertFalse(folder.exists())
        self.assertIn('permanentemente',self.messages[-1])

    def test_renaming_invalidates_selected_folder(self):
        asyncio.run(self.ui.project_dialog());controls=self.dialogs[-1].content.content.controls
        controls[-1].controls[0].value=str(self.root/'Study')
        controls[0].value='Different';controls[0].on_change(None)
        self.assertEqual(controls[-1].controls[0].value,'')

    def test_desktop_folder_browser_uses_the_selected_folder(self):
        destination=self.root/'pasta escolhida'
        self.ui.page.web=False
        destination.mkdir()
        field,folder=self.ui.folder_controls(str(destination))
        field.error='Informe a pasta do projeto.'
        self.assertTrue(field.read_only)
        asyncio.run(field.suffix_icon.on_tap(None))
        asyncio.run(self.dialogs[-1].actions[-1].on_click(None))
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

    def test_folder_browser_hides_dot_folders_and_creates_under_current_directory(self):
        from biomolexplorer.ui.folder_browser import FolderBrowser
        base=self.root/'Downloads';base.mkdir();(base/'.hidden').mkdir();(base/'visible').mkdir()
        selected=[];browser=FolderBrowser(self.ui,selected.append,str(base))
        asyncio.run(browser.open())
        self.assertEqual([c.content.value for c in browser.entries.controls],['visible'])
        asyncio.run(browser.prompt_create());dialog=self.dialogs[-1]
        self.assertEqual(dialog.content.controls[0].value,str(base))
        dialog.content.controls[1].value='teste'
        asyncio.run(dialog.actions[-1].on_click(None))
        self.assertTrue((base/'teste').is_dir());self.assertEqual(selected,[str(base/'teste')])
        self.assertFalse(browser.dialog.open);self.assertFalse(dialog.open)

    def test_project_palette_saves_selected_color_without_tags_field(self):
        from biomolexplorer.workspace import COLORS
        asyncio.run(self.ui.project_dialog())
        dialog=self.dialogs[-1];controls=dialog.content.content.controls
        self.assertEqual(len(controls),4)
        self.assertFalse(any(getattr(c,'label',None)=='Tags' for c in controls))
        swatches=controls[2].controls[1].controls
        self.assertEqual(len(swatches),len(COLORS))
        swatches[3].on_click(None)
        destination=self.root/'palette-project';destination.mkdir()
        controls[0].value='Palette'
        field=controls[-1].controls[0];field.value=str(destination/'Palette')
        asyncio.run(field.suffix_icon.on_tap(None))
        asyncio.run(self.dialogs[-1].actions[-1].on_click(None))
        asyncio.run(dialog.actions[-1].on_click(None))
        self.assertEqual(self.ui.current['color'],COLORS[3]);self.assertEqual(self.ui.current['tags'],[])

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
