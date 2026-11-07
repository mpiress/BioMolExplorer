"""Dictionary coverage, splash lifecycle and adding a locale without UI code."""
import json
import shutil
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import flet as ft
from biomolexplorer import i18n
from biomolexplorer.ui import localization
from biomolexplorer.ui.configuration_loading import ProcessingSplash
from biomolexplorer.ui.processing import MESSAGES,DEFAULT
from test_localization import Page
from scripts.audit_localization import audit


def clear_catalogs():
    i18n.locales.cache_clear();i18n.dictionary.cache_clear();i18n.localized_dictionary.cache_clear()
    localization.catalogs.cache_clear();localization.default_stage_names.cache_clear()


class LocalizationCatalogTests(unittest.TestCase):
    def test_source_and_all_dictionary_namespaces_pass_audit(self):
        self.assertEqual(audit(Path(__file__).resolve().parents[1])['errors'],[])

    def test_all_processing_messages_translate_in_reused_modal_and_overlay(self):
        page=localization.LocalizedPage(Page(),'en')
        for modal in (True,False):
            splash=ProcessingSplash(page)
            for index,(title,description) in enumerate(set(MESSAGES.values())|{DEFAULT}):
                splash.show(index,title,description,modal=modal)
                self.assertNotEqual(splash.title.value,title,title)
                self.assertNotEqual(splash.description.value,description,description)
                self.assertEqual(splash.ring.semantics_label,splash.title.value)
                if index:splash.hide(index-1)
            page.set_language('pt');page.update()
            self.assertEqual(splash.title.value,title)
            self.assertEqual(splash.description.value,description)
            splash.close();page.set_language('en')

    def test_tooltip_objects_and_overlay_updates_follow_session_language(self):
        raw=Page();page=localization.LocalizedPage(raw,'en')
        button=ft.IconButton(icon=ft.Icons.INFO,tooltip=ft.Tooltip(message='Carregando projeto…'))
        raw.overlay=[ft.Text('Preparando execução…')];page.add(button);page.update()
        self.assertEqual(button.tooltip.message,'Loading project…')
        self.assertEqual(raw.overlay[0].value,'Preparing execution…')
        raw.overlay[0].value='Salvando projeto…';page.update()
        self.assertEqual(raw.overlay[0].value,'Saving project…')
        page.set_language('pt');page.update()
        self.assertEqual(button.tooltip.message,'Carregando projeto…')
        self.assertEqual(raw.overlay[0].value,'Salvando projeto…')

    def test_new_locale_loads_from_dictionaries_and_preserves_reordered_arguments(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for source in i18n.ROOT.glob('*.json'):shutil.copyfile(source,root/source.name)
            manifest=json.loads((root/'locales.json').read_text())
            manifest['zz']={'name':'Test locale','language_code':'zz'}
            (root/'locales.json').write_text(json.dumps(manifest))
            (root/'zz.json').write_text(json.dumps({'Salvar':'ZZ save','Processando dados…':'ZZ processing'}))
            (root/'patterns_zz.json').write_text(json.dumps({'Arquivo {0} ({1}): {2}':'ZZ {2} / {1} / {0}'}))
            (root/'runtime_zz.json').write_text(json.dumps({'log.details':'ZZ details'}))
            (root/'viewer_zz.json').write_text(json.dumps({'viewer':'ZZ structure'}))
            with patch.object(i18n,'ROOT',root):
                clear_catalogs()
                try:
                    translate=localization.Translator('zz')
                    self.assertEqual(translate('Salvar'),'ZZ save')
                    self.assertEqual(translate('Save'),'ZZ save')
                    self.assertEqual(translate('Arquivo Novo projeto.csv (A): C42'),'ZZ C42 / A / Novo projeto.csv')
                    self.assertEqual(translate('Carregando projeto…'),'Loading project…')
                    self.assertEqual(i18n.RuntimeMessages('zz')('log.details'),'ZZ details')
                    from biomolexplorer.pdb_view import viewer_document
                    self.assertIn('ZZ structure',viewer_document('Novo projeto.pdb','zz'))
                    page=localization.LocalizedPage(Page(),'zz')
                    self.assertEqual(page.locale_configuration.current_locale.language_code,'zz')
                    from biomolexplorer.ui.app import WorkspaceUI
                    ui=WorkspaceUI(Page(),SimpleNamespace(),SimpleNamespace(),language='zz')
                    self.addCleanup(ui.structure_viewer.close)
                    _,choice,_=ui.top_bar()
                    self.assertIn('zz',[item.data for item in choice.items])
                finally:clear_catalogs()
        clear_catalogs()


if __name__=='__main__':unittest.main()
