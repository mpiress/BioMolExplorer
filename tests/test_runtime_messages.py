"""Run dialogs and worker messages keep the viewing session's language."""
import ast
import asyncio
import copy
import json
import os
from pathlib import Path
import string
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import AsyncMock,patch

import flet as ft
import test_localization as fixtures
from biomolexplorer.i18n import RuntimeMessages,runtime_catalog
from biomolexplorer.progress import report_progress,read_progress
from biomolexplorer.ui.app import WorkspaceUI
from biomolexplorer.ui.localization import LocalizedPage,Translator,catalogs
from biomolexplorer.ui.run_progress import RunProgress
from biomolexplorer.ui.runtime_messages import stage_log


class RuntimeMessageTests(unittest.TestCase):
    def run_fixture(self,language):
        ui=object.__new__(WorkspaceUI);ui.language=language;ui.token='session'
        ui.current={'id':'project','role':'editor'};ui.current_run='run'
        ui.page=LocalizedPage(fixtures.Page(),language)
        run=dict(id='run',project_id='project',status='running',created=100,updated=120,error=None,
            stages=[dict(id='s',name='Recuperar ChEMBL',operation='retrieve_compounds',status='running',
                started_at=105,configuration={'enabled':True},progress=dict(message='Worker legado',
                    message_key='chembl.target',message_args={'target':'Meu arquivo.csv'},completed=2,total=5))])
        ui.run_progress=RunProgress(ui,run);ui.run_progress.show()
        return ui,run

    def test_catalogs_have_identical_keys_and_interpolation_contracts(self):
        pt,en=runtime_catalog('pt'),runtime_catalog('en');self.assertEqual(set(pt),set(en))
        fields=lambda v:{field for _,field,_,_ in string.Formatter().parse(v) if field is not None}
        for key in pt:self.assertEqual(fields(pt[key]),fields(en[key]),key)
        with self.assertRaises(KeyError):RuntimeMessages('pt')('missing.key')

    def test_portuguese_catalog_is_explicit_and_does_not_reverse_english_at_runtime(self):
        catalogs.cache_clear()
        original=Path.read_text
        def read(path,*args,**kwargs):
            if path.name=='en.json':raise AssertionError('English catalog used for Portuguese')
            return original(path,*args,**kwargs)
        with patch.object(Path,'read_text',read):
            self.assertEqual(Translator('pt')('Your pipeline is running'),'Seu pipeline está em execução')
        catalogs.cache_clear()

    def test_worker_progress_stores_key_arguments_and_legacy_message(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'progress.json'
            with patch.dict(os.environ,BIOMOL_PROGRESS_FILE=str(path)):
                report_progress(None,2,5,key='chembl.target',arguments={'target':'Novo projeto.csv'})
            phase=read_progress(path)
            self.assertEqual(phase['message_key'],'chembl.target')
            self.assertIn('Novo projeto.csv',phase['message'])
            self.assertIn('Consultando',phase['message'])
            self.assertEqual(phase['completed'],2)
            self.assertEqual(RuntimeMessages('en')(phase['message_key'],**phase['message_args']),
                'Querying target Novo projeto.csv in ChEMBL…')

    def test_large_progress_arguments_stay_within_the_reader_boundary(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'progress.json'
            with patch.dict(os.environ,BIOMOL_PROGRESS_FILE=str(path)):
                report_progress(key='chembl.target',arguments={'target':'C'*12000})
            self.assertIsNotNone(read_progress(path))
            self.assertLess(len(path.read_text()),8192)

    def test_invalid_monitoring_metadata_cannot_fail_the_scientific_work(self):
        with tempfile.TemporaryDirectory() as folder,patch('biomolexplorer.diagnostics.get_logger') as logger:
            with patch.dict(os.environ,BIOMOL_PROGRESS_FILE=str(Path(folder)/'progress.json')):
                report_progress(key='missing.key')
            logger.return_value.warning.assert_called_once()

    def test_dock6_preparation_step_is_localized_and_compound_code_is_preserved(self):
        ui,run=self.run_fixture('en')
        run['stages'][0]['progress'].update(message_key='preparation.step',
            message_args=dict(step='Minimização DOCK6',code='Novo composto'))
        ui.run_progress.update(run)
        self.assertIn('DOCK6 minimization: Novo composto',ui.run_progress.message.value)

    def test_popup_and_accessible_progress_do_not_flip_during_polling(self):
        for language,heading,processed in [('pt','Seu pipeline está em execução','Registros processados'),
                                          ('en','Your pipeline is running','Records processed')]:
            ui,run=self.run_fixture(language);p=ui.run_progress
            for _ in range(20):
                p.update(run,now=130);ui.page.update()
                self.assertEqual(p.heading.value,heading)
                self.assertIn(processed,p.message.value)
                self.assertEqual(p.message_semantics.label,p.message.value)
                self.assertIn('Meu arquivo.csv',p.message.value)
                self.assertEqual(p.messages.language,language)
            self.assertEqual(run['stages'][0]['progress']['message'],'Worker legado')

    def test_language_does_not_leak_across_sessions_and_a_run_keeps_its_snapshot(self):
        pt,run=self.run_fixture('pt');en,other=self.run_fixture('en')
        pt.language='en'  # A background assignment must not change an existing run.
        pt.run_progress.update(run);en.run_progress.update(other)
        self.assertEqual(pt.run_progress.heading.value,'Seu pipeline está em execução')
        self.assertEqual(en.run_progress.heading.value,'Your pipeline is running')

    def test_all_run_states_are_loaded_from_the_selected_catalog(self):
        for lang in ('pt','en'):
            ui,run=self.run_fixture(lang)
            for status in ('queued','awaiting_input','succeeded','failed','cancelled','interrupted'):
                run['status']=status;run['stages'][0]['status']=status
                ui.run_progress.update(run);ui.page.update()
                self.assertEqual(ui.run_progress.heading.value,RuntimeMessages(lang)('run.heading.'+status))

    def test_structured_progress_preserves_user_titles_and_filenames(self):
        ui,run=self.run_fixture('en');run['stages'][0]['name']='Seu pipeline está em execução'
        run['stages'][0]['progress'].update(message_key='graphs.analysis',message_args=dict(index=1,total=2,file='Executando.csv'))
        ui.run_progress.update(run);ui.page.update()
        self.assertIn('Current stage: Seu pipeline está em execução',ui.run_progress.message.value)
        self.assertIn('Graph analysis 1/2 · Executando.csv',ui.run_progress.message.value)

    def test_legacy_progress_is_localized_without_rewriting_saved_data(self):
        ui,run=self.run_fixture('en');phase=dict(message='Baixando estruturas PDB…')
        run['stages'][0]['progress']=phase;ui.run_progress.update(run)
        self.assertIn('Downloading PDB structures',ui.run_progress.message.value)
        self.assertEqual(phase,dict(message='Baixando estruturas PDB…'))

    def test_language_selector_cannot_change_an_active_run(self):
        ui,_=self.run_fixture('pt');bar,selector,_=ui.top_bar()
        self.assertTrue(selector.disabled)
        english=next(item for item in selector.items if item.data=='en')
        asyncio.run(english.on_click(SimpleNamespace(control=english)))
        self.assertEqual(ui.language,'pt');self.assertEqual(ui.page.translator.language,'pt')

    def test_event_log_uses_locale_and_preserves_context_and_native_tool_output(self):
        prefix='2026-10-26T00:00:00 INFO worker event=worker.succeeded [stage_id=Novo projeto] '
        original=prefix+'Scientific worker completed\nVINA RESULT: -7.5  0.000\n'
        rendered=stage_log(original,'pt')
        self.assertEqual(rendered,prefix+'Executor científico concluído\nVINA RESULT: -7.5  0.000\n')
        self.assertEqual(stage_log(original,'en'),original)

    def test_all_literal_progress_ids_exist_in_both_catalogs(self):
        root=Path(__file__).parents[1]/'src'
        for file in root.rglob('*.py'):
            for node in ast.walk(ast.parse(file.read_text())):
                if not isinstance(node,ast.Call) or not isinstance(node.func,ast.Name) or node.func.id!='report_progress':continue
                key=next((k.value for k in node.keywords if k.arg=='key'),node.args[0] if node.args else None)
                if isinstance(key,ast.Constant) and isinstance(key.value,str):
                    for language in ('pt','en'):self.assertIn(key.value,runtime_catalog(language),str(file))


class LogDialogTests(unittest.IsolatedAsyncioTestCase):
    async def test_log_dialog_localizes_view_and_downloads_original_bytes(self):
        ui=object.__new__(WorkspaceUI);ui.language='pt';ui.token='t';ui.current={'id':'p'}
        ui.page=LocalizedPage(fixtures.Page(),'pt')
        original=b'INFO backend event=worker.succeeded [stage_id=s] Scientific worker completed\n'
        ui.store=SimpleNamespace(read_file=lambda *args:original);ui.picker=SimpleNamespace(save_file=AsyncMock())
        async def call(fn,*args):return fn(*args)
        ui.call=call
        await ui.open_stage_log('p','/some/stage.log')
        dialog=ui.page._page.dialogs[-1];self.assertEqual(dialog.title.value,'Log da etapa')
        self.assertIn('Executor científico concluído',dialog.content.content.value)
        await dialog.actions[0].on_click(None)
        ui.picker.save_file.assert_awaited_once_with(file_name='stage.log',src_bytes=original)
