"""Language belongs to the UI session; scientific and user data stay unchanged."""
import asyncio
import copy
import unittest
from types import SimpleNamespace

import flet as ft
import flet.canvas as canvas
from biomolexplorer.catalog import LABELS,TITLES
from biomolexplorer.ui.app import WorkspaceUI
from biomolexplorer.ui.localization import LocalizedPage,Translator,verbatim


class Page:
    def __init__(self):
        self.controls=[];self.dialogs=[];self.width=1200;self.height=900;self.web=False
    def add(self,*controls):self.controls.extend(controls)
    def update(self,*controls):pass
    def show_dialog(self,dialog):dialog.open=True;self.dialogs.append(dialog)
    def pop_dialog(self):return self.dialogs.pop() if self.dialogs else None


def walk(root):
    if isinstance(root,(list,tuple)):
        for item in root:yield from walk(item)
        return
    if root is None or isinstance(root,(str,int,float,bool)):return
    yield root
    for name in LocalizedPage._children:yield from walk(getattr(root,name,None))


class LocalizationTests(unittest.TestCase):
    def test_stage_catalog_and_parameter_labels_have_english_translations(self):
        translate=Translator('en')
        for title,description,category in TITLES.values():
            self.assertNotEqual(translate(title),title)
            self.assertNotEqual(translate(description),description)
        for name,label in LABELS.items():
            if name not in ('morgan','maccs'):
                self.assertNotEqual(translate(label),label,(name,label))

    def test_dynamic_messages_preserve_identifiers_and_filenames(self):
        translate=Translator('en')
        self.assertEqual(translate('Convite para Novo projeto'),'Invitation to Novo projeto')
        self.assertEqual(translate('Baixar compounds com espaço.csv'),'Download compounds com espaço.csv')
        self.assertEqual(translate('3 arquivos · Página 2'),'3 files · Page 2')
        self.assertEqual(translate('Etapa atual: Meu alvo\nProcessando os dados…'),
                         'Current stage: Meu alvo\nProcessing data…')

    def test_file_validation_translates_contract_and_retains_filename(self):
        from biomolexplorer.input_validation import contract
        message='Novo projeto.csv: Selecione um arquivo CSV.\nPadrão esperado: '+contract('compounds')
        translated=Translator('en')(message)
        self.assertTrue(translated.startswith('Novo projeto.csv: Select a CSV file.'))
        self.assertIn('Expected format: UTF-8 CSV',translated)
        self.assertNotIn('Também aceitamos',translated)

    def test_controls_translate_and_updates_do_not_change_parameters_or_user_content(self):
        page=LocalizedPage(Page(),'en')
        field=ft.TextField(label='Nome do projeto',value='Novo projeto',error='Informe a pasta do projeto.')
        dropdown=ft.Dropdown(label='Como processar os arquivos selecionados?',value='individual',
            options=[ft.DropdownOption(key='individual',text='Processar individualmente')])
        name=verbatim(ft.Text('Novo projeto'))
        status=ft.Text('Aguardando')
        picture=canvas.Canvas(shapes=[canvas.Text(1,2,'Componente 2')])
        original=copy.deepcopy((field.value,dropdown.value,dropdown.options[0].key))
        page.add(ft.Column([field,dropdown,name,status,picture]))
        self.assertEqual(field.label,'Project name');self.assertEqual(field.error,'Enter the project folder.')
        self.assertEqual(name.value,'Novo projeto');self.assertEqual(status.value,'Queued')
        self.assertEqual(picture.shapes[0].value,'Component 2')
        self.assertEqual((field.value,dropdown.value,dropdown.options[0].key),original)
        status.value='Concluído';page.update()
        self.assertEqual(status.value,'Completed')
        page.set_language('pt');page.update()
        self.assertEqual(field.label,'Nome do projeto');self.assertEqual(status.value,'Concluído')
        self.assertEqual(dropdown.options[0].text,'Processar individualmente')

    def test_language_is_independent_between_simultaneous_sessions_and_dialogs(self):
        english,portuguese=LocalizedPage(Page(),'en'),LocalizedPage(Page(),'pt')
        a,b=ft.AlertDialog(title=ft.Text('Novo projeto')),ft.AlertDialog(title=ft.Text('Novo projeto'))
        english.show_dialog(a);portuguese.show_dialog(b)
        self.assertEqual(a.title.value,'New project');self.assertEqual(b.title.value,'Novo projeto')
        portuguese.set_language('en');portuguese.update()
        english.set_language('pt');english.update()
        self.assertEqual(a.title.value,'Novo projeto');self.assertEqual(b.title.value,'New project')

    def test_login_language_switch_keeps_registration_fields(self):
        page=Page();ui=WorkspaceUI(page,SimpleNamespace(),SimpleNamespace(),language='pt')
        ui.show_login(True)
        controls=list(walk(page.controls))
        for label,value in [('Seu nome','Novo projeto'),('E-mail','someone@example.org'),('Senha','secret-password'),('Confirme a senha','secret-password')]:
            next(c for c in controls if isinstance(c,ft.TextField) and c.label==label).value=value
        language=next(c for c in controls if isinstance(c,ft.PopupMenuButton))
        english=next(item for item in language.items if item.data=='en')
        asyncio.run(english.on_click(SimpleNamespace(control=english)))
        controls=list(walk(page.controls))
        for label,value in [('Your name','Novo projeto'),('E-mail','someone@example.org'),('Password','secret-password'),('Confirm password','secret-password')]:
            self.assertEqual(next(c for c in controls if isinstance(c,ft.TextField) and c.label==label).value,value)
        self.assertTrue(any(isinstance(c,ft.Button) and c.content=='Create account' for c in controls))
        self.assertEqual(ui.language,'en');self.assertIsNone(ui.token)

    def test_authentication_errors_are_translated_without_changing_backend(self):
        def login(*args):raise ValueError('E-mail ou senha inválidos.')
        page=Page();ui=WorkspaceUI(page,SimpleNamespace(login=login),SimpleNamespace(),language='en')
        async def call(function,*args,**kwargs):return function(*args,**kwargs)
        ui.call=call
        ui.show_login()
        button=next(c for c in walk(page.controls) if isinstance(c,ft.Button) and c.content=='Sign in to workspace')
        asyncio.run(button.on_click(None))
        self.assertTrue(any(isinstance(c,ft.Text) and c.value=='Invalid email or password.' for c in walk(page.controls)))


if __name__=='__main__':unittest.main()
