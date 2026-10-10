"""Typed import table, immediate validation and mixed-file pipeline handoff."""
import asyncio
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import AsyncMock

import test_pipeline_execution as execution
from biomolexplorer.catalog import new_stage
from biomolexplorer.flow import output_types,compatible
from biomolexplorer.pipeline import PipelineService,validate_pipeline
from biomolexplorer.ui.guided import GuidedForm

ATOM='ATOM      1  C   ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n'


class ImportFilesTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    upload=execution.PipelineExecutionTests.upload

    def form(self,stage,assets=None,writable=True):
        async def call(function,*args):return function(*args)
        async def guard(function):await function()
        ui=SimpleNamespace(token=self.token,store=self.store,current={'id':self.project_id,'pipeline':[stage]},
            page=SimpleNamespace(update=lambda:None),call=call,guard=guard)
        return GuidedForm(ui,stage,assets if assets is not None else self.store.assets(self.token,self.project_id),writable)

    def test_add_disk_files_validate_table_remove_and_reopen(self):
        csv=self.upload();pdb=self.upload('1ABC.pdb',ATOM,'structures')
        stage=new_stage('import_results');form=self.form(stage);panel=form.import_files
        form.ui.pick_uploads=AsyncMock(return_value=[csv,pdb])
        asyncio.run(panel.upload(None))
        self.assertEqual(len(panel.table.rows),2)
        self.assertEqual([c.label.value for c in panel.table.columns],['Tipo','Arquivo','Ações'])
        saved=form.read()
        self.assertEqual(saved['parameters']['asset_types'],{csv:'compounds',pdb:'structures'})
        self.assertEqual(output_types(saved),{'compounds','structures'})
        self.assertTrue(compatible(saved,new_stage('admet'),'base_input_path'))
        self.assertTrue(compatible(saved,new_stage('prepare_structures'),'base_input_path'))
        reopened=self.form(saved).import_files
        reopened.remove(pdb)
        self.assertEqual(list(reopened.entries),[csv])
        self.assertTrue(self.store.asset_path(self.token,self.project_id,pdb).exists())

    def test_remove_table_action_updates_saved_selection_without_deleting_files(self):
        first=self.upload('first.csv');second=self.upload('second.csv')
        stage=new_stage('import_results');stage['parameters']['asset_ids']=[first,second]
        form=self.form(stage);panel=form.import_files
        button=panel.table.rows[0].cells[2].content.content
        self.assertEqual(button.content,'Remover')
        button.on_click(None)
        saved=form.read()
        self.assertEqual(saved['parameters']['asset_ids'],[second])
        self.assertEqual(saved['parameters']['asset_types'],{second:'compounds'})
        self.assertTrue(self.store.asset_path(self.token,self.project_id,first).exists())
        reopened=self.form(saved).import_files
        self.assertEqual(len(reopened.table.rows),1)
        reopened.table.rows[0].cells[2].content.content.on_click(None)
        self.assertTrue(reopened.empty.visible)
        self.assertFalse(reopened.table_view.visible)

    def test_empty_import_configuration_explains_how_to_add_files(self):
        form=self.form(new_stage('import_results'));panel=form.import_files
        self.assertTrue(panel.empty.visible)
        self.assertFalse(panel.table_view.visible)
        self.assertIn('Adicione arquivos',panel.empty.value)
        form.layout()

    def test_type_change_rejects_wrong_contents_and_restores_selection(self):
        csv=self.upload();stage=new_stage('import_results');stage['parameters']['asset_ids']=[csv]
        panel=self.form(stage).import_files;control=panel.type_controls[csv];control.value='structures'
        with self.assertRaisesRegex(ValueError,'Padrão esperado'):
            asyncio.run(control.on_select(None))
        self.assertEqual(control.value,'compounds')
        self.assertEqual(panel.entries[csv]['kind'],'compounds')

    def test_mixed_import_normalizes_compounds_and_preserves_receptor_layout(self):
        csv=self.upload('zinc.csv','zinc_id,smiles\nZINC1,CCO\n')
        pdb=self.upload('1ABC.pdb',ATOM,'structures')
        stage=new_stage('import_results');stage['parameters'].update(kind='other',asset_ids=[csv,pdb],
            asset_types={csv:'compounds',pdb:'structures'})
        service=PipelineService(self.store)
        try:
            output=self.store.project_dir(self.project_id)/'mixed'
            paths=service._import(self.project_id,stage['parameters'],output)
            self.assertIn('molecule_chembl_id', (output/'zinc.csv').read_text())
            self.assertIn('ZINC1',(output/'zinc.csv').read_text())
            self.assertTrue((output/'MeuAlvo/1ABC.pdb').exists())
            receptor,_,_=service._materialize_inputs(self.project_id,new_stage('prepare_structures'),'base_input_path',
                [{'stage':stage['id']}],{stage['id']:paths})
            self.assertTrue((receptor/'MeuAlvo/1ABC.pdb').exists())
        finally:service.close()

    def test_types_must_cover_exactly_selected_files(self):
        stage=new_stage('import_results');stage['parameters'].update(asset_ids=['a'*32],asset_types={'b'*32:'compounds'})
        with self.assertRaisesRegex(ValueError,'cada arquivo'):validate_pipeline([stage])
        stage['parameters']['asset_types']={'a'*32:['compounds']}
        with self.assertRaises(ValueError):validate_pipeline([stage])

    def test_upload_with_declared_type_rejects_invalid_file(self):
        with self.assertRaisesRegex(ValueError,'Faltam as colunas'):
            self.upload('invalid.csv','wrong,header\n1,2\n','compounds')
        self.assertEqual(self.store.assets(self.token,self.project_id),[])

    def test_readonly_table_cannot_remove_files(self):
        asset=self.upload();stage=new_stage('import_results');stage['parameters']['asset_ids']=[asset]
        panel=self.form(stage,writable=False).import_files
        panel.remove(asset)
        self.assertIn(asset,panel.entries)
        self.assertTrue(panel.type_controls[asset].disabled)
        self.assertTrue(panel.table.rows[0].cells[2].content.content.disabled)

    def test_raw_and_prepared_receptor_choices_in_mixed_import_keep_settings_in_sync(self):
        origin=new_stage('import_results');origin['parameters'].update(kind='other',asset_ids=['a'*32,'b'*32],
            asset_types={'a'*32:'structures','b'*32:'prepared_structures'})
        stage=new_stage('prepare_structures')
        stage['bindings']={'base_input_path':{'stage':origin['id'],'selector':'1ABC.pdb'}}
        form=self.form(stage,[]);form.ui.current['pipeline'].insert(0,origin)
        form.ui.artifact_choices={origin['id']:{'1ABC.pdb','1ABC_A.dockprep.pdbqt','1ABC_LIG_1A.lig.pdbqt'}}
        # Rebuild after loading producer metadata, as the application does.
        form=GuidedForm(form.ui,stage,[],True);editor=form.input_editors['base_input_path']
        self.assertEqual({o.key for o in editor.rows[0][1].options},{'auto','1ABC.pdb','1ABC_A.dockprep.pdbqt'})
        self.assertFalse(form.preparation_settings.receptor_inactive)
        editor.rows[0][1].value='1ABC_A.dockprep.pdbqt';editor.changed()
        self.assertTrue(form.preparation_settings.receptor_inactive)
