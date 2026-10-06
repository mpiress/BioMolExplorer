import copy
import unittest
from types import SimpleNamespace
try:
    from biomolexplorer.ui.guided import GuidedForm
    from biomolexplorer.ui.flow_canvas import FlowCanvas
except ImportError:
    GuidedForm=FlowCanvas=None
from biomolexplorer.catalog import TITLES,new_stage
from biomolexplorer.templates import validate_templates

@unittest.skipIf(GuidedForm is None,'Install the ui extra to test Flet controls')
class GuidedUITests(unittest.TestCase):
    def setUp(self):
        self.ui=SimpleNamespace(current={'pipeline':[new_stage(k) for k in TITLES]},selected=None,service=SimpleNamespace(dock6_path=None),
            page=SimpleNamespace(update=lambda:None),event=lambda *args:None,save_draft=None,preset_dialog=None,run_pipeline=None)
    def test_all_operations_construct_and_generate_valid_templates(self):
        FlowCanvas(self.ui,True).build()
        for stage in self.ui.current['pipeline']:
            with self.subTest(stage=stage['operation']):
                result=GuidedForm(self.ui,stage,[],True).read()
                validate_templates(result['templates'])
                self.assertEqual(stage['name'],result['name'])
    def test_complex_records_preserve_resolution(self):
        stage=new_stage('redocking'); stage['parameters']['pdb_codes']=[['1ABC','LIG',123,'A',1.7]]
        result=GuidedForm(self.ui,stage,[],True).read()
        self.assertEqual(result['parameters']['pdb_codes'],stage['parameters']['pdb_codes'])
    def test_retrieval_form_always_uses_file_selection(self):
        for operation in ('retrieve_compounds','retrieve_structures','retrieve_zinc'):
            stage=new_stage(operation)
            form=GuidedForm(self.ui,stage,[],True)
            self.assertFalse(hasattr(form,'process_all'))
            self.assertFalse(form.read()['process_all'])
            saved=form.read()
            saved['process_all']=True
            self.assertFalse(GuidedForm(self.ui,saved,[],False).read()['process_all'])
    def test_vina_single_reference_opens_as_a_record(self):
        stage=new_stage('docking_vina'); stage['parameters']['pdb_code']=['1ABC','LIG',123,'A']
        result=GuidedForm(self.ui,stage,[],True).read()
        self.assertEqual(result['parameters']['pdb_code'],[['1ABC','LIG',123,'A']])
    def test_ic50_only_selection_survives_layout(self):
        import json
        import flet as ft
        form=GuidedForm(self.ui,new_stage('retrieve_compounds'),[],True)
        def walk(control):
            yield control
            for child in getattr(control,'controls',[]) or []: yield from walk(child)
            child=getattr(control,'content',None)
            if isinstance(child,ft.Control): yield from walk(child)
        layout=form.layout()
        for control in walk(layout):
            if isinstance(control,ft.Checkbox) and control.label=='Ki': control.value=False
        result=form.read()
        self.assertEqual(json.loads(result['templates']['crawlers/bioactivity.json'])['standard_type__in'],['IC50'])
        self.assertGreaterEqual(layout.spacing,24)
    def test_login_and_registration_use_bundled_logos(self):
        from biomolexplorer.ui.app import WorkspaceUI
        from biomolexplorer.ui.branding import image_bytes
        ui=object.__new__(WorkspaceUI)
        ui.page=SimpleNamespace(width=1440,height=1000,controls=[],add=lambda c:None);ui.token=None
        ui.show_login();ui.show_login(True)
        self.assertTrue(image_bytes('logo.png').startswith(b'\x89PNG'))
        ui.page.width=390;ui.page.height=844;ui.show_login()
    def test_direct_retrieval_mode_shows_relevant_controls_and_preserves_query(self):
        stage=new_stage('retrieve_compounds')
        stage['parameters'].update(search_mode='substructure',search_term='F/C=C/F')
        form=GuidedForm(self.ui,stage,[],True);form.layout()
        self.assertFalse(form.field_controls['max_targets'].visible)
        self.assertFalse(form.filter_tiles['target'].visible)
        self.assertFalse(form.filter_tiles['bioactivity'].visible)
        self.assertTrue(form.filter_tiles['molecules'].visible)
        self.assertEqual(form.read()['parameters']['search_term'],'F/C=C/F')
        form.field_controls['search_mode'].value='similarity';form.sync_retrieval()
        self.assertTrue(form.field_controls['similarity_threshold'].visible)
        form.field_controls['search_mode'].value='target';form.sync_retrieval()
        self.assertTrue(form.filter_tiles['bioactivity'].visible)
        self.assertTrue(form.conditional_cells['max_targets'].visible)

    def test_default_filters_have_no_redundant_overrides(self):
        result=GuidedForm(self.ui,new_stage('retrieve_compounds'),[],True).read()
        self.assertEqual(result['templates'],{})
    def test_undo_restores_removed_block_and_connection(self):
        from biomolexplorer import flow
        import asyncio
        a,b=new_stage('retrieve_compounds'),new_stage('admet')
        flow.connect([a,b],a['id'],b['id'],'base_input_path')
        self.ui.current['pipeline']=[a,b]; self.ui.selected=a['id']
        editor=FlowCanvas(self.ui,True); original=copy.deepcopy(editor.stages)
        asyncio.run(editor.delete(None)); asyncio.run(editor.undo(None))
        self.assertEqual(editor.stages,original)
