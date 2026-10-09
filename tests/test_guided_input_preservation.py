"""The visual editor preserves configured inputs until their source is changed."""
from types import SimpleNamespace
import unittest

try:
    import flet as ft
    from biomolexplorer.ui.guided import GuidedForm
except ImportError:
    GuidedForm = None

from biomolexplorer.catalog import LABELS, PATH_FIELDS, TITLES, new_stage, operation_fields


@unittest.skipIf(GuidedForm is None, 'Install the ui extra to test Flet controls')
class GuidedInputPreservationTests(unittest.TestCase):
    def test_zinc_download_workers_default_and_edit_survive_form(self):
        stage=new_stage('retrieve_zinc')
        self.assertEqual(stage['parameters']['download_workers'],4)
        form=GuidedForm(self.ui([stage]),stage,[],True)
        field=form.field_controls['download_workers']
        self.assertEqual(field.label,'Downloads simultâneos')
        self.assertIn('1 a 16',field.helper)
        field.value='8'
        self.assertEqual(form.read()['parameters']['download_workers'],8)

    def controls(self, form):
        def walk(control):
            yield control
            for child in getattr(control,'controls',[]) or []:yield from walk(child)
        return [c for root in form.controls for c in walk(root)]

    def ui(self, stages):
        return SimpleNamespace(current={'pipeline':stages},
            service=SimpleNamespace(dock6_path=None), page=SimpleNamespace(update=lambda:None))

    def test_all_direct_inputs_survive_visual_configuration(self):
        for operation in TITLES:
            if operation == 'import_results':
                continue
            stage = new_stage(operation)
            if operation=='redocking':stage['parameters']['pdb_codes']=[['1ABC','LIG',123,'A']]
            paths = {f['name']: '/tmp/project/prepared data/' + f['name']
                for f in operation_fields(operation) if f['name'] in PATH_FIELDS}
            if not paths:
                continue
            with self.subTest(operation=operation):
                stage['parameters'].update(paths)
                for writable in (True, False):
                    form = GuidedForm(self.ui([stage]),stage,[],writable)
                    form.layout()
                    result = form.read()
                    self.assertEqual({key:result['parameters'][key] for key in paths},paths)
                    self.assertEqual(result['bindings'],{})
                    self.assertEqual({key:stage['parameters'][key] for key in paths},paths)

    def test_selected_files_and_processing_mode_survive_configuration_in_all_blocks(self):
        from biomolexplorer.flow import INPUTS
        for operation,fields in INPUTS.items():
            stage=new_stage(operation)
            if operation=='redocking':stage['parameters']['pdb_codes']=[['1ABC','LIG',123,'A']]
            origins=[]
            expected={}
            for field,kinds in fields.items():
                source=new_stage('similarity' if operation=='graphs' else 'import_results')
                if source['operation']=='import_results':source['parameters']['kind']=next(iter(sorted(kinds)))
                origins.append(source)
                suffix='.dockprep.pdbqt' if operation in ('docking_vina','docking_dock6','prepare_structures') and field=='base_input_path' else '.csv'
                expected[field]={'sources':[{'stage':source['id'],'selector':'first/selected'+suffix},
                                           {'stage':source['id'],'selector':'second/selected'+suffix}]}
            stage['bindings']=expected
            for mode in ('individual','merge'):
                with self.subTest(operation=operation,mode=mode):
                    stage['input_processing']=mode
                    form=GuidedForm(self.ui(origins+[stage]),stage,[],True)
                    result=form.read()
                    self.assertEqual(result['bindings'],expected)
                    self.assertEqual(result['input_processing'],mode)

    def test_new_connection_replaces_configured_folder(self):
        source, stage = new_stage('retrieve_compounds'), new_stage('admet')
        stage['parameters']['base_input_path'] = '/tmp/project/prepared'
        form = GuidedForm(self.ui([source,stage]),stage,[],True)
        dropdown = next(c for c in self.controls(form) if isinstance(c,ft.Dropdown)
            and c.label == LABELS['base_input_path'])
        self.assertEqual(dropdown.value,'configured-path')
        dropdown.value = 'stage:' + source['id']
        dropdown.on_select(None)
        result = form.read()
        self.assertNotIn('base_input_path',result['parameters'])
        self.assertEqual(result['bindings']['base_input_path'],{'stage':source['id'],'selector':'auto'})

    def test_clearing_source_removes_the_configured_folder(self):
        stage = new_stage('admet')
        stage['parameters']['base_input_path'] = '/tmp/project/prepared'
        form = GuidedForm(self.ui([stage]),stage,[],True)
        dropdown = next(c for c in self.controls(form) if isinstance(c,ft.Dropdown)
            and c.label == LABELS['base_input_path'])
        dropdown.value = ''
        result = form.read()
        self.assertNotIn('base_input_path',result['parameters'])
        self.assertNotIn('base_input_path',result['bindings'])

    def test_changing_consensus_vina_removes_the_old_derived_dependency(self):
        from biomolexplorer.flow import connect
        original, replacement = new_stage('docking_vina'), new_stage('docking_vina')
        consensus = new_stage('consensus')
        connect([original,replacement,consensus],original['id'],consensus['id'],'base_vina_path')
        form = GuidedForm(self.ui([original,replacement,consensus]),consensus,[],True)
        self.assertFalse(any(isinstance(c,ft.Dropdown) and c.label==LABELS['base_input_path']
            for c in self.controls(form)))
        dropdown = next(c for c in self.controls(form) if isinstance(c,ft.Dropdown)
            and c.label==LABELS['base_vina_path'])
        dropdown.value = 'stage:' + replacement['id']
        dropdown.on_select(None)
        result = form.read()
        self.assertNotIn('base_input_path',result['bindings'])
        self.assertEqual(result['bindings']['base_vina_path']['stage'],replacement['id'])
        self.assertFalse(any(b.get('stage')==original['id'] for b in result['bindings'].values()))

    def test_new_consensus_uses_only_the_two_docking_inputs(self):
        consensus = new_stage('consensus')
        form = GuidedForm(self.ui([consensus]),consensus,[],True)
        labels = [c.label for c in self.controls(form) if isinstance(c,ft.Dropdown) and c.label in (LABELS['base_input_path'],LABELS['base_vina_path'],LABELS['base_dock6_path'])]
        self.assertNotIn(LABELS['base_input_path'],labels)
        self.assertIn(LABELS['base_vina_path'],labels)
        self.assertIn(LABELS['base_dock6_path'],labels)


if __name__ == '__main__':
    unittest.main()
