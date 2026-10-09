"""Dependent controls block edits without discarding configured values."""
import json
import unittest
from types import SimpleNamespace

from biomolexplorer.catalog import new_stage
try:
    import flet as ft
    from biomolexplorer.ui.guided import GuidedForm
except ImportError:
    ft=GuidedForm=None


@unittest.skipIf(GuidedForm is None,'Install the ui extra to test Flet controls')
class GuidedDependenciesTests(unittest.TestCase):
    def form(self,operation='retrieve_compounds',writable=True,stage=None):
        ui=SimpleNamespace(page=SimpleNamespace(update=lambda:None),current={'pipeline':[]},
                           service=SimpleNamespace(dock6_path=None))
        return GuidedForm(ui,stage or new_stage(operation),[],writable)

    def toggle(self,control,value):
        control.value=value
        control.on_change(SimpleNamespace(control=control))

    def test_molecule_types_are_closed_and_roundtrip_lowercase(self):
        stage=new_stage('retrieve_compounds')
        stage['parameters']['chembl_filters']={'molecules':{'molecule_type':'Small molecule'},
                                               'similars':{'molecule_type':'Protein'}}
        form=self.form(stage=stage)
        for group,value in (('molecules','small molecule'),('similarmols','protein')):
            control=form.template_controls[f'crawlers/{group}.json']['molecule_type']
            self.assertIsInstance(control,ft.Dropdown)
            self.assertFalse(control.editable)
            self.assertEqual(control.value,value)
            self.assertIn('antibody drug conjugate',{o.key for o in control.options})
            self.assertEqual(json.loads(form.read()['templates'][f'crawlers/{group}.json'])['molecule_type'],value)
            control.value=''
            self.assertNotIn('molecule_type',json.loads(form.template_readers[f'crawlers/{group}.json']()))

    def test_invalid_legacy_molecule_type_is_not_added_as_an_option(self):
        stage=new_stage('retrieve_compounds')
        stage['parameters']['chembl_filters']={'molecules':{'molecule_type':'typo'}}
        with self.assertRaisesRegex(ValueError,'tipo de molécula válido'):self.form(stage=stage)

    def test_chembl_form_has_no_pubchem_expansion_settings(self):
        form=self.form();form.layout()
        self.assertNotIn('include_pubchem',form.field_controls)
        self.assertNotIn('pubchem_threshold',form.field_controls)
        self.assertFalse(form.read()['parameters']['include_pubchem'])

    def test_chembl_toggle_blocks_filters_after_layout_and_preserves_selection(self):
        form=self.form();form.layout()
        control=form.template_controls['crawlers/similarmols.json']['molecule_type']
        control.value='protein'
        self.toggle(form.field_controls['expand_chembl'],False)
        self.assertTrue(form.filter_tiles['similarmols'].visible)
        self.assertTrue(control.disabled)
        reopened=self.form(stage=form.read())
        self.assertEqual(reopened.template_controls['crawlers/similarmols.json']['molecule_type'].value,'protein')
        self.toggle(form.field_controls['expand_chembl'],True)
        self.assertFalse(control.disabled)
        form.field_controls['search_mode'].value='molecule_id'
        form.field_controls['search_mode'].on_select(SimpleNamespace(control=form.field_controls['search_mode']))
        self.assertTrue(control.disabled)
        self.assertFalse(form.filter_cells['similarmols'].visible)

    def test_readonly_controls_stay_disabled_when_switches_change(self):
        form=self.form(writable=False)
        self.toggle(form.field_controls['expand_chembl'],True)
        self.assertTrue(form.template_controls['crawlers/similarmols.json']['molecule_type'].disabled)

    def test_automatic_complex_detection_blocks_manual_entry(self):
        form=self.form('docking_vina')
        tile=next(c for c in form.controls if isinstance(c,ft.ExpansionTile) and any(isinstance(x,ft.Checkbox) for x in c.controls))
        automatic,records,append=tile.controls
        self.assertTrue(records.disabled);self.assertTrue(append.disabled)
        self.toggle(automatic,False)
        self.assertFalse(records.disabled);self.assertFalse(append.disabled)

    def test_nested_dock6_dependencies_preserve_edits_and_templates(self):
        form=self.form('docking_dock6')
        controls=form.template_controls['dock6/docking.template']
        controls['internal_energy_cutoff'].value='91'
        self.toggle(controls['use_internal_energy'],False)
        self.assertTrue(controls['internal_energy_cutoff'].disabled)
        self.assertIn('91',form.read()['templates']['dock6/docking.template'])
        self.toggle(controls['minimize_ligand'],False)
        self.assertTrue(controls['minimize_flexible_growth'].disabled)
        self.assertTrue(controls['minimize_flexible_growth_ramp'].disabled)
        self.assertTrue(controls['simplex_max_cycles'].disabled)
        self.toggle(controls['minimize_ligand'],True)
        self.assertFalse(controls['simplex_max_cycles'].disabled)
        self.toggle(controls['score_molecules'],False)
        self.assertTrue(controls['grid_score_primary'].disabled)
        self.assertTrue(controls['grid_score_vdw_scale'].disabled)
        self.toggle(controls['score_molecules'],True)
        self.assertFalse(controls['grid_score_vdw_scale'].disabled)

    def test_redocking_preparation_and_cofactors_follow_switches(self):
        stage=new_stage('redocking');stage['parameters']['pdb_codes']=[['1ABC','LIG',1,'A']]
        form=self.form(stage=stage)
        _,_,(has_cofactors,cofactors),settings,_=form.redocking_pairs.rows[0]
        self.assertTrue(cofactors.disabled)
        self.toggle(has_cofactors,True);self.assertFalse(cofactors.disabled)
        cofactors.value='ATP'
        self.toggle(form.field_controls['prepare_complex'],False)
        self.assertTrue(cofactors.disabled);self.assertTrue(has_cofactors.disabled)
        self.assertTrue(settings['ligand']['charge_type'].disabled)
        self.toggle(form.field_controls['prepare_complex'],True)
        self.assertFalse(settings['ligand']['charge_type'].disabled)
        self.assertEqual(cofactors.value,'ATP');self.assertFalse(cofactors.disabled)
