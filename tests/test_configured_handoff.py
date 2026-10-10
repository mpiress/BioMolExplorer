"""Saved input choices must not require repeated execution-time curation."""
import unittest
from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import requires_curation,configured_inputs

class ConfiguredHandoffTests(unittest.TestCase):
    def test_explicit_bindings_skip_review_for_each_supported_operation(self):
        for op in ('admet','fingerprints','similarity','graphs','retrieve_pubchem','retrieve_zinc','expand_similar_compounds','docking_vina'):
            with self.subTest(operation=op):
                stage=new_stage(op);stage['bindings']['base_input_path']={'stage':'a'*32,'selector':'chosen.csv'}
                self.assertFalse(requires_curation(stage))
                stage['parameters'].pop('filename',None)
                stage['parameters'].pop('input_file',None)
                stage['bindings']['base_input_path']['selector']='auto'
                self.assertTrue(requires_curation(stage))

    def test_saved_filenames_materialize_explicit_selectors_without_mutating_draft(self):
        for op,key,value in [('admet','input_file','chosen.csv'),('retrieve_pubchem','input_file','chosen.csv'),
                             ('fingerprints','files',['chosen.csv','second.csv']),('similarity','filename','chosen.csv')]:
            with self.subTest(operation=op):
                stage=new_stage(op);stage['parameters'][key]=value
                stage['bindings']['base_input_path']={'stage':'a'*32,'selector':'auto'}
                result=configured_inputs(stage)
                self.assertFalse(requires_curation(result))
                self.assertEqual(stage['bindings']['base_input_path']['selector'],'auto')

    def test_any_unconfigured_input_still_requires_selection(self):
        stage=new_stage('docking_vina');stage['bindings']={
            'base_input_path':{'stage':'a'*32,'selector':'receptor.pdb'},
            'base_selected_mols':{'stage':'b'*32,'selector':'auto'}}
        self.assertTrue(requires_curation(stage))
