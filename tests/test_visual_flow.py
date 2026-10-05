import copy
import unittest
from biomolexplorer.catalog import new_stage
from biomolexplorer import flow

class VisualFlowTests(unittest.TestCase):
    def test_graph_block_only_has_similarity_port_and_mcs_configuration(self):
        stage=new_stage('graphs')
        self.assertNotIn('metric',stage['parameters'])
        self.assertNotIn('fingerprint',stage['parameters'])
        self.assertNotIn('threshold',stage['parameters'])
        self.assertEqual([p['field'] for p in flow.input_ports(stage)],['similarity_path'])
    def test_compatible_link_and_target(self):
        source=new_stage('retrieve_compounds'); target=new_stage('admet')
        flow.connect([source,target],source['id'],target['id'],'base_input_path')
        self.assertEqual(target['bindings']['base_input_path']['stage'],source['id'])
    def test_incompatible_link_is_transactional(self):
        stages=[new_stage('retrieve_structures'),new_stage('admet')]; before=copy.deepcopy(stages)
        with self.assertRaises(ValueError): flow.connect(stages,stages[0]['id'],stages[1]['id'],'base_input_path')
        self.assertEqual(stages,before)
    def test_cycle_is_transactional(self):
        a,b=new_stage('admet'),new_stage('admet'); stages=[a,b]
        flow.connect(stages,a['id'],b['id'],'base_input_path'); before=copy.deepcopy(stages)
        with self.assertRaises(ValueError): flow.connect(stages,b['id'],a['id'],'base_input_path')
        self.assertEqual(stages,before)
    def test_import_types_and_removal(self):
        source=new_stage('import_results'); target=new_stage('docking_vina')
        source['parameters']['kind']='prepared_structures'
        self.assertTrue(flow.compatible(source,target,'base_input_path'))
        self.assertFalse(flow.compatible(source,target,'base_selected_mols'))
        stages=[source,target]; flow.connect(stages,source['id'],target['id'],'base_input_path')
        flow.remove(stages,source['id']); self.assertEqual(target['bindings'],{})
    def test_large_layout_has_no_overlapping_positions(self):
        stages=[new_stage('admet') for _ in range(100)]
        for i in range(1,100): flow.connect(stages,stages[i-1]['id'],stages[i]['id'],'base_input_path')
        flow.arrange(stages)
        self.assertEqual(len({tuple(s['position'].values()) for s in stages}),100)
        parallel=[new_stage('retrieve_compounds') for _ in range(100)]; flow.arrange(parallel)
        self.assertEqual(len({tuple(s['position'].values()) for s in parallel}),100)
    def test_new_block_finds_free_position_after_deletion(self):
        remaining=new_stage('admet'); remaining['position']={'x':410,'y':64}
        candidate=flow.available_position([remaining],new_stage('import_results'))
        self.assertNotEqual(candidate,remaining['position'])
    def test_consensus_alias_removed_together(self):
        a,b=new_stage('docking_vina'),new_stage('consensus')
        flow.connect([a,b],a['id'],b['id'],'base_vina_path')
        self.assertEqual(b['bindings']['base_input_path'],b['bindings']['base_vina_path'])
        flow.disconnect(b,'base_vina_path'); self.assertEqual(b['bindings'],{})

if __name__=='__main__': unittest.main()
