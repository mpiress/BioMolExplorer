"""Exhaustive typed connections, independent of canvas contract implementation."""
import copy
import unittest
from biomolexplorer.catalog import TITLES,new_stage
from biomolexplorer.flow import compatible,connect,disconnect,input_ports
from biomolexplorer.input_validation import output_kind
from biomolexplorer.pipeline import validate_pipeline

# Independent reference for the public stage contract (not imported from flow).
OUTPUTS={'retrieve_compounds':{'compounds','chembl'},'expand_similar_compounds':{'compounds'},
    'retrieve_pubchem':{'compounds'},
    'retrieve_structures':{'structures'},'retrieve_zinc':{'compounds'},'prepare_structures':{'prepared_structures','compounds'},
    'admet':{'compounds'},'fingerprints':{'fingerprints','compounds'},'similarity':{'similarity'},'graphs':{'compounds'},
    'redocking':{'structures','prepared_structures'},'docking_vina':{'vina','compounds'},'docking_dock6':{'dock6','compounds'},'consensus':{'scores','compounds'}}
INPUTS={'expand_similar_compounds':{'base_input_path':{'chembl'}},'retrieve_zinc':{'base_input_path':{'other','zinc_urls'}},
    'retrieve_pubchem':{'base_input_path':{'compounds'}},
    'prepare_structures':{'base_input_path':{'structures','prepared_structures'},'base_selected_mols':{'compounds'}},'admet':{'base_input_path':{'compounds'}},
    'fingerprints':{'base_input_path':{'compounds'}},'similarity':{'base_input_path':{'fingerprints'}},
    'graphs':{'similarity_path':{'similarity'}},'redocking':{'base_input_path':{'structures'}},
    'docking_vina':{'base_input_path':{'structures','prepared_structures'},'base_selected_mols':{'compounds','vina','dock6'}},
    'docking_dock6':{'base_input_path':{'structures','prepared_structures'},'base_selected_mols':{'compounds','vina','dock6'}},
    'consensus':{'base_vina_path':{'vina'},'base_dock6_path':{'dock6'}}}
KINDS=('compounds','structures','prepared_structures','fingerprints','similarity','vina','dock6','scores','other','zinc_urls')


class ConnectionMatrixTests(unittest.TestCase):
    def test_every_native_imported_and_supplied_result_against_every_port(self):
        producers=[]
        for operation,kinds in OUTPUTS.items():
            producers.append((new_stage(operation),kinds))
            provided=new_stage(operation)
            provided['provided_results']={'kind':output_kind(operation),'asset_ids':['f'*32]}
            producers.append((provided,{output_kind(operation)}))
        for kind in KINDS:
            producer=new_stage('import_results');producer['parameters']['kind']=kind
            producers.append((producer,{kind}))
        accepted=rejected=0
        for producer,kinds in producers:
            for operation,ports in INPUTS.items():
                variants=(True,False) if operation=='redocking' else (None,)
                for prepare in variants:
                    for field,expected in ports.items():
                        target=new_stage(operation)
                        if prepare is not None:
                            target['parameters']['prepare_complex']=prepare
                            expected={'structures'} if prepare else {'prepared_structures'}
                        allowed=bool(kinds&expected) and (operation!='graphs' or producer['operation']=='similarity')
                        with self.subTest(source=producer['operation'],kinds=kinds,target=operation,field=field,prepare=prepare):
                            self.assertEqual(compatible(producer,target,field),allowed)
                            stages=copy.deepcopy([producer,target]);before=copy.deepcopy(stages)
                            if allowed:
                                connect(stages,producer['id'],target['id'],field)
                                self.assertEqual(validate_pipeline(stages),[producer['id'],target['id']])
                                accepted+=1
                            else:
                                with self.assertRaises(ValueError):connect(stages,producer['id'],target['id'],field)
                                self.assertEqual(stages,before);rejected+=1
        self.assertGreater(accepted,70);self.assertGreater(rejected,400)

    def test_consensus_multiple_sources_preserve_one_alias_and_disconnect(self):
        a,b,target=new_stage('docking_vina'),new_stage('docking_vina'),new_stage('consensus')
        stages=[a,b,target]
        for source in (a,b):connect(stages,source['id'],target['id'],'base_vina_path')
        self.assertEqual(target['bindings']['base_input_path'],target['bindings']['base_vina_path'])
        self.assertIsNot(target['bindings']['base_input_path'],target['bindings']['base_vina_path'])
        disconnect(target,'base_vina_path');self.assertEqual(target['bindings'],{})

    def test_prepared_redocking_rejects_raw_inputs_before_execution(self):
        raw,prepared,target=new_stage('retrieve_structures'),new_stage('prepare_structures'),new_stage('redocking')
        self.assertEqual(input_ports(target)[0]['types'],{'structures'})
        self.assertTrue(compatible(raw,target,'base_input_path'));self.assertFalse(compatible(prepared,target,'base_input_path'))
        target['parameters']['prepare_complex']=False
        self.assertEqual(input_ports(target)[0]['types'],{'prepared_structures'})
        self.assertFalse(compatible(raw,target,'base_input_path'));self.assertTrue(compatible(prepared,target,'base_input_path'))

    def test_self_links_and_disabled_sources_never_mutate_pipeline(self):
        for operation in TITLES:
            target=new_stage(operation)
            for port in input_ports(target):self.assertFalse(compatible(target,target,port['field']))
        source,target=new_stage('admet'),new_stage('fingerprints');source['enabled']=False
        stages=[source,target];before=copy.deepcopy(stages)
        with self.assertRaises(ValueError):connect(stages,source['id'],target['id'],'base_input_path')
        self.assertEqual(stages,before)

    def test_rejected_cycle_preserves_receptor_preparation_mode(self):
        source,target=new_stage('redocking'),new_stage('prepare_structures')
        source['depends_on']=[target['id']]
        stages=[source,target];before=copy.deepcopy(stages)
        with self.assertRaises(ValueError):connect(stages,source['id'],target['id'],'base_input_path')
        self.assertEqual(stages,before)


if __name__=='__main__':unittest.main()
