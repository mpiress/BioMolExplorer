"""Exercise the actual Amber Score pipeline boundary with mixed legacy ports."""
import copy
from pathlib import Path
import tempfile
import unittest

from biomolexplorer.workspace import WorkspaceStore
from biomolexplorer.catalog import new_stage
from biomolexplorer.artifacts import records,reference
from biomolexplorer.artifact_registry import SemanticArtifactResolver
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.md_inputs import MDInputs,pair_id


class MDSemanticInputsTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.store=WorkspaceStore(self.tmp.name)
        self.token=self.store.register('Owner','md-port@example.org','md-port-password')
        self.pid=self.store.create_project(self.token,'Mixed target ports')['id']
        root=self.store.project_dir(self.pid)
        self.target=root/'4M0E_A.noH.pdb'
        self.target.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n')
        pose=root/'pose.pdbqt';pose.write_text('REMARK VINA RESULT: -5 0 0\n')
        self.producer=new_stage('redocking')
        item=dict(self.producer,status='succeeded',artifacts=[str(self.target),str(pose)])
        descriptors=records(item)
        target=next(r for r in descriptors if r['path']==str(self.target))
        self.ref=dict(reference(target),output_port='structures')
        target.pop('output_port',None);target['output']='4m0e-a-noh'
        other=next(r for r in descriptors if r['path']==str(pose));other['output_port']='redocking_results'
        self.context=[dict(self.producer,_execution_artifact_records=descriptors)]
        self.results={self.producer['id']:item['artifacts']}
        ligand=root/'CHEMBL92663.mol2'
        ligand.write_text('@<TRIPOS>MOLECULE\nLIG\n1 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 LIG 0\n@<TRIPOS>BOND\n')
        self.md=new_stage('dynamics_dock6')
        self.md['bindings']={'base_input_path':copy.deepcopy(self.ref)}
        self.md['parameters']['base_selected_mols']=str(ligand)
        inputs=MDInputs(self.store,self.pid,self.results)
        r=inputs.files(self.md,'base_input_path')[0];l=inputs.files(self.md,'base_selected_mols')[0]
        self.md['parameters']['selected_pairs']=[pair_id(r,l)]
        self.service=object.__new__(PipelineService);self.service.store=self.store;self.service.dock6_path=root

    def resolve(self):
        return self.service._resolve(self.pid,self.store.user(self.token)['id'],self.md,self.results,pipeline=self.context)

    def test_mixed_legacy_target_ports_reach_actual_md_input_resolution(self):
        before=copy.deepcopy(self.context)
        with self.assertRaisesRegex(ValueError,'porta semântica'):
            SemanticArtifactResolver(self.store,self.pid,self.context,self.results).resolve(self.ref)
        result=self.resolve()
        self.assertEqual(result['md_inputs'][0]['receptor']['path'],str(self.target))
        self.assertEqual(len(result['md_inputs']),1)
        self.assertEqual(self.context,before)
        self.assertEqual(self.md['bindings']['base_input_path'],self.ref)

    def test_wrong_identity_and_wrong_port_are_not_accepted(self):
        for changes in ({'artifact_id':'0'*64},{'output_port':'redocking_results'}):
            with self.subTest(changes=changes):
                self.md['bindings']['base_input_path']=dict(self.ref,**changes)
                with self.assertRaisesRegex(ValueError,'porta semântica'):self.resolve()

    def test_unavailable_target_and_declared_conflicting_port_are_not_repaired(self):
        original=copy.deepcopy(self.results)
        self.results[self.producer['id']]=[p for p in original[self.producer['id']] if p!=str(self.target)]
        with self.assertRaisesRegex(ValueError,'porta semântica'):self.resolve()
        self.results=original
        target=next(r for r in self.context[0]['_execution_artifact_records'] if r['path']==str(self.target))
        target['output_port']='redocking_results'
        with self.assertRaisesRegex(ValueError,'porta semântica'):self.resolve()

    def test_pinned_reference_never_uses_current_run_compatibility(self):
        self.md['bindings']['base_input_path']=dict(self.ref,mode='pinned',run_id='a'*32)
        with self.assertRaisesRegex(ValueError,'histórica'):self.resolve()

    def test_each_declared_target_producer_uses_its_own_typed_port(self):
        from biomolexplorer.artifact_contracts import infer_output_port
        for operation in ('retrieve_structures','prepare_structures','docking_vina','docking_dock6','dynamics_dock6'):
            with self.subTest(operation=operation):
                self.context[0]['operation']=operation
                self.md['bindings']['base_input_path']=dict(self.ref,output_port=infer_output_port(operation,'receptor'))
                self.assertEqual(self.resolve()['md_inputs'][0]['receptor']['path'],str(self.target))


    def test_ligand_only_producers_do_not_trigger_target_file_inspection(self):
        from unittest.mock import patch
        from biomolexplorer.md_semantic import semantic_context
        consumer=dict(self.md,bindings={'base_selected_mols':dict(self.ref)})
        with patch('biomolexplorer.md_inputs.receptor_input_files',side_effect=AssertionError('Unrelated target inspected')):
            context=semantic_context(self.context,self.results,consumer)
        self.assertIs(context[0],self.context[0])



class MDTargetLabelsTests(unittest.TestCase):
    import test_md_pair_selection as fixtures
    setUp=fixtures.AmberPairSelectionTests.setUp
    form=fixtures.AmberPairSelectionTests.form

    def test_canvas_form_pair_table_and_tooltips_use_targets(self):
        from biomolexplorer.flow import input_ports
        form=self.form('en');self.ui.page.add(form.layout());settings=form.md_settings
        self.assertEqual(input_ports(form.stage)[0]['label'],'Targets')
        self.assertEqual(form.input_editors['base_input_path'].control.controls[0].value,'Target PDB files for Amber Score')
        self.assertEqual(settings.table.columns[1].label.value,'PDB / target')
        self.assertEqual(settings.search.label,'Search targets or ligands')
        self.assertEqual(settings.controls['amber_score_receptor_movable_atom_expr'].label,'Target NAB expression')
        self.assertTrue(settings.controls['amber_score_receptor_movable_atom_expr'].tooltip.startswith('Target NAB expression'))
        self.ui.page.set_language('pt');self.ui.page.update()
        self.assertEqual(form.input_editors['base_input_path'].control.controls[0].value,'Targets PDB para Amber Score')
        self.assertEqual(settings.search.label,'Buscar target ou ligante')
