"""Typed graph contracts, MD molecular content and deliberate receptor choices."""
import copy
import tempfile
from pathlib import Path
import unittest
from biomolexplorer import flow
from biomolexplorer.catalog import new_stage
from biomolexplorer.md_inputs import validate_structure,receptor_input_files
from biomolexplorer.pipeline import validate_pipeline
from test_md_pair_selection import AmberPairSelectionTests


class StrictPortTests(unittest.TestCase):
    def test_generic_outputs_are_rejected_on_every_input_port_transactionally(self):
        for operation in flow.INPUTS:
            stage=new_stage(operation)
            origin=new_stage('import_results');origin['parameters']['kind']='other'
            for port in flow.input_ports(stage):
                with self.subTest(operation=operation,field=port['field']):
                    stages=[origin,stage];before=copy.deepcopy(stages)
                    self.assertFalse(flow.compatible(origin,stage,port['field']))
                    with self.assertRaisesRegex(ValueError,'Dados incompatíveis'):
                        flow.connect(stages,origin['id'],stage['id'],port['field'])
                    self.assertEqual(stages,before)
                    forged=copy.deepcopy(stages)
                    forged[-1]['bindings'][port['field']]={'stage':origin['id']}
                    with self.assertRaisesRegex(ValueError,'Dados incompatíveis'):validate_pipeline(forged)

    def test_md_roles_are_distinct_but_native_multi_role_producers_remain_usable(self):
        target=new_stage('dynamics_dock6')
        for operation in ('retrieve_structures','redocking','prepare_structures','docking_vina','docking_dock6','dynamics_dock6'):
            source=new_stage(operation)
            self.assertTrue(flow.compatible(source,target,'base_input_path'),operation)
        for operation in ('retrieve_compounds','retrieve_pubchem','retrieve_zinc','admet','consensus'):
            source=new_stage(operation)
            self.assertFalse(flow.compatible(source,target,'base_input_path'),operation)
            self.assertTrue(flow.compatible(source,target,'base_selected_mols'),operation)
        source=new_stage('import_results')
        for kind,receivers,ligands in [('structures',True,False),('receptor_pdb',True,False),
                ('compounds',False,True),('molecular_poses',False,True),('vina',False,True),
                ('dock6',False,True),('scores',False,False),('dynamics',False,False),('other',False,False)]:
            source['parameters']['kind']=kind
            self.assertEqual(flow.compatible(source,target,'base_input_path'),receivers,kind)
            self.assertEqual(flow.compatible(source,target,'base_selected_mols'),ligands,kind)

    def test_hidden_ports_and_manual_reference_mode_cannot_receive_edges(self):
        source=new_stage('retrieve_compounds');target=new_stage('retrieve_pubchem')
        target['parameters']['reference_source']='manual'
        self.assertFalse(flow.compatible(source,target,'base_input_path'))
        target=new_stage('admet');target['provided_results']={'kind':'compounds','asset_ids':['a'*32]}
        self.assertFalse(flow.compatible(source,target,'base_input_path'))
        self.assertFalse(flow.compatible(source,new_stage('graphs'),'base_input_path'))

    def test_pdb_ligand_and_polymer_mol2_are_not_interchangeable(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);ligand=root/'ligand.pdb';receptor=root/'receptor.pdb'
            ligand.write_text('HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n')
            receptor.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n')
            validate_structure(receptor,'receptor')
            broken=root/'Prepared'/'receptor_A.noH.pdb';broken.parent.mkdir();broken.write_text(ligand.read_text())
            self.assertEqual(receptor_input_files([broken,receptor]),[receptor])
            with self.assertRaisesRegex(ValueError,'macromolécula'):validate_structure(ligand,'receptor')
            self.assertEqual(receptor_input_files([ligand,receptor]),[receptor])
            mol2=root/'protein.mol2';mol2.write_text('@<TRIPOS>MOLECULE\nProtein\n1 0 1 0 0\nPROTEIN\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n1 CA 0 0 0 C.3 1 ALA 0\n@<TRIPOS>BOND\n')
            with self.assertRaisesRegex(ValueError,'molécula pequena'):validate_structure(mol2,'ligand')
            # Cache is keyed by stat signature: changing content cannot retain
            # the earlier receptor role after a file replacement.
            receptor.write_text(ligand.read_text())
            self.assertEqual(receptor_input_files([receptor]),[])

    def test_invalid_canvas_drop_reports_feedback_without_adding_an_edge(self):
        import asyncio
        from types import SimpleNamespace
        from unittest.mock import Mock
        from biomolexplorer.ui.flow_canvas import FlowCanvas
        from test_frontend_navigation import Page
        source=new_stage('import_results');source['parameters']['kind']='other'
        target=new_stage('dynamics_dock6');stages=[source,target]
        ui=SimpleNamespace(current={'pipeline':stages},selected=None,page=Page(),notify=Mock(),service=SimpleNamespace(dock6_path=None))
        canvas=FlowCanvas(ui,True);before=copy.deepcopy(stages)
        asyncio.run(canvas.link(source['id'],target['id'],'base_input_path'))
        self.assertEqual(stages,before)
        self.assertIn('incompatíveis',canvas.status.value)
        ui.notify.assert_called_once()



class DeliberateReceptorSelectionTests(AmberPairSelectionTests):
    def test_auto_connection_starts_empty_and_adds_only_chosen_files(self):
        form=self.form(choose_receptors=False);editor=form.input_editors['base_input_path']
        picker=editor.selection.pickers['base_input_path']
        self.assertEqual(len(picker.available.options),34)
        self.assertFalse(any(flag.value for flag,_ in picker.rows))
        self.assertEqual(picker.pair_table.table.rows,[])
        self.assertEqual(form.md_settings.pairs,[])
        picker.set_selected([7],True)
        self.assertEqual(len(picker.pair_table.table.rows),1)
        self.assertEqual(len(form.md_settings.pairs),1)
        form.md_settings.all(None)
        saved=form.read();previous_pairs=saved['parameters']['selected_pairs']
        self.stage=saved
        restored=self.form(choose_receptors=False)
        picker=restored.input_editors['base_input_path'].selection.pickers['base_input_path']
        self.assertEqual(sum(bool(flag.value) for flag,_ in picker.rows),1)
        self.assertEqual(restored.read()['parameters']['selected_pairs'],previous_pairs)

    def test_old_auto_binding_restores_only_receptors_with_saved_pair_ids(self):
        form=self.form();settings=form.md_settings
        settings.search.value='A007';settings.filter(None);settings.all(None)
        saved=form.read();saved['bindings']['base_input_path']={'stage':self.producer['id'],'selector':'auto'}
        self.stage=saved
        restored=self.form(choose_receptors=False)
        picker=restored.input_editors['base_input_path'].selection.pickers['base_input_path']
        self.assertEqual(sum(bool(flag.value) for flag,_ in picker.rows),1)
        self.assertEqual(len(restored.md_settings.pairs),1)
        self.assertEqual(restored.read()['parameters']['selected_pairs'],saved['parameters']['selected_pairs'])

    def test_uploaded_receptors_keep_their_typed_role_without_offering_other_formats(self):
        form=self.form(choose_receptors=False)
        assets=[dict(id='r',name='receptor.pdb',kind='structures'),dict(id='l',name='ligand.mol2',kind='molecular_poses'),
            dict(id='c',name='molecules.csv',kind='compounds'),dict(id='x',name='picture.png',kind='other'),
            dict(id='legacy',name='legacy.mol2',kind='other')]
        for editor in form.input_editors.values():editor.assets[:]=assets
        receptor_options={o.key for o in form.input_editors['base_input_path'].source_options()}
        ligand_options={o.key for o in form.input_editors['base_selected_mols'].source_options()}
        self.assertIn('asset:r',receptor_options);self.assertNotIn('asset:l',receptor_options)
        self.assertIn('asset:l',ligand_options);self.assertIn('asset:c',ligand_options)
        self.assertNotIn('asset:r',ligand_options)
        self.assertNotIn('asset:x',receptor_options|ligand_options)
        self.assertNotIn('asset:legacy',receptor_options|ligand_options)
        editor=form.input_editors['base_selected_mols']
        editor.initial_references.append({'asset':'legacy'})
        self.assertIn('asset:legacy',{o.key for o in editor.source_options()})
