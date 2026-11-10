import tempfile
import unittest
from pathlib import Path

import flet as ft

from biomolexplorer.artifact_naming import output_role
from biomolexplorer.artifacts import index_run
from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.ui.docking_files import DockingFilePicker
from biomolexplorer.ui.file_selection import FileSelection
from biomolexplorer.ui.input_editor import InputEditor
from biomolexplorer.catalog import operation_fields
from biomolexplorer.bindings import pack,sources
from types import SimpleNamespace
from biomolexplorer.workspace import WorkspaceStore


class DockingReferenceRegressionTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.store=WorkspaceStore(self.tmp.name)
        self.token=self.store.register('Owner','docking-reference@example.org','owner-password')
        self.project=self.store.create_project(self.token,'Docking references');self.pid=self.project['id']

    def tearDown(self):self.tmp.cleanup()

    def test_engine_receptors_publish_receptor_port_instead_of_pose_port(self):
        for operation,suffix in (('docking_vina','pdbqt'),('docking_dock6','mol2'),
                ('dynamics_dock6','pdbqt'),('redocking','pdbqt')):
            with self.subTest(operation=operation):
                stage=new_stage(operation);path=self.store.project_dir(self.pid)/stage['id']/f'4M0E_A.dockprep.{suffix}'
                path.parent.mkdir();path.write_text('ATOM\n' if suffix=='pdbqt' else '@<TRIPOS>MOLECULE\nREC\n')
                stage.update(status='succeeded',artifacts=[str(path)],configuration=dict(stage))
                index_run([stage],'a'*32)
                self.assertEqual(output_role(operation,path,stage['parameters']),'receptor')
                self.assertEqual(stage['artifact_records'][0]['output_port'],
                    'structures' if operation=='redocking' else 'receptors')

    def test_reference_metadata_does_not_create_unavailable_duplicate(self):
        producer=new_stage('prepare_structures');consumer=new_stage('docking_vina')
        root=self.store.project_dir(self.pid)/producer['id']/'Target'/'Prepared';root.mkdir(parents=True)
        receptor=root/'4M0E_A.dockprep.pdbqt';receptor.write_text('ATOM\n')
        item=dict(producer,status='succeeded',artifacts=[str(receptor)])
        index_run([item],'b'*32);record=item['artifact_records'][0]
        consumer['bindings']['base_input_path']={'stage':producer['id'],'selector':record['selector'],
            'artifact_id':record['artifact_id'],'output_port':'prepared_receptors','mode':'current'}
        form=FileSelection({'stages':[item]},{'configuration':consumer,'name':consumer['name']})
        rows=form.rows['base_input_path']
        self.assertEqual(len(rows),1);self.assertTrue(rows[0][0].available);self.assertTrue(rows[0][0].value)
        saved=form.read()['bindings']['base_input_path']
        self.assertEqual(saved['artifact_id'],record['artifact_id'])
        self.assertIsInstance(form.pickers['base_input_path'],DockingFilePicker)
        self.assertTrue(any(isinstance(control,ft.DataTable) for control in
            form.pickers['base_input_path'].table_host.content.controls))

    def test_semantic_receptor_keeps_companion_metadata_for_redocking(self):
        producer=new_stage('prepare_structures');stage=new_stage('redocking')
        stage['parameters'].update(target='Target',prepare_complex=False)
        root=self.store.project_dir(self.pid)/'source'/'Target';prepared=root/'Prepared';prepared.mkdir(parents=True)
        receptor=prepared/'4M0E_A.dockprep.pdbqt';receptor.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n')
        mol2=prepared/'4M0E_A.dockprep.mol2';mol2.write_text('@<TRIPOS>MOLECULE\nREC\n1 0 1 0 0\nPROTEIN\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n1 CA 0 0 0 C.3 1 ALA 0\n')
        ligand=prepared/'4M0E_LIG_1A.lig.pdbqt';ligand.write_text('HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00           C\n')
        metadata=root/'pdb_codes.csv';metadata.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN\n4M0E,LIG,1,A\n')
        ref={'stage':producer['id'],'selector':receptor.name,'_resolved_paths':[str(receptor)],
            '_bundle_paths':list(map(str,[receptor,mol2,ligand,metadata]))}
        service=PipelineService(self.store)
        try:
            destination,_,selected=service._materialize_inputs(self.pid,stage,'base_input_path',[ref],{})
        finally:service.close()
        self.assertTrue((destination/'Target'/'pdb_codes.csv').is_file())
        self.assertTrue((destination/'Target'/'Prepared'/ligand.name).is_file())
        self.assertIn(metadata,list(map(Path,selected)))

    def test_three_receptors_remain_three_after_editor_save_and_reopen(self):
        for operation in ('docking_vina','docking_dock6'):
            producer=new_stage('prepare_structures');stage=new_stage(operation)
            root=self.store.project_dir(self.pid)/producer['id']/'Prepared';root.mkdir(parents=True)
            paths=[]
            for code in ('4M0E','4EY7','1B41'):
                path=root/(code+'_A.dockprep.pdbqt');path.write_text('ATOM\n');paths.append(str(path))
            item=dict(producer,status='succeeded',artifacts=paths)
            index_run([item],'a'*32)
            from biomolexplorer.artifacts import reference
            refs=[dict(reference(record),mode='current',selection={'pdb_id':Path(record['path']).name[:4]})
                for record in item['artifact_records']]
            stage['bindings']['base_input_path']=pack(refs)
            ui=SimpleNamespace(current={'pipeline':[producer,stage]},artifact_items={producer['id']:item},
                artifact_choices={},page=SimpleNamespace(update=lambda:None))
            field=next(field for field in operation_fields(operation) if field['name']=='base_input_path')
            for iteration in range(3):
                editor=InputEditor(ui,stage,field,[],True)
                picker=editor.selection.pickers['base_input_path']
                self.assertEqual(len(picker.rows),3)
                self.assertEqual(len(picker.table.rows),3)
                self.assertFalse(editor.box.visible)
                selected=sources(editor.read())
                self.assertEqual(len(selected),3)
                self.assertTrue(all(flag.available for flag,_ in picker.rows))
                stage['bindings']['base_input_path']=pack(selected)

    def test_old_manifest_receptor_role_is_repaired_before_selection(self):
        producer=new_stage('redocking');stage=new_stage('docking_vina')
        root=self.store.project_dir(self.pid)/producer['id']/'Prepared';root.mkdir(parents=True)
        path=root/'4M0E_A.dockprep.pdbqt';path.write_text('ATOM\n')
        item=dict(producer,status='succeeded',artifacts=[str(path)],
            artifact_records=[{'path':str(path),'output':'4m0e-a-dockprep','artifact_id':'b'*64}])
        stage['bindings']['base_input_path']={'stage':producer['id'],'selector':path.name,
            'artifact_id':'b'*64,'output_port':'redocking_results','mode':'current'}
        form=FileSelection({'stages':[item]},{'configuration':stage,'name':stage['name']})
        self.assertEqual(len(form.rows['base_input_path']),1)
        self.assertEqual(form.read()['bindings']['base_input_path']['output_port'],'structures')

    def test_ligand_actions_open_selection_and_remove_only_the_chosen_file(self):
        producer=new_stage('retrieve_compounds');stage=new_stage('docking_vina')
        root=self.store.project_dir(self.pid)/'ligands';root.mkdir()
        path=root/'compounds.csv';path.write_text('molecule_chembl_id,canonical_smiles\nM1,CCO\nM2,CCN\n')
        stage['bindings']['base_selected_mols']={'stage':producer['id'],'selector':path.name}
        dialogs=[]
        page=SimpleNamespace(update=lambda:None,show_dialog=dialogs.append,pop_dialog=lambda:dialogs.pop())
        item=dict(producer,status='succeeded',artifacts=[str(path)])
        form=FileSelection({'stages':[item]},{'configuration':stage,'name':stage['name']},page=page)
        picker=form.pickers['base_selected_mols']
        actions=picker.table.rows[0].cells[-1].content.controls
        select=next(button for button in actions if button.tooltip=='Selecionar compostos')
        select.on_click(None)
        self.assertEqual(len(dialogs),1)
        choice=form.compound_fields[id(form.rows['base_selected_mols'][0][0])]
        choice.value='M2'
        dialogs[0].actions[0].on_click(None)
        self.assertEqual(form.read()['bindings']['base_selected_mols']['compound_id'],'M2')
        remove=next(button for button in picker.table.rows[0].cells[-1].content.controls if button.tooltip=='Remover arquivo')
        remove.on_click(None)
        self.assertEqual(picker.table.rows,[])
        self.assertEqual(form.selected_references(form.rows['base_selected_mols']),[])


if __name__=='__main__':unittest.main()
