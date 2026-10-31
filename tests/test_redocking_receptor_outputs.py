"""Redocking reference ligands remain internal and never feed downstream inputs."""
from pathlib import Path
import tempfile
import unittest
from biomolexplorer.catalog import new_stage
from biomolexplorer import flow
from biomolexplorer.pipeline import validate_pipeline
from biomolexplorer.artifact_choices import choices
from biomolexplorer.ui.file_selection import FileSelection


class RedockingReceptorOutputsTests(unittest.TestCase):
    def test_redocking_cannot_connect_to_any_compound_or_score_input(self):
        source=new_stage('redocking')
        for operation in ('dynamics_dock6','docking_vina','docking_dock6','prepare_structures','admet','fingerprints','retrieve_pubchem','consensus'):
            target=new_stage(operation)
            for port in flow.input_ports(target):
                if port['types'] & {'structures','prepared_structures','receptor_pdb'}:continue
                with self.subTest(operation=operation,field=port['field']):
                    self.assertFalse(flow.compatible(source,target,port['field']))
                    with self.assertRaises(ValueError):flow.connect([source,target],source['id'],target['id'],port['field'])
                    target['bindings'][port['field']]={'stage':source['id']}
                    with self.assertRaises(ValueError):validate_pipeline([source,target])
                    target['bindings'].clear()
        for operation in ('dynamics_dock6','docking_vina','docking_dock6','prepare_structures'):
            self.assertTrue(flow.compatible(source,new_stage(operation),'base_input_path'))

    def test_public_choices_contain_only_receptors_and_hide_crystal_ligands(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);prepared=root/'Prepared';prepared.mkdir()
            receptor=prepared/'4M0E_A.noH.pdb'
            receptor.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n')
            pdbqt=prepared/'4M0E_A.dockprep.pdbqt';pdbqt.write_text('ATOM\n')
            ligand=prepared/'4M0E_LIG_1A.lig.mol2';ligand.write_text('@<TRIPOS>MOLECULE\nLIG\n1 0 1 0 0\nSMALL\nUSER_CHARGES\n')
            metadata=prepared/'pdb_codes.csv';metadata.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN\n4M0E,LIG,1,A\n')
            source=new_stage('redocking');item=dict(source,status='succeeded',artifacts=list(map(str,[receptor,pdbqt,ligand,metadata])))
            self.assertEqual(choices(item),{receptor.name,pdbqt.name})
            target=new_stage('dynamics_dock6');target['bindings']={field:{'stage':source['id']} for field in ('base_input_path','base_selected_mols')}
            form=FileSelection({'stages':[item]},{'configuration':target,'name':target['name']})
            self.assertEqual(len(form.rows['base_input_path']),1)
            self.assertEqual(form.rows['base_selected_mols'],[])
            self.assertTrue(ligand.is_file());self.assertTrue(metadata.is_file())
