"""The curated file selection is the executable pair plan, without reselection."""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from biomolexplorer.catalog import new_stage
from biomolexplorer.redocking_config import validate_pairs,pair_key
from biomolexplorer.ui.file_selection import FileSelection
from biomolexplorer.ui.guided import GuidedForm

class CuratedHandoffTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        root=Path(self.temp.name)
        self.producer=new_stage('retrieve_structures');self.stage=new_stage('redocking')
        paths=[]
        for code in ('1ABC','2XYZ'):
            path=root/(code+'.pdb');path.write_text('END\n');paths.append(str(path))
        metadata=root/'pdb_codes.csv'
        metadata.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RESOLUTION\n1ABC,ATP,102,A,1.7\n2XYZ,LIG,23,B,2.0\n')
        self.item=dict(self.producer,configuration=self.producer,status='succeeded',artifacts=paths+[str(metadata)])
        self.stage['bindings']['base_input_path']={'stage':self.producer['id']}
        self.page=SimpleNamespace(update=lambda:None)
        self.form=FileSelection({'stages':[self.item]},{'configuration':self.stage,'name':self.stage['name']},page=self.page)
        self.paths=paths

    def guided(self,configuration):
        ui=SimpleNamespace(current={'pipeline':[self.producer,configuration]},service=SimpleNamespace(dock6_path=None),
            page=self.page,artifact_items={self.producer['id']:self.item},redocking_records={})
        return GuidedForm(ui,configuration,[],True)

    def test_reviewed_pairs_survive_stale_empty_run_snapshot_and_main_confirmation(self):
        self.form.pickers['base_input_path'].add_all_sources()
        self.form.apply_pdb_reviews({self.paths[0]:dict(reviewed=True,ligands=[dict(LIGAND='ATP',RESNUM='102',CHAIN='A',RESOLUTION='1.7')])})
        selected=self.form.read()
        self.assertEqual(selected['parameters']['pdb_codes'],[['1ABC','ATP',102,'A',1.7],['2XYZ','LIG',23,'B',2.]])
        main=self.guided(selected)
        confirmed=main.read()['parameters']
        self.assertEqual(confirmed['pdb_codes'],selected['parameters']['pdb_codes'])
        validate_pairs(confirmed['pdb_codes'],confirmed['preparation_pairs'])
        self.assertFalse(main.redocking_pairs.available.visible)
        self.assertFalse(main.redocking_pairs.select_all_button.visible)

    def test_general_receptor_and_ligand_options_apply_to_every_pair_and_reopen(self):
        self.form.pickers['base_input_path'].add_all_sources()
        main=self.guided(self.form.read());preparation=main.redocking_pairs.preparation
        preparation.settings['receptor']['remove_hydrogens'].value=False
        preparation.settings['ligand']['minimize'].value=False
        preparation.settings['ligand']['charge_type'].value='am1'
        saved=main.read()
        for record in saved['parameters']['pdb_codes']:
            config=saved['parameters']['preparation_pairs'][pair_key(record)]
            self.assertFalse(config['receptor']['remove_hydrogens'])
            self.assertFalse(config['ligand']['minimize'])
            self.assertEqual(config['ligand']['charge_type'],'am1')
        self.assertEqual(self.guided(saved).read(),saved)

    def test_empty_curated_structure_does_not_restore_old_metadata_or_pairs(self):
        self.form.pickers['base_input_path'].add_all_sources()
        self.form.apply_pdb_reviews({self.paths[0]:dict(reviewed=True,ligands=[])})
        self.assertEqual([r[0] for r in self.form.read()['parameters']['pdb_codes']],['2XYZ'])

    def test_unchecked_structure_prunes_pairs_and_preparation(self):
        self.form.pickers['base_input_path'].add_all_sources()
        configuration=self.guided(self.form.read()).read()
        second=FileSelection({'stages':[self.item]},{'configuration':configuration,'name':configuration['name']},page=self.page)
        second.rows['base_input_path'][0][0].value=False
        selected=second.read()['parameters']
        self.assertEqual([r[0] for r in selected['pdb_codes']],['2XYZ'])
        self.assertEqual(set(selected['preparation_pairs']),{'2XYZ|LIG|23|B'})

    def test_execution_prunes_only_orphan_settings_without_changing_saved_project(self):
        from biomolexplorer.pipeline import execution_stage
        stage=new_stage('redocking')
        stage['parameters']['pdb_codes']=[['1ABC','ATP',102,'A']]
        stage['parameters']['preparation_pairs']={
            '1ABC|ATP|102|A':{'ligand':{'minimize':False}},
            '2XYZ|LIG|23|B':{'ligand':{'minimize':True}}}
        effective=execution_stage(stage,'curated')
        self.assertEqual(set(effective['parameters']['preparation_pairs']),{'1ABC|ATP|102|A'})
        self.assertFalse(effective['parameters']['preparation_pairs']['1ABC|ATP|102|A']['ligand']['minimize'])
        self.assertEqual(len(stage['parameters']['preparation_pairs']),2)
