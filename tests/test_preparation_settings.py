"""Preparation labels, independent molecule settings and scientific handoff."""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import flet as ft
from biomolexplorer.catalog import new_stage
from biomolexplorer.operations import validate_operation
from biomolexplorer.redocking_config import DEFAULTS, pair_key, configure_template
from biomolexplorer.templates import RESOURCE_ROOT
from biomolexplorer.ui.guided import GuidedForm
from biomolexplorer.ui.localization import Translator


class PreparationSettingsTests(unittest.TestCase):
    def form(self,stage=None,writable=True):
        stage=stage or new_stage('prepare_structures')
        ui=SimpleNamespace(current={'pipeline':[stage]},service=SimpleNamespace(dock6_path=None),page=SimpleNamespace(update=lambda:None))
        return GuidedForm(ui,stage,[],writable)

    def test_title_and_legacy_default_names_use_prepare_for_docking(self):
        stage=new_stage('prepare_structures')
        self.assertEqual(Translator('en')(stage['name']),'Prepare for docking')
        for old in ('Preparar meus complexos','Prepare my complexes'):
            stage['name']=old
            self.assertEqual(self.form(stage).read()['name'],'Preparar para docking')
            self.assertEqual(Translator('en')(old),'Prepare for docking')
        stage['name']='My custom study'
        self.assertEqual(self.form(stage).read()['name'],'My custom study')

    def test_only_two_scoped_preparation_groups_are_shown(self):
        form=self.form();panel=form.preparation_settings
        self.assertEqual([g.title.value for g in panel.control.controls],
            ['Preparação do receptor','Preparação e conformação do ligante'])
        self.assertFalse(any(n.startswith('chimera/') for n in form.template_readers))
        self.assertNotIn('charge_type',form.field_controls)
        self.assertIn('receptor e ligante',form.field_controls['pH'].label)
        receptor=panel.control.controls[0].controls
        self.assertIn(panel.has_cofactors,receptor);self.assertIn(panel.cofactors,receptor)
        for role in ('receptor','ligand'):
            self.assertEqual(set(panel.settings[role]),set(DEFAULTS))
            self.assertNotIn('Chimera',panel.settings[role]['charge_type'].label)

    def test_independent_settings_survive_save_and_reopen(self):
        form=self.form();panel=form.preparation_settings
        panel.settings['receptor']['add_hydrogens'].value=False
        panel.settings['ligand']['minimize'].value=False
        panel.settings['ligand']['charge_type'].value='am1'
        panel.has_cofactors.value=True;panel.cofactors.value='fad, mg'
        result=form.read();settings=result['parameters']['preparation_options']
        self.assertFalse(settings['receptor']['add_hydrogens'])
        self.assertTrue(settings['ligand']['add_hydrogens'])
        self.assertTrue(settings['receptor']['minimize'])
        self.assertFalse(settings['ligand']['minimize'])
        self.assertEqual(settings['cofactors'],['FAD','MG'])
        self.assertNotIn('charge_type',result['parameters'])
        reopened=self.form(result)
        self.assertEqual(reopened.read()['parameters']['preparation_options'],settings)
        self.assertTrue(all(c.disabled for o in self.form(result,False).preparation_settings.settings.values() for c in o.values()))

    def test_legacy_template_options_are_preserved_in_scoped_settings(self):
        stage=new_stage('prepare_structures');stage['parameters']['charge_type']='am1'
        stage['templates']={'chimera/prepare_receptor.template':
            (RESOURCE_ROOT/'chimera/prepare_receptor.template').read_text().replace('addh\n','')}
        result=self.form(stage).read()
        settings=result['parameters']['preparation_options']
        self.assertFalse(settings['receptor']['add_hydrogens'])
        self.assertTrue(settings['ligand']['add_hydrogens'])
        self.assertEqual(settings['ligand']['charge_type'],'am1')
        self.assertEqual(result['templates'],stage['templates'])

    def test_complex_separation_survives_disabled_preparation_options(self):
        settings={role:{**DEFAULTS,'remove_solvent':False,'remove_hydrogens':False,'add_hydrogens':False,'minimize':False}
                  for role in ('receptor','ligand')}
        source=(RESOURCE_ROOT/'chimera/prepare_complex.template').read_text()
        result=configure_template(source,'prepare_complex.template',settings)
        self.assertIn('select #0:.{chain}',result)
        self.assertIn('select invert\ndelete selected',result)
        self.assertNotIn('delete solvent',result)

    def test_backend_applies_options_to_metadata_selected_complexes(self):
        from wrappers.redocking import prepare_structures
        config=self.form().read()['parameters']['preparation_options']
        config['ligand']['charge_type']='am1';config['receptor']['minimize']=False
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);source=root/'input'/'Target';source.mkdir(parents=True)
            (source/'pdb_codes.csv').write_text('PDB_CODE,LIGAND,RESNUM,CHAIN\n1ABC,LIG,1,A\n2ABC,LIG,2,B\n')
            records=[['1ABC','LIG',1,'A'],['2ABC','LIG',2,'B']]
            with patch('wrappers.redocking.Docking') as engine:
                engine.return_value.prepare_for_docking.return_value=records
                prepare_structures(str(root/'input'),'Target',str(root/'output'),preparation_options=config)
                settings=engine.return_value.prepare_for_docking.call_args.kwargs['preparation_pairs']
                self.assertEqual(settings,{pair_key(r):config for r in records})
                self.assertTrue((root/'output/Target/pdb_codes.csv').exists())

    def test_invalid_scoped_options_are_rejected_before_execution(self):
        params={'base_input_path':'input','target':'Target','preparation_options':{'ligand':{'minimize':'false'}}}
        with self.assertRaisesRegex(ValueError,'booleana'):validate_operation('prepare_structures',params)



class CandidatePreparationTests(unittest.TestCase):
    form=PreparationSettingsTests.form

    def test_redocking_selector_hides_reference_ligands_and_disables_only_receptor(self):
        stage=new_stage('prepare_structures');origin=new_stage('redocking')
        stage['bindings']={'base_input_path':{'stage':origin['id'],'selector':'1ABC_A.dockprep.pdbqt'}}
        ui=SimpleNamespace(current={'pipeline':[origin,stage]},service=SimpleNamespace(dock6_path=None),
            page=SimpleNamespace(update=lambda:None),artifact_choices={origin['id']:{
                '1ABC_A.dockprep.pdbqt','1ABC_A.dockprep.mol2','1ABC_LIG_1A.pdbqt','1ABC.pdb','centers.csv'}})
        form=GuidedForm(ui,stage,[],True);panel=form.preparation_settings
        editor=form.input_editors['base_input_path']
        self.assertEqual({o.key for o in editor.rows[0][1].options},{'auto','1ABC_A.dockprep.pdbqt'})
        self.assertTrue(all(c.disabled for c in panel.settings['receptor'].values()))
        self.assertTrue(all(not c.disabled for c in panel.settings['ligand'].values()))
        self.assertTrue(form.read()['parameters']['receptor_prepared'])
        raw=new_stage('retrieve_structures');ui.current['pipeline'].append(raw)
        editor.rows[0][0].value='stage:'+raw['id'];editor.rows[0][1].value='auto';editor.changed()
        self.assertTrue(all(not c.disabled for c in panel.settings['receptor'].values()))
        self.assertFalse(form.read()['parameters']['receptor_prepared'])

    def test_multiple_compound_sources_and_engine_output_contract(self):
        from biomolexplorer.flow import connect,compatible
        from biomolexplorer.bindings import sources
        stage=new_stage('prepare_structures');origins=[new_stage(op) for op in ('retrieve_compounds','retrieve_pubchem','retrieve_zinc')]
        stages=origins+[stage]
        for origin in origins:connect(stages,origin['id'],stage['id'],'base_selected_mols')
        self.assertEqual(len(sources(stage['bindings']['base_selected_mols'])),3)
        stage['parameters']['docking_engines']='vina'
        self.assertTrue(compatible(stage,new_stage('docking_vina'),'base_selected_mols'))
        self.assertFalse(compatible(stage,new_stage('docking_dock6'),'base_input_path'))
        stage['parameters']['docking_engines']='both'
        for engine in ('vina','dock6'):
            for field in ('base_input_path','base_selected_mols'):
                self.assertTrue(compatible(stage,new_stage('docking_'+engine),field))

    def test_mixed_raw_and_prepared_receptors_cannot_be_saved(self):
        raw,prepared,stage=new_stage('retrieve_structures'),new_stage('redocking'),new_stage('prepare_structures')
        stage['bindings']={'base_input_path':{'sources':[{'stage':raw['id']},{'stage':prepared['id']}]}}
        ui=SimpleNamespace(current={'pipeline':[raw,prepared,stage]},service=SimpleNamespace(dock6_path=None),page=SimpleNamespace(update=lambda:None))
        form=GuidedForm(ui,stage,[],True)
        with self.assertRaisesRegex(ValueError,'brutos ou preparados'):
            form.read()

    def test_prepared_receptor_is_copied_without_running_preparation(self):
        from wrappers.redocking import prepare_structures
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);source=root/'input/Target';prepared=source/'Prepared';prepared.mkdir(parents=True)
            (source/'pdb_codes.csv').write_text('PDB_CODE,LIGAND,RESNUM,CHAIN\n1ABC,LIG,1,A\n')
            for name in ('1ABC_A.dockprep.pdbqt','1ABC_A.dockprep.mol2','1ABC_A.noH.pdb','1ABC_LIG_1A.pdbqt','centers.csv'):
                (prepared/name).write_text('unchanged '+name)
            with patch('wrappers.redocking.Docking') as engine:
                prepare_structures(str(root/'input'),'Target',str(root/'output'),receptor_prepared=True)
                engine.assert_not_called()
            files={p.name for p in (root/'output/Target/Prepared').iterdir()}
            self.assertNotIn('1ABC_LIG_1A.pdbqt',files)
            self.assertIn('centers.csv',files)
            self.assertEqual((root/'output/Target/Prepared/1ABC_A.dockprep.pdbqt').read_text(),'unchanged 1ABC_A.dockprep.pdbqt')

    def test_exported_candidates_are_reused_without_conversion(self):
        from biomolexplorer.docking_data import write_csv,prepare_ligands,copy_input_records,input_records
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);a=root/'A.pdbqt';b=root/'A.mol2';a.write_text('vina prepared');b.write_text('dock6 prepared')
            table=root/'compounds.csv'
            write_csv(table,[dict(molecule_chembl_id='A',canonical_smiles='CCO',conformer_file=str(b),
                prepared_pdbqt=str(a),prepared_mol2=str(b),docking_engines='both',prepared_origin='pose')])
            copied=copy_input_records(input_records([table],[table,a,b]),root/'materialized')
            with patch('biomolexplorer.docking_data.convert_structure') as convert:
                for format,expected in (('pdbqt','vina prepared'),('mol2','dock6 prepared')):
                    outputs=prepare_ligands(copied,root/format,format)
                    self.assertEqual(outputs[0].read_text(),expected)
                convert.assert_not_called()
            write_csv(table,[dict(molecule_chembl_id='A',canonical_smiles='CCO',prepared_pdbqt=str(a),docking_engines='vina')])
            with self.assertRaisesRegex(ValueError,'formato solicitado'):prepare_ligands(table,root/'invalid','mol2')

    def test_charge_fields_have_vertical_separation(self):
        panel=self.form().preparation_settings
        for group in panel.control.controls:
            charge=group.controls[-1]
            self.assertIsInstance(charge,ft.Container)
            self.assertGreaterEqual(charge.padding.top,24)

    def test_preparation_errors_keep_receptor_filenames_in_english(self):
        translate=Translator('en')
        self.assertEqual(translate('Arquivo do receptor necessário para both: 1ABC_A.noH.pdb'),
            'Receptor file required for both: 1ABC_A.noH.pdb')
        self.assertEqual(translate('Centro do sítio de docking ausente ou inválido para 1ABC_A.'),
            'The docking site center is missing or invalid for 1ABC_A.')

if __name__=='__main__':unittest.main()
