"""Redocking contracts: shared options, table curation and engine-specific forms."""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch
import flet as ft
from biomolexplorer.catalog import new_stage
from biomolexplorer.redocking_config import effective_preparation, pair_key, shared_preparation, validate_redocking_parameters
from biomolexplorer.ui.guided import GuidedForm
from biomolexplorer.ui.localization import LocalizedPage
import test_localization as localization


class RedockingMethodsTests(unittest.TestCase):
    def form(self,language='en',records=None):
        stage=new_stage('redocking');stage['parameters']['pdb_codes']=records or [['1ABC','LIG',12,'A',1.7],['2ABC','ATP',13,'B',2.0]]
        page=SimpleNamespace(update=Mock(),show_dialog=Mock(),width=900)
        ui=SimpleNamespace(language=language,page=page,current={'pipeline':[stage]},redocking_records={},service=SimpleNamespace(dock6_path=None))
        return GuidedForm(ui,stage,[],True)

    def test_common_preparation_is_authoritative_and_copied_per_pair(self):
        records=[['1ABC','LIG',12,'A'],['1ABC','ATP',13,'A']]
        common={'ligand':{'minimize':False},'receptor':{'charge_type':'am1'},'cofactors':['FAD']}
        params={'pdb_codes':records,'preparation_options':common,'preparation_pairs':{pair_key(records[0]):{'cofactors':['MG']}}}
        result=effective_preparation(params)
        self.assertEqual(result[pair_key(records[1])],common)
        result[pair_key(records[0])]['cofactors'].clear()
        self.assertEqual(common['cofactors'],['FAD'])
        self.assertEqual(result[pair_key(records[1])]['cofactors'],['FAD'])
        _,conflict=shared_preparation({'preparation_pairs':{'a':{},'b':{'ligand':{'minimize':False}}}})
        self.assertTrue(conflict)

    def test_table_has_read_only_identity_and_only_edit_remove_actions(self):
        form=self.form();pairs=form.redocking_pairs
        self.assertEqual(len(pairs.table.columns),5)
        for row in pairs.table.rows:
            self.assertTrue(all(isinstance(c.content,ft.Text) for c in row.cells[:4]))
            self.assertEqual(len(row.cells[4].content.controls),2)
        pairs.resize(SimpleNamespace(width=1120));self.assertEqual(pairs.table.width,1120)
        pairs.resize(SimpleNamespace(width=360));self.assertEqual(pairs.table.width,680)
        pairs.table.rows[1].cells[4].content.controls[1].on_click(None)
        self.assertEqual(len(pairs.read()[0]),1)
        self.assertEqual(pairs.summary.value,'Pares selecionados: 1')

    def editor(self,form,index=0):
        pairs=form.redocking_pairs;pairs.open_editor(pairs.rows[index])
        dialog=form.ui.page.show_dialog.call_args.args[0]
        fields=[cell.content for cell in dialog.content.content.controls[0].controls]
        return dialog,fields

    def test_edit_is_transactional_keeps_resolution_and_rejects_duplicates(self):
        form=self.form();pairs=form.redocking_pairs
        dialog,fields=self.editor(form);fields[0].value='3ABC'
        dialog.actions[0].on_click(None)
        self.assertEqual(pairs.read()[0][0],['1ABC','LIG',12,'A',1.7])
        dialog,fields=self.editor(form);fields[0].value='3ABC';dialog.actions[1].on_click(None)
        self.assertEqual(pairs.read()[0][0],['3ABC','LIG',12,'A'])
        dialog,fields=self.editor(form)
        for field,value in zip(fields,['2ABC','ATP','13','B']):field.value=value
        dialog.actions[1].on_click(None)
        self.assertTrue(dialog.content.content.controls[1].visible)
        self.assertEqual(pairs.read()[0][0][0],'3ABC')

    def test_engine_switch_ignores_inactive_values_and_roundtrips_shared_options(self):
        form=self.form();layout=form.layout()
        form.field_controls['density'].value='not a number'
        result=form.read();self.assertEqual(result['parameters']['redocking_engine'],'vina')
        form.field_controls['density'].value='.75'
        form.field_controls['redocking_engine'].value='dock6';form.redocking_settings.changed()
        form.field_controls['exhaustiveness'].value='not a number'
        form.preparation_settings.settings['ligand']['minimize'].value=False
        result=form.read()
        self.assertEqual(result['parameters']['density'],.75)
        self.assertEqual(result['parameters']['max_orientations'],1000)
        self.assertFalse(form.redocking_settings.panels['vina'].visible)
        self.assertTrue(form.redocking_settings.panels['dock6'].visible)
        for options in result['parameters']['preparation_pairs'].values():self.assertFalse(options['ligand']['minimize'])
        reopened=GuidedForm(form.ui,result,[],True)
        self.assertEqual(reopened.read(),result)
        self.assertIn('dock6/redocking.template',result['templates'])

    def test_cofactors_are_optional_and_role_panels_have_matching_geometry(self):
        form=self.form();layout=form.layout();preparation=form.preparation_settings
        receptor,ligand=(preparation.groups[k] for k in ('receptor','ligand'))
        self.assertEqual([c.height for c in receptor.controls],[c.height for c in ligand.controls])
        self.assertFalse(preparation.cofactors.visible)
        preparation.has_cofactors.value=True;preparation.has_cofactors.on_change(None)
        self.assertTrue(preparation.cofactors.visible)
        preparation.cofactors.value='FAD';self.assertEqual(form.read()['parameters']['preparation_options']['cofactors'],['FAD'])
        preparation.has_cofactors.value=False;preparation.has_cofactors.on_change(None)
        self.assertFalse(preparation.cofactors.visible)
        self.assertEqual(form.read()['parameters']['preparation_options']['cofactors'],[])

    def test_entire_shared_section_uses_session_language(self):
        for language in ('en','pt'):
            form=self.form(language);raw=localization.Page();page=LocalizedPage(raw,language)
            form.ui.page=page;layout=form.layout();page.add(layout)
            shared=layout.controls[2]
            self.assertEqual(shared.content.controls[0].value,'Shared preparation for the collection' if language=='en' else 'Preparação única do conjunto')
            role=form.preparation_settings.settings['receptor']
            self.assertEqual(role['remove_solvent'].label,'Remove solvent' if language=='en' else 'Remover solvente')

    def test_engine_parameters_validate_only_current_engine(self):
        validate_redocking_parameters({'redocking_engine':'dock6','sizeof_box':[],'exhaustiveness':0})
        validate_redocking_parameters({'redocking_engine':'vina','density':-1})
        for params in ({'redocking_engine':'bad'},{'redocking_engine':'dock6','dock6_seed':0},{'sizeof_box':[0,24,24]}):
            with self.assertRaises(ValueError):validate_redocking_parameters(params)


class ReceptorFormatTests(unittest.TestCase):
    def test_selectors_accept_dock6_only_and_prefer_existing_pdbqt(self):
        from biomolexplorer.docking_inputs import preferred_receptors,prepared_receptor
        files=['A.dockprep.mol2','A.dockprep.pdbqt','B.dockprep.mol2','B.noH.pdb']
        self.assertTrue(prepared_receptor(files[2]))
        self.assertEqual(preferred_receptors(files),files[1:])

    def test_vina_format_export_retains_geometry_and_charges_without_chimera(self):
        from biomolexplorer.docking_preparation import export_receptor_pdbqt
        with tempfile.TemporaryDirectory() as folder:
            prepared=Path(folder);(prepared/'1ABC_A.dockprep.mol2').write_text('prepared')
            with patch('biomolexplorer.processes.run_command') as command:
                export_receptor_pdbqt(prepared,[['1ABC','LIG',1,'A'],['1ABC','ATP',2,'A']])
            args=command.call_args.args[0]
            self.assertEqual(command.call_count,1)
            self.assertIn('-xr',args);self.assertNotIn('-p',args);self.assertNotIn('--partialcharge',args)


class RedockingScientificContractTests(unittest.TestCase):
    MOL2='@<TRIPOS>MOLECULE\nligand\n1 0 1 0 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 C 0 0 0 C.3 1 LIG 0.0\n@<TRIPOS>BOND\n@<TRIPOS>SUBSTRUCTURE\n1 LIG 1\n'

    def test_rmsd_retains_translation_and_checks_pose_identity(self):
        from biomolexplorer.redocking_dock6 import heavy_atom_rmsd
        with tempfile.TemporaryDirectory() as folder:
            reference=Path(folder)/'reference.mol2';pose=Path(folder)/'pose.mol2'
            reference.write_text(self.MOL2)
            pose.write_text('########## Grid_Score: -5.0\n'+self.MOL2.replace('C 0 0 0','C 3 0 0'))
            self.assertAlmostEqual(heavy_atom_rmsd(reference,pose),3.)
            pose.write_text('########## Grid_Score: -5.0\n'+self.MOL2.replace('C.3','N.3'))
            with self.assertRaisesRegex(ValueError,'identidade química'):heavy_atom_rmsd(reference,pose)

    def test_crystal_reference_overrides_preparation_coordinates(self):
        from biomolexplorer.redocking_dock6 import crystal_reference
        with tempfile.TemporaryDirectory() as folder:
            reference=Path(folder)/'reference.mol2';prepared=Path(folder)/'prepared.mol2';raw=Path(folder)/'1ABC.pdb'
            prepared.write_text(self.MOL2)
            raw.write_text('HETATM    1  C   LIG A   1       5.000   0.000   0.000  1.00 20.00           C\n')
            self.assertEqual(crystal_reference(prepared,raw,['1ABC','LIG',1,'A'],reference),'crystal')
            self.assertIn('5.000000',reference.read_text())
            self.assertEqual(crystal_reference(prepared,Path(folder)/'absent.pdb',['1ABC','LIG',1,'A'],reference),'prepared')

    def test_dock6_keys_partition_engines_and_ignore_vina_only_knobs(self):
        from biomolexplorer.redocking_cache import PairCache
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);(root/'1ABC.pdb').write_text('raw')
            tools=root/'tools'
            for name in ('bin/dock6','bin/sphgen','bin/sphere_selector','bin/showbox','bin/grid','parameters/vdw_AMBER_parm99.defn','parameters/flex.defn','parameters/flex_drive.tbl'):
                p=tools/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(name)
            record=['1ABC','LIG',1,'A'];params={'redocking_engine':'dock6','dock6_app_path':str(tools)}
            a=PairCache(root,root/'output',[record],{},params)
            b=PairCache(root,root/'output',[record],{},dict(params,exhaustiveness=999,sizeof_box=[50]*3))
            c=PairCache(root,root/'output',[record],{},dict(params,dock6_seed=74))
            self.assertEqual(a.keys(record),b.keys(record))
            self.assertNotEqual(a.keys(record)[1],c.keys(record)[1])
            vina=PairCache(root,root/'output',[record],{},{})
            self.assertNotEqual(a.keys(record),vina.keys(record))
            (tools/'bin/dock6').write_text('changed binary')
            updated=PairCache(root,root/'output',[record],{},params)
            self.assertNotEqual(a.keys(record),updated.keys(record))
