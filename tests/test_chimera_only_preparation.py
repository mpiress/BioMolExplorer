"""Chimera is the only preparation route, including historical native settings."""
import asyncio
import copy
import unittest
from types import SimpleNamespace
from unittest.mock import Mock,AsyncMock,patch
from biomolexplorer.catalog import new_stage
from biomolexplorer.redocking_config import chimera_preparation,chimera_pairs
from biomolexplorer.pipeline import execution_stage
from biomolexplorer.ui.guided import GuidedForm
from biomolexplorer.ui.app import WorkspaceUI
from biomolexplorer.ui.updates import UpdateNotice
from test_frontend_navigation import Page,walk
import flet as ft

class ChimeraOnlyConfigurationTests(unittest.TestCase):
    def test_saved_native_configuration_migrates_without_losing_chemistry_options(self):
        config={'preparation_backend':'python','ambertools_env':'old','ambertools_path':'/old',
                'vina_charge_policy':'fast','receptor':{'minimize':False,'charge_type':'gas'},
                'ligand':{'charge_type':'am1'},'cofactors':['FAD'],'ligand_chain':'B'}
        before=copy.deepcopy(config);actual=chimera_preparation(config)
        self.assertEqual(config,before)
        self.assertEqual(actual['preparation_backend'],'legacy_chimera')
        self.assertEqual(actual['ligand']['charge_type'],'am1')
        self.assertEqual(actual['cofactors'],['FAD']);self.assertEqual(actual['ligand_chain'],'B')
        self.assertFalse({'ambertools_env','ambertools_path','vina_charge_policy'} & actual.keys())
        for operation in ('redocking','prepare_structures','docking_vina','docking_dock6'):
            stage=new_stage(operation)
            if operation=='redocking':
                stage['parameters'].update(pdb_codes=[['1ABC','ATP',1,'A']],preparation_pairs={'1ABC|ATP|1|A':config})
                effective=execution_stage(stage,'curated')['parameters']['preparation_pairs']['1ABC|ATP|1|A']
            else:
                stage['parameters']['preparation_options']=config
                effective=execution_stage(stage,'curated')['parameters']['preparation_options']
            self.assertEqual(effective,actual)

    def test_guided_forms_expose_only_chimera_and_gas_am1(self):
        for operation in ('prepare_structures','docking_vina','docking_dock6','redocking'):
            with self.subTest(operation=operation):
                stage=new_stage(operation)
                ui=SimpleNamespace(current={'pipeline':[stage]},service=SimpleNamespace(dock6_path=None),page=Page())
                form=GuidedForm(ui,stage,[],True)
                preparation=form.redocking_pairs.preparation if operation=='redocking' else form.preparation_settings
                self.assertEqual([o.key for o in preparation.backend.options],['legacy_chimera'])
                controls=list(walk(preparation.control))
                self.assertFalse(any(hasattr(preparation,name) for name in ('ambertools','ambertools_env','vina_policy')))
                self.assertFalse(any('AmberTools' in str(getattr(c,'label','')) for c in controls))
                for role in ('receptor','ligand'):
                    self.assertEqual([o.key for o in preparation.settings[role]['charge_type'].options],['gas','am1'])
                self.assertEqual(preparation.read()['preparation_backend'],'legacy_chimera')

class IdleFeedbackTests(unittest.IsolatedAsyncioTestCase):
    async def test_periodic_update_checks_never_show_splash_over_idle_configuration(self):
        page=Page();ui=SimpleNamespace(page=page,token=None,tr=lambda s:s)
        async def call(function,*args,**kwargs):return await WorkspaceUI.call(ui,function,*args,**kwargs)
        ui.call=call
        configuration=ft.AlertDialog(title=ft.Text('Configuração'));page.show_dialog(configuration)
        checker=SimpleNamespace(check=Mock(return_value=None))
        notice=UpdateNotice(ui,checker)
        sleeps=AsyncMock(side_effect=[None,None,asyncio.CancelledError()])
        with patch('biomolexplorer.ui.updates.asyncio.sleep',sleeps):
            with self.assertRaises(asyncio.CancelledError):await notice.poll()
        self.assertEqual(checker.check.call_count,3)
        self.assertTrue(configuration.open)
        self.assertFalse(hasattr(ui,'_configuration_loading'))
        self.assertEqual(len(page.dialogs),1)

class ChimeraCacheRoutingTests(unittest.TestCase):
    def test_preparation_keys_never_request_native_topology_or_ambertools(self):
        import tempfile
        from pathlib import Path
        from biomolexplorer.structure_cache import StructureCache
        from biomolexplorer.redocking_cache import PairCache
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);source=root/'source';source.mkdir()
            (source/'1ABC.pdb').write_text('END\n')
            row=['1ABC','LIG',1,'A'];settings={'1ABC|LIG|1|A':{'preparation_backend':'python'}}
            with patch('biomolexplorer.amber_preparation.parameterization_identity',side_effect=AssertionError('AmberTools called')),\
                 patch('biomolexplorer.structure_selection.component_path',side_effect=AssertionError('Native topology called')),\
                 patch('biomolexplorer.structure_cache.preparation_identity',return_value={'chimera':'test'}):
                StructureCache(source,root/'prepared',settings,'gas',7.4).keys(row)
                PairCache(source,root/'output',[row],settings,dict(prepare_complex=True,pH=7.4,charge_type='gas',sizeof_box=[24]*3,exhaustiveness=1,num_modes=1)).keys(row)

    def test_native_history_is_not_adopted_as_chimera_preparation(self):
        import tempfile
        from pathlib import Path
        from biomolexplorer.redocking_cache import legacy_indexes
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'structures'/'Target'/'pdb_codes.csv';path.parent.mkdir(parents=True)
            path.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RESOLUTION,RMSD\n1ABC,LIG,1,A,2.0,0.1\n')
            config=new_stage('redocking');config['parameters']['preparation_pairs']={'1ABC|LIG|1|A':{'preparation_backend':'python'}}
            self.assertEqual(list(legacy_indexes(dict(configuration=config,artifacts=[str(path)]))),[])

class SeparatePreparationFolderTests(unittest.TestCase):
    def test_prepared_files_are_validated_in_the_output_folder_not_the_raw_source(self):
        import tempfile
        from pathlib import Path
        from contextlib import nullcontext
        from biomolexplorer.redocking_preparation import prepare_pairs
        atom='ATOM      1  CA  ALA A   2       0.000   0.000   0.000  1.00  0.00           C\n'
        ligand='HETATM    2  C   LIG A   1       1.000   0.000   0.000  1.00  0.00           C\n'
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);source=root/'raw';source.mkdir();destination=root/'results'/'Prepared';destination.mkdir(parents=True)
            (source/'1ABC.pdb').write_text(atom+ligand)
            row=['1ABC','LIG',1,'A']
            def prepare(*args,**kwargs):
                (destination/'1ABC_A.dockprep.pdbqt').write_text(atom)
                (destination/'1ABC_LIG_1A.lig.pdbqt').write_text(ligand)
                (destination/'centers.csv').write_text('1ABC_LIG_1A\n1\n0\n0\n')
                return [row]
            dock=SimpleNamespace(complexpath=str(source),outputpath=str(destination),logger=Mock(),prepare_for_docking=Mock(side_effect=prepare))
            cache=Mock();cache.locked.return_value=nullcontext();cache.restore.return_value=(False,False)
            with patch('biomolexplorer.structure_cache.StructureCache',return_value=cache):
                completed,failures=prepare_pairs(dock,[row],{},'gas',7.4)
            self.assertEqual(completed,[row]);self.assertEqual(failures,[])
            self.assertFalse((source/'Prepared').exists())

    def test_minimization_does_not_run_implicit_dock_prep_or_replace_selected_charges(self):
        from biomolexplorer.redocking_config import configure_template
        from biomolexplorer.templates import RESOURCE_ROOT
        for role in ('receptor','ligand'):
            name='prepare_'+role+'.template'
            for method in ('gas','am1'):
                script=configure_template((RESOURCE_ROOT/'chimera'/name).read_text(),name,{role:{'charge_type':method,'minimize':True}})
                self.assertIn('method '+method,script)
                self.assertIn('minimize spec #0 prep false',script)

    def test_docking_helper_without_options_still_routes_am1_through_chimera(self):
        import tempfile
        from pathlib import Path
        from biomolexplorer.docking_data import prepare_ligands,write_csv
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);prepared=root/'prepared.mol2'
            prepared.write_text('@<TRIPOS>MOLECULE\nA\n1 0 0 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n1 C 0 0 0 C.3 1 LIG 0.0\n@<TRIPOS>BOND\n')
            table=root/'prepared.csv';write_csv(table,[dict(molecule_chembl_id='A',canonical_smiles='C',prepared_mol2=str(prepared),prepared_origin='pose')])
            write_csv(root/'input.csv',[dict(molecule_chembl_id='A',canonical_smiles='C')])
            with patch('biomolexplorer.docking_preparation.prepare_candidates',return_value=table) as prepare,\
                 patch('biomolexplorer.molecular_adapters.AM1BCCEngine',side_effect=AssertionError('External AmberTools called')):
                outputs=prepare_ligands(root/'input.csv',root/'out','mol2',charge_type='am1')
            self.assertEqual(len(outputs),1)
            config=prepare.call_args.args[2]
            self.assertEqual(config['preparation_backend'],'legacy_chimera')
            self.assertEqual(config['ligand']['charge_type'],'am1')
