"""Scientific mode boundaries, coordinate preservation and mode-scoped caches."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from biomolexplorer.catalog import new_stage,operation_fields
from biomolexplorer.dock6_protocol import mode_for_row,apply_policy
from caad.docking import Dock6
import test_dock6_incremental as incremental


class ModeTests(unittest.TestCase):
    def test_auto_uses_bound_provenance_and_explicit_references_require_coordinates(self):
        self.assertEqual(mode_for_row('auto',dict(canonical_smiles='CCO')),'denovo')
        self.assertEqual(mode_for_row('auto',dict(conformer_file='library.mol2',prepared_origin='library')),'denovo')
        self.assertEqual(mode_for_row('auto',dict(conformer_file='vina.pdbqt',engine='vina',receptor_id='1ABC_A')),'guided')
        self.assertEqual(mode_for_row('restrained',dict(conformer_file='crystal.mol2')),'restrained')
        with self.assertRaises(ValueError):mode_for_row('score_only',dict(canonical_smiles='CCO'))
        with self.assertRaises(ValueError):mode_for_row('invalid',dict(conformer_file='pose.mol2'))

    def test_restricted_modes_override_unsafe_template_switches_without_duplicates(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'run.in'
            for mode in ('guided','score_only','restrained'):
                path.write_text('conformer_search_type flex\norient_ligand yes\norient_ligand yes\nminimize_ligand yes\n')
                apply_policy(path,mode)
                lines=[line.split() for line in path.read_text().splitlines()]
                self.assertEqual(sum(parts[0]=='orient_ligand' for parts in lines),1)
                values=dict(lines)
                self.assertEqual(values['conformer_search_type'],'rigid')
                self.assertEqual(values['orient_ligand'],'no')
                self.assertEqual(values['minimize_ligand'],'no' if mode=='score_only' else 'yes')
                self.assertEqual(values['simplex_restraint_min'],'yes' if mode=='restrained' else 'no')

    def test_configuration_defaults_preserve_independent_redocking_and_auto_docking(self):
        self.assertEqual(new_stage('docking_dock6')['parameters']['dock6_mode'],'auto')
        self.assertEqual(new_stage('redocking')['parameters']['dock6_mode'],'redock')
        choices=next(field['choices'] for field in operation_fields('redocking') if field['name']=='dock6_mode')
        self.assertEqual(choices,['redock','score_only','restrained'])

    def test_smiles_cannot_start_reference_evaluation_or_receptor_preparation(self):
        from wrappers.docking import perform_docking_dock6
        from biomolexplorer.docking_data import write_csv
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            write_csv(root/'compounds.csv',[dict(molecule_chembl_id='M',canonical_smiles='CCO')])
            with patch('wrappers.docking.prepare_docking_receptor',side_effect=AssertionError('Preflight must run first')):
                with self.assertRaisesRegex(ValueError,'requer uma pose'):
                    perform_docking_dock6(str(root),'Target',str(root/'out'),str(root),'/dock6','gas','compounds',dock6_mode='score_only')


class CacheModeTests(unittest.TestCase):
    setUp=incremental.Dock6IncrementalTests.setUp

    def test_same_grid_is_shared_but_search_and_refinement_results_are_distinct(self):
        # Use the existing fixture's real file preparation and cache machinery.
        root=self.root
        ligands=root/'inputs';ligands.mkdir(exist_ok=True)
        from test_dock6_incremental import MOL2
        (ligands/'A.lig.mol2').write_text(MOL2.format(code='A'))
        receptor=root/'rec'
        if not receptor.exists():
            receptor.mkdir()
            (receptor/'1ABC_A.dockprep.mol2').write_text(MOL2.format(code='R'))
            (receptor/'1ABC_A.noH.pdb').write_text('ATOM\n')
        engine=Dock6(dock6_path=str(self.install),ligand_input_path=str(ligands),
            receptor_input_path=str(receptor),base_output_path=str(root/'mode-out'),pdb_code='1ABC_A',binding_site_center=[0.,0.,0.])
        try:
            checkpoint=engine._checkpoints();keys={};sites=[]
            for mode in ('guided','score_only','restrained','denovo'):
                engine.compound_modes={'A':mode}
                keys[mode]=checkpoint.compound_key('perform_parallel_docking','A')
                sites.append(checkpoint.site_key())
            self.assertEqual(len(set(keys.values())),4)
            self.assertEqual(len(set(sites)),1)
            engine.compound_modes={'A':'guided'}
            with patch.object(engine,'_run_compound_jobs') as run:
                engine.prepare_minimization()
                self.assertEqual(run.call_args.args[1],[])
        finally:engine.export_results()


if __name__=='__main__':
    unittest.main()
