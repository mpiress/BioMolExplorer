"""Pair curation reaches scientific preparation without repeated options."""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock
from biomolexplorer.catalog import new_stage
from biomolexplorer.redocking_config import configure_template, metadata_records, pair_key, validate_pairs
from biomolexplorer.templates import RESOURCE_ROOT
from biomolexplorer.ui.guided import GuidedForm
from caad.docking import Docking


class RedockingConfigurationTests(unittest.TestCase):
    def ui(self, stages, records=None):
        return SimpleNamespace(current={'pipeline': stages}, redocking_records=records or {},
            service=SimpleNamespace(dock6_path=None), page=SimpleNamespace(update=lambda: None))

    def test_selection_required_even_when_metadata_is_available(self):
        source, stage = new_stage('retrieve_structures'), new_stage('redocking')
        stage['bindings']['base_input_path'] = {'stage': source['id']}
        ui = self.ui([source, stage], {'stage:'+source['id']: [['1ABC','LIG',12,'A',1.7]]})
        form = GuidedForm(ui, stage, [], True)
        self.assertEqual(len(form.redocking_pairs.available.options), 1)
        with self.assertRaisesRegex(ValueError, 'Selecione pelo menos um par'): form.read()
        form.redocking_pairs.available.value = '0'
        form.redocking_pairs.choose(None)
        result = form.read()
        self.assertEqual(result['parameters']['pdb_codes'], [['1ABC','LIG',12,'A',1.7]])
        self.assertNotIn('charge_type', result['parameters'])
        self.assertFalse(any(n.startswith('chimera/') for n in form.template_readers))

    def test_pair_settings_survive_reopening_and_chain_change(self):
        stage = new_stage('redocking')
        stage['parameters']['pdb_codes'] = [['1ABC','LIG',12,'A',1.7]]
        form = GuidedForm(self.ui([stage]), stage, [], True)
        fields, _, (has_cofactors, cofactors), options, _ = form.redocking_pairs.rows[0]
        fields[3].value = 'B'; has_cofactors.value = True; cofactors.value = 'FAD, MG'
        options['ligand']['add_hydrogens'].value = False
        options['receptor']['charge_type'].value = 'am1'
        result = form.read()
        self.assertEqual(result['parameters']['pdb_codes'], [['1ABC','LIG',12,'B']])
        config = result['parameters']['preparation_pairs']['1ABC|LIG|12|B']
        self.assertNotIn('ligand_chain', config)
        self.assertEqual(config['cofactors'], ['FAD', 'MG'])
        self.assertEqual(GuidedForm(self.ui([result]), result, [], True).read(), result)

    def test_cofactor_codes_and_conflicting_receptors_are_rejected(self):
        records = [['1ABC','LIG',12,'A'], ['1ABC','ATP',13,'A']]
        with self.assertRaisesRegex(ValueError, 'mesmos cofatores'):
            validate_pairs(records, {pair_key(records[0]): {'cofactors': ['FAD']}})
        with self.assertRaisesRegex(ValueError, 'cofatores válidos'):
            validate_pairs(records[:1], {pair_key(records[0]): {'cofactors': ['FAD;delete']}})

    def test_generated_scripts_keep_cofactors_and_use_exact_ligand_chain(self):
        record = ['1ABC','LIG',12,'BC']
        config = {'cofactors': ['FAD'], 'ligand_chain': 'D',
                  'ligand': {'add_hydrogens': False, 'minimize': False, 'charge_type': 'am1'},
                  'receptor': {'remove_solvent': False}}
        with tempfile.TemporaryDirectory() as folder:
            dock = object.__new__(Docking)
            dock.outputpath = folder; dock.complexpath = folder; dock.logger = MagicMock()
            dock.process_in_parallel = MagicMock()
            dock.prepare_on_obabel = MagicMock()
            dock.calculate_ligand_centerofmass = MagicMock(return_value=[1,2,3])
            dock.prepare_for_docking([record], 'gas', 7.4, False, {pair_key(record): config})
            complex_script = Path(folder, 'prepare_complex_1ABC_BC.com').read_text()
            self.assertIn('select #0:.BC | :FAD', complex_script)
            self.assertNotIn('delete solvent', complex_script)
            receptor_script = Path(folder, 'prepare_receptor_1ABC_BC.com').read_text()
            self.assertIn('select protein | :FAD | solvent', receptor_script)
            self.assertNotIn('delete ligand', receptor_script)
            self.assertNotIn('delete solvent', receptor_script)
            ligand_script = Path(folder, 'prepare_ligand_1ABC_LIG_12BC.com').read_text()
            self.assertIn('select :12.BC', ligand_script)
            self.assertIn('method am1', ligand_script)
            self.assertNotIn('\naddh\n', ligand_script)
            self.assertNotIn('minimize', ligand_script)
            ligand_call = dock.prepare_on_obabel.call_args_list[-1]
            self.assertEqual(ligand_call.args[2], [])
        conformation = configure_template((RESOURCE_ROOT/'chimera/prepare_better_conform.template').read_text(),
                                          'prepare_better_conform.template', config)
        self.assertIn('method am1', conformation)
        self.assertNotIn('\naddh\n', conformation)
        self.assertNotIn('minimize', conformation)

    def test_selection_can_be_deferred_until_retrieval_but_not_execution(self):
        from biomolexplorer.operations import validate_operation
        params={'base_input_path':'pending','target':'MeuAlvo'}
        validate_operation('redocking',params,defer_redocking_selection=True)
        with self.assertRaisesRegex(ValueError,'Selecione pelo menos um par'):
            validate_operation('redocking',params)

    def test_no_cofactor_and_metadata_resolution(self):
        script = configure_template((RESOURCE_ROOT/'chimera/prepare_complex.template').read_text(), 'prepare_complex.template', {})
        self.assertIn('select #0:.{chain}\n', script)
        self.assertNotIn(':FAD', script)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder, 'pdb_codes.csv')
            path.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RESOLUTION\n1ABC,LIG,12,A,1.7\n')
            self.assertEqual(metadata_records([path]), [['1ABC','LIG',12,'A',1.7]])

    def test_pandas_records_reach_script_generation_and_engine_errors_propagate(self):
        import pandas as pd
        from unittest.mock import patch
        records = pd.DataFrame([['4M0E','1YL',604,'A',2.0]],
            columns=['PDB_CODE','LIGAND','RESNUM','CHAIN','RESOLUTION']).to_records(index=False)
        self.assertEqual(pair_key(records[0]), '4M0E|1YL|604|A')
        with tempfile.TemporaryDirectory(prefix='prepared structures ') as folder:
            dock = object.__new__(Docking)
            dock.outputpath = folder; dock.complexpath = folder; dock.logger = MagicMock()
            dock.process_in_parallel = MagicMock(side_effect=RuntimeError('Chimera failed'))
            with self.assertRaisesRegex(RuntimeError, 'Chimera failed'):
                dock.prepare_for_docking(records, 'am1', 7.4, True, {})
            source = Path(folder,'prepare_ligand_4M0E_1YL_604A.com').read_text()
            self.assertIn('select :604.A', source)
            self.assertIn('method am1', source)
            self.assertIn('open '+folder+'/', source)
            self.assertNotIn('open "', source)
            self.assertTrue(dock.logger.exception.called)
            dock.process_in_parallel.side_effect = None
            dock.prepare_on_obabel = MagicMock(return_value=True)
            dock.calculate_ligand_centerofmass = MagicMock(return_value=[1,2,3])
            prepared = dock.prepare_for_docking(records, 'am1', 7.4, True, {})
            self.assertEqual(prepared, [('4M0E','1YL',604,'A',2.0)])
            self.assertTrue(Path(folder,'centers.csv').is_file())

    def test_preflight_checks_residue_chain_cofactors_and_tools(self):
        from biomolexplorer.redocking_config import validate_structure_pairs, validate_redocking_tools
        from unittest.mock import patch
        record = ['1ABC','LIG',1,'A']
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder,'1ABC.pdb')
            path.write_text('ATOM      1  C   ALA A   2       0.000   0.000   0.000  1.00  0.00           C\n'
                            'HETATM    2  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n')
            validate_structure_pairs(folder,[record],{})
            with self.assertRaisesRegex(ValueError,'Cadeia do receptor'):
                validate_structure_pairs(folder,[['1ABC','LIG',1,'B']],{})
            with self.assertRaisesRegex(ValueError,'Ligante ou resíduo'):
                validate_structure_pairs(folder,[['1ABC','LIG',99,'A']],{})
            with self.assertRaisesRegex(ValueError,'Cofator não encontrado'):
                validate_structure_pairs(folder,[record],{pair_key(record):{'cofactors':['FAD']}})
        with patch('shutil.which',return_value=None):
            with self.assertRaisesRegex(ValueError,'chimera'):
                validate_redocking_tools()

    def test_verbose_is_hidden_and_overrides_cannot_enable_vina_verbosity(self):
        from biomolexplorer.catalog import operation_fields
        from biomolexplorer.templates import materialize_templates, quiet_template
        normalized = quiet_template('vina/config.template', '  Verbosity = 9\nsize_x = 24\nverbose = 2\n')
        self.assertEqual(normalized.count('verbosity = 0'), 1)
        self.assertIn('size_x = 24', normalized)
        self.assertNotIn('= 9', normalized)
        self.assertNotIn('= 2\n', normalized)
        for operation in ('redocking','docking_vina','admet','retrieve_zinc'):
            self.assertNotIn('verbose',[field['name'] for field in operation_fields(operation)])
        stage = new_stage('redocking');stage['parameters']['pdb_codes']=[['1ABC','LIG',1,'A']]
        stage['parameters']['verbose']=True
        stage['templates']['vina/config.template']=(RESOURCE_ROOT/'vina/config.template').read_text().replace('verbosity = 0','verbosity = 9')
        form=GuidedForm(self.ui([stage]),stage,[],True)
        self.assertNotIn('verbose',form.read()['parameters'])
        self.assertIn('verbosity = 0',form.template_readers['vina/config.template']())
        with tempfile.TemporaryDirectory() as folder:
            resource=materialize_templates(Path(folder)/'resources',stage['templates'])
            self.assertIn('verbosity = 0',(resource/'vina/config.template').read_text())

    def test_vina_requires_output_and_finite_rmsd_before_updating_metadata(self):
        import pandas as pd
        from unittest.mock import patch
        from caad.docking import DockVina
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);prepared=root/'Prepared';prepared.mkdir();output=root/'Vina';output.mkdir()
            records=pd.DataFrame([['1ABC','LIG',1,'A',2.0]],columns=['PDB_CODE','LIGAND','RESNUM','CHAIN','RESOLUTION'])
            dock=DockVina(ligand_input_path=str(prepared),receptor_input_path=str(prepared),
                complex_input_path=str(root),output_path=str(output),pdb_codes=records)
            dock.retrieve_centerofmass_dataset=MagicMock(return_value=[0,0,0])
            dock.generate_docking_script=MagicMock();dock.perform_vina_evaluation=MagicMock()
            (prepared/'1ABC_LIG_1A.lig.pdbqt').write_text('reference')
            with self.assertRaisesRegex(ValueError,'Saída de redocking ausente'):
                dock.redocking(7.4)
            self.assertFalse((root/'pdb_codes.csv').exists())
            (output/'1ABC_LIG_1A.lig.pdbqt').write_text('pose')
            with patch('caad.docking.Descriptors.calcRMSD',return_value=float('nan')):
                with self.assertRaisesRegex(ValueError,'RMSD inválido'):dock.redocking(7.4)
            with patch('caad.docking.Descriptors.calcRMSD',return_value=1.25):dock.redocking(7.4)
            self.assertEqual(pd.read_csv(root/'pdb_codes.csv')['RMSD'].tolist(),[1.25])
