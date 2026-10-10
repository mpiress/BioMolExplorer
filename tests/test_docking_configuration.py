"""Docking configuration follows prepared-input state and shares role controls."""
import unittest
from types import SimpleNamespace
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from biomolexplorer.catalog import new_stage, operation_fields
from biomolexplorer.flow import input_ports, compatible
from biomolexplorer.ui.guided import GuidedForm
from biomolexplorer.docking_data import prepare_ligands, write_csv
from wrappers.docking import prepare_docking_receptor


class DockingConfigurationTests(unittest.TestCase):
    def form(self, operation, origin='redocking'):
        source, stage=new_stage(origin),new_stage(operation)
        stage['bindings']['base_input_path']={'stage':source['id']}
        ui=SimpleNamespace(current={'pipeline':[source,stage]},service=SimpleNamespace(dock6_path=None),
                           page=SimpleNamespace(update=lambda:None))
        return GuidedForm(ui,stage,[],True),source,stage

    def test_prepared_receptor_hides_preparation_and_complex_for_both_engines(self):
        for operation in ('docking_vina','docking_dock6'):
            with self.subTest(operation=operation):
                form,source,stage=self.form(operation)
                form.layout()
                self.assertFalse(form.preparation_settings.groups['receptor'].visible)
                self.assertTrue(form.preparation_settings.groups['ligand'].visible)
                self.assertFalse(form.field_controls['pdb_code'].visible)
                self.assertFalse(any(n.startswith('chimera/') for n in form.template_readers))
                form.preparation_settings.settings['ligand']['charge_type'].value='am1'
                result=form.read()
                self.assertTrue(result['parameters']['receptor_prepared'])
                self.assertNotIn('pdb_code',result['parameters'])
                self.assertEqual(result['parameters']['preparation_options']['ligand']['charge_type'],'am1')
                reopened=GuidedForm(form.ui,result,[],True)
                self.assertEqual(reopened.read(),result)

    def test_switching_to_raw_pdb_restores_receptor_options_and_layout(self):
        for operation in ('docking_vina','docking_dock6'):
            form,source,stage=self.form(operation)
            raw=new_stage('retrieve_structures');form.ui.current['pipeline'].append(raw)
            form.layout()
            editor=form.input_editors['base_input_path']
            editor.rows[0][0].value='stage:'+raw['id'];editor.changed()
            self.assertTrue(form.preparation_settings.groups['receptor'].visible)
            self.assertTrue(form.field_controls['pdb_code'].visible)
            self.assertFalse(form.read()['parameters']['receptor_prepared'])
            editor.rows[0][0].value='stage:'+source['id'];editor.changed()
            self.assertFalse(form.preparation_settings.groups['receptor'].visible)
            self.assertFalse(form.field_controls['pdb_code'].visible)

    def test_dock6_has_one_candidate_port_accepting_vina(self):
        stage=new_stage('docking_dock6');vina=new_stage('docking_vina')
        self.assertNotIn('base_vina_path',{f['name'] for f in operation_fields('docking_dock6')})
        self.assertEqual({p['field'] for p in input_ports(stage)},{'base_input_path','base_selected_mols'})
        self.assertTrue(compatible(vina,stage,'base_selected_mols'))
        stage['bindings']['base_vina_path']={'stage':vina['id'],'selector':'docking_results.csv','compound_id':'M1'}
        ui=SimpleNamespace(current={'pipeline':[vina,stage]},service=SimpleNamespace(dock6_path=None),page=SimpleNamespace(update=lambda:None))
        form=GuidedForm(ui,stage,[],True);result=form.read()
        self.assertNotIn('base_vina_path',form.input_editors)
        self.assertNotIn('base_vina_path',result['bindings'])
        self.assertEqual(result['bindings']['base_selected_mols']['compound_id'],'M1')

    def test_receptor_preparation_runs_only_for_raw_inputs(self):
        options={'receptor':{'charge_type':'am1'},'ligand':{'charge_type':'gas'}}
        with patch('wrappers.redocking.prepare_structures') as prepare:
            self.assertEqual(prepare_docking_receptor('input','Target','out',None,7.4,options,True),'input')
            prepare.assert_not_called()
            result=prepare_docking_receptor('input','Target','out',['1ABC','LIG',1,'A'],7.4,options,False)
            self.assertEqual(result,str(Path('out/Receptor')))
            self.assertEqual(prepare.call_args.kwargs['pdb_codes'],[['1ABC','LIG',1,'A']])
            self.assertEqual(prepare.call_args.kwargs['preparation_options'],options)

    def test_ligand_settings_reach_preparation_and_prepared_inputs_are_reused(self):
        with TemporaryDirectory() as folder:
            root=Path(folder);source=root/'M1.pdbqt';source.write_text('prepared ligand')
            dataset=root/'input.csv';write_csv(dataset,[dict(molecule_chembl_id='M1',canonical_smiles='CCO',prepared_pdbqt=str(source))])
            options={'ligand':{'charge_type':'am1','minimize':False}}
            with patch('biomolexplorer.docking_preparation.subprocess.run') as runner:
                paths=prepare_ligands(dataset,root/'out','pdbqt',preparation_options=options)
                runner.assert_not_called()
                self.assertEqual(paths[0].read_text(),'prepared ligand')

import test_pipeline_execution as execution
from biomolexplorer.pipeline import PipelineService


class DockingRawInputTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    upload=execution.PipelineExecutionTests.upload

    def test_raw_pdb_is_materialized_and_marked_for_preparation(self):
        from test_docking_handoffs import ATOM
        for operation in ('docking_vina','docking_dock6'):
            stage=new_stage(operation)
            receptor=self.upload('1ABC.pdb',ATOM,'structures')
            compounds=self.upload('compounds.csv','name,smiles\nM1,CCO\n')
            stage['bindings']={'base_input_path':{'asset':receptor},'base_selected_mols':{'asset':compounds}}
            stage['parameters']['pdb_code']=['1ABC','LIG',1,'A']
            if operation=='docking_dock6':stage['parameters']['charge_type']='gas'
            service=PipelineService(self.store,dock6_path=self.store.project_dir(self.project_id))
            try:
                params=service._resolve(self.project_id,self.store.user(self.token)['id'],stage,{})
                self.assertFalse(params['receptor_prepared'])
                self.assertTrue((Path(params['base_input_path'])/'MeuAlvo/1ABC.pdb').is_file())
                self.assertTrue((Path(params['base_selected_mols'])/(params['mol_filename']+'.csv')).is_file())
            finally:service.close()

if __name__=='__main__':unittest.main()
