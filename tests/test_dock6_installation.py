"""DOCK6 uses one server installation during preflight and worker execution."""
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import test_pipeline_execution as execution
from biomolexplorer.catalog import new_stage
from biomolexplorer.docking_tools import dock6_root
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.ui.app import WorkspaceUI


class Dock6DiscoveryTests(unittest.TestCase):
    def installation(self,root):
        for name in ('dock6','grid','sphgen','sphere_selector','showbox'):
            file=root/'bin'/name;file.parent.mkdir(parents=True,exist_ok=True)
            file.write_text('#!/bin/sh\nexit 0\n');file.chmod(0o755)
        file=root/'parameters/vdw_AMBER_parm99.defn';file.parent.mkdir(parents=True);file.write_text('parameters')
        return root.resolve()

    def test_root_from_path_and_environment_and_configured_precedence(self):
        with tempfile.TemporaryDirectory() as folder:
            root=self.installation(Path(folder)/'dock6')
            with patch.dict(os.environ,{'BIOMOL_DOCK6_ROOT':'','PATH':str(root/'bin')}):
                self.assertEqual(dock6_root(),root)
            with patch.dict(os.environ,{'BIOMOL_DOCK6_ROOT':str(root),'PATH':''}):
                self.assertEqual(dock6_root(),root)
                self.assertEqual(dock6_root(Path(folder)/'configured'),Path(folder)/'configured')
            (root/'parameters/vdw_AMBER_parm99.defn').unlink()
            with patch.dict(os.environ,{'BIOMOL_DOCK6_ROOT':str(root),'PATH':''}),patch('pathlib.Path.home',return_value=Path(folder)):
                self.assertIsNone(dock6_root())

    def test_new_block_uses_vina_at_the_single_candidate_port(self):
        receptor,compounds,vina,stage=[new_stage(op) for op in ('redocking','retrieve_compounds','docking_vina','docking_dock6')]
        ui=SimpleNamespace(service=SimpleNamespace(dock6_path=Path('/engine/dock6')))
        WorkspaceUI.auto_bind(ui,stage,[receptor,compounds,vina])
        self.assertEqual(stage['bindings']['base_selected_mols']['stage'],vina['id'])
        self.assertEqual(stage['bindings']['base_input_path']['stage'],receptor['id'])
        self.assertNotIn('base_vina_path',stage['bindings'])


class Dock6PreflightTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown

    def test_saved_block_without_installation_field_uses_server_configuration(self):
        stage=new_stage('docking_dock6')
        stage['parameters'].update(target='Estruturas',mol_filename='molecules',
            preparation_options={'ligand':{'charge_type':'am1'}})
        stage['bindings']={'base_input_path':{'stage':new_stage('redocking')['id']},
                           'base_selected_mols':{'stage':new_stage('docking_vina')['id']}}
        engine=self.store.root/'dock6'
        service=PipelineService(self.store,dock6_path=engine)
        try:
            service._validate_parameters(self.project_id,stage,partial=True)
            params=service._runtime_parameters(stage)
            self.assertEqual(params['dock6_app_path'],str(engine))
            self.assertEqual(params['charge_type'],'am1')
            self.assertNotIn('dock6_app_path',stage['parameters'])
        finally:service.close()

    def test_missing_installation_reports_setup_instead_of_input_files(self):
        with patch('biomolexplorer.docking_tools.dock6_root',return_value=None):
            service=PipelineService(self.store)
        try:
            with self.assertRaisesRegex(ValueError,'Instalação DOCK6 não encontrada'):
                service._validate_parameters(self.project_id,new_stage('docking_dock6'),partial=True)
        finally:service.close()

if __name__=='__main__':unittest.main()
