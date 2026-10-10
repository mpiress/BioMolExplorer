import tempfile
import unittest
from pathlib import Path
from biomolexplorer.artifact_choices import choices, matches_selector
from biomolexplorer.pipeline import select_input
from biomolexplorer.ui.file_selection import file_selector
from biomolexplorer.docking_inputs import dock6_variant_matches


class StableSelectorsTests(unittest.TestCase):
    def test_legacy_pdb_selection_survives_new_worker_job(self):
        old='individual/001_PDB/.biomolexplorer/jobs/old/artifacts/PDB/Estruturas/4M0E.pdb'
        new=Path('/project/runs/new/stage/individual/001_PDB/.biomolexplorer/jobs/new/artifacts/PDB/Estruturas/4M0E.pdb')
        self.assertEqual(select_input([new], 'base_input_path', old), new.parent)
        self.assertTrue(matches_selector(new, old))
        self.assertFalse(matches_selector(new, old.replace('4M0E', '4EY4')))

    def test_choices_preserve_batches_without_worker_identity(self):
        files=[Path(f'/project/stage/individual/{i}/.biomolexplorer/jobs/job{i}/artifacts/results/docking_results.csv') for i in ('001', '002')]
        names=choices({'id':'stage','operation':'docking_vina','batches':[{}],'artifacts':files})
        self.assertEqual(names, {f'individual/{i}/results/docking_results.csv' for i in ('001', '002')})
        for name in names:
            selected=[p for p in files if matches_selector(p,name)]
            self.assertEqual(len(selected), 1)
            self.assertEqual(select_input(files,'base_selected_mols',name),selected[0].parent)
        selector=file_selector(files[0],files)
        self.assertNotIn('job001',selector)
        self.assertTrue(matches_selector(str(files[0]).replace('job001','new-job'),selector))
        with self.assertRaisesRegex(ValueError,'mais de um'):
            select_input(files, 'base_selected_mols', 'docking_results.csv')

    def test_removed_batch_does_not_fall_back_to_another_same_named_file(self):
        file=Path('/project/individual/002/.biomolexplorer/jobs/new/artifacts/results/4M0E.pdb')
        with self.assertRaisesRegex(ValueError,'não está'):
            select_input([file],'base_input_path','individual/001/.biomolexplorer/jobs/old/artifacts/results/4M0E.pdb')

    def test_vina_and_dock6_result_materialization_accepts_old_worker_paths(self):
        from biomolexplorer.workspace import WorkspaceStore
        from biomolexplorer.catalog import new_stage
        from biomolexplorer.pipeline import PipelineService
        from biomolexplorer.docking_data import write_csv, read_compounds
        with tempfile.TemporaryDirectory() as temporary:
            store=WorkspaceStore(temporary)
            token=store.register('Owner','owner@example.com','test-password-123')
            project=store.create_project(token,'Docking review')['id']
            service=PipelineService(store)
            try:
                for engine in ('vina','dock6'):
                    with self.subTest(engine=engine):
                        producer=new_stage('docking_'+engine)
                        stage=new_stage('docking_'+('dock6' if engine=='vina' else 'vina'))
                        folder=store.project_dir(project)/'runs/new'/producer['id']/'individual/001/.biomolexplorer/jobs/new/artifacts/results'
                        folder.mkdir(parents=True)
                        pose=folder/('M1.pdbqt' if engine=='vina' else 'M1_scored.mol2')
                        pose.write_text('pose')
                        table=folder/'docking_results.csv'
                        write_csv(table,[dict(molecule_chembl_id='M1',canonical_smiles='CCO',receptor_id='1ABC_A',
                            engine=engine,score=-1,conformer_file=pose.name)])
                        old='individual/001/.biomolexplorer/jobs/old/artifacts/results/docking_results.csv'
                        path,_,_=service._materialize_inputs(project,stage,'base_selected_mols',
                            [{'stage':producer['id'],'selector':old}],{producer['id']:[str(table),str(pose)]})
                        rows=read_compounds(path/'compounds.csv')
                        self.assertEqual(rows[0]['molecule_chembl_id'],'M1')
                        self.assertTrue(Path(rows[0]['conformer_file']).is_file())
            finally:service.close()
