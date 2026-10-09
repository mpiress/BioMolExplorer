"""Docking identity, batch intersection, pose handoffs and curated result access."""
import csv
import json
import pickle
import time
import unittest
from pathlib import Path
from unittest.mock import patch

import test_pipeline_execution as execution
from biomolexplorer.catalog import new_stage
from biomolexplorer.docking_data import read_results, write_csv, consensus_rows
from biomolexplorer.docking_results import DockingResults
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.stage_cache import artifact_manifest
from biomolexplorer.workspace import AccessDenied


class DockingInteroperabilityTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    upload=execution.PipelineExecutionTests.upload

    def results(self, folder, engine, codes=('A','B'), receptor='1ABC_A'):
        folder.mkdir(parents=True,exist_ok=True)
        rows=[]
        for i,code in enumerate(codes):
            pose=folder/(code+('.pdbqt' if engine=='vina' else '_scored.mol2'))
            pose.write_text('pose '+code+'\n')
            rows.append(dict(molecule_chembl_id=code,canonical_smiles='CCO',receptor_id=receptor,
                             engine=engine,score=-7-i,conformer_file=pose.name))
        write_csv(folder/'docking_results.csv',rows)
        return list(folder.iterdir())

    def test_all_batches_feed_consensus_without_nested_duplicates(self):
        root=self.store.project_dir(self.project_id)
        files=self.results(root/'batch1','vina',('A',))+self.results(root/'batch2','vina',('B',))
        nested=self.results(root/'batch1'/'receptor','vina',('REMOVED',))
        source=new_stage('docking_vina');stage=new_stage('consensus')
        service=PipelineService(self.store)
        try:
            path,_,_=service._materialize_inputs(self.project_id,stage,'base_vina_path',
                [{'stage':source['id'],'selector':'auto'}],{source['id']:[str(p) for p in files+nested]})
            vina=read_results(path,'vina')
            dock6=self.results(root/'dock6','dock6',('A','B','C'))
            rows=consensus_rows(vina,read_results(root/'dock6','dock6'))
            self.assertEqual([r['molecule_chembl_id'] for r in rows],['A','B'])
            self.assertEqual(len(vina),2)
        finally:service.close()

    def test_both_engines_accept_actual_poses_and_keep_compound_identity(self):
        root=self.store.project_dir(self.project_id)
        for source_engine,target_engine in (('vina','dock6'),('dock6','vina')):
            files=self.results(root/source_engine,source_engine)
            source=new_stage('docking_'+source_engine);target=new_stage('docking_'+target_engine)
            service=PipelineService(self.store)
            try:
                path,stem,_=service._materialize_inputs(self.project_id,target,'base_selected_mols',
                    [{'stage':source['id'],'selector':'auto'}],{source['id']:[str(p) for p in files]})
                with (path/(stem+'.csv')).open() as stream:rows=list(csv.DictReader(stream))
                self.assertEqual([r['molecule_chembl_id'] for r in rows],['A','B'])
                for row in rows:
                    self.assertEqual(Path(row['conformer_file']).read_text(),'pose '+row['molecule_chembl_id']+'\n')
                from biomolexplorer.docking_data import prepare_ligands
                with patch('biomolexplorer.docking_data.convert_structure') as convert,patch('biomolexplorer.visualizations.molecule_sdf') as generate:
                    prepare_ligands(path/(stem+'.csv'),root/('prepared_'+target_engine),'mol2' if target_engine=='dock6' else 'pdbqt')
                    self.assertEqual(convert.call_count,2);generate.assert_not_called()
                    self.assertEqual(convert.call_args_list[0].args[0],Path(rows[0]['conformer_file']))
            finally:service.close()

    def test_independent_user_compounds_accept_aliases_for_both_engines(self):
        asset=self.upload('isolated.csv','name,smiles\nCUSTOM1,CCO\nCUSTOM2,CCC\n')
        service=PipelineService(self.store)
        try:
            for engine in ('vina','dock6'):
                path,stem,_=service._materialize_inputs(self.project_id,new_stage('docking_'+engine),
                    'base_selected_mols',[{'asset':asset}],{})
                with (path/(stem+'.csv')).open() as stream:rows=list(csv.DictReader(stream))
                self.assertEqual([r['molecule_chembl_id'] for r in rows],['CUSTOM1','CUSTOM2'])
                self.assertEqual([r['canonical_smiles'] for r in rows],['CCO','CCC'])
        finally:service.close()

    def test_consensus_requires_matching_receptors_and_structures(self):
        root=self.store.project_dir(self.project_id)
        self.results(root/'vina','vina',('A',));self.results(root/'dock6','dock6',('A',),'2DEF_A')
        a,b=read_results(root/'vina','vina'),read_results(root/'dock6','dock6')
        self.assertEqual(consensus_rows(a,b),[])
        b[0]['receptor_id']=a[0]['receptor_id'];b[0]['canonical_smiles']='CCC'
        with self.assertRaisesRegex(ValueError,'estruturas diferentes'):consensus_rows(a,b)

    def test_generated_ligand_is_centered_without_changing_mol2_charges(self):
        from biomolexplorer.docking_data import center_mol2
        path=self.store.project_dir(self.project_id)/'generated.mol2'
        path.write_text('@<TRIPOS>MOLECULE\nM1\n2 1 0 0 0\nSMALL\nUSER_CHARGES\n'
            '@<TRIPOS>ATOM\n1 C 0 0 0 C.3 1 LIG -0.2\n2 O 2 0 0 O.3 1 LIG 0.2\n'
            '@<TRIPOS>BOND\n1 1 2 1\n')
        center_mol2(path,[10.,20.,30.])
        self.assertIn('1 C 9.0000 20.0000 30.0000 C.3 1 LIG -0.2',path.read_text())
        self.assertIn('2 O 11.0000 20.0000 30.0000 O.3 1 LIG 0.2',path.read_text())
        self.assertIn('@<TRIPOS>BOND\n1 1 2 1',path.read_text())

    def test_parallel_tool_failure_preserves_diagnostic(self):
        from biomolexplorer.processes import ScientificToolError
        error=pickle.loads(pickle.dumps(ScientificToolError('grid failed','TOOL_EXIT_FAILED','Inspect grid output')))
        self.assertEqual(str(error),'grid failed')
        self.assertEqual(error.error_code,'TOOL_EXIT_FAILED')
        self.assertEqual(error.action,'Inspect grid output')

    def test_dock6_stages_long_paths_and_exports_artifacts(self):
        from caad.docking import Dock6
        root=self.store.project_dir(self.project_id)/('long path '*15)
        ligands,receptor,output=root/'ligands',root/'receptor',root/'results'
        ligands.mkdir(parents=True);receptor.mkdir()
        (ligands/'A.lig.mol2').write_text('input ligand')
        (receptor/'1ABC_A.noH.pdb').write_text('input receptor')
        engine=Dock6(dock6_path=str(root/'engine'),ligand_input_path=str(ligands),
                     receptor_input_path=str(receptor),base_output_path=str(output),pdb_code='1ABC_A')
        runtime=Path(engine.runtime_path)
        try:
            self.assertLess(len(engine.outputpath),50)
            self.assertEqual((Path(engine.receptorpath)/'1ABC_A.noH.pdb').read_text(),'input receptor')
            self.assertTrue(engine._Dock6__dock6_path.endswith('/'))
            result=Path(engine.outputpath)/'rigid'/'A_scored.mol2';result.parent.mkdir()
            result.write_text('docked geometry')
            engine.export_results()
            self.assertEqual((output/'rigid'/'A_scored.mol2').read_text(),'docked geometry')
            self.assertFalse(runtime.exists())
        finally:
            import shutil
            if runtime.exists():shutil.rmtree(runtime)

    def test_individual_consensus_pose_keeps_identity_and_selected_geometry(self):
        from biomolexplorer.docking_data import input_records
        root=self.store.project_dir(self.project_id)
        a,b=root/'0_vina_pose.pdbqt',root/'0_dock6_pose.mol2'
        a.write_text('Vina geometry');b.write_text('DOCK6 geometry')
        table=root/'consensus.csv'
        write_csv(table,[dict(molecule_chembl_id='CUSTOM1',canonical_smiles='CCO',
                             conformer_file=b.name,vina_pose=a.name,dock6_pose=b.name,vina=-7,dock6=-5)])
        with patch('biomolexplorer.docking_data.subprocess.run') as convert:
            for pose in (a,b):
                rows=input_records([pose],[table,a,b])
                self.assertEqual(rows[0]['molecule_chembl_id'],'CUSTOM1')
                self.assertEqual(rows[0]['conformer_file'],str(pose))
            convert.assert_not_called()

    def test_curated_summaries_are_authoritative_and_preview_is_scoped(self):
        root=self.store.project_dir(self.project_id)
        files=self.results(root/'batch1','dock6',('A',))+self.results(root/'batch2','dock6',('B',))
        files+=self.results(root/'batch1'/'receptor','dock6',('A',))
        stage=new_stage('docking_dock6');rid='b'*32
        item=dict(id=stage['id'],operation=stage['operation'],configuration=stage,status='succeeded',
                  artifacts=[str(p) for p in files],artifact_manifest=artifact_manifest([str(p) for p in files]))
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',(rid,self.project_id,self.store.user(self.token)['id'],
                'succeeded',json.dumps([item]),time.time(),time.time(),None))
        service=DockingResults(self.store);args=(self.token,self.project_id,rid,stage['id'])
        tables=service.tables(*args);self.assertEqual(len(tables),2)
        table=str(root/'batch1'/'docking_results.csv');page=service.page(*args,table)
        self.assertEqual(Path(service.pose(*args,table,0,page['version'])).read_text(),'pose A\n')
        with self.assertRaises(AccessDenied):service.pose(*args,table,0,page['version'],'unauthorized')
        service.remove(*args,table,0,page['version'])
        self.assertEqual(service.page(*args,table)['total'],0)
        # Neither the nested receptor CSV nor its pose can resurrect A.
        self.assertEqual([r['molecule_chembl_id'] for r in read_results(root,'dock6')],['B'])
        consumer=new_stage('consensus');pipeline=PipelineService(self.store)
        try:
            path,_,_=pipeline._materialize_inputs(self.project_id,consumer,'base_dock6_path',
                [{'stage':stage['id'],'selector':'auto'}],{stage['id']:item['artifacts']})
            self.assertEqual([r['molecule_chembl_id'] for r in read_results(path,'dock6')],['B'])
        finally:pipeline.close()


if __name__=='__main__':unittest.main()
