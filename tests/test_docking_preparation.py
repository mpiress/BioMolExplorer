"""Independent compound sources and prepared receptors at pipeline boundaries."""
import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch

import test_pipeline_execution as execution
import test_docking_handoffs as handoffs
from biomolexplorer.catalog import new_stage
from biomolexplorer.docking_data import read_compounds,write_csv
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.ui.file_selection import FileSelection
from biomolexplorer.docking_preparation import validate_receptor_outputs,prepare_candidates


class ReceptorOutputTests(unittest.TestCase):
    def test_engine_formats_and_finite_centers_are_required(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);records=[['1ABC','LIG',1,'A']]
            (root/'1ABC_A.dockprep.pdbqt').write_text('prepared receptor')
            (root/'centers.csv').write_text('1ABC_LIG_1A\n1\n2\n3\n')
            validate_receptor_outputs(root,records,'vina')
            for engine in ('dock6','both'):
                with self.assertRaisesRegex(ValueError,'dockprep.mol2'):
                    validate_receptor_outputs(root,records,engine)
            (root/'1ABC_A.dockprep.mol2').write_text('charged receptor')
            with self.assertRaisesRegex(ValueError,'noH.pdb'):
                validate_receptor_outputs(root,records,'both')
            (root/'1ABC_A.noH.pdb').write_text('surface receptor')
            validate_receptor_outputs(root,records,'both')
            (root/'centers.csv').write_text('1ABC_LIG_1A\n1\nNaN\n3\n')
            with self.assertRaisesRegex(ValueError,'Centro.*1ABC_A'):
                validate_receptor_outputs(root,records,'both')

    def test_incomplete_redocking_bundle_fails_before_candidate_preparation(self):
        from wrappers.redocking import prepare_structures
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);source=root/'input/Target/Prepared';source.mkdir(parents=True)
            (source/'1ABC_A.dockprep.pdbqt').write_text('prepared receptor')
            with patch('wrappers.redocking.Docking') as receptor,patch('biomolexplorer.docking_preparation.prepare_candidates') as candidates:
                with self.assertRaisesRegex(ValueError,'dockprep.mol2'):
                    prepare_structures(str(root/'input'),'Target',str(root/'output'),
                        pdb_codes=[['1ABC','LIG',1,'A']],receptor_prepared=True,
                        base_selected_mols=str(root/'compounds'),docking_engines='both')
                receptor.assert_not_called();candidates.assert_not_called()

    def test_invalid_engine_does_not_read_or_prepare_candidates(self):
        with patch('biomolexplorer.docking_preparation.read_compounds') as read:
            with self.assertRaisesRegex(ValueError,'Vina, DOCK6'):
                prepare_candidates('missing.csv','unused',engines='invalid')
            read.assert_not_called()

    def test_each_engine_choice_exports_from_one_candidate_preparation(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);dataset=root/'compounds.csv'
            write_csv(dataset,[dict(molecule_chembl_id='CID1',canonical_smiles='CCO')])
            def convert(source,destination,**options):
                Path(destination).write_text('converted candidate')
            def prepare(command,**options):
                work=Path(options['cwd'])
                commands=(work/'prepare.com').read_text()
                self.assertIn('method am1',commands)
                self.assertNotIn('minimize spec #0',commands)
                (work/'prepared.mol2').write_text('one prepared conformation')
                return SimpleNamespace(returncode=0)
            settings={'ligand':{'charge_type':'am1','minimize':False}}
            for engine,formats in (('vina',{'pdbqt'}),('dock6',{'mol2'}),('both',{'pdbqt','mol2'})):
                with self.subTest(engine=engine),patch('biomolexplorer.docking_preparation.convert_structure',side_effect=convert),patch('biomolexplorer.docking_preparation.subprocess.run',side_effect=prepare) as runner:
                    output=prepare_candidates(dataset,root/engine,settings,engines=engine)
                    row=read_compounds(output)[0]
                    runner.assert_called_once()
                    self.assertEqual(row['docking_engines'],engine)
                    self.assertEqual({f for f in ('pdbqt','mol2') if row.get('prepared_'+f)},formats)
                    for format in formats:self.assertTrue(Path(row['prepared_'+format]).is_file())


class PreparationPipelineTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    upload=execution.PipelineExecutionTests.upload
    prepared=handoffs.DockingHandoffTests.prepared

    def test_multiple_candidate_sources_merge_by_identifier(self):
        ids=[self.upload(name,content) for name,content in (
            ('chembl.csv','name,smiles\nCHEMBL1,CCO\n'),
            ('pubchem.csv','name,smiles\nCID1,CCC\n'),
            ('zinc.csv','name,smiles\nZINC1,CCN\n'))]
        service=PipelineService(self.store)
        try:
            path,stem,_=service._materialize_inputs(self.project_id,new_stage('prepare_structures'),'base_selected_mols',
                [{'asset':identifier} for identifier in ids],{})
            self.assertEqual({r['molecule_chembl_id'] for r in read_compounds(path/(stem+'.csv'))},{'CHEMBL1','CID1','ZINC1'})
        finally:service.close()

    def test_redocking_receiver_infers_prepared_mode_and_selected_metadata(self):
        origin,stage=new_stage('redocking'),new_stage('prepare_structures')
        files=self.prepared(origin)+self.prepared(origin,'2ABC')
        asset=self.upload()
        stage['bindings']={'base_input_path':{'stage':origin['id'],'selector':'2ABC_A.dockprep.pdbqt'},
                           'base_selected_mols':{'asset':asset}}
        service=PipelineService(self.store)
        try:
            params=service._resolve(self.project_id,self.store.user(self.token)['id'],stage,{origin['id']:files})
            self.assertTrue(params['receptor_prepared'])
            receptor=Path(params['base_input_path'])/'MeuAlvo/Prepared'
            self.assertTrue((receptor/'2ABC_A.dockprep.mol2').exists())
            self.assertFalse((receptor/'1ABC_A.dockprep.pdbqt').exists())
            self.assertFalse((receptor/'2ABC_LIG_1A.lig.pdbqt').exists())
        finally:service.close()

    def test_redocking_popup_shows_only_prepared_receptors(self):
        origin,stage=new_stage('redocking'),new_stage('prepare_structures')
        files=self.prepared(origin)
        raw=Path(files[0]).parent/'1ABC.pdb';raw.write_text('raw PDB')
        stage['bindings']={'base_input_path':{'stage':origin['id'],'selector':'auto'}}
        popup=FileSelection({'stages':[dict(origin,status='succeeded',artifacts=files+[str(raw)])]},
            dict(stage,name=stage['name'],configuration=stage))
        names=[ref['selector'] for check,ref in popup.rows['base_input_path']]
        self.assertEqual(names,['1ABC_A.dockprep.pdbqt'])

    def test_every_preparation_batch_reaches_both_engines(self):
        root=self.store.project_dir(self.project_id);origin=new_stage('prepare_structures');files=[]
        for code in ('A','B'):
            folder=root/code;folder.mkdir()
            pdbqt=folder/(code+'.lig.pdbqt');mol2=folder/(code+'.lig.mol2')
            pdbqt.write_text('prepared vina');mol2.write_text('prepared dock6')
            write_csv(folder/'compounds.csv',[dict(molecule_chembl_id=code,canonical_smiles='CCO',
                conformer_file=str(mol2),prepared_pdbqt=str(pdbqt),prepared_mol2=str(mol2),docking_engines='both')])
            files.extend(map(str,folder.iterdir()))
        service=PipelineService(self.store)
        try:
            for engine in ('vina','dock6'):
                path,stem,_=service._materialize_inputs(self.project_id,new_stage('docking_'+engine),'base_selected_mols',
                    [{'stage':origin['id'],'selector':'auto'}],{origin['id']:files})
                rows=read_compounds(path/(stem+'.csv'))
                self.assertEqual({r['molecule_chembl_id'] for r in rows},{'A','B'})
                self.assertTrue(all(Path(r['prepared_pdbqt']).is_file() and Path(r['prepared_mol2']).is_file() for r in rows))
        finally:service.close()
