"""Regression contracts for preparation, pose pairing and docking score handoffs.

External engines are replaced at their command boundary; file/score processing is real.
"""
import csv
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock,patch
import pandas as pd
import test_pipeline_execution as execution
from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService,select_input
from biomolexplorer.input_validation import validate_file
from biomolexplorer.ui.file_selection import FileSelection
from biomolexplorer.docking_inputs import matching_poses
from caad.docking import Docking,DockVina,Dock6
from wrappers.docking import generate_consensus,perform_docking_dock6
from wrappers.redocking import perform_redocking

ATOM='ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n'
MOL2='@<TRIPOS>MOLECULE\nreference\n1 0 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C 0.0 0.0 0.0 C.3 1 LIG 0.0\n'

class DockingHandoffTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    upload=execution.PipelineExecutionTests.upload

    def prepared(self,source,code='1ABC',chain='A'):
        root=self.store.project_dir(self.project_id)/f'{code}_{chain}'/'MeuAlvo'
        folder=root/'Prepared';folder.mkdir(parents=True)
        (root/'pdb_codes.csv').write_text(f'PDB_CODE,LIGAND,RESNUM,CHAIN\n{code},LIG,1,{chain}\n')
        (folder/'centers.csv').write_text(f'{code}_LIG_1{chain}\n1\n2\n3\n')
        for suffix in ('.dockprep.pdbqt','.dockprep.pdb','.noH.pdb'):(folder/f'{code}_{chain}{suffix}').write_text(ATOM)
        (folder/f'{code}_{chain}.dockprep.mol2').write_text(MOL2)
        for suffix in ('.lig.pdb','.lig.pdbqt'):(folder/f'{code}_LIG_1{chain}{suffix}').write_text(ATOM)
        return [str(p) for p in root.rglob('*') if p.is_file()]

    def test_dock6_surface_uses_native_dms_and_keeps_sphgen_handoff(self):
        root=self.store.project_dir(self.project_id)
        receptor=root/'receptor';receptor.mkdir()
        ligands=root/'ligands';ligands.mkdir()
        (receptor/'1ABC.noH.pdb').write_text(ATOM.replace('LIG','ALA'))
        (ligands/'M1.mol2').write_text(MOL2)
        obj=object.__new__(Dock6);obj.logger=MagicMock()
        obj.receptorpath=str(receptor);obj.ligandpath=str(ligands)
        obj._Dock6__base_output_path=str(root/'out');obj._Dock6__pdb_code='1ABC'
        obj._Dock6__density=.5;obj._Dock6__radius=1.4;obj._Dock6__distance=10.
        obj.generate_docking_script=MagicMock()
        def external(command, cwd=None):
            if isinstance(command,list) and command[0]=='sphere_selector':
                (Path(cwd)/'selected_spheres.sph').write_text('spheres')
        obj.perform_subprocess=MagicMock(side_effect=external)
        obj.prepare_surface()
        surface=root/'out/surface/1ABC.dms'
        self.assertIn('SC0',surface.read_text())
        commands=[call.args[0] for call in obj.perform_subprocess.call_args_list]
        self.assertEqual(commands[0],'sphgen -i INSPH -o OUTSPH')
        self.assertEqual(commands[1][0],'sphere_selector')
        self.assertEqual(len(commands),2)
        self.assertEqual((root/'out/surface/Molecules/M1.sph').read_text(),'spheres')

    def test_preparation_copies_native_centers_ligands_and_dock6_companions(self):
        source,stage=new_stage('prepare_structures'),new_stage('docking_dock6')
        files=self.prepared(source)+self.prepared(source,'2ABC')
        service=PipelineService(self.store)
        try:
            path,_,selected=service._materialize_inputs(self.project_id,stage,'base_input_path',
                [{'stage':source['id'],'selector':'1ABC_A.dockprep.pdbqt'}],{source['id']:files})
            folder=path/'MeuAlvo'/'Prepared'
            self.assertTrue((folder/'1ABC_A.dockprep.mol2').is_file())
            self.assertTrue((folder/'1ABC_A.noH.pdb').is_file())
            self.assertTrue((folder/'1ABC_LIG_1A.lig.pdb').is_file())
            self.assertFalse(any('2ABC' in p.name for p in folder.iterdir()))
            self.assertEqual(pd.read_csv(folder/'centers.csv').columns.tolist(),['1ABC_LIG_1A'])
            self.assertEqual(pd.read_csv(path/'MeuAlvo'/'pdb_codes.csv')['PDB_CODE'].tolist(),['1ABC'])
            validate_file(folder/'1ABC_A.dockprep.mol2','prepared_structures')
        finally:service.close()

    def test_center_lookup_accepts_native_and_legacy_columns(self):
        obj=object.__new__(Docking);obj.logger=MagicMock()
        for key in ('1ABC_LIG_1A','1ABC_LIG_1_A'):
            obj.centers=pd.DataFrame({key:[1.,2.,3.]})
            self.assertEqual(list(obj.retrieve_centerofmass_dataset('unused','1ABC','LIG','1','A')),[1.,2.,3.])

    def test_individual_consensus_pairs_inputs_instead_of_cartesian_product(self):
        stage=new_stage('consensus')
        vina=[self.upload(f'{name}.lig.pdbqt',ATOM+'REMARK VINA RESULT: -7.0 0 0\n','vina') for name in ('A','B')]
        dock6=[self.upload(f'{name}_scored.mol2',MOL2+'########## Grid_Score: -5.0\n','dock6') for name in ('A','B')]
        stage['bindings']={'base_vina_path':{'sources':[{'asset':a} for a in vina]},
                           'base_dock6_path':{'sources':[{'asset':a} for a in dock6]}}
        stage['bindings']['base_input_path']=stage['bindings']['base_vina_path'].copy()
        service=PipelineService(self.store)
        try:
            variants=list(service._variants(self.project_id,stage));self.assertEqual(len(variants),2)
            for variant,_ in variants:self.assertEqual(variant['bindings']['base_input_path'],variant['bindings']['base_vina_path'])
            form=FileSelection({'stages':[]},{'configuration':stage,'name':stage['name']})
            self.assertEqual(set(form.rows),{'base_vina_path','base_dock6_path'})
            stage['bindings']['base_vina_path']={'asset':vina[0]};stage['bindings']['base_input_path']=stage['bindings']['base_vina_path'].copy()
            stage['bindings']['base_dock6_path']={'asset':dock6[1]}
            with self.assertRaisesRegex(ValueError,'mesmos identificadores'):list(service._variants(self.project_id,stage))
        finally:service.close()

    def test_dock6_batches_pair_receptor_compounds_and_poses(self):
        source,stage=new_stage('prepare_structures'),new_stage('docking_dock6')
        structures=self.prepared(source)+self.prepared(source,'2ABC')
        refs=[{'stage':source['id'],'selector':f'{code}_A.dockprep.pdbqt'} for code in ('1ABC','2ABC')]
        compounds=[self.upload(f'{code}.csv',f'name,smiles\n{code},CCO\n') for code in ('M1','M2')]
        poses=[self.upload(f'{code}_LIG_1A_{mol}.lig.pdbqt',ATOM+'REMARK VINA RESULT: -7 0 0\n','vina')
               for code,mol in (('1ABC','M1'),('2ABC','M2'))]
        stage['bindings']={'base_input_path':{'sources':refs},'base_selected_mols':{'sources':[{'asset':a} for a in compounds]},
            'base_vina_path':{'sources':[{'asset':a} for a in poses]}}
        service=PipelineService(self.store)
        try:self.assertEqual(len(list(service._variants(self.project_id,stage,{source['id']:structures}))),2)
        finally:service.close()

    def test_ambiguous_filename_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'mais de um resultado'):
            select_input(['/some/a/file.csv','/some/b/file.csv'],'base_input_path','file.csv')

    def test_pose_identity_supports_native_and_previous_files(self):
        poses=['1ABC_LIG_1A_M1.lig.pdbqt','1ABC_A_M1.lig.pdbqt','1ABC_M1.lig.pdbqt','M1.pdbqt',
               '1ABC_LIG_1B_M1.lig.pdbqt','1ABC_LIG_1A_M2.lig.pdbqt']
        self.assertEqual(len(matching_poses(poses,[['1ABC','LIG','1','A']],{'M1'})),4)

    def test_dock6_uses_selected_compounds_and_a_string_receptor_identifier(self):
        root=self.store.project_dir(self.project_id);poses=root/'poses';poses.mkdir()
        (root/'compounds.csv').write_text('name,smiles\nM1,CCO\n')
        for name in ('1ABC_LIG_1A_M1','1ABC_LIG_1A_M2','1ABC_LIG_1B_M1'):
            (poses/(name+'.lig.pdbqt')).write_text(ATOM)
        with patch('wrappers.docking.Dock6') as engine:
            perform_docking_dock6(str(root),'MeuAlvo',str(root/'out'),str(root),
                '/engines/dock6','gas','compounds',['1ABC','LIG',1,'A'],base_vina_path=str(poses))
            self.assertEqual(engine.call_args.kwargs['pdb_code'],'1ABC_A')
            engine.return_value.recover_better_conforms_of_vina.assert_called_once_with(
                charge_type='gas',filename=['1ABC_LIG_1A_M1.lig.pdbqt'])

    def test_vina_reference_names_do_not_collide_or_skip_other_chains(self):
        root=self.store.project_dir(self.project_id);ligands=root/'ligands';output=root/'poses'
        ligands.mkdir();output.mkdir()
        for mol in ('M1','M2'):(ligands/f'{mol}.lig.pdbqt').write_text(ATOM)
        (output/'1ABC_LIG_1A_M1.lig.pdbqt').write_text(ATOM)
        obj=object.__new__(DockVina);obj.ligandpath=str(ligands);obj.outputpath=str(output);obj.receptorpath=str(root)
        obj.logger=MagicMock();obj._DockVina__pdb_codes=[['1ABC','LIG',1,'A'],['1ABC','LIG',2,'B']]
        obj._DockVina__centerofmasspath=str(root);obj._DockVina__sizeof_box=[24]*3
        obj._DockVina__exhaustiveness=20;obj._DockVina__num_modes=10
        obj.retrieve_centerofmass_dataset=MagicMock(return_value=[1,2,3])
        obj.generate_docking_script=MagicMock();obj.perform_subprocess=MagicMock(return_value=True)
        obj.docking(str(root))
        names=[c.kwargs['out'] for c in obj.generate_docking_script.call_args_list]
        self.assertEqual(set(names),{'1ABC_LIG_1A_M2.lig.pdbqt','1ABC_LIG_2B_M1.lig.pdbqt','1ABC_LIG_2B_M2.lig.pdbqt'})

    def test_prepared_redocking_normalizes_manual_four_column_records(self):
        root=self.store.project_dir(self.project_id);source=new_stage('prepare_structures');files=self.prepared(source)
        base=Path(files[0]).parent
        if base.name=='Prepared':base=base.parent
        def redock(**kwargs):
            frame=pd.read_csv(base/'pdb_codes.csv');frame['RMSD']=0.8;frame.to_csv(base/'pdb_codes.csv',index=False)
        with patch('wrappers.redocking.DockVina') as engine:
            engine.return_value.redocking.side_effect=redock
            perform_redocking(str(base.parent),'MeuAlvo',str(root/'redocking'),
                              pdb_codes=[['1ABC','LIG',1,'A']],prepare_complex=False)
            frame=engine.call_args.kwargs['pdb_codes']
            self.assertEqual(frame.columns.tolist(),['PDB_CODE','LIGAND','RESNUM','CHAIN','RESOLUTION'])
            self.assertEqual(frame.iloc[0]['PDB_CODE'],'1ABC')

    def test_consensus_singleton_and_equal_scores_are_finite_and_importable(self):
        root=self.store.project_dir(self.project_id);vina=root/'vina';dock6=root/'dock6';out=root/'consensus'
        vina.mkdir();dock6.mkdir();out.mkdir()
        for size in (1,2):
            name=f'M{size}'
            (vina/f'{name}.pdbqt').write_text(ATOM+'REMARK VINA RESULT: -7 0 0\n')
            (dock6/f'{name}_scored.mol2').write_text(MOL2+'########## Grid_Score: -5e0\n')
            generate_consensus(str(root),str(out),'example',base_vina_path=str(vina),base_dock6_path=str(dock6))
            path=out/'example.csv';validate_file(path,'scores')
            frame=pd.read_csv(path)
            self.assertEqual(len(frame),size)
            self.assertEqual(frame['vina'].tolist(),[-7.]*size);self.assertEqual(frame['dock6'].tolist(),[-5.]*size)
            self.assertTrue((frame[['z-score','min-max']]==0).all().all())

    def test_consensus_missing_pose_and_nonfinite_scores_fail_clearly(self):
        root=self.store.project_dir(self.project_id);vina=root/'vina';dock6=root/'dock6';out=root/'out'
        vina.mkdir();dock6.mkdir();out.mkdir()
        (dock6/'M1_scored.mol2').write_text(MOL2+'########## Grid_Score: -5.0\n')
        args=(str(root),str(out),'example')
        kwargs={'base_vina_path':str(vina),'base_dock6_path':str(dock6)}
        with self.assertRaisesRegex(ValueError,'pose Vina correspondente'):generate_consensus(*args,**kwargs)
        (vina/'M1.lig.pdbqt').write_text(ATOM+'REMARK VINA RESULT: nan 0 0\n')
        with self.assertRaisesRegex(ValueError,'não finito'):generate_consensus(*args,**kwargs)

if __name__=='__main__':unittest.main()
