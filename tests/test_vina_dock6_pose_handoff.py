"""Keep target-specific Vina minima through materialization and DOCK6 preparation."""
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import MagicMock,patch

from biomolexplorer.catalog import new_stage
from biomolexplorer.docking_data import best_vina_model,read_compounds,write_csv
from biomolexplorer.pipeline import PipelineService
from test_preparation_failures import MOL2,PDBQT


class PoseHandoffTests(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup)
        self.root=Path(temporary.name)

    def pose(self,name,values):
        path=self.root/name
        path.write_text(''.join(f'MODEL {i}\nREMARK VINA RESULT: {score} 0 0\n'+PDBQT+'ENDMDL\n'
            for i,score in enumerate(values,1)))
        return path

    def test_lowest_model_is_selected_even_when_file_is_not_sorted(self):
        self.assertEqual(best_vina_model(self.pose('models.pdbqt',[-6,-9,-7])),(2,-9.))

    def test_semantic_table_keeps_authorized_pose_companions_and_each_target(self):
        a=self.pose('a.pdbqt',[-7]);b=self.pose('b.pdbqt',[-9])
        table=self.root/'docking_results.csv'
        write_csv(table,[dict(molecule_chembl_id='M',canonical_smiles='CCO',receptor_id=rid,
            engine='vina',score=value,conformer_file=pose.name) for rid,value,pose in [('1ABC_A',-7,a),('2ABC_A',-9,b)]])
        service=object.__new__(PipelineService)
        service.store=SimpleNamespace(scoped_path=lambda pid,path:Path(path),project_dir=lambda pid:self.root)
        ref=dict(stage='vina',selector='auto',_resolved_paths=[str(table)])
        selected,available=service._docking_sources('project',[ref],{'vina':list(map(str,[table,a,b]))})
        destination=PipelineService._write_docking_compounds([(ref,selected,available)],self.root/'materialized',by_receptor=True)
        rows=read_compounds(destination)
        self.assertEqual({row['receptor_id'] for row in rows},{'1ABC_A','2ABC_A'})
        self.assertEqual(len(rows),2)
        self.assertEqual([Path(row['conformer_file']).read_text() for row in rows],[a.read_text(),b.read_text()])

    def test_dock6_prepares_the_best_model_for_the_matching_target_only(self):
        from wrappers.docking import perform_docking_dock6
        a=self.pose('a.pdbqt',[-7,-8]);b=self.pose('b.pdbqt',[-9])
        table=self.root/'compounds.csv'
        write_csv(table,[dict(molecule_chembl_id='M',canonical_smiles='CCO',receptor_id=rid,
            engine='vina',score=value,conformer_file=str(pose)) for rid,value,pose in [('1ABC_A',-7,a),('2ABC_A',-9,b)]])
        prepared=[]
        def prepare(dataset,destination,*args,**kwargs):
            rows=read_compounds(dataset);prepared.extend(rows)
            destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
            path=destination/'M.lig.mol2';path.write_text(MOL2)
            return [path]
        def factory(**options):
            engine=MagicMock();engine.excluded_compounds=set()
            def dock():
                folder=Path(options['base_output_path'])/'flex';folder.mkdir(parents=True,exist_ok=True)
                (folder/'M_scored.mol2').write_text(MOL2+'\nGrid_Score: -5\n')
            engine.perform_dock6_evaluation.side_effect=dock
            return engine
        with patch('wrappers.docking.prepare_docking_receptor',return_value=str(self.root)), \
             patch('biomolexplorer.docking_data.prepare_ligands',side_effect=prepare), \
             patch('biomolexplorer.docking_data.binding_center',return_value=[0.,0.,0.]), \
             patch('wrappers.docking.Dock6',side_effect=factory):
            perform_docking_dock6(str(self.root),'Target',str(self.root/'out'),str(self.root),'/dock6','gas','compounds',
                pdb_code=[['1ABC','LIG',1,'A'],['2ABC','LIG',1,'A']])
        self.assertEqual([(row['receptor_id'],row['conformer_index'],row['score']) for row in prepared],
                         [('1ABC_A','2','-8.0'),('2ABC_A','1','-9.0')])
        self.assertEqual([Path(row['conformer_file']) for row in prepared],[a,b])
        outputs=read_compounds(self.root/'out/docking_results.csv')
        self.assertEqual(len(outputs),2)
        self.assertEqual([(row['vina_input_model'],row['vina_input_score']) for row in outputs],[('2','-8.0'),('1','-9.0')])


if __name__=='__main__':
    unittest.main()
