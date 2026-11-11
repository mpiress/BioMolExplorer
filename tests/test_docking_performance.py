"""Bounded publication and scientifically scoped preparation reuse."""
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch,MagicMock

from biomolexplorer.artifacts import annotate_files,records
from biomolexplorer.catalog import new_stage
from biomolexplorer.docking_data import write_csv,write_results,read_compounds
from biomolexplorer.docking_preparation import prepare_candidates
import test_preparation_failures as fixtures
MOL2,PDBQT=fixtures.MOL2,fixtures.PDBQT


class DockingPerformanceTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)

    def test_table_annotation_never_decorates_thousands_of_pose_files(self):
        stage=new_stage('docking_vina')
        table=self.root/'docking_results.csv';table.write_text('molecule_chembl_id,score\nA,-7\n')
        stage['artifacts']=[str(self.root/(str(i)+'.lig.pdbqt')) for i in range(10000)]+[str(table)]
        from biomolexplorer.artifact_naming import name_records
        sizes=[]
        def bounded(entries):
            entries=list(entries);sizes.append(len(entries));return name_records(entries)
        with patch('biomolexplorer.artifact_naming.name_records',side_effect=bounded):
            result=annotate_files(stage,[dict(path=str(table),name=table.name)])
        self.assertEqual(sizes,[1]);self.assertEqual(len(result),1)

    def test_descriptor_uses_its_batch_inputs_instead_of_aggregate_manifest(self):
        stage=new_stage('docking_vina');paths=[str(self.root/(str(i)+'.pdbqt')) for i in range(2)]
        stage.update(artifacts=paths,input_manifest={'unrelated':'aggregate'},batches=[
            dict(artifacts=[p],input_manifest={'own':i},processing_runtime={'native':i}) for i,p in enumerate(paths)])
        result=records(stage,decorate=False)
        self.assertEqual([r['dependency_manifest'] for r in result],[{'own':0},{'own':1}])
        self.assertNotEqual(result[0]['runtime_id'],result[1]['runtime_id'])

    def test_receptor_context_is_copied_once_per_receptor_not_per_compound(self):
        prepared=self.root/'Prepared';prepared.mkdir()
        receptor=prepared/'4M0E_A.noH.pdb';receptor.write_text('ATOM\n')
        table=self.root/'compounds.csv'
        write_csv(table,[dict(molecule_chembl_id='M'+str(i),canonical_smiles='CCO') for i in range(3)])
        output=self.root/'out';output.mkdir();poses=[]
        for i in range(3):
            pose=output/('4M0E_LIG_1A_M'+str(i)+'.lig.pdbqt')
            pose.write_text('REMARK VINA RESULT: -7 0 0\n'+PDBQT);poses.append(pose)
        import shutil
        with patch('shutil.copy2',wraps=shutil.copy2) as copier:
            rows=write_results('vina',poses,table,[['4M0E','LIG',1,'A']],output,prepared)
        self.assertEqual(copier.call_count,1);self.assertEqual(len(rows),3)
        self.assertTrue(all(row['score']==-7 for row in rows))

    def test_vina_and_dock6_reuse_preparation_across_receptors_and_keep_original_geometry(self):
        table=self.root/'compounds.csv';write_csv(table,[dict(molecule_chembl_id='A',canonical_smiles='CCO')])
        cache=self.root/'shared-cache'
        def convert(source,destination,**kwargs):
            Path(destination).write_text(PDBQT if Path(destination).suffix=='.pdbqt' else MOL2)
        with patch.dict(os.environ,{'BIOMOL_PREPARATION_CACHE':str(cache),'BIOMOL_PREPARATION_WORKERS':'1'}), \
             patch('biomolexplorer.ligand_cache.preparation_identity',return_value={'tool':'unchanged'}), \
             patch('biomolexplorer.docking_preparation.convert_structure',side_effect=convert), \
             patch('biomolexplorer.docking_preparation.subprocess.run',side_effect=fixtures.PreparationFailureTests.chimera) as chimera:
            vina=prepare_candidates(table,self.root/'receptor-one',engines='vina')
            dock6=prepare_candidates(table,self.root/'receptor-two',engines='dock6')
            again=prepare_candidates(table,self.root/'receptor-three',engines='vina')
            self.assertEqual(chimera.call_count,1)
            self.assertEqual((vina.parent/'A.lig.mol2').read_bytes(),(dock6.parent/'A.lig.mol2').read_bytes())
            self.assertEqual((vina.parent/'A.lig.pdbqt').read_bytes(),(again.parent/'A.lig.pdbqt').read_bytes())
            prepare_candidates(table,self.root/'different-ph',ph=6.0,engines='dock6')
            self.assertEqual(chimera.call_count,2)

    def test_dock6_prepares_once_and_places_independent_copies_for_each_receptor(self):
        from wrappers.docking import perform_docking_dock6
        table=self.root/'compounds.csv';write_csv(table,[dict(molecule_chembl_id='A',canonical_smiles='CCO')])
        output=self.root/'out'
        def prepare(dataset,destination,*args,**kwargs):
            destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
            ligand=destination/'A.lig.mol2';ligand.write_text(MOL2)
            write_csv(destination/'compounds.csv',[dict(molecule_chembl_id='A',canonical_smiles='CCO',prepared_origin='smiles')])
            return [ligand]
        def engine(**kwargs):
            instance=MagicMock()
            instance.excluded_compounds=set()
            def dock():
                folder=Path(kwargs['base_output_path'])/'flex';folder.mkdir(parents=True,exist_ok=True)
                (folder/'A_scored.mol2').write_text(MOL2+'\n########## Grid_Score: -5\n')
            instance.perform_dock6_evaluation.side_effect=dock
            # An optional footprint failure must not erase a valid docking score.
            instance.prepare_footprint.side_effect=lambda **options:instance.excluded_compounds.add('A')
            return instance
        with patch('wrappers.docking.prepare_docking_receptor',return_value=str(self.root)), \
             patch('biomolexplorer.docking_data.prepare_ligands',side_effect=prepare) as preparation, \
             patch('biomolexplorer.docking_data.binding_center',side_effect=lambda path,record:[0.,0.,0.] if record[3]=='A' else [10.,0.,0.]), \
             patch('wrappers.docking.Dock6',side_effect=engine):
            perform_docking_dock6(str(self.root),'Target',str(output),str(self.root),'/dock6','gas','compounds',
                pdb_code=[['4M0E','LIG',1,'A'],['4M0E','LIG',1,'B']])
        self.assertEqual(preparation.call_count,1)
        self.assertEqual((output/'Molecules'/'A.lig.mol2').read_text(),MOL2)
        self.assertIn('10.0000',(output/'Target'/'4M0E_B'/'Molecules'/'A.lig.mol2').read_text())
        self.assertEqual(len(read_compounds(output/'docking_results.csv')),2)

    def test_completed_docking_batches_are_readable_during_organization_and_after_interruption(self):
        from biomolexplorer.workspace import WorkspaceStore
        from biomolexplorer.docking_results import DockingResults
        store=WorkspaceStore(self.root/'workspace')
        token=store.register('Owner','checkpoint@local.test','temporary-password')
        pid=store.create_project(token,'Docking')['id'];root=store.project_dir(pid)
        table=root/'docking_results.csv';pending=root/'pending.csv'
        for path in (table,pending):
            write_csv(path,[dict(molecule_chembl_id='A',canonical_smiles='CCO',engine='vina',score=-7)])
        stage=new_stage('docking_vina');stage.update(status='running',artifacts=[str(table),str(pending)],
            batches=[dict(status='succeeded',artifacts=[str(table)])])
        with store.connect() as db:
            db.execute('INSERT INTO runs VALUES(?,?,?,?,?,?,?,?)',('checkpoint',pid,store.user(token)['id'],
                'running',json.dumps([stage]),1.,1.,None))
        service=DockingResults(store)
        with patch.object(store,'get_result_run',wraps=store.get_result_run) as reader:
            self.assertEqual([row['path'] for row in service.tables(token,pid,'checkpoint',stage['id'])],[str(table)])
            page=service.page(token,pid,'checkpoint',stage['id'],str(table))
            self.assertEqual(page['total'],1);self.assertFalse(page['can_edit'])
            self.assertEqual(reader.call_count,1)
        stage['status']='interrupted'
        with store.connect() as db:
            db.execute('UPDATE runs SET status=?,stages=?,updated=? WHERE id=?',
                ('interrupted',json.dumps([stage]),2.,'checkpoint'))
        self.assertEqual(service.page(token,pid,'checkpoint',stage['id'],str(table))['total'],1)


if __name__=='__main__':unittest.main()
