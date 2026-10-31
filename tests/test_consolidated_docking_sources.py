"""Primary docking inputs expose final score tables, retaining saved poses."""
import tempfile
from pathlib import Path
import unittest
from biomolexplorer.catalog import new_stage
from biomolexplorer.artifact_choices import consolidated_docking_tables,choices
from biomolexplorer.ui.file_selection import FileSelection

class ConsolidatedSourcesTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
    def file(self,name,content):
        path=self.root/name;path.parent.mkdir(exist_ok=True,parents=True);path.write_text(content);return path
    def test_native_docking_hides_poses_context_and_logs(self):
        final=self.file('Target/docking_results.csv','molecule_chembl_id,canonical_smiles,engine,score,conformer_file\nM1,CCO,vina,-7,pose.pdbqt\n')
        pose=self.file('Target/Vina/pose.pdbqt','ATOM\n');context=self.file('Target/context/context.csv',final.read_text());log=self.file('run.out','log')
        producer=new_stage('docking_vina');item=dict(producer,status='succeeded',artifacts=list(map(str,[final,pose,context,log])))
        self.assertEqual(consolidated_docking_tables(item['artifacts']),[final])
        self.assertEqual(choices(item),{'docking_results.csv'})
        stage=new_stage('dynamics_dock6');stage['bindings']['base_selected_mols']={'stage':producer['id']}
        form=FileSelection({'stages':[item]},{'configuration':stage,'name':stage['name']})
        self.assertEqual([Path(flag.physical_file) for flag,_ in form.rows['base_selected_mols']],[final])
        stage['bindings']['base_selected_mols']={'stage':producer['id'],'selector':'pose.pdbqt'}
        form=FileSelection({'stages':[item]},{'configuration':stage,'name':stage['name']})
        flags=[flag for flag,_ in form.rows['base_selected_mols'] if flag.value]
        self.assertEqual(len(flags),1);self.assertTrue(flags[0].available);self.assertTrue(flags[0].legacy_selection)
        form.pickers['base_selected_mols'].source.value='stage:'+producer['id']
        self.assertNotIn(flags[0],[flag for _,flag,_ in form.pickers['base_selected_mols'].available_rows()])
    def test_consensus_table_is_a_compatible_md_source(self):
        final=self.file('MeuAlvo_consensus.csv','molecule_chembl_id,canonical_smiles,vina,dock6,conformer_file\nM1,CCO,-7,-20,pose.mol2\n')
        producer=new_stage('consensus');item=dict(producer,status='succeeded',artifacts=[str(final)])
        stage=new_stage('dynamics_dock6');stage['bindings']['base_selected_mols']={'stage':producer['id']}
        form=FileSelection({'stages':[item]},{'configuration':stage,'name':stage['name']})
        self.assertEqual([Path(flag.physical_file) for flag,_ in form.rows['base_selected_mols']],[final])

    def test_consensus_expands_both_engine_poses_with_distinct_identities(self):
        from biomolexplorer.md_inputs import MDInputs
        from test_md_pair_selection import Store
        vina=self.file('poses/vina.pdbqt','ATOM\n');dock6=self.file('poses/dock6.mol2','@<TRIPOS>MOLECULE\n')
        table=self.file('MeuAlvo_consensus.csv','molecule_chembl_id,canonical_smiles,receptor_id,vina,dock6,vina_pose,dock6_pose\nM1,CCO,1ABC_A,-7,-20,poses/vina.pdbqt,poses/dock6.mol2\n')
        stage=new_stage('dynamics_dock6');stage['parameters']['base_selected_mols']=str(table)
        choices=MDInputs(Store(self.root),'project',{}).files(stage,'base_selected_mols')
        self.assertEqual(len(choices),2)
        self.assertEqual({Path(c['compound']['conformer_file']) for c in choices},{vina,dock6})
        self.assertEqual(len({c['identity'] for c in choices}),2)
