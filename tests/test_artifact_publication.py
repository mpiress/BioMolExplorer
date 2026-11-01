"""Publication, persistent IDs, reference handoffs and independent block reuse."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from biomolexplorer.artifact_naming import assign_block_ids,publish,canonical_filename,parse_filename
from biomolexplorer.artifact_choices import choices
from biomolexplorer.artifacts import index_run,records,reference,validate_reference
from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.fingerprint_selection import filename_kind
import test_pipeline_execution as execution


class PublicationTests(unittest.TestCase):
    def test_mol_sim_and_downstream_admet_keep_meaningful_roles(self):
        from biomolexplorer.artifact_naming import human_name
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=new_stage('retrieve_compounds');admet=new_stage('admet')
            assign_block_ids([source,admet]);files=[]
            for role in ('MOLS','SIMS'):
                path=root/('CHEMBL220_'+role+'.csv')
                path.write_text('molecule_chembl_id,canonical_smiles\nA,CCO\n');files.append(str(path))
            outputs=publish(source,files,root)
            self.assertEqual({Path(p).name for p in outputs},{'B01__chembl__mol.csv','B01__chembl__sim.csv'})
            self.assertEqual(len(choices(dict(source,artifacts=outputs))),2)
            generated=[]
            for role in ('BBB+','BBB-','HIA+','HIA-'):
                path=root/('B01__chembl__mol_'+role+'.csv')
                path.write_text('molecule_chembl_id,canonical_smiles\nA,CCO\n');generated.append(str(path))
            result=publish(admet,generated,root)
            self.assertEqual({Path(p).name for p in result},{'B02__admet__'+role+'.csv' for role in ('bbb-pos','bbb-neg','hia-pos','hia-neg')})
            self.assertEqual({parse_filename(human_name(r))['result_type'] for r in records(dict(admet,artifacts=result))},{'BBB+','BBB-','HIA+','HIA-'})

    def test_short_names_counts_and_graph_mapping_preserve_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);stage=new_stage('fingerprints');assign_block_ids([stage])
            original=root/'morgan_CHEMBL220_ADMET_BBB+_VERY_LONG_INPUT.csv'
            original.write_text('molecule_chembl_id,canonical_smiles,fingerprint\nA,CCO,"[1,0]"\n')
            before=original.read_bytes();result=publish(stage,[str(original)],root)
            self.assertEqual(Path(result[0]).name,'B01__fingerprint__morgan.csv')
            self.assertEqual(Path(result[0]).read_bytes(),before);self.assertFalse(original.exists())
            self.assertFalse(any(p.is_symlink() for p in root.iterdir()))
            item=dict(stage,configuration=stage,artifacts=result)
            record=records(item,'run')[0]
            self.assertEqual(record['record_count'],1);self.assertEqual(record['block_id'],'B01')
            self.assertEqual(filename_kind(result[0]),'morgan')
            self.assertEqual(validate_reference(reference(record),result),result[0])
            self.assertEqual(publish(stage,result,root),result)

    def test_lineage_is_artifact_references_not_filenames(self):
        a,b=new_stage('retrieve_compounds'),new_stage('admet');assign_block_ids([a,b])
        a=dict(a,configuration=copy.deepcopy(a),artifacts=['/one/compounds.csv'])
        b['bindings']={'base_input_path':{'stage':a['id']}}
        b=dict(b,configuration=copy.deepcopy(b),artifacts=['/two/filtered.csv'],
            input_files={'base_input_path':['/one/compounds.csv']})
        index_run([a,b],'run')
        child=b['artifact_records'][0]
        self.assertEqual(child['parent_artifact_ids'],[a['artifact_records'][0]['artifact_id']])
        self.assertEqual(child['lineage'],child['parent_artifact_ids'])

    def test_collision_and_reference_only_inputs_are_handled(self):
        with tempfile.TemporaryDirectory() as directory,tempfile.TemporaryDirectory() as elsewhere:
            root=Path(directory);stage=new_stage('admet');assign_block_ids([stage])
            files=[]
            for name in ('A','B'):
                path=root/(name+'_ADMET_BBB+.csv');path.write_text('id\nA\n');files.append(str(path))
            external=Path(elsewhere)/'legacy.csv';external.write_text('id\nA\n');files.append(str(external))
            result=publish(stage,files,root)
            self.assertEqual(len(set(result)),3);self.assertEqual(result[-1],str(external))
            self.assertTrue(all(len(Path(p).name)<80 for p in result))


class PersistentPublicationTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    finish_selected=execution.PipelineExecutionTests.finish_selected
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def test_ids_survive_reorder_delete_new_block_and_duplicate(self):
        a,b=new_stage('admet'),new_stage('fingerprints');self.save([a,b])
        self.assertEqual([s['id'] for s in self.store.project(self.token,self.project_id)['pipeline']],[a['id'],b['id']])
        self.save([b,a]);self.assertEqual([s['id'] for s in self.store.project(self.token,self.project_id)['pipeline']],[b['id'],a['id']])
        self.save([a]);c=new_stage('similarity');self.save([a,c])
        self.assertEqual(self.store.project(self.token,self.project_id)['pipeline'][0]['id'],a['id'])
        self.assertNotIn(c['id'],{a['id'],b['id']})

    def test_equivalent_blocks_reference_prior_files_and_tamper_invalidates(self):
        a,b=new_stage('retrieve_compounds'),new_stage('retrieve_compounds');self.save([a,b])
        manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store)
            try:
                run=self.finish_selected(service,service.submit(self.token,self.project_id))
                self.assertEqual(run['status'],'succeeded',run.get('error'))
                self.assertEqual(len(calls),1)
                first,second=run['stages']
                self.assertEqual(first['artifacts'],second['artifacts'])
                self.assertTrue(second['reused'])
                self.assertEqual(first['artifact_records'][0]['artifact_id'],second['artifact_records'][0]['artifact_id'])
                validate_reference(reference(second['artifact_records'][0]),second['artifacts'])
                self.assertEqual(parse_filename(first['artifact_records'][0]['public_name'])['operation'],'ChEMBL')
                Path(first['artifacts'][0]).write_text('tampered')
                again=self.finish_selected(service,service.submit(self.token,self.project_id))
                self.assertEqual(again['status'],'succeeded',again.get('error'))
                self.assertEqual(len(calls),2)
            finally:service.close()


if __name__=='__main__':unittest.main()
