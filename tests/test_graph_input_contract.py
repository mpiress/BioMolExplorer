"""Only similarity stages or validated external edge tables feed graph blocks."""
import copy
import json
import tempfile
import unittest
from pathlib import Path

from biomolexplorer.catalog import new_stage
from biomolexplorer.flow import compatible, connect, input_ports, stage_issues
from biomolexplorer.graph_contract import normalize_legacy_graphs
from biomolexplorer.graph_inputs import GraphInputs
from biomolexplorer.operations import execute_operation, validate_operation
from biomolexplorer.pipeline import validate_pipeline
from biomolexplorer.visualizations import SUFFIX, load_view
from biomolexplorer.workspace import WorkspaceStore


class GraphInputContractTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.store=WorkspaceStore(self.temp.name)
        self.token=self.store.register('Owner','owner@graph-contract.org','graph-password')
        self.pid=self.store.create_project(self.token,'Graph contract')['id']
        self.root=self.store.project_dir(self.pid)

    def tearDown(self): self.temp.cleanup()

    def upload(self,name,content,kind):
        ticket=self.store.prepare_upload(self.token,self.pid,name,kind)
        (self.store.staging/ticket).write_text(content)
        return self.store.finish_upload(self.token,ticket)

    def test_multiple_similarity_stages_share_one_port_and_other_operations_are_rejected(self):
        graph=new_stage('graphs');left=new_stage('similarity');right=new_stage('similarity')
        stages=[left,right,graph]
        connect(stages,left['id'],graph['id'],'similarity_path')
        connect(stages,right['id'],graph['id'],'similarity_path')
        self.assertEqual(len(graph['bindings']['similarity_path']['sources']),2)
        self.assertEqual([p['field'] for p in input_ports(graph)],['similarity_path'])
        self.assertEqual(stage_issues(graph),[])
        for operation,field in [('fingerprints','fingerprints_path'),('fingerprints','similarity_path'),
                                ('retrieve_compounds','base_input_path'),('import_results','similarity_path')]:
            with self.subTest(operation=operation):
                source=new_stage(operation)
                if operation=='import_results':source['parameters']['kind']='similarity'
                before=copy.deepcopy(graph)
                self.assertFalse(compatible(source,graph,field))
                with self.assertRaises(ValueError):connect(stages+[source],source['id'],graph['id'],field)
                self.assertEqual(graph,before)
                bad=copy.deepcopy(graph);bad['bindings'][field]={'stage':source['id']}
                with self.assertRaises(ValueError):validate_pipeline([left,right,source,bad])

    def test_external_edges_without_smiles_produce_a_custom_graph(self):
        asset=self.upload('custom.csv','source,target,value\n001,002,0.85\n002,003,0.7\n','similarity')
        graph=new_stage('graphs');graph['bindings']['similarity_path']={'asset':asset}
        params=GraphInputs(self.store,self.pid,[graph],{}).resolve(graph)
        result=execute_operation('graphs',params,self.root/'out')
        model=load_view(Path(next(p for p in result.artifacts if p.endswith(SUFFIX))).read_bytes())
        self.assertEqual({n['id'] for n in model['nodes']},{'001','002','003'})
        self.assertEqual(len(model['edges']),2)
        self.assertEqual(model['fragment']['status'],'missing_structures')
        self.assertTrue(all(n['properties']['canonical_smiles'] is None for n in model['nodes']))

    def test_optional_external_structures_preserve_isolates_and_allow_mcs(self):
        edge=self.upload('custom.csv','source,target,value\nA,B,0.85\n','similarity')
        compounds=self.upload('compounds.csv','molecule_chembl_id,canonical_smiles\nA,CCO\nB,CCCO\nC,c1ccccc1\n','compounds')
        graph=new_stage('graphs');graph['bindings']={'similarity_path':{'asset':edge},'base_input_path':{'asset':compounds}}
        validate_pipeline([graph])
        params=GraphInputs(self.store,self.pid,[graph],{}).resolve(graph)
        result=execute_operation('graphs',params,self.root/'out')
        model=load_view(Path(next(p for p in result.artifacts if p.endswith(SUFFIX))).read_bytes())
        self.assertEqual({n['id'] for n in model['nodes']},{'A','B','C'})
        self.assertEqual(model['mcc'],['A','B'])
        self.assertEqual(model['fragment']['status'],'complete')

    def test_invalid_external_files_report_expected_format_before_execution(self):
        for content in ('molecule_chembl_id,fingerprint\nA,"[1,0]"\n',
                        'source,target\nA,B\n','source,target,value,value\nA,B,0.8,0.9\n'):
            file=self.root/'external.csv';file.write_text(content)
            graph=new_stage('graphs');graph['parameters']['similarity_path']='external.csv'
            with self.subTest(content=content),self.assertRaisesRegex(ValueError,'[Pp]adrão esperado'):
                GraphInputs(self.store,self.pid,[graph],{}).resolve(graph)

    def test_bad_similarity_records_are_removed_during_execution_with_a_report(self):
        file=self.root/'external.csv'
        file.write_text('source,target,value\nA,B,0.8\nA,B,NaN\nA,B,1.7\nA,B,bad\n')
        graph=new_stage('graphs');graph['parameters']['similarity_path']='external.csv'
        params=GraphInputs(self.store,self.pid,[graph],{}).resolve(graph)
        result=execute_operation('graphs',params,self.root/'clean')
        model=load_view(Path(next(p for p in result.artifacts if p.endswith(SUFFIX))).read_bytes())
        self.assertEqual(len(model['edges']),1)
        self.assertEqual(result.details['excluded_records'],3)

    def test_worker_contract_disallows_all_fingerprint_entry_points(self):
        for params in ({'fingerprints_path':'anything'},{'metric':'Dice'},{'threshold':70},{'fingerprint':'morgan'},
                       {'graph_inputs':[{'kind':'fingerprints','file':'x.csv','compound_files':[]}]}):
            with self.subTest(params=params),self.assertRaises(ValueError):validate_operation('graphs',params)

    def test_legacy_upgrade_does_not_hide_invalid_references(self):
        for binding in ({'sources':None},{'stage':{}},{'stage':'f'*32}):
            graph=new_stage('graphs');graph['bindings']['similarity_path']=binding
            with self.subTest(binding=binding),self.assertRaises(ValueError):
                validate_pipeline(normalize_legacy_graphs([graph]))

    def test_legacy_project_upgrade_preserves_similarity_links_and_experiment_files(self):
        fp=new_stage('fingerprints');sim=new_stage('similarity');graph=new_stage('graphs')
        graph['parameters'].update(metric='Dice',fingerprint='morgan',threshold=70)
        graph['bindings']={'fingerprints_path':{'stage':fp['id']},'similarity_path':{'stage':sim['id']}}
        before=copy.deepcopy([fp,sim,graph])
        normalized=normalize_legacy_graphs(before)
        self.assertEqual(before[2]['parameters']['metric'],'Dice')
        validate_pipeline(normalized)
        self.assertEqual(normalized[2]['bindings'],{'similarity_path':{'stage':sim['id']}})
        self.assertNotIn('metric',normalized[2]['parameters'])
        self.assertEqual(normalize_legacy_graphs(normalized),normalized)
        file=self.root/'old-result.csv';file.write_text('experiment retained')
        with self.store.connect() as db:
            db.execute('UPDATE projects SET pipeline=? WHERE id=?',(json.dumps(before),self.pid))
        project=self.store.project(self.token,self.pid)
        self.assertEqual(project['pipeline'],normalized)
        self.store.save_pipeline(self.token,self.pid,project['pipeline'],project['revision'])
        self.assertEqual(file.read_text(),'experiment retained')


if __name__=='__main__': unittest.main()
