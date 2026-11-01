"""The user's Alzheimer topology retains six separate ADMET/PubChem graphs."""
import copy
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import pandas as pd
from biomolexplorer.analysis_selection import expand_selection
from biomolexplorer.bindings import sources
from biomolexplorer.operations import execute_operation
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.result_files import ResultFiles
from biomolexplorer.result_origins import ResultOrigins
from biomolexplorer.visualizations import SUFFIX
import test_pipeline_execution as execution

FIXTURE=Path(__file__).parent/'fixtures/alzheimer_analysis_pipeline.json'


class AnalysisOriginsIntegrationTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    finish_selected=execution.PipelineExecutionTests.finish_selected
    wait=execution.PipelineExecutionTests.wait

    def test_three_old_pubchem_selections_expand_to_six_real_graphs_with_exact_origins(self):
        pipeline=json.loads(FIXTURE.read_text())
        by_op={s['operation']:s for s in pipeline}
        fp=by_op['fingerprints'];sim=by_op['similarity'];graph=by_op['graphs']
        direct=sources(fp['bindings']['base_input_path'])[3:]
        fp['bindings']['base_input_path']['sources']=sources(fp['bindings']['base_input_path'])[:3]
        self.save(pipeline);calls=[]
        class Manager:
            def __init__(self,config):self.config=config
            def submit(self,operation,params):
                calls.append(operation);output=self.config.workspace/'.biomolexplorer/jobs/worker/artifacts';output.mkdir(parents=True,exist_ok=True)
                if operation=='retrieve_compounds':
                    query=params['search_term'];target=output/'ChEMBL/DrugBank';target.mkdir(parents=True)
                    pd.DataFrame({'molecule_chembl_id':[query+'_A',query+'_B'],
                        'canonical_smiles':['CCO','CCO']}).to_csv(target/(query+'_FULL.csv'),index=False)
                elif operation=='admet':
                    source=Path(params['base_input_path'])/params['input_file'];frame=pd.read_csv(source)
                    frame=frame.assign(TPSA=20.,WLOGP=1.,BBB='BBB+',HIA='HIA+')
                    frame.to_csv(output/source.name,index=False)
                    frame.to_csv(output/(source.stem+'_BBB+.csv'),index=False)
                elif operation=='retrieve_pubchem':
                    frame=pd.read_csv(Path(params['base_input_path'])/params['input_file'])
                    frame['molecule_chembl_id']='PUBCHEM_'+frame['molecule_chembl_id']
                    target=output/'compounds/PubChem';target.mkdir(parents=True)
                    frame.to_csv(target/'compounds.csv',index=False)
                else:
                    result=execute_operation(operation,params,output).to_dict()
                if operation in ('retrieve_compounds','admet','retrieve_pubchem'):
                    result=dict(artifacts=[str(p) for p in output.rglob('*') if p.is_file()],details={})
                self.job=dict(id='worker',status='succeeded',log_path=str(output/'execution.log'),progress=None,error=None,result=result)
                return self.job
            def get(self,job):return self.job
            def cancel(self,job):self.job['status']='cancelled'
            def close(self):pass
        # This fixture represents outputs saved before canonical publication.
        # Keep those historical engine filenames while exercising selection
        # expansion and exact lineage against the real analytical workers.
        with patch('biomolexplorer.pipeline.JobManager',Manager),patch(
                'biomolexplorer.artifact_naming.publish',side_effect=lambda stage,files,root,**context:files):
            service=PipelineService(self.store);self.addCleanup(service.close)
            first=self.finish_selected(service,service.submit(self.token,self.project_id))
            self.assertEqual(first['status'],'succeeded',first.get('error'))
            fp['bindings']['base_input_path']['sources']+=direct
            self.save(pipeline)
            second=self.finish_selected(service,service.submit(self.token,self.project_id))
        self.assertEqual(second['status'],'succeeded',second.get('error'))
        stages={s['operation']:s for s in second['stages']}
        self.assertEqual(len(stages['fingerprints']['batches']),6)
        self.assertEqual(len(stages['similarity']['batches']),6)
        self.assertEqual(len(stages['graphs']['batches']),6)
        self.assertEqual(stages['similarity']['input_expansions'][0]['previous'],3)
        self.assertEqual(stages['graphs']['input_expansions'][0]['current'],6)
        saved={s['operation']:s for s in self.store.project(self.token,self.project_id)['pipeline']}
        self.assertEqual(len(sources(saved['similarity']['bindings']['base_input_path'])),6)
        self.assertEqual(len(sources(saved['graphs']['bindings']['similarity_path'])),6)
        results=ResultFiles(self.store)
        datasets=results.graph_datasets(self.token,self.project_id,second['id'],graph['id'])
        self.assertEqual(len(datasets),6)
        self.assertEqual(len({d['label'] for d in datasets}),6)
        branch_counts={}
        for item in datasets:
            model=results.graph_model(self.token,self.project_id,second['id'],graph['id'],item['path'])
            origin=model['dataset_origin'];expanded='Retrieve PubChem' in origin['origin']
            branch_counts[(origin['source_label'],expanded)]=True
            self.assertTrue(all(str(n['id']).startswith('PUBCHEM_')==expanded for n in model['nodes']))
            target=origin['source_label'].split(' · ')[-1]
            self.assertTrue(all(target in str(n['id']) for n in model['nodes']))
            self.assertEqual(len(model['nodes']),2);self.assertEqual(len(model['edges']),1)
            self.assertIn('CHEMBL',model['title'])
            # Companion files with identical names are scoped to this exact batch.
            same=[p for p in results.files(self.token,self.project_id,second['id'],graph['id'])
                  if Path(p['path']).is_relative_to(Path(item['path']).parent.parent)]
            self.assertTrue(same)
            self.assertTrue(all(p['origin']==origin['origin'] for p in same))
        self.assertEqual(len(branch_counts),6)
        for name,target in [('ACHE','CHEMBL220'),('BCHE','CHEMBL1914'),('BACE1','CHEMBL4822')]:
            self.assertIn((name+' · '+target,False),branch_counts)
            self.assertIn((name+' · '+target,True),branch_counts)
        fpfiles=results.files(self.token,self.project_id,second['id'],fp['id'])
        self.assertEqual(len({f['origin'] for f in fpfiles}),6)
        from biomolexplorer.ui.file_selection import FileSelection
        chooser=FileSelection(second,dict(configuration=saved['similarity'],name=sim['name']))
        picker=chooser.pickers['base_input_path'];picker.search.value='CHEMBL220'
        self.assertEqual(len(picker.available_rows()),2)
        self.assertTrue(all('CHEMBL220' in picker.label(flag) for _,flag,_ in picker.available_rows()))
        self.assertEqual(calls.count('retrieve_compounds'),3)
        self.assertEqual(calls.count('admet'),3)
        self.assertEqual(calls.count('retrieve_pubchem'),3)


class OriginRulesTests(unittest.TestCase):
    def test_partial_selection_is_not_expanded_without_complete_prior_coverage(self):
        from biomolexplorer.catalog import new_stage
        parent=new_stage('fingerprints');child=new_stage('similarity')
        files=['/old/1/morgan_compounds.csv','/old/2/morgan_compounds.csv']
        parent.update(status='succeeded',artifacts=files)
        child['bindings']['base_input_path']={'stage':parent['id'],'selector':'1/morgan_compounds.csv'}
        before=copy.deepcopy(child)
        old=dict(child,status='succeeded',configuration=copy.deepcopy(child))
        result,changes=expand_selection(child,[parent],[[parent,old]])
        self.assertEqual(result,before);self.assertEqual(changes,[])
        child['follow_inputs']=False
        self.assertEqual(expand_selection(child,[parent],[[parent,old]])[1],[])

    def test_cycles_terminate_and_do_not_invent_target_names(self):
        from biomolexplorer.catalog import new_stage
        a=new_stage('fingerprints');b=new_stage('similarity')
        for stage,parent,path in [(a,b,'/a.csv'),(b,a,'/b.csv')]:
            stage.update(artifacts=[path],artifact_records=[dict(path=path,
                inputs={'base_input_path':[dict(path='/b.csv' if stage is a else '/a.csv',stage_id=parent['id'])]})])
        description=ResultOrigins([a,b],'run').describe('/a.csv',a['id'],'run')
        self.assertNotIn('CHEMBL220',description['origin'])
        self.assertLess(len(description['source_routes']),4)

    def test_cached_lineage_loads_the_original_run_instead_of_a_new_target(self):
        from biomolexplorer.catalog import new_stage
        root=new_stage('retrieve_compounds');root['name']='ACHE';root['parameters']['search_term']='CHEMBL220'
        root.update(artifacts=['/old/compounds.csv'])
        old=dict(id='old',stages=[copy.deepcopy(root)])
        root['name']='BCHE';root['parameters']['search_term']='CHEMBL1914';root['artifacts']=['/new/compounds.csv']
        fp=new_stage('fingerprints');fp.update(artifacts=['/old/morgan_compounds.csv'],
            artifact_records=[dict(path='/old/morgan_compounds.csv',run_id='old',
                inputs={'base_input_path':[dict(path='/old/compounds.csv',name='compounds.csv',stage_id=root['id'],run_id='old')]})])
        resolver=ResultOrigins([root,fp],'new',load_run=lambda run:old if run=='old' else None)
        description=resolver.describe('/old/morgan_compounds.csv',fp['id'],'new')
        self.assertEqual(description['source_label'],'ACHE · CHEMBL220')
        self.assertNotIn('CHEMBL1914',description['origin'])


if __name__=='__main__':unittest.main()
