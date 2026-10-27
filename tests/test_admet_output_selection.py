"""All consumers see generated ADMET datasets with their original target."""
import copy
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from biomolexplorer.catalog import new_stage
from biomolexplorer.artifact_choices import produced_files,choices
from biomolexplorer.artifacts import index_run
from biomolexplorer.ui.file_selection import FileSelection

class ADMETOutputSelectionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.sources=[];self.generated=[];self.raw=[]
        self.admet=new_stage('admet');self.admet['block_id']='B04'
        batches=[]
        for index,target in enumerate(('ACHE','BCHE','BACE1'),1):
            source=new_stage('retrieve_compounds');source['parameters']['target']=target
            source['parameters']['search_term']=('CHEMBL220','CHEMBL1914','CHEMBL4822')[index-1]
            source['block_id']=f'B{index:02d}';paths=[]
            for kind in ('mol','sim'):
                folder=self.root/source['id']/target;folder.mkdir(parents=True,exist_ok=True)
                raw=folder/f'B{index:02d}__chembl__{kind}.csv';raw.write_text('molecule_chembl_id,canonical_smiles\nA,CCO\n')
                paths.append(str(raw));self.raw.append(str(raw))
                output=self.root/self.admet['id']/'individual'/f'{target}_{kind}';output.mkdir(parents=True)
                artifacts=[]
                for role in ('compounds','bbb-pos','bbb-neg','hia-pos','hia-neg','evaluated','selected','rejected'):
                    p=output/f'B04__admet__{role}.csv'
                    p.write_text('molecule_chembl_id,canonical_smiles,BBB,HIA,ADMET_selected\nA,CCO,BBB+,HIA+,True\n')
                    artifacts.append(str(p));self.generated.append(str(p))
                config=copy.deepcopy(self.admet);config['bindings']={'base_input_path':{'stage':source['id']}}
                batches.append(dict(label=f'B{index:02d}__chembl__{kind}',artifacts=artifacts,
                    configuration=config,input_files={'base_input_path':[str(raw)]}))
            self.sources.append(dict(source,configuration=copy.deepcopy(source),status='succeeded',artifacts=paths))
        self.item=dict(self.admet,configuration=copy.deepcopy(self.admet),status='succeeded',
            artifacts=self.generated+self.raw,batches=batches)
        self.run={'id':'run','stages':self.sources+[self.item]}
        index_run(self.run['stages'],'run')

    def form(self,operation):
        stage=new_stage(operation);field='base_selected_mols' if operation in ('docking_vina','docking_dock6','prepare_structures') else 'base_input_path'
        stage['bindings']={field:{'stage':self.admet['id']}}
        return FileSelection(self.run,dict(configuration=stage,name=stage['name']),page=SimpleNamespace(update=lambda:None)),field

    def test_every_compound_consumer_lists_outputs_only_and_can_find_bace1(self):
        for operation in ('retrieve_pubchem','fingerprints','graphs','admet','prepare_structures','docking_vina','docking_dock6'):
            with self.subTest(operation=operation):
                form,field=self.form(operation);rows=form.rows[field]
                self.assertEqual(len(rows),48)
                self.assertEqual({flag.physical_file for flag,_ in rows},set(self.generated))
                picker=form.pickers[field];self.assertEqual(len(picker.available.options),48)
                self.assertIn('de 48',picker.available_count.value)
                picker.search.value='BACE1';picker.search_changed(None)
                self.assertEqual(len(picker.available.options),16)
                self.assertTrue(all('BACE1' in o.text for o in picker.available.options))
                picker.add_all_sources();self.assertEqual(len(form.read()['bindings'][field]['sources']),48)

    def test_choices_and_legacy_outputs_exclude_raw_input_copies(self):
        self.assertEqual(set(produced_files(self.item)),set(map(Path,self.generated)))
        self.assertEqual(len(choices(self.item)),48)
        folder=self.root/'legacy';folder.mkdir()
        processed=folder/'ACHE.csv';processed.write_text('molecule_chembl_id,canonical_smiles,BBB,HIA\nA,CCO,BBB+,HIA+\n')
        raw=folder/'BCHE.csv';raw.write_text('molecule_chembl_id,canonical_smiles\nA,CCO\n')
        item=dict(self.item,artifacts=[str(processed),str(raw)])
        self.assertEqual(produced_files(item),[processed])
        imported=copy.deepcopy(item);imported['configuration']['provided_results']={'kind':'compounds'}
        self.assertEqual(produced_files(imported),[processed,raw])

class ADMETResolverTests(unittest.TestCase):
    def test_runtime_accepts_processed_output_and_rejects_input_from_the_admet_port(self):
        import test_pipeline_execution as execution
        from biomolexplorer.pipeline import PipelineService
        fixture=execution.PipelineExecutionTests();fixture.setUp();self.addCleanup(fixture.tearDown)
        producer=new_stage('admet');consumer=new_stage('retrieve_pubchem')
        root=fixture.store.project_dir(fixture.project_id)
        raw=root/'raw.csv';raw.write_text('molecule_chembl_id,canonical_smiles\nA,CCO\n')
        output=root/'B01__admet__bbb-pos.csv';output.write_text('molecule_chembl_id,canonical_smiles,BBB,HIA\nA,CCO,BBB+,HIA+\n')
        service=PipelineService(fixture.store);self.addCleanup(service.close)
        uid=fixture.store.user(fixture.token)['id'];results={producer['id']:[str(raw),str(output)]}
        consumer['bindings']={'base_input_path':{'stage':producer['id'],'selector':raw.name}}
        with self.assertRaisesRegex(ValueError,'arquivo escolhido|arquivo compatível|resultado'):
            service._resolve(fixture.project_id,uid,consumer,results,cache_only=True,pipeline=[producer,consumer])
        consumer['bindings']['base_input_path']['selector']=output.name
        params=service._resolve(fixture.project_id,uid,consumer,results,cache_only=True,pipeline=[producer,consumer])
        self.assertEqual(params['input_file'],output.name)
        self.assertTrue(raw.is_file());self.assertTrue(output.is_file())
