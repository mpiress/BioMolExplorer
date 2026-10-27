"""Producer/run identities, exact lineage and compact multi-source selection."""
import copy
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from biomolexplorer.artifacts import records, reference, validate_reference, index_run, locator
from biomolexplorer.catalog import new_stage, operation_fields
from biomolexplorer.bindings import pack, sources
from biomolexplorer.pipeline import PipelineService, validate_pipeline
import test_pipeline_execution as execution


class ArtifactRegistryTests(unittest.TestCase):
    def test_effective_batch_context_survives_catalog_rediscovery_and_reuse(self):
        stage=new_stage('docking_vina')
        path='/project/runs/run/'+stage['id']+'/batch/docking_results.csv'
        stage.update(artifacts=[path],batches=[{'artifacts':[path],
            'artifact_context':{'target':'effective-target'}}])
        index_run([stage],'run')
        self.assertEqual(records(stage,'run')[0]['context'],{'target':'effective-target'})
        stage.update(reused=True,reused_from_run='run')
        index_run([stage],'next-run')
        self.assertEqual(records(stage,'next-run')[0]['context'],{'target':'effective-target'})
        self.assertEqual(records(stage,'next-run')[0]['run_id'],'run')

    def test_identity_separates_producers_and_occurrences_and_survives_relocation(self):
        a=new_stage('retrieve_compounds');b=new_stage('retrieve_compounds')
        a.update(artifacts=['/project/runs/old/'+a['id']+'/compounds.csv'])
        b.update(artifacts=['/project/runs/old/'+b['id']+'/compounds.csv'])
        old=records(a,'old')[0];other=records(b,'old')[0]
        self.assertNotEqual(old['artifact_id'],other['artifact_id'])
        a['artifacts']=['/migrated/runs/new/'+a['id']+'/compounds.csv']
        new=records(a,'new')[0]
        self.assertEqual(old['artifact_id'],new['artifact_id'])
        self.assertNotEqual(old['occurrence_id'],new['occurrence_id'])
        a.update(artifact_records=[new],reused_from_run='new')
        self.assertEqual(records(a,'reused')[0]['run_id'],'new')

    def test_redocking_cache_materialization_keeps_downstream_identity(self):
        stage='a'*32
        direct=f'/project/runs/old/{stage}/individual/001_4M0E/structures/Prepared/4M0E_A.dockprep.pdbqt'
        cached=f'/project/runs/new/{stage}/individual/001_4M0E/structures/cache/redocking-0/structures/Prepared/4M0E_A.dockprep.pdbqt'
        self.assertEqual(locator(stage,direct),locator(stage,cached))

    def test_five_thousand_files_render_only_one_page_and_bulk_update_once(self):
        from biomolexplorer.ui.artifact_picker import ArtifactPicker
        updates=[]
        stage_id='a'*32
        rows=[(SimpleNamespace(value=False,label=f'{i:04d}.csv',available=True,on_change=None),
               {'stage':stage_id,'selector':f'{i:04d}.csv'}) for i in range(5000)]
        picker=ArtifactPicker(rows,{'stage:'+stage_id:'Retrieval'},SimpleNamespace(update=lambda:updates.append(1)))
        self.assertEqual(len(picker.available.options),25)
        picker.add_all_sources()
        self.assertEqual(sum(flag.value for flag,_ in rows),5000)
        self.assertEqual(len(picker.selected.controls),25)
        self.assertEqual(len(updates),1)

    def test_adding_a_duplicate_name_preserves_identity_and_resolves_by_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            stage=new_stage('retrieve_compounds');root=Path(tmp)/stage['id']
            first=root/'one'/'compounds.csv';first.parent.mkdir(parents=True);first.write_text('data')
            stage['artifacts']=[str(first)]
            original=records(stage,'first')[0];ref=reference(original)
            second=root/'two'/'compounds.csv';second.parent.mkdir();second.write_text('data')
            stage['artifacts'].append(str(second))
            updated=records(stage,'second')
            self.assertEqual(original['artifact_id'],updated[0]['artifact_id'])
            self.assertNotEqual(original['selector'],updated[0]['selector'])
            self.assertEqual(str(validate_reference(ref,stage['artifacts'])),str(first))

    def test_stale_ambiguous_and_mismatched_references_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);files=[]
            for folder in ('one','two'):
                p=root/folder/'compounds.csv';p.parent.mkdir();p.write_text('data');files.append(str(p))
            stage=new_stage('retrieve_compounds');stage['artifacts']=files
            ref=reference(records(stage,'run')[0]);validate_reference(ref,files)
            with self.assertRaisesRegex(ValueError,'ambíguo'):validate_reference(dict(stage=stage['id'],selector='compounds.csv'),files)
            with self.assertRaisesRegex(ValueError,'identidade'):validate_reference(dict(ref,artifact_id='f'*64),files)
            Path(files[0]).unlink()
            with self.assertRaisesRegex(ValueError,'ausente'):validate_reference(ref,files)


class ArtifactSelectionIntegrationTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def test_two_equal_filenames_remain_distinct_in_individual_and_merge_runs(self):
        from biomolexplorer.ui.file_selection import FileSelection
        from biomolexplorer.result_files import ResultFiles
        for mode in ('individual','merge'):
            with self.subTest(mode=mode):
                a,b=new_stage('retrieve_compounds'),new_stage('retrieve_compounds')
                a['name']='Retrieval AChE';b['name']='Retrieval BChE'
                b['parameters']['search_term']='CHEMBL221'
                consumer=new_stage('admet');consumer['bindings']['base_input_path']=pack([{'stage':a['id']},{'stage':b['id']}])
                self.save([a,b,consumer]);manager,calls=self.immediate_manager()
                with patch('biomolexplorer.pipeline.JobManager',manager):
                    service=PipelineService(self.store)
                    try:
                        run=self.wait(service.submit(self.token,self.project_id,reuse_results=False))
                        pending=run['stages'][-1];form=FileSelection(run,pending)
                        picker=form.pickers['base_input_path']
                        picker.add_all_sources();picker.add_all_sources()
                        self.assertEqual(sum(flag.value for flag,_ in picker.rows),2)
                        form.processing.value=mode;chosen=form.read();validate_pipeline([a,b,chosen])
                        restored=FileSelection(run,dict(pending,configuration=chosen))
                        self.assertEqual(restored.read()['bindings'],chosen['bindings'])
                        finished=self.wait(service.resume(self.token,run['id'],chosen))
                        self.assertEqual(finished['status'],'succeeded',finished.get('error'))
                        output=finished['stages'][-1]
                        self.assertEqual(len(output['artifacts']),2 if mode=='individual' else 1)
                        all_origins=set()
                        for artifact in output['artifact_records']:
                            inputs=artifact['inputs']['base_input_path']
                            self.assertEqual(len(inputs),1 if mode=='individual' else 2)
                            all_origins.update(r['stage_id'] for r in inputs)
                            self.assertTrue(all(r['run_id']==finished['id'] for r in inputs))
                        self.assertEqual(all_origins,{a['id'],b['id']})
                        inspected=ResultFiles(self.store).files(self.token,self.project_id,finished['id'],consumer['id'])
                        self.assertTrue(all('Retrieval' in r['description'] for r in inspected))
                        self.assertEqual({r['artifact_id'] for r in inspected},{r['artifact_id'] for r in output['artifact_records']})
                    finally:service.close()

    def test_large_picker_search_batch_selection_removal_and_missing_file(self):
        from biomolexplorer.ui.file_selection import FileSelection
        a,b=new_stage('retrieve_compounds'),new_stage('retrieve_compounds')
        consumer=new_stage('fingerprints');consumer['bindings']['base_input_path']=pack([{'stage':a['id']},{'stage':b['id']}])
        root=self.store.project_dir(self.project_id)
        completed=[]
        for source in (a,b):
            paths=[]
            for i in range(101):
                p=root/source['id']/f'{i:03d}.csv';p.parent.mkdir(exist_ok=True)
                p.write_text('molecule_chembl_id,canonical_smiles\nM1,CCO\n');paths.append(str(p))
            completed.append(dict(source,status='succeeded',artifacts=paths))
        pending={'name':consumer['name'],'configuration':consumer}
        form=FileSelection({'id':'run','stages':completed},pending)
        picker=form.pickers['base_input_path'];self.assertLessEqual(len(picker.available.options),25)
        picker.available.value='0';picker.add_file();picker.add_file()
        self.assertEqual(sum(f.value for f,_ in picker.rows),1)
        picker.source.value='stage:'+b['id'];picker.source_changed(None);picker.add_all()
        self.assertEqual(sum(f.value for f,_ in picker.rows),102)
        self.assertEqual(len(picker.selected.controls),25)
        picker.search.value='100.csv';picker.search_changed(None)
        self.assertEqual(len(picker.available.options),1)
        picker.set_selected([0],False);self.assertEqual(sum(f.value for f,_ in picker.rows),101)
        picker.add_all_sources();self.assertEqual(sum(f.value for f,_ in picker.rows),202)
        saved=form.read();self.assertEqual(len(sources(saved['bindings']['base_input_path'])),202)
        Path(completed[0]['artifacts'][0]).unlink()
        reopened=FileSelection({'stages':completed},dict(pending,configuration=saved))
        with self.assertRaisesRegex(ValueError,'indisponível'):reopened.read()
        picker.clear();self.assertFalse(any(f.value for f,_ in picker.rows))

    def test_docking_batch_tables_are_all_offered_without_redundant_nested_tables(self):
        from biomolexplorer.ui.file_selection import FileSelection
        producer=new_stage('docking_vina');consumer=new_stage('docking_dock6')
        consumer['bindings']['base_selected_mols']={'stage':producer['id']}
        paths=[]
        for batch in ('one','two'):
            for nested in ('','receptor'):
                path=self.store.project_dir(self.project_id)/producer['id']/batch/nested/'docking_results.csv'
                path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text('molecule_chembl_id,canonical_smiles,engine,score,conformer_file\nM1,CCO,vina,-1,pose.pdbqt\n')
                paths.append(str(path))
        run={'stages':[dict(producer,status='succeeded',artifacts=paths)]}
        form=FileSelection(run,{'name':consumer['name'],'configuration':consumer},compound_references={})
        self.assertEqual(len(form.rows['base_selected_mols']),2)
        self.assertTrue(all('/receptor/' not in r['selector'] for _,r in form.rows['base_selected_mols']))

    def test_guided_editor_shares_picker_and_preserves_other_origin_when_connection_changes(self):
        from biomolexplorer.ui.input_editor import InputEditor
        a,b,c=[new_stage('retrieve_compounds') for _ in range(3)]
        stage=new_stage('admet');stage['bindings']['base_input_path']=pack([{'stage':a['id']},{'stage':b['id']}])
        items={}
        for origin in (a,b,c):
            p=self.store.project_dir(self.project_id)/origin['id']/'compounds.csv';p.parent.mkdir();p.write_text('molecule_chembl_id,canonical_smiles\nM1,CCO\n')
            items[origin['id']]=dict(origin,status='succeeded',artifacts=[str(p)])
        ui=SimpleNamespace(current={'pipeline':[a,b,c,stage]},artifact_items=items,artifact_choices={},page=SimpleNamespace(update=lambda:None))
        field=next(f for f in operation_fields('admet') if f['name']=='base_input_path')
        editor=InputEditor(ui,stage,field,[],True)
        editor.selection.pickers['base_input_path'].add_all_sources()
        selected=editor.read();self.assertEqual(len(sources(selected)),2)
        editor.rows[1][0].value='stage:'+c['id'];editor.rows[1][1].value='auto';editor.changed()
        remaining=sources(editor.read());self.assertEqual({r['stage'] for r in remaining},{a['id']})
        self.assertTrue(all(r.get('artifact_id') for r in remaining))
