"""Typed ports, entity selections, registry persistence and reporting."""
import io
import json
import tempfile
import time
import unittest
import zipfile
from pathlib import Path
from types import SimpleNamespace
from uuid import uuid4

from biomolexplorer.artifact_contracts import INPUT_PORTS,OUTPUT_PORTS,operation_contract
from biomolexplorer.artifact_registry import SemanticArtifactResolver,index_entities,list_artifacts,list_entities,register_run
from biomolexplorer.artifacts import index_run
from biomolexplorer.catalog import TITLES,new_stage
from biomolexplorer.experimental_report import build_report,report_docx,report_html,report_pdf
from biomolexplorer.flow import connect
from biomolexplorer.pipeline import validate_pipeline
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.selection_views import csv_view,path_matches
from biomolexplorer.semantic_references import migrate_pipeline_references,migrate_project_references,parse_selection
from biomolexplorer.stage_cache import artifact_manifest
from biomolexplorer.workspace import WorkspaceStore


class SemanticArtifactTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.store=WorkspaceStore(self.tmp.name)
        self.token=self.store.register('Owner','semantic@example.org','owner-password')
        self.project=self.store.create_project(self.token,'Alzheimer semantic pipeline');self.pid=self.project['id']

    def tearDown(self):self.tmp.cleanup()

    def test_every_scientific_operation_has_one_shared_typed_contract(self):
        for operation in TITLES:
            with self.subTest(operation=operation):
                contract=operation_contract(operation)
                if operation!='import_results':self.assertTrue(contract['outputs'])
                self.assertEqual(set(contract['inputs']),set(INPUT_PORTS.get(operation,{})))
                self.assertEqual(set(contract['outputs']),set(OUTPUT_PORTS.get(operation,{})))

    def test_future_pdb_entity_selection_is_valid_before_execution(self):
        source=new_stage('retrieve_structures');source['name']='PDB ACHE'
        target=new_stage('prepare_structures')
        connect([source,target],source['id'],target['id'],'base_input_path')
        reference=target['bindings']['base_input_path']
        self.assertEqual(reference['output_port'],'receptor_structures')
        reference['selection']=parse_selection('pdb:4M0E, chain:A')
        self.assertEqual(validate_pipeline([source,target]),[source['id'],target['id']])
        self.assertEqual(reference['selection'],{'pdb_id':'4M0E','chain_id':'A'})

    def test_editor_preserves_future_port_selection_and_reference_mode(self):
        from biomolexplorer.ui.input_editor import InputEditor
        source=new_stage('retrieve_structures');target=new_stage('prepare_structures')
        target['bindings']['base_input_path']={'stage':source['id'],'output_port':'receptor_structures',
            'selector':'auto','mode':'current','selection':{'pdb_id':'4M0E','chain_id':'A'}}
        ui=SimpleNamespace(current={'pipeline':[source,target]},page=SimpleNamespace(update=lambda:None))
        field=next(field for field in __import__('biomolexplorer.catalog',fromlist=['operation_fields']).operation_fields('prepare_structures')
                   if field['name']=='base_input_path')
        editor=InputEditor(ui,target,field,[],True)
        self.assertFalse(editor.rows[0][1].visible)
        labels=[control.value for control in editor.rows[0][2].controls if hasattr(control,'value') and isinstance(control.value,str)]
        self.assertTrue(any('receptor_structures' in label for label in labels))
        saved=editor.read()
        self.assertEqual(saved['output_port'],'receptor_structures')
        self.assertEqual(saved['mode'],'current')
        self.assertEqual(saved['selection'],{'pdb_id':'4M0E','chain_id':'A'})

    def test_compound_and_structure_views_preserve_sources(self):
        root=Path(self.tmp.name);source=root/'compounds.csv'
        source.write_text('molecule_chembl_id,canonical_smiles,TPSA\nCHEMBL1,CCO,20\nCHEMBL2,CCN,90\n')
        selected=csv_view(source,{'entity_ids':['CHEMBL2'],'properties':{'TPSA':{'min':80}}},root/'views')
        self.assertEqual(source.read_text().count('\n'),3)
        self.assertIn('CHEMBL2',selected.read_text());self.assertNotIn('CHEMBL1',selected.read_text())
        self.assertTrue(path_matches(root/'4M0E_A.pdb',{'pdb_id':'4M0E','chain_id':'A'}))
        self.assertFalse(path_matches(root/'1ABC_A.pdb',{'pdb_id':'4M0E'}))
        with self.store.connect() as db:
            self.assertEqual(index_entities(db,self.pid,'a'*64,source,'compounds'),2)
            self.assertEqual([row['entity_id'] for row in list_entities(db,self.pid,'a'*64)],['CHEMBL1','CHEMBL2'])

    def test_executor_resolves_entity_selection_at_the_adapter_boundary(self):
        producer=new_stage('retrieve_compounds');consumer=new_stage('fingerprints')
        source=self.store.project_dir(self.pid)/'source'/'compounds.csv';source.parent.mkdir()
        source.write_text('molecule_chembl_id,canonical_smiles\nCHEMBL1,CCO\nCHEMBL2,CCN\n')
        reference={'stage':producer['id'],'output_port':'compounds','selector':'auto',
                   'selection':{'entity_ids':['CHEMBL2']}}
        consumer['bindings']['base_input_path']=reference
        service=PipelineService(self.store)
        try:
            folder,filename,selected=service._materialize_inputs(self.pid,consumer,'base_input_path',[reference],
                {producer['id']:[str(source)]},self.store.project_dir(self.pid)/'adapter')
            materialized=folder/filename
            self.assertIn('CHEMBL2',materialized.read_text());self.assertNotIn('CHEMBL1',materialized.read_text())
            self.assertEqual(len(selected),1)
        finally:service.close()

    def test_registry_manifest_and_report_are_persistent_and_path_free(self):
        stage=new_stage('retrieve_structures');stage['block_id']='PDB_ACHE'
        path=self.store.project_dir(self.pid)/'runs'/'run'/stage['id']/'artifacts'/'4M0E.pdb'
        path.parent.mkdir(parents=True);path.write_text('HEADER    4M0E\n')
        stage.update(status='succeeded',artifacts=[str(path)],artifact_manifest=artifact_manifest([str(path)]),
            configuration=dict(stage),processing_version='test-v1',producer_software='test-tool:1')
        run_id=uuid4().hex;index_run([stage],run_id)
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',(run_id,self.pid,self.store.user(self.token)['id'],'succeeded',json.dumps([stage]),time.time(),time.time(),None))
            register_run(db,self.pid,run_id,[stage])
        with self.store.connect() as db:
            artifacts=list_artifacts(db,self.pid,stage_id=stage['id'],run_id=run_id)
            manifest=json.loads(db.execute('SELECT manifest FROM execution_manifests WHERE run_id=?',(run_id,)).fetchone()[0])
        self.assertEqual(len(artifacts),1);self.assertEqual(artifacts[0]['state'],'available')
        self.assertIn(artifacts[0]['occurrence_id'],sum(manifest['outputs'].values(),[]))
        report=build_report(self.store,self.token,self.pid)
        self.assertEqual(report['schema_version'],1)
        self.assertNotIn(str(self.store.root),json.dumps(report['pipeline']))
        self.assertEqual(report['artifacts'][0]['content_sha256'],stage['artifact_manifest'][str(path.resolve())])
        html=report_html(report,'en')
        self.assertIn(b'Reproducible Experimental Report',html)
        self.assertIn(b'data:image/png;base64,',html)
        docx=report_docx(report,'en')
        with zipfile.ZipFile(io.BytesIO(docx)) as document:
            self.assertIn('word/media/logo.png',document.namelist())
            self.assertIn(b'Reproducible Experimental Report',document.read('word/document.xml'))
        pdf=report_pdf(report,'en')
        self.assertTrue(pdf.startswith(b'%PDF-'))
        self.assertIn(b'/Image',pdf)
        with self.store.connect() as db:
            self.assertEqual(list_entities(db,self.pid,stage['artifact_records'][0]['artifact_id'])[0]['entity_id'],'4M0E')

    def test_report_summarizes_tables_embeds_figures_and_keeps_json_lightweight(self):
        from PIL import Image
        from biomolexplorer.experimental_report import public_report
        stage=new_stage('admet');root=self.store.project_dir(self.pid)/'report-evidence';root.mkdir()
        table=root/'metrics.csv';table.write_text('compound,score,rmsd\nA,1.5,2.0\nB,2.5,4.0\n')
        figure=root/'result.png';Image.new('RGB',(320,180),'white').save(figure)
        stage.update(status='succeeded',artifacts=[str(table),str(figure)])
        run_id=uuid4().hex
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',(run_id,self.pid,self.store.user(self.token)['id'],'succeeded',json.dumps([stage]),time.time(),time.time(),None))
        report=build_report(self.store,self.token,self.pid)
        self.assertEqual(report['results']['tables'][0]['numeric_summary']['score']['mean'],2)
        self.assertTrue(report['results']['figures'][0]['data'])
        self.assertNotIn('data',public_report(report)['results']['figures'][0])
        docx=report_docx(report,'en')
        with zipfile.ZipFile(io.BytesIO(docx)) as document:
            self.assertIn('word/header1.xml',document.namelist())
            self.assertTrue(any(name.startswith('word/media/figure1.') for name in document.namelist()))
            self.assertIn(b'w:orient="landscape"',document.read('word/document.xml'))
        pdf=report_pdf(report,'en')
        self.assertTrue(pdf.startswith(b'%PDF-'))
        self.assertIn(b'/MediaBox [ 0 0 841.8898 595.2756 ]',pdf)

    def test_report_bounds_database_history_before_materializing_it(self):
        now=time.time()
        with self.store.connect() as db:
            db.executemany('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',[(uuid4().hex,self.pid,
                self.store.user(self.token)['id'],'succeeded','[]',now+index,now+index,None) for index in range(55)])
            rows=[]
            for index in range(275):
                identity=f'{index:064x}'
                rows.append((identity,identity,self.pid,'run',None,'stage',None,None,'file',1,
                    str(self.store.project_dir(self.pid)/('missing-'+str(index))),None,'available','{}','[]','{}',now+index))
            db.executemany('INSERT INTO artifact_instances VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',rows)
        report=build_report(self.store,self.token,self.pid)
        self.assertEqual(len(report['executions']),50)
        self.assertEqual(len(report['artifacts']),250)
        self.assertEqual(report['document_scope']['execution_count'],55)
        self.assertEqual(report['document_scope']['artifact_count'],275)
        self.assertTrue(report['document_scope']['bounded'])

    def test_current_manifest_port_selects_records_without_filename_inference(self):
        producer=new_stage('prepare_structures');root=self.store.project_dir(self.pid)/'manifest'
        receptor=root/'opaque-a';compound=root/'opaque-b';root.mkdir()
        receptor.write_text('receptor');compound.write_text('compound')
        records=[{'path':str(receptor),'output_port':'prepared_receptors','artifact_id':'a'*64},
                 {'path':str(compound),'output_port':'compounds','artifact_id':'b'*64}]
        pipeline=[dict(producer,_execution_artifact_records=records)]
        resolver=SemanticArtifactResolver(self.store,self.pid,pipeline,{producer['id']:[str(receptor),str(compound)]})
        self.assertEqual(resolver.resolve({'stage':producer['id'],'output_port':'prepared_receptors'}),[str(receptor)])

    def test_pinned_manifest_resolves_shared_occurrence_through_logical_consumer(self):
        producer=new_stage('retrieve_compounds');consumer=new_stage('fingerprints')
        path=self.store.project_dir(self.pid)/'shared.csv';path.write_text('molecule_chembl_id,canonical_smiles\nCHEMBL1,CCO\n')
        record={'occurrence_id':'c'*64,'artifact_id':'a'*64,'path':str(path),'output_port':'compounds',
            'artifact_type':'compounds','run_id':'d'*32,'stage_id':producer['id'],'content_sha256':'e'*64}
        producer['artifact_records']=[record];consumer['artifact_records']=[dict(record)]
        run_id='d'*32
        with self.store.connect() as db:register_run(db,self.pid,run_id,[producer,consumer])
        resolver=SemanticArtifactResolver(self.store,self.pid,[],{})
        reference={'stage':consumer['id'],'output_port':'compounds','mode':'pinned','run_id':run_id,'artifact_id':'a'*64}
        self.assertEqual(resolver.resolve(reference),[str(path)])

    def test_legacy_reference_migration_is_idempotent_and_records_ambiguity(self):
        source=new_stage('retrieve_structures');target=new_stage('prepare_structures')
        target['bindings']['base_input_path']={'stage':source['id'],'selector':'4M0E.pdb'}
        migrated,issues,changed=migrate_pipeline_references([source,target])
        self.assertTrue(changed);self.assertEqual(issues,[])
        self.assertEqual(migrated[1]['bindings']['base_input_path']['output_port'],'receptor_structures')
        again,issues,changed=migrate_pipeline_references(migrated)
        self.assertFalse(changed);self.assertEqual(again,migrated)

        incompatible=new_stage('fingerprints')
        incompatible['bindings']['base_input_path']={'stage':source['id'],'selector':'legacy.csv'}
        with self.store.connect() as db:
            db.execute('UPDATE projects SET pipeline=? WHERE id=?',(json.dumps([source,incompatible]),self.pid))
            issues,changed=migrate_project_references(db,self.pid,[source,incompatible])
            stored=db.execute('SELECT reason,details FROM artifact_migration_issues WHERE project_id=?',(self.pid,)).fetchall()
        self.assertFalse(changed);self.assertEqual(issues[0]['reason'],'incompatible')
        self.assertEqual(len(stored),1);self.assertIn(source['id'],stored[0]['details'])


if __name__=='__main__':unittest.main()
