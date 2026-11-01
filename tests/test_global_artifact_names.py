"""Global aliases preserve science bytes, identities and exact input lineage."""
import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from biomolexplorer.artifact_naming import (public_filename,parse_filename,name_records,
    publish,publication_metadata,download_name,bundle_manifest,native_name,PUBLIC_OPERATIONS)
from biomolexplorer.artifacts import index_run,records,reference,validate_reference
from biomolexplorer.catalog import new_stage,TITLES
from biomolexplorer.bindings import pack
from biomolexplorer.stage_cache import stage_key
from biomolexplorer.workspace import WorkspaceStore
from biomolexplorer.pipeline import PipelineService
import test_pipeline_execution as execution


class GlobalNamesTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)

    def output(self,stage,name,text='data',inputs=None):
        path=self.root/stage['id']/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text)
        return dict(stage,configuration=copy.deepcopy(stage),status='succeeded',artifacts=[str(path)],input_files=inputs or {})

    def flow(self):
        source=new_stage('retrieve_compounds');source['name']='ACHE'
        mol=self.output(source,'ACHE_MOLS.csv','molecule_chembl_id,canonical_smiles\nM1,CCO\n')
        admet=new_stage('admet');admet['bindings']['base_input_path']={'stage':source['id']}
        bbb=self.output(admet,'ACHE_MOL_BBB+.csv','molecule_chembl_id,canonical_smiles,BBB,HIA\nM1,CCO,True,True\n',{'base_input_path':mol['artifacts']})
        fps=new_stage('fingerprints');fps['bindings']['base_input_path']={'stage':admet['id']}
        fp=self.output(fps,'morgan_very_long_old_hierarchy.csv','molecule_chembl_id,canonical_smiles,fingerprint\nM1,CCO,"[1,0]"\n',{'base_input_path':bbb['artifacts']})
        similarity=new_stage('similarity');similarity['bindings']['base_input_path']={'stage':fps['id']}
        sim=self.output(similarity,'matrix.csv','source,target,value\nM1,M1,1\n',{'base_input_path':fp['artifacts']})
        return [mol,bbb,fp,sim]

    def test_signs_extensions_and_parser_roundtrip(self):
        for kind in ('BBB+','BBB-','HIA+','HIA-','Fingerprints-Morgan','RMSD'):
            for ext in ('csv','sdf','pdb','pdbqt','mol2','csv.gz','biomol-view.json','nc'):
                with self.subTest(kind=kind,ext=ext):
                    name=public_filename('ACHE','MOL','ADMET',kind,ext,'abcdef12')
                    parsed=parse_filename(name)
                    self.assertEqual(parsed['result_type'],kind);self.assertEqual(parsed['extension'],ext)
                    self.assertEqual(parsed['target'],'ACHE');self.assertEqual(parsed['identifier'],'abcdef12')
        self.assertEqual(parse_filename('B01__admet__bbb-neg.csv')['scheme'],0)
        with self.assertRaises(ValueError):public_filename('ACHE','MOL','ADMET','BBB+','../csv')

    def test_immediate_source_not_ancestral_hierarchy_and_real_targets(self):
        stages=self.flow();before={p:Path(p).read_bytes() for s in stages for p in s['artifacts']}
        index_run(stages,'run');tags=[parse_filename(s['artifact_records'][0]['public_name']) for s in stages]
        self.assertEqual([p['source'] for p in tags],['ChEMBL','MOL','BBB+','Fingerprints-Morgan'])
        self.assertEqual({p['target'] for p in tags},{'ACHE'})
        self.assertEqual(tags[1]['result_type'],'BBB+');self.assertEqual(tags[2]['result_type'],'Fingerprints-Morgan')
        self.assertTrue(tags[3]['result_type'].startswith('Similarity-'))
        for s in stages:
            r=s['artifact_records'][0];self.assertLess(len(r['public_name']),130)
            self.assertEqual(validate_reference(reference(r),s['artifacts']),r['path'])
        self.assertEqual({p:Path(p).read_bytes() for p in before},before)
        parent=stages[1]['artifact_records'][0];child=stages[2]['artifact_records'][0]
        self.assertEqual(child['naming']['source_artifact_ids'],[parent['artifact_id']])
        self.assertEqual(child['parent_artifact_ids'],[parent['artifact_id']])

    def test_collision_suffix_is_stable_across_order_and_restart(self):
        a,b=new_stage('retrieve_compounds'),new_stage('retrieve_compounds')
        for s in (a,b):s['name']='ACHE'
        stages=[self.output(s,'ACHE_MOLS.csv','molecule_chembl_id,canonical_smiles\nM1,CCO\n') for s in (a,b)]
        index_run(stages,'run');names={r['artifact_id']:r['public_name'] for s in stages for r in s['artifact_records']}
        self.assertEqual(len(set(names.values())),2)
        restored=json.loads(json.dumps(stages));index_run(restored[::-1],'run')
        self.assertEqual(names,{r['artifact_id']:r['public_name'] for s in restored for r in s['artifact_records']})
        relocated=json.loads(json.dumps(stages).replace(str(self.root),str(self.root/'new-location')))
        index_run(relocated,'run')
        self.assertEqual(names,{r['artifact_id']:r['public_name'] for s in relocated for r in s['artifact_records']})

    def test_all_blocks_have_a_public_operation_and_truthful_content_type(self):
        self.assertEqual(set(TITLES)-set(PUBLIC_OPERATIONS),set())
        examples={'import_results':('user.csv','molecule_chembl_id,canonical_smiles\nM1,CCO\n','MOL'),
            'retrieve_compounds':('x_MOLS.csv','molecule_chembl_id,canonical_smiles\nM1,CCO\n','MOL'),
            'expand_similar_compounds':('x_SIMS.csv','molecule_chembl_id,canonical_smiles\nM1,CCO\n','SIM'),
            'retrieve_pubchem':('compounds.csv','molecule_chembl_id,canonical_smiles\nM1,CCO\n','SIM'),
            'retrieve_zinc':('molecules.sdf','SDF','MOL'),
            'retrieve_structures':('4M0E.pdb','PDB','Complex-4M0E'),
            'prepare_structures':('4M0E_A.dockprep.pdbqt','PDBQT','Receptor-4M0E_A'),
            'admet':('x_BBB-.csv','molecule_chembl_id,BBB,HIA\nM1,False,True\n','BBB-'),
            'fingerprints':('maccs_x.csv','molecule_chembl_id,fingerprint\nM1,"[1,0]"\n','Fingerprints-MACCS'),
            'similarity':('data.csv','source,target,value\nM1,M1,1\n','Similarity-'),
            'graphs':('full_graph.csv','source,target,value\nM1,M1,1\n','Graph-Full'),
            'redocking':('pdb_codes.csv','PDB_CODE,LIGAND,RESNUM,CHAIN,RMSD\n4M0E,LIG,1,A,0.2\n','RMSD'),
            'docking_vina':('docking_results.csv','molecule_chembl_id,engine,score\nM1,vina,-1\n','Scores'),
            'docking_dock6':('docking_results.csv','molecule_chembl_id,engine,score\nM1,dock6,-1\n','Scores'),
            'consensus':('scores.csv','molecule_chembl_id,vina,dock6\nM1,-1,-2\n','Scores'),
            'dynamics_dock6':('energies.json','{"score":-1}','AmberScore')}
        self.assertEqual(set(examples),set(TITLES))
        for op,(name,text,kind) in examples.items():
            with self.subTest(op=op):
                stage=self.output(new_stage(op),name,text);index_run([stage],'run')
                parsed=parse_filename(stage['artifact_records'][0]['public_name']);self.assertTrue(parsed['result_type'].startswith(kind),parsed)

    def test_no_artificial_target_for_independent_molecule_library(self):
        source=self.output(new_stage('retrieve_zinc'),'molecules.sdf')
        fp=new_stage('fingerprints');fp['parameters']['target']='MeuAlvo';fp['bindings']['base_input_path']={'stage':source['id']}
        child=self.output(fp,'morgan.csv','molecule_chembl_id,fingerprint\nA,"[1,0]"\n',{'base_input_path':source['artifacts']})
        index_run([source,child],'run')
        self.assertEqual(child['artifact_records'][0]['naming']['targets'],[])
        self.assertEqual(parse_filename(child['artifact_records'][0]['public_name'])['target'],'NA')
        provided=new_stage('fingerprints');provided['provided_results']={'target':'BCHE'}
        imported=self.output(provided,'morgan.csv','molecule_chembl_id,fingerprint\nA,"[1,0]"\n')
        index_run([imported],'run')
        self.assertEqual(imported['artifact_records'][0]['naming']['target'],'BCHE')

    def test_docking_target_comes_from_receptor_not_compound_collection(self):
        mol,bbb,*_=self.flow();pdb=new_stage('retrieve_structures');pdb['parameters']['target']='BCHE'
        receptor=self.output(pdb,'6ZWI.pdb')
        docking=new_stage('docking_vina');docking['bindings']={'base_input_path':{'stage':pdb['id']},'base_selected_mols':{'stage':bbb['id']}}
        output=self.output(docking,'docking_results.csv','molecule_chembl_id,engine,score\nM1,vina,-1\n',{'base_input_path':receptor['artifacts'],'base_selected_mols':bbb['artifacts']})
        index_run([mol,bbb,receptor,output],'run');name=parse_filename(output['artifact_records'][0]['public_name'])
        self.assertEqual(name['target'],'BCHE');self.assertEqual(name['source'],'BBB+')

    def test_publication_preserves_native_bundles_and_legacy_files(self):
        stage=new_stage('redocking');native=self.output(stage,'Prepared/4M0E_A.dockprep.pdbqt')
        path=Path(native['artifacts'][0]);before=path.read_bytes()
        names=publish(stage,native['artifacts'],self.root/stage['id'],run_id='run')
        self.assertEqual(names,native['artifacts']);self.assertEqual(path.read_bytes(),before)
        meta=publication_metadata(path);self.assertTrue(parse_filename(meta['public_name']))
        self.assertEqual(validate_reference(reference(records(dict(native,artifact_records=[]),'run')[0]),names),str(path))
        self.assertEqual(publish(stage,names,self.root/stage['id'],run_id='run'),names)

    def test_export_names_and_bundle_metadata_keep_full_references(self):
        stages=self.flow();index_run(stages,'run');record=stages[1]['artifact_records'][0]
        self.assertEqual(download_name(record),record['public_name'])
        figure=download_name(record,result_type='Figure-BBB+',extension='png')
        self.assertEqual(parse_filename(figure)['result_type'],'Figure-BBB+')
        self.assertNotEqual(figure,record['public_name'])
        manifest=bundle_manifest([record],stages)
        self.assertEqual(manifest['artifacts'][0]['artifact_id'],record['artifact_id'])
        self.assertEqual(manifest['producers'][stages[1]['id']]['parameters'],stages[1]['parameters'])

    def test_canonical_prepared_bundle_routes_to_native_names_and_validates_content(self):
        from biomolexplorer.input_validation import validate_bundle
        from biomolexplorer.docking_inputs import prepared_receptor
        from biomolexplorer.fingerprint_selection import filename_kind
        samples=[('Receptor-4M0E_A','pdbqt','ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n'),
            ('Pairs','csv','PDB_CODE,LIGAND,RESNUM,CHAIN\n4M0E,LIG,1,A\n'),('Centers','csv','4M0E_LIG_1A\n0\n0\n0\n')]
        files=[]
        for i,(kind,ext,text) in enumerate(samples):
            p=self.root/public_filename('ACHE','PDB','Preparation',kind,ext,hashlib.sha256(str(i).encode()).hexdigest()[:8]);p.write_text(text);files.append(p)
        validate_bundle(files,'prepared_structures')
        self.assertTrue(prepared_receptor(files[0]));self.assertEqual(native_name(files[0]),'4M0E_A.dockprep.pdbqt')
        self.assertEqual(filename_kind(public_filename('ACHE','MOL','Fingerprints','Fingerprints-Morgan','csv')),'morgan')
        files[0].write_text('invalid')
        with self.assertRaises(ValueError):validate_bundle(files,'prepared_structures')


class NamingPersistenceTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    finish_selected=execution.PipelineExecutionTests.finish_selected
    immediate_manager=execution.PipelineExecutionTests.immediate_manager

    def test_restart_incremental_reuse_and_archive_restore_keep_aliases_and_ids(self):
        a,b=new_stage('retrieve_compounds'),new_stage('admet');a['name']='ACHE'
        b['bindings']['base_input_path']={'stage':a['id']};self.save([a,b]);manager,calls=self.immediate_manager()
        with patch('biomolexplorer.pipeline.JobManager',manager):
            service=PipelineService(self.store)
            try:run=self.finish_selected(service,service.submit(self.token,self.project_id))
            finally:service.close()
            self.assertEqual(run['status'],'succeeded',run.get('error'))
            saved=[(r['artifact_id'],r['public_name']) for s in run['stages'] for r in s['artifact_records']]
            hashes={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for s in run['stages'] for p in s['artifacts']}
            restored=WorkspaceStore(self.store.root);service=PipelineService(restored)
            try:again=self.finish_selected(service,service.submit(self.token,self.project_id))
            finally:service.close()
            self.assertEqual(len(calls),2)
            self.assertEqual(saved,[(r['artifact_id'],r['public_name']) for s in again['stages'] for r in s['artifact_records']])
            self.assertEqual(hashes,{p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in hashes})
            # A legacy project can lack both persisted descriptors and naming
            # sidecars. Export must carry the lazy index across run-ID remaps.
            with restored.connect() as db:
                for row in db.execute('SELECT id,stages FROM runs WHERE project_id=?',(self.project_id,)).fetchall():
                    stages=json.loads(row['stages'])
                    for stage in stages:stage.pop('artifact_records',None)
                    db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),row['id']))
            for path in hashes:Path(path).parent.joinpath('.artifact-publication.json').unlink(missing_ok=True)
            archive=restored.export_project(self.token,self.project_id)
            import zipfile
            with zipfile.ZipFile(archive) as z:
                data=json.loads(z.read('manifest.json'));self.assertIn('artifact_names',data)
                entries=[r for run in data['artifact_names'].values() for r in run['artifacts']]
                self.assertTrue(all(parse_filename(r['public_name']) for r in entries))
            imported=restored.import_project(self.token,archive,str(self.store.root/'Restored names'))
            new_runs=restored.list_runs(self.token,imported['id'])
            self.assertTrue(new_runs)
            for item in new_runs[0]['stages']:
                for r in records(item,new_runs[0]['id']):self.assertIn(r['public_name'],[name for identity,name in saved])

    def test_naming_metadata_does_not_change_scientific_cache_key(self):
        stage=new_stage('fingerprints');stage['parameters']['target']='MeuAlvo'
        before=stage_key(stage,stage['parameters'],{},'test')
        stage['artifact_records']=[dict(public_name='TG[ACHE]-SRC[MOL]-OP[Fingerprints]-Fingerprints.csv',naming={'scheme':1})]
        self.assertEqual(before,stage_key(stage,stage['parameters'],{},'test'))

    def test_publication_sidecar_is_not_scientific_input(self):
        from biomolexplorer.stage_cache import path_manifest
        root=self.store.root/'inputs';root.mkdir();(root/'molecules.csv').write_text('id,smiles\nA,CCO\n')
        before=path_manifest(root)
        (root/'.artifact-publication.json').write_text(json.dumps({'molecules.csv':{'public_name':'new alias'}}))
        self.assertEqual(before,path_manifest(root))
        (root/'molecules.csv').write_text('id,smiles\nA,CCC\n')
        self.assertNotEqual(before,path_manifest(root))
