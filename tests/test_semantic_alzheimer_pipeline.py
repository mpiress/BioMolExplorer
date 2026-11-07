"""Architecture-level Alzheimer DAG coverage with small coherent fixtures.

These tests validate references, manifests and persistence.  They deliberately
do not claim scientific validation of external docking or MD engines.
"""
import tempfile
import unittest
from pathlib import Path

from biomolexplorer.artifact_registry import SemanticArtifactResolver,register_run
from biomolexplorer.catalog import new_stage
from biomolexplorer.flow import connect
from biomolexplorer.pipeline import validate_pipeline
from biomolexplorer.selection_views import csv_view,path_matches
from biomolexplorer.workspace import WorkspaceStore


class SemanticAlzheimerPipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.store=WorkspaceStore(self.tmp.name)
        self.token=self.store.register('Owner','alzheimer-semantic@example.org','owner-password')
        self.project=self.store.create_project(self.token,'Alzheimer semantic integration')

    def tearDown(self):self.tmp.cleanup()

    @staticmethod
    def link(stages,source,target,field):
        connect(stages,source['id'],target['id'],field)

    def test_branched_analysis_and_molecular_dags_persist_semantic_references(self):
        stages=[];analyses={}
        for name,target in (('AChE','CHEMBL220'),('BChE','CHEMBL1914'),('BACE1','CHEMBL4822')):
            retrieval=new_stage('retrieve_compounds');retrieval.update(name=name)
            retrieval['parameters']['search_term']=target
            admet=new_stage('admet');fingerprints=new_stage('fingerprints')
            similarity=new_stage('similarity');graphs=new_stage('graphs')
            branch=[retrieval,admet,fingerprints,similarity,graphs];stages.extend(branch)
            self.link(stages,retrieval,admet,'base_input_path')
            self.link(stages,admet,fingerprints,'base_input_path')
            self.link(stages,fingerprints,similarity,'base_input_path')
            self.link(stages,similarity,graphs,'similarity_path')
            analyses[name]=branch

        pdb=new_stage('retrieve_structures');pdb['name']='PDB AChE'
        prepare=new_stage('prepare_structures');redocking=new_stage('redocking')
        vina=new_stage('docking_vina');dock6=new_stage('docking_dock6');consensus=new_stage('consensus')
        amber=new_stage('dynamics_dock6')
        stages.extend([pdb,prepare,redocking,vina,dock6,consensus,amber])
        self.link(stages,pdb,prepare,'base_input_path')
        prepare['bindings']['base_input_path']['selection']={'pdb_id':'4M0E','chain_id':'A'}
        self.link(stages,analyses['AChE'][1],prepare,'base_selected_mols')
        self.link(stages,pdb,redocking,'base_input_path')
        self.link(stages,prepare,vina,'base_input_path');self.link(stages,prepare,vina,'base_selected_mols')
        self.link(stages,prepare,dock6,'base_input_path');self.link(stages,prepare,dock6,'base_selected_mols')
        self.link(stages,vina,consensus,'base_vina_path');self.link(stages,dock6,consensus,'base_dock6_path')
        self.link(stages,prepare,amber,'base_input_path');self.link(stages,consensus,amber,'base_selected_mols')

        self.assertEqual(len(validate_pipeline(stages)),len(stages))
        for stage in stages:
            for binding in stage.get('bindings',{}).values():
                refs=binding.get('sources',[]) if 'sources' in binding else [binding]
                self.assertTrue(all('output_port' in ref for ref in refs if 'stage' in ref))
        saved=self.store.save_pipeline(self.token,self.project['id'],stages,self.project['revision'])
        reopened=self.store.project(self.token,self.project['id'])
        self.assertEqual(reopened['pipeline'],saved['pipeline'])
        persisted=next(s for s in reopened['pipeline'] if s['id']==prepare['id'])
        self.assertEqual(persisted['bindings']['base_input_path']['selection'],{'pdb_id':'4M0E','chain_id':'A'})

    def test_manifest_selection_resolves_only_4m0e_and_multiple_compounds(self):
        producer=new_stage('retrieve_structures');consumer=new_stage('prepare_structures')
        root=self.store.project_dir(self.project['id'])/'fixtures';root.mkdir()
        atom='ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n'
        chosen=root/'opaque-a.pdb';other=root/'opaque-b.pdb'
        chosen.write_text('HEADER'.ljust(62)+'4M0E\n'+atom);other.write_text('HEADER'.ljust(62)+'1ABC\n'+atom)
        records=[dict(path=str(chosen),output_port='receptor_structures',artifact_id='a'*64,occurrence_id='c'*64,
                      artifact_type='structures',stage_id=producer['id'],run_id='e'*32),
                 dict(path=str(other),output_port='receptor_structures',artifact_id='b'*64,occurrence_id='d'*64,
                      artifact_type='structures',stage_id=producer['id'],run_id='e'*32)]
        producer['artifact_records']=records
        with self.store.connect() as db:register_run(db,self.project['id'],'e'*32,[producer])
        resolver=SemanticArtifactResolver(self.store,self.project['id'],
            [dict(producer,_execution_artifact_records=records)],{producer['id']:[str(chosen),str(other)]})
        selected=resolver.resolve({'stage':producer['id'],'output_port':'receptor_structures',
                                   'selection':{'pdb_id':'4M0E','chain_id':'A'}})
        # Opaque storage names demonstrate that the registry entity, not a
        # display filename, drives selection.
        self.assertEqual(selected,[str(chosen)])
        named=root/'4M0E.pdb';named.write_text(atom)
        self.assertTrue(path_matches(named,{'pdb_id':'4M0E','chain_id':'A'}))

        compounds=root/'compounds.csv'
        compounds.write_text('molecule_chembl_id,canonical_smiles\nCHEMBL25,CCO\nCHEMBL120,CCN\nCHEMBL999,CCC\n')
        view=csv_view(compounds,{'entity_ids':['CHEMBL25','CHEMBL120']},root/'views')
        text=view.read_text();self.assertIn('CHEMBL25',text);self.assertIn('CHEMBL120',text);self.assertNotIn('CHEMBL999',text)


if __name__=='__main__':unittest.main()
