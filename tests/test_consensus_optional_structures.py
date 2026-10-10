"""Consensus scores remain usable without receptor, reference or pose files."""
import csv
import json
import unittest
from pathlib import Path
from unittest.mock import patch
import test_pipeline_execution as execution
from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import PipelineService
from biomolexplorer.operations import execute_operation
from biomolexplorer.docking_data import write_csv
from biomolexplorer.docking_results import DockingResults


class ConsensusOptionalStructureTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown

    def inputs(self,stale=False):
        root=self.store.project_dir(self.project_id)
        results={};stage=new_stage('consensus');producers=[]
        for engine in ('vina','dock6'):
            producer=new_stage('docking_'+engine);producers.append(producer)
            folder=root/engine;folder.mkdir();rows=[]
            for i in (1,2):
                row=dict(molecule_chembl_id='M'+str(i),score=-i*(10 if engine=='dock6' else 1))
                if stale:row.update(canonical_smiles='CCO',engine=engine,receptor_id='1ABC_A',conformer_file='missing/pose.mol2',
                    receptor_file='receptors/1ABC_A.noH.pdb',reference_file='missing/reference.pdb',footprint_file='missing/footprint.pdf')
                rows.append(row)
            write_csv(folder/'docking_results.csv',rows)
            from biomolexplorer.input_validation import validate_file
            validate_file(folder/'docking_results.csv',engine)
            results[producer['id']]=[str(folder/'docking_results.csv')]
            stage['bindings']['base_'+engine+'_path']={'stage':producer['id'],'selector':'auto'}
        return stage,results,producers

    def run_consensus(self,stale):
        stage,results,producers=self.inputs(stale)
        service=PipelineService(self.store)
        try:
            params=service._resolve(self.project_id,self.store.user(self.token)['id'],stage,results,pipeline=producers+[stage])
            result=execute_operation('consensus',params,str(self.store.project_dir(self.project_id)/'out'))
            table=next(Path(p) for p in result.artifacts if Path(p).name=='MeuAlvo.csv')
            with table.open() as stream:rows=list(csv.DictReader(stream))
            self.assertEqual([r['molecule_chembl_id'] for r in rows],['M2','M1'])
            self.assertEqual([float(r['dock6']) for r in rows],[-20,-10])
            self.assertEqual([float(r['normalized_score']) for r in rows],[1,0])
            self.assertTrue(all(not r['conformer_file'] for r in rows))
        finally:service.close()

    def test_codes_and_scores_alone_reach_worker_and_produce_consensus(self):self.run_consensus(False)
    def test_stale_structural_metadata_does_not_block_consensus(self):self.run_consensus(True)

    def test_direct_score_tables_need_no_structural_metadata(self):
        from wrappers.docking import generate_consensus
        stage,results,producers=self.inputs(False)
        root=self.store.project_dir(self.project_id)
        frame=generate_consensus(str(root),str(root/'direct'),'scores',
            base_vina_path=str(root/'vina'),base_dock6_path=str(root/'dock6'))
        self.assertEqual(frame['molecule_chembl_id'].tolist(),['M2','M1'])
        self.assertEqual(frame['dock6'].tolist(),[-20,-10])

    def test_legacy_dock6_pdf_reaches_consensus_without_a_pose(self):
        stage,results,producers=self.inputs(False)
        root=self.store.project_dir(self.project_id)
        pdf=root/'dock6/footprint/plots/M1.pdf';pdf.parent.mkdir(parents=True)
        pdf.write_bytes(b'%PDF-1.4\nfixture')
        table=root/'dock6/docking_results.csv'
        write_csv(table,[dict(molecule_chembl_id='M1',score=-10,footprint_file='footprint/plots/M1.pdf',footprint_origin='docked_pose')])
        # The legacy producer manifest contains only the table, not the PDF.
        service=PipelineService(self.store)
        try:
            params=service._resolve(self.project_id,self.store.user(self.token)['id'],stage,results,pipeline=producers+[stage])
            result=execute_operation('consensus',params,str(root/'out'))
            output=next(Path(p) for p in result.artifacts if Path(p).name=='MeuAlvo.csv')
            with output.open() as stream:row=next(csv.DictReader(stream))
            attached=output.parent/row['footprint_file']
            self.assertEqual(attached.read_bytes(),pdf.read_bytes())
            self.assertIn(str(attached),result.artifacts)
            self.assertFalse(row['conformer_file'])
        finally:service.close()

    def test_optional_context_is_resolved_at_its_source_and_raw_score_is_preserved(self):
        from biomolexplorer.docking_data import score_input_records,consensus_rows
        root=self.store.project_dir(self.project_id);folder=root/'source';folder.mkdir()
        receptor=folder/'receptor.pdb';receptor.write_text('receptor')
        pose=folder/'pose.mol2';pose.write_text('Internal_energy_repulsive: 10\n')
        table=folder/'docking_results.csv'
        write_csv(table,[dict(molecule_chembl_id='M1',canonical_smiles='CCO',score=-20,conformer_file=pose.name,receptor_file=receptor.name)])
        records=score_input_records([table],[table,pose,receptor],'dock6')
        self.assertEqual(records[0]['receptor_file'],str(receptor))
        other=dict(records[0],engine='vina',score=-7)
        self.assertEqual(consensus_rows([other],records)[0]['dock6'],-20)

if __name__=='__main__':unittest.main()
