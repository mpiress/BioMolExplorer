"""Execution quality policy across retrieval, analysis, similarity and graph output."""
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from biomolexplorer.molecule_quality import QualityReport, merge_clean_csv, sanitize_csv
from biomolexplorer.operations import execute_operation
from biomolexplorer.visualizations import SUFFIX, load_view


class MoleculeQualityTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)

    def tearDown(self):self.temp.cleanup()

    def source(self,name,contents):
        path=self.root/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(contents)
        return path

    def test_invalid_smiles_and_ambiguous_identifiers_are_excluded_without_editing_original(self):
        path=self.source('raw.csv','name,smiles,source\n001,CCO,own\nBAD,invalid,own\nEMPTY,,own\nDUP,CCC,own\nDUP,CCN,own\n')
        original=path.read_bytes();report=QualityReport('fingerprints')
        clean=sanitize_csv(path,self.root/'clean.csv',report=report)
        rows=pd.read_csv(clean,dtype={'molecule_chembl_id':str})
        self.assertEqual(rows['molecule_chembl_id'].tolist(),['001'])
        self.assertEqual(rows['source'].tolist(),['own'])
        self.assertEqual(len(report.records),4)
        self.assertEqual(path.read_bytes(),original)

    def test_bad_fingerprints_and_outlier_width_are_excluded(self):
        path=self.source('fp.csv','molecule_chembl_id,fingerprint\nA,"[1,0]"\nB,"[0,1]"\nBAD,nope\nWIDTH,"[1,0,1]"\nBITS,"[1,2]"\n')
        report=QualityReport()
        clean=sanitize_csv(path,self.root/'clean.csv',report=report)
        self.assertEqual(pd.read_csv(clean)['molecule_chembl_id'].tolist(),['A','B'])
        self.assertEqual(len(report.records),3)

    def test_merge_removes_cross_file_conflicts_and_keeps_distinct_codes_for_same_structure(self):
        left=self.source('left.csv','molecule_chembl_id,canonical_smiles\nA,CCO\nB,CCC\nSHARED,CCN\n')
        right=self.source('right.csv','molecule_chembl_id,canonical_smiles\nALIAS,CCO\nSHARED,CO\nBAD,invalid\n')
        report=QualityReport();out=self.root/'merged.csv'
        merge_clean_csv([left,right],out,'compounds',report)
        self.assertEqual(set(pd.read_csv(out)['molecule_chembl_id']),{'A','B','ALIAS'})
        self.assertEqual(len(report.records),3)

    def test_broken_schema_still_fails(self):
        path=self.source('bad.csv','name,wrong_column\nA,CCO\n')
        with self.assertRaisesRegex(ValueError,'colunas'):
            sanitize_csv(path,self.root/'clean.csv','compounds')

    def test_bad_records_do_not_reappear_after_real_fingerprint_similarity_and_graph_stages(self):
        source=self.source('input/compounds.csv','molecule_chembl_id,canonical_smiles\nA,CCO\nB,CCCO\nBAD,invalid\nBLANK,\n')
        original=source.read_bytes()
        with patch.dict(os.environ,{'BIOMOL_CPU_WORKERS':'1'}):
            fp=execute_operation('fingerprints',{'base_input_path':str(source.parent),
                'maccs':False,'pharmacophore':False,'chunk_size':1},self.root/'fingerprints')
            self.assertEqual(fp.details['excluded_records'],2)
            fp_file=next(Path(p) for p in fp.artifacts if p.endswith('.csv'))
            self.assertEqual(set(pd.read_csv(fp_file)['molecule_chembl_id']),{'A','B'})
            sim=execute_operation('similarity',{'base_input_path':str(fp_file.parent),
                'threshold':1,'approximate':False},self.root/'similarity')
            sim_file=next(p for p in sim.artifacts if p.endswith('.csv'))
            graph=execute_operation('graphs',{'graph_inputs':[{'kind':'similarity','file':sim_file,
                'compound_files':[str(fp_file)]}]},self.root/'graphs')
        model=load_view(Path(next(p for p in graph.artifacts if p.endswith(SUFFIX))).read_bytes())
        self.assertEqual({n['id'] for n in model['nodes']},{'A','B'})
        self.assertEqual(source.read_bytes(),original)
        # Reexecution applies the same policy to the original input, independent of caches.
        with patch.dict(os.environ,{'BIOMOL_CPU_WORKERS':'1'}):
            repeated=execute_operation('fingerprints',{'base_input_path':str(source.parent),
                'maccs':False,'pharmacophore':False},self.root/'repeated')
        self.assertEqual(repeated.details['excluded_records'],2)

    def test_invalid_admet_input_is_removed_before_property_calculation(self):
        path=self.source('input/compounds.csv','molecule_chembl_id,canonical_smiles\nA,CCO\nBAD,invalid\n')
        result=execute_operation('admet',{'base_input_path':str(path.parent),'input_file':path.name},self.root/'admet')
        report=json.loads(Path(result.details['exclusion_report']).read_text())
        self.assertEqual(report['excluded_records'],1)
        dataset=pd.read_csv(next(p for p in result.artifacts if Path(p).name=='compounds.csv'))
        self.assertNotIn('BAD',set(dataset['molecule_chembl_id']))

    def test_orphan_edges_and_invalid_smiles_are_removed_together(self):
        compounds=self.source('compounds.csv','molecule_chembl_id,canonical_smiles\nA,CCO\nB,CCCO\nBAD,invalid\n')
        edges=self.source('similarity.csv','source,target,value\nA,B,0.8\nA,BAD,0.9\nA,MISSING,0.7\nA,B,NaN\n')
        result=execute_operation('graphs',{'graph_inputs':[{'kind':'similarity','file':str(edges),
            'compound_files':[str(compounds)]}]},self.root/'graphs')
        model=load_view(Path(next(p for p in result.artifacts if p.endswith(SUFFIX))).read_bytes())
        self.assertEqual({n['id'] for n in model['nodes']},{'A','B'})
        self.assertEqual(len(model['edges']),1)
        self.assertEqual(result.details['excluded_records'],4)
        self.assertEqual(set(pd.read_csv(self.root/'graphs/Molecules/molecules.csv')['molecule_chembl_id']),{'A','B'})

    def test_no_valid_records_produces_an_empty_graph_and_an_explicit_report(self):
        compounds=self.source('compounds.csv','molecule_chembl_id,canonical_smiles\nBAD,invalid\n')
        edges=self.source('similarity.csv','source,target,value\nBAD,MISSING,0.8\n')
        result=execute_operation('graphs',{'graph_inputs':[{'kind':'similarity','file':str(edges),
            'compound_files':[str(compounds)]}]},self.root/'empty')
        model=load_view(Path(next(p for p in result.artifacts if p.endswith(SUFFIX))).read_bytes())
        self.assertEqual(model['nodes'],[]);self.assertEqual(model['edges'],[])
        self.assertEqual(result.details['excluded_records'],2)

    def test_adapter_failures_are_not_suppressed_by_quality_cleanup(self):
        path=self.source('input/compounds.csv','molecule_chembl_id,canonical_smiles\nA,CCO\nBAD,invalid\n')
        with patch('wrappers.molecular_analyzer.generate_fingerprints',side_effect=RuntimeError('engine unavailable')):
            with self.assertRaisesRegex(RuntimeError,'engine unavailable'):
                execute_operation('fingerprints',{'base_input_path':str(path.parent)},self.root/'failed')
        self.assertTrue((self.root/'failed/molecule_exclusions.json').exists())
