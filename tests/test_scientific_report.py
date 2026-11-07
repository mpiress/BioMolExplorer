"""Scientific report statistics use persisted values without overclaiming."""
import base64
import io
import unittest

from biomolexplorer.scientific_report import _analyze_table,pipeline_diagram,relevant_parameters
from biomolexplorer.catalog import new_stage
from biomolexplorer.experimental_report import report_content,report_pdf


def metric(values,name):return next(item['value'] for item in values if item['name']==name)


class ScientificReportTests(unittest.TestCase):
    def test_pipeline_overview_uses_readable_block_cards(self):
        first,second=new_stage('retrieve_compounds'),new_stage('fingerprints')
        first['block_id']='B01';second['block_id']='B02'
        second['bindings']['base_input_path']={'stage':first['id']}
        figure=pipeline_diagram([first,second])
        from PIL import Image
        image=Image.open(io.BytesIO(base64.b64decode(figure['data'])))
        self.assertEqual(figure['mime'],'image/png')
        self.assertGreater(image.width,image.height)
        self.assertIn('Canvas snapshot',figure['caption'])

    def test_pipeline_connection_map_disambiguates_duplicate_block_ids(self):
        first,second=new_stage('retrieve_compounds'),new_stage('fingerprints')
        first['block_id']=second['block_id']='B17'
        second['bindings']['base_input_path']={'stage':first['id'],'output_port':'compounds'}
        base={'schema':'biomolexplorer.experimental-report','schema_version':1,'generated_at':'now',
            'application':{'version':'test'},'project':{'id':'p','name':'Study','description':'','revision':0},
            'pipeline':[first,second],'executions':[],'artifacts':[],'analyses':[],'results':{'figures':[]},
            'pipeline_figure':pipeline_diagram([first,second]),'notes':{'en':'note'}}
        _,sections=report_content(base,'en');rendered=str(sections)
        self.assertIn('Connection map',rendered);self.assertIn('C01',rendered)
        self.assertIn('B01 ·',rendered);self.assertIn('B02 ·',rendered)

    def test_retrieval_keeps_incompatible_activity_measures_separate(self):
        fields=['molecule_chembl_id','standard_type','standard_value','standard_units']
        rows=[dict(molecule_chembl_id='A',standard_type='IC50',standard_value='10',standard_units='nM'),
              dict(molecule_chembl_id='B',standard_type='IC50',standard_value='20',standard_units='nM'),
              dict(molecule_chembl_id='A',standard_type='Ki',standard_value='2',standard_units='uM')]
        metrics,tables,_=_analyze_table('retrieve_compounds',fields,rows,{})
        self.assertEqual(metric(metrics,'retrieved_records'),3);self.assertEqual(metric(metrics,'unique_compounds'),2)
        self.assertEqual({table['title'] for table in tables},{'Activity IC50','Activity Ki'})

    def test_pdb_admet_and_fingerprint_statistics_are_exact(self):
        metrics,tables,_=_analyze_table('retrieve_structures',['PDB_CODE','RESOLUTION','LIGAND'],[
            {'PDB_CODE':'4M0E','RESOLUTION':'1.8','LIGAND':'FAD'},{'PDB_CODE':'4M0E','RESOLUTION':'1.8','LIGAND':'FAD'}],{})
        self.assertEqual(metric(metrics,'retrieved_structures'),1);self.assertEqual(metric(metrics,'resolution')['unit'],'Å')
        metrics,_,_=_analyze_table('admet',['compound','ADMET_selected','MW'],[
            {'compound':'A','ADMET_selected':'true','MW':'100'},{'compound':'B','ADMET_selected':'false','MW':'200'}],{})
        self.assertEqual(metric(metrics,'accepted_compounds'),1);self.assertEqual(metric(metrics,'retention'),50)
        metrics,_,_=_analyze_table('fingerprints',['compound','fingerprint'],[
            {'compound':'A','fingerprint':'[1, 0, 1, 0]'},{'compound':'B','fingerprint':'invalid'}],{})
        self.assertEqual(metric(metrics,'valid_fingerprints'),1);self.assertEqual(metric(metrics,'invalid_fingerprints'),1)

    def test_similarity_and_graph_metrics_preserve_real_denominators(self):
        rows=[{'source':'A','target':'B','value':'0.8'},{'source':'B','target':'C','value':'0.6'}]
        metrics,_,_=_analyze_table('similarity',['source','target','value'],rows,{})
        self.assertEqual(metric(metrics,'molecules'),3);self.assertAlmostEqual(metric(metrics,'similarity_coefficients')['mean'],.7)
        metrics,_,_=_analyze_table('graphs',['source','target','value'],rows,{})
        self.assertEqual(metric(metrics,'connected_components'),1);self.assertEqual(metric(metrics,'largest_component'),3)
        self.assertAlmostEqual(metric(metrics,'network_density'),2/3)

    def test_redocking_docking_consensus_and_amber_do_not_invent_results(self):
        rows=[{'PDB_CODE':'A','LIGAND':'L1','RMSD':'1.0'},{'PDB_CODE':'B','LIGAND':'L2','RMSD':'3.0'}]
        metrics,_,_=_analyze_table('redocking',list(rows[0]),rows,{'rmsd_threshold':2})
        success=next(item for item in metrics if item['name']=='success_rate')
        self.assertEqual(success['value'],50);self.assertEqual(success['denominator'],2);self.assertIn('2.0 Å',success['criterion'])
        docking=[{'receptor':'R','ligand':'A','affinity':'-8.0'},{'receptor':'R','ligand':'B','affinity':'-6.0'}]
        for operation in ('docking_vina','docking_dock6','consensus','dynamics_dock6'):
            metrics,tables,_=_analyze_table(operation,list(docking[0]),docking,{})
            self.assertEqual(metric(metrics,'reported_results'),2);self.assertEqual(tables[0]['rows'][0]['ligand'],'A')
        empty,_,charts=_analyze_table('docking_vina',['affinity'],[],{})
        self.assertEqual(metric(empty,'reported_results'),0);self.assertEqual(charts,[])

    def test_only_nondefault_or_methodologically_essential_parameters_are_shown(self):
        stage=new_stage('fingerprints');stage['parameters']['radius']=3;stage['parameters']['chunk_size']=stage['parameters'].get('chunk_size',500)
        relevant=relevant_parameters(stage)
        self.assertEqual(relevant['radius'],3)
        if 'chunk_size' in stage['parameters']:self.assertNotIn('chunk_size',relevant)

    def test_configured_and_partial_pipelines_never_receive_expected_results(self):
        stage=new_stage('docking_vina');stage['contract']={'outputs':{'docking_results':{}}}
        base={'schema':'biomolexplorer.experimental-report','schema_version':1,'generated_at':'now',
            'application':{'version':'test'},'project':{'id':'p','name':'Study','description':'','revision':0},
            'pipeline':[stage],'executions':[],'artifacts':[],'analyses':[],'notes':{'en':'note'}}
        _,sections=report_content(base,'en');rendered=str(sections)
        self.assertIn('configured, not executed',rendered);self.assertNotIn('Reported results',rendered)
        partial=dict(base,executions=[{'id':'run','status':'failed','created':0,'stages':[{'id':stage['id'],'status':'failed'}]}],
            analyses=[{'stage_id':stage['id'],'name':'Docking','operation':'docking_vina','status':'failed','method':'Docking method.',
                'parameters':{},'metrics':[],'tables':[],'charts':[],'warnings':['tool failed']}])
        _,sections=report_content(partial,'en');rendered=str(sections)
        self.assertIn('tool failed',rendered);self.assertNotIn('reported_results',rendered)

    def test_target_fields_do_not_break_metrics_and_long_parameters_fit_pdf(self):
        stage=new_stage('graphs');stage['parameters']['target']='ACHE'
        report={'schema':'biomolexplorer.experimental-report','schema_version':1,'generated_at':'now',
            'application':{'version':'test'},'project':{'id':'p','name':'Large study','description':'','revision':0},
            'pipeline':[stage],'executions':[{'id':'run','status':'succeeded','created':0,
                'stages':[{'id':stage['id'],'status':'succeeded'}]}],'artifacts':[],'results':{'figures':[]},
            'analyses':[{'stage_id':stage['id'],'name':'Graph','operation':'graphs','status':'succeeded',
                'method':'Network analysis.','parameters':{'graph_inputs':[{'name':'dataset-'+str(i),'files':['x'*80]} for i in range(200)]},
                'metrics':[{'name':'molecules','value':10,'kind':'count'}],'tables':[],'charts':[],'warnings':[]}],
            'notes':{'en':'note'},'pipeline_figure':None}
        _,sections=report_content(report,'en')
        self.assertIn('Molecular targets: ACHE',str(sections))
        self.assertTrue(report_pdf(report,'en').startswith(b'%PDF-'))


if __name__=='__main__':unittest.main()
