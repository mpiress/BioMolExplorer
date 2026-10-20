"""Explicit source joins, molecular statistics and linked visual selection."""
import copy
import json
import unittest
import asyncio
from pathlib import Path
from unittest.mock import patch
from biomolexplorer.molecular_network import enrich,numeric_properties,histogram,association
from biomolexplorer.result_files import ResultFiles
from biomolexplorer.workspace import AccessDenied
from test_network_analysis import AnalyticalViewerTests,model
from test_graph_explorer import GraphProvenanceTests
import networkx as nx


class MolecularDataTests(unittest.TestCase):
    def test_cached_source_keeps_physical_producing_run_and_survives_root_relocation(self):
        from biomolexplorer.artifact_sources import records
        stage={'id':'a'*32,'artifacts':['/old/projects/p/runs/'+('1'*32)+'/'+('a'*32)+'/data.csv'],
               'reused_from_run':'2'*32,'artifact_records':None}
        original=next(records(stage,'3'*32))
        self.assertEqual(original['run_id'],'1'*32)
        stage['artifacts']=[stage['artifacts'][0].replace('/old/','/migrated/')]
        self.assertEqual(next(records(stage,'3'*32))['artifact_id'],original['artifact_id'])
    def test_nested_chembl_properties_from_real_csv_representation_are_available(self):
        nodes=[{'id':'CHEMBL1','properties':{'molecule_properties':"{'full_mwt': '399.54', 'alogp': '5.14', 'full_molformula': 'C26H29N3O'}"}}]
        properties,values=numeric_properties(nodes)
        self.assertEqual(values['molecule_properties.full_mwt']['CHEMBL1'],399.54)
        self.assertEqual(values['molecule_properties.alogp']['CHEMBL1'],5.14)
        self.assertEqual(properties['CHEMBL1']['molecule_properties.full_molformula'],'C26H29N3O')
    def test_join_rejects_structure_conflicts_and_preserves_graph_values(self):
        data=model(nx.path_graph(3));data['nodes'][0]['properties']={'canonical_smiles':'CCO','MW':46.07}
        data['nodes'][1]['properties']={'canonical_smiles':'CCN'}
        original=copy.deepcopy(data)
        source={'artifact_id':'a'*64,'stage_name':'ADMET A','run_id':'run','context':{'target':'A'}}
        result=enrich(data,'molecule_chembl_id,canonical_smiles,MW,TPSA\n0,OCC,99,20\n1,CCC,44,0\n2,CCO,46,20\n2,CCN,45,22\n',source)
        self.assertEqual(result['matched'],1);self.assertEqual(set(result['conflicts']),{'1','2'})
        properties,numeric=numeric_properties(data['nodes'],result)
        self.assertEqual(properties['0']['MW'],46.07);self.assertEqual(numeric['TPSA']['0'],20)
        self.assertEqual(data,original)

    def test_histograms_and_associations_do_not_impute_missing_values(self):
        x={'A':1,'B':2,'C':3,'D':4};y={'A':2,'B':4,'C':6}
        pairs,r=association(x,y);self.assertEqual(len(pairs),3);self.assertAlmostEqual(r,1)
        self.assertEqual(sum(count for _,_,count in histogram(y.values())),3)
        self.assertIsNone(association({'A':1},{'A':2})[1])
        _,fields=numeric_properties([{'id':'A','properties':{'MW':float('nan'),'LogP':'missing','TPSA':'20.5'}}])
        self.assertEqual(fields,{'TPSA':{'A':20.5}})


class MolecularDashboardTests(AnalyticalViewerTests):
    def test_lazy_analysis_and_linked_interval_preserve_topology(self):
        view=self.viewer(nx.path_graph(5));original=copy.deepcopy(view.data)
        for i,node in enumerate(view.data['nodes']):node['properties']['MW']=100+i*10
        molecular=view.explorer.molecular;molecular.refresh_options()
        self.assertFalse(molecular.drawing.shapes)
        asyncio.run(molecular.analyze());self.assertEqual(len(molecular.points),5)
        molecular.lower.value='120';molecular.upper.value='130'
        asyncio.run(molecular.filter_range());self.assertEqual(set(view.points),{'2','3'})
        self.assertEqual(len(view.graph),5);self.assertEqual(view.graph.number_of_edges(),4)
        asyncio.run(view.explorer.reset_filters());self.assertEqual(len(view.points),5)
        view.close();self.assertFalse(view.drawing.shapes);self.assertFalse(molecular.drawing.shapes)


class MolecularSourceAuthorizationTests(unittest.TestCase):
    setUp=GraphProvenanceTests.setUp
    tearDown=GraphProvenanceTests.tearDown
    resolve=GraphProvenanceTests.resolve
    fixture_run=GraphProvenanceTests.fixture_run

    def test_initialization_selects_exact_molecular_csv_and_excludes_edge_csv(self):
        from biomolexplorer.ui.results_viewer import ResultsViewer
        stage=self.fixture_run();service=ResultFiles(self.store)
        filename=service.graph_datasets(self.token,self.pid,self.rid,self.graph['id'])[0]['path']
        view=ResultsViewer(GraphProvenanceTests.ui(self),self.pid,service.graph_model(self.token,self.pid,self.rid,self.graph['id'],filename))
        view.build();asyncio.run(view.initialize())
        molecular=view.explorer.molecular;selected=molecular.sources[molecular.source.value]
        self.assertIn('Molecules',Path(selected['path']).parts)
        self.assertEqual(Path(selected['path']).stem,Path(filename).name.removesuffix('.biomol-view.json'))
        self.assertGreater(view.molecular_enrichment['matched'],0)
        self.assertEqual(molecular.scope.value,'mcc')

    def test_persisted_source_identity_and_authorized_join(self):
        stage=self.fixture_run();service=ResultFiles(self.store)
        source=Path(stage['artifacts'][0]).parent/'molecular_metadata.csv'
        source.write_text('molecule_chembl_id,MW,TPSA\nA,46.07,20.23\nB,60.1,20.23\n')
        stage['artifacts'].append(str(source))
        with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps([stage]),self.rid))
        selected=next(s for s in service.graph_sources(self.token,self.pid,self.rid) if s['path']==str(source))
        filename=service.graph_datasets(self.token,self.pid,self.rid,self.graph['id'])[0]['path']
        before=Path(filename).read_bytes()
        result=service.graph_enrichment(self.token,self.pid,self.rid,self.graph['id'],filename,selected['artifact_id'])
        self.assertEqual(result['source']['run_id'],self.rid);self.assertEqual(result['matched'],2)
        self.assertEqual(len(result['source']['content_sha256']),64)
        self.assertEqual(Path(filename).read_bytes(),before)
        with self.assertRaises(AccessDenied):service.graph_enrichment(self.token,self.pid,self.rid,self.graph['id'],filename,'f'*64)
        digest=service.graph_model(self.token,self.pid,self.rid,self.graph['id'],filename)['dashboard_context']['sha256']
        changed=json.loads(before);changed['title']='Updated experiment'
        Path(filename).write_text(json.dumps(changed))
        with self.assertRaisesRegex(ValueError,'mudou'):
            service.graph_enrichment(self.token,self.pid,self.rid,self.graph['id'],filename,selected['artifact_id'],digest)

    def test_oversized_graph_does_not_block_other_saved_analyses(self):
        stage=self.fixture_run();service=ResultFiles(self.store)
        from biomolexplorer.visualizations import MAX_GRAPH_BYTES,SUFFIX
        path=Path(stage['artifacts'][0]).parent/('large'+SUFFIX)
        with path.open('wb') as file:file.truncate(MAX_GRAPH_BYTES+1)
        stage['artifacts'].insert(0,str(path))
        with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps([stage]),self.rid))
        choices=service.graph_datasets(self.token,self.pid,self.rid,self.graph['id'])
        self.assertEqual(len(choices),3);self.assertTrue(choices[0]['unavailable'])
        data=service.graph_model(self.token,self.pid,self.rid,self.graph['id'],choices[1]['path'])
        self.assertEqual(data['kind'],'graph')

    def test_large_fingerprint_source_streams_only_matching_molecular_metadata(self):
        import hashlib
        stage=self.fixture_run();service=ResultFiles(self.store)
        source=Path(stage['artifacts'][0]).parent/'fingerprint_metadata.csv'
        source.write_text('molecule_chembl_id,canonical_smiles,MW,fingerprint\nA,CCO,46.07,"'+('0,'*2048)+'0"\nB,CCCO,60.1,"0,1"\n')
        stage['artifacts'].append(str(source))
        with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps([stage]),self.rid))
        selected=next(s for s in service.graph_sources(self.token,self.pid,self.rid) if s['path']==str(source))
        filename=service.graph_datasets(self.token,self.pid,self.rid,self.graph['id'])[0]['path']
        with patch('biomolexplorer.result_files.MAX_VIEW_BYTES',64):
            result=service.graph_enrichment(self.token,self.pid,self.rid,self.graph['id'],filename,selected['artifact_id'])
        self.assertEqual(result['matched'],2)
        self.assertNotIn('fingerprint',result['matches']['A']['properties'])
        self.assertEqual(result['source']['content_sha256'],hashlib.sha256(source.read_bytes()).hexdigest())

