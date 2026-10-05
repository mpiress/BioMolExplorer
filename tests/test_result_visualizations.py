import asyncio
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import networkx as nx
import pandas as pd

from biomolexplorer.visualizations import SUFFIX, PointIndex, egg_view, graph_view, load_view, molecule_image
from caad.complex_network import GraphAnalysis
from kernel.descriptors import fingerprints, similarityFunctions

try:
    from biomolexplorer.ui.results_viewer import ResultsViewer
except ImportError:
    ResultsViewer = None


class VisualizationTests(unittest.TestCase):
    def setUp(self):
        self.compounds=pd.DataFrame([
            {'molecule_chembl_id':'A','canonical_smiles':'CCO','source':'Imported','TPSA':20.23,'WLOGP':-.001,'BBB':'BBB+'},
            {'molecule_chembl_id':'B','canonical_smiles':'CCC','source':'ChEMBL','TPSA':0.,'WLOGP':1.42,'BBB':'BBB+'},
            {'molecule_chembl_id':'C','canonical_smiles':'O','source':'PubChem','TPSA':31.5,'WLOGP':-1.,'BBB':'BBB-'}])
        self.graph=nx.Graph();self.graph.add_nodes_from(['A','B','C']);self.graph.add_edge('A','B',value=.8)

    def test_full_graph_includes_isolates_and_mcc_preserves_compound_metadata(self):
        model=load_view(json.dumps(graph_view(self.graph,self.compounds,'Graph')).encode())
        self.assertEqual({n['id'] for n in model['nodes']},{'A','B','C'})
        self.assertEqual(model['mcc'],['A','B'])
        self.assertEqual(model['edges'][0]['value'],.8)
        isolated=next(n for n in model['nodes'] if n['id']=='C')
        self.assertEqual(isolated['properties']['degree'],0)
        self.assertEqual(isolated['properties']['source'],'PubChem')
        self.assertEqual(model,graph_view(self.graph,self.compounds,'Graph'))

    def test_graph_worker_exports_full_mcc_and_png_even_without_edges(self):
        empty=pd.DataFrame(columns=['source','target','value'])
        for dataset,edges in ((self.compounds,pd.DataFrame([{'source':'A','target':'B','value':.8}])),
                              (self.compounds,empty),(self.compounds.iloc[:0],empty)):
            with self.subTest(empty=edges.empty),tempfile.TemporaryDirectory() as temp:
                root=Path(temp);source=root/'input';similarity=root/'similarity';output=root/'output'
                source.mkdir();similarity.mkdir()
                dataset.to_csv(source/'compounds.csv',index=False)
                # Use the actual enum prefixes rather than assuming capitalization.
                filename=similarityFunctions.TanimotoSimilarity.value+'_'+fingerprints.Morgan.value+'_compounds'
                edges.to_csv(similarity/(filename+'.csv'),index=False)
                analyzer=GraphAnalysis(str(source)+'/',str(output)+'/')
                analyzer.prepare_graph_analysis(similarity_path=str(similarity))
                model=load_view((output/'plots'/(filename+SUFFIX)).read_bytes())
                self.assertEqual(len(model['nodes']),len(dataset))
                self.assertEqual(model['mcc'],[] if dataset.empty else ['A'] if edges.empty else ['A','B'])
                self.assertTrue((output/'plots'/(filename+'.png')).read_bytes().startswith(b'\x89PNG'))
                selected=pd.read_csv(output/'Molecules'/'molecules.csv')
                self.assertEqual(set(selected['molecule_chembl_id']),set(model['mcc']))

    def test_admet_exports_interactive_egg_png_and_original_metadata(self):
        from wrappers.admet import ADMETWrapper
        for dataset in (self.compounds,self.compounds.iloc[:0]):
            with self.subTest(empty=dataset.empty),tempfile.TemporaryDirectory() as temp:
                source=Path(temp)/'source';output=Path(temp)/'output';source.mkdir()
                dataset.to_csv(source/'compounds.csv',index=False)
                ADMETWrapper(str(output),str(source),input_file='compounds.csv').run_pipeline()
                model=load_view((output/('compounds_egg'+SUFFIX)).read_bytes())
                self.assertEqual(len(model['nodes']),len(dataset))
                if not dataset.empty:
                    self.assertEqual(model['nodes'][0]['properties']['source'],'Imported')
                    self.assertIn('MW',model['nodes'][0]['properties'])
                self.assertTrue((output/'compounds_egg.png').read_bytes().startswith(b'\x89PNG'))

    def test_empty_graph_and_equal_components_have_consistent_mcc(self):
        empty=nx.Graph()
        self.assertEqual(graph_view(empty,self.compounds,'Empty')['nodes'],[])
        edges=pd.DataFrame([{'source':'C','target':'D','value':.7},{'source':'A','target':'B','value':.7}])
        result=GraphAnalysis().max_conected_component(edges)
        self.assertEqual(set(result[0]),{'A','B'})
        self.assertEqual(graph_view(result[3],self.compounds,'Tie')['mcc'],['A','B'])
        self.assertEqual(GraphAnalysis().get_statisticals(empty)['Coeficiente de Cluster Medio'],0)

    def test_egg_preserves_values_and_handles_empty_duplicates_and_missing_values(self):
        model=load_view(json.dumps(egg_view(self.compounds,'EGG')).encode())
        self.assertEqual(len(model['nodes']),3)
        self.assertEqual(model['nodes'][0]['properties']['BBB'],'BBB+')
        self.assertEqual(len(egg_view(pd.concat([self.compounds,self.compounds]),'EGG')['nodes']),3)
        self.assertEqual(egg_view(self.compounds.iloc[:0],'Empty')['nodes'],[])
        invalid=self.compounds.copy();invalid.loc[0,'TPSA']=float('nan')
        self.assertEqual(len(egg_view(invalid,'Missing')['nodes']),2)
        conflict=self.compounds.iloc[:1].copy();conflict.loc[0,'canonical_smiles']='CCC'
        with self.assertRaises(ValueError):egg_view(pd.concat([self.compounds,conflict]),'Conflict')

    def test_rejects_invalid_artifact_and_spatial_hit_testing(self):
        model=graph_view(self.graph,self.compounds,'Graph')
        for modification in ({'version':2},{'kind':'html'},{'nodes':[{'id':'A','x':float('nan'),'y':0,'properties':{}}]},
                             {'edges':[{'source':'A','target':'OUTSIDE'}]},{'mcc':['OUTSIDE']}):
            with self.subTest(modification=modification),self.assertRaises(ValueError):
                load_view(json.dumps(dict(model,**modification)).encode())
        index=PointIndex({'A':(23,23),'B':(48,48)})
        self.assertEqual(index.nearest(24,24),'A')
        self.assertEqual(index.nearest(49,47),'B')
        self.assertIsNone(index.nearest(300,300))
        self.assertEqual(PointIndex({'A':(5,5),'B':(5,5),'C':(50,50)}).nearby(5,5),['A','B'])
        self.assertTrue(molecule_image('CCO').startswith(b'\x89PNG'))
        self.assertIsNone(molecule_image(None))

    def test_malformed_metadata_produces_clear_errors_before_building_controls(self):
        model=graph_view(self.graph,self.compounds,'Graph')
        bad_node=json.loads(json.dumps(model['nodes'][0]));bad_node['properties']['component']=None
        bad_degree=json.loads(json.dumps(model['nodes'][0]));bad_degree['properties']['degree']=-1
        for modification in ({'version':True},{'nodes':[bad_node]},{'nodes':[bad_degree]},
                             {'edges':[{'source':[],'target':'B'}]},
                             {'edges':[{'source':'A','target':'B','value':float('inf')}]},
                             {'mcc':['A','A']},{'kind':'egg'}):
            with self.subTest(modification=modification),self.assertRaises(ValueError):
                load_view(json.dumps(dict(model,**modification)).encode())

    @unittest.skipIf(ResultsViewer is None,'Install the ui extra to test Flet controls')
    def test_hover_click_modes_and_egg_outliers_use_correct_compound(self):
        checks=[]
        async def call(function,*args):return function(*args)
        async def guard(action):await action()
        ui=SimpleNamespace(page=SimpleNamespace(update=lambda:None),token='test-token',call=call,guard=guard,
                           store=SimpleNamespace(project=lambda token,pid:checks.append((token,pid))))
        viewer=ResultsViewer(ui,'project',graph_view(self.graph,self.compounds,'Graph'))
        viewer.build()
        x,y=viewer.points['A'];event=SimpleNamespace(local_position=SimpleNamespace(x=x,y=y))
        viewer.hover(event)
        self.assertEqual(viewer.tip_text.value,'A')
        asyncio.run(viewer.click(event))
        self.assertEqual(viewer.selected,'A')
        self.assertTrue(checks)
        self.assertTrue(all(check==('test-token','project') for check in checks))
        self.assertTrue(any(getattr(c,'value',None)=='A' for c in viewer.details.controls))
        viewer.mode='mcc';viewer.redraw()
        self.assertEqual(set(viewer.points),{'A','B'})
        data=self.compounds.copy();data.loc[0,'TPSA']=350;data.loc[0,'WLOGP']=12
        egg=ResultsViewer(ui,'project',egg_view(data,'EGG'));egg.build()
        for x,y in egg.points.values():
            self.assertTrue(0<x<egg.WIDTH)
            self.assertTrue(0<y<egg.HEIGHT)

    @unittest.skipIf(ResultsViewer is None,'Install the ui extra to test Flet controls')
    def test_coincident_compounds_offer_an_authorized_choice_and_revocation_blocks_details(self):
        from biomolexplorer.workspace import AccessDenied
        dialogs=[]
        authorized=[True]
        def project(token,pid):
            if not authorized[0]:raise AccessDenied('Acesso revogado')
        async def call(function,*args):return function(*args)
        async def guard(action):await action()
        ui=SimpleNamespace(page=SimpleNamespace(update=lambda:None,show_dialog=dialogs.append,pop_dialog=lambda:None),
                           token='test',call=call,guard=guard,store=SimpleNamespace(project=project))
        data=self.compounds.copy();data.loc[1,['TPSA','WLOGP']]=data.loc[0,['TPSA','WLOGP']]
        viewer=ResultsViewer(ui,'project',egg_view(data,'EGG'))
        asyncio.run(viewer.pick('A'))
        choices=dialogs[0].content.content.controls
        self.assertEqual([choice.content for choice in choices],['A','B'])
        asyncio.run(choices[1].on_click(None))
        self.assertEqual(viewer.selected,'B')
        authorized[0]=False
        with self.assertRaises(AccessDenied):asyncio.run(viewer.select('A'))
        self.assertEqual(viewer.selected,'B')

    @unittest.skipIf(ResultsViewer is None,'Install the ui extra to test Flet controls')
    def test_later_selection_wins_when_structure_rendering_finishes_out_of_order(self):
        async def scenario():
            started,released=asyncio.Event(),asyncio.Event()
            async def call(function,*args):
                if function is molecule_image and args==('CCO',):
                    started.set();await released.wait()
                return function(*args)
            ui=SimpleNamespace(page=SimpleNamespace(update=lambda:None),token='test',call=call,
                               store=SimpleNamespace(project=lambda token,pid:None))
            viewer=ResultsViewer(ui,'project',egg_view(self.compounds,'EGG'))
            slow=asyncio.create_task(viewer.select('A'))
            await started.wait()
            await viewer.select('B')
            released.set();await slow
            self.assertEqual(viewer.selected,'B')
        asyncio.run(scenario())

    @unittest.skipIf(ResultsViewer is None,'Install the ui extra to test Flet controls')
    def test_revocation_during_structure_rendering_does_not_publish_details(self):
        from biomolexplorer.workspace import AccessDenied
        async def scenario():
            authorized=[True]
            def project(token,pid):
                if not authorized[0]:raise AccessDenied('Acesso revogado')
            async def call(function,*args):
                result=function(*args)
                if function is molecule_image:authorized[0]=False
                return result
            ui=SimpleNamespace(page=SimpleNamespace(update=lambda:None),token='test',call=call,
                               store=SimpleNamespace(project=project))
            viewer=ResultsViewer(ui,'project',egg_view(self.compounds,'EGG'))
            original=list(viewer.details.controls)
            with self.assertRaises(AccessDenied):await viewer.select('A')
            self.assertIsNone(viewer.selected)
            self.assertEqual(viewer.details.controls,original)
        asyncio.run(scenario())
