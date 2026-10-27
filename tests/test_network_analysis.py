"""Scientific scopes, non-destructive exploration and bounded analytical work."""
import asyncio
import copy
import csv
import io
import unittest
from types import SimpleNamespace
from unittest.mock import Mock,patch
import networkx as nx
import flet as ft
from biomolexplorer.network_analysis import NetworkAnalysis
from biomolexplorer.ui.results_viewer import ResultsViewer


def model(graph):
    graph=nx.relabel_nodes(graph,str)
    return dict(version=1,kind='graph',title='Saved relationship network',
        nodes=[dict(id=n,x=i,y=0,properties={}) for i,n in enumerate(graph)],
        edges=[dict(source=a,target=b,**attrs) for a,b,attrs in graph.edges(data=True)],
        mcc=sorted(max(nx.connected_components(graph),key=len,default=set())))


class NetworkAnalysisTests(unittest.TestCase):
    def test_universes_metrics_and_isolates_are_correct(self):
        data=model(nx.disjoint_union_all([nx.path_graph(3),nx.empty_graph(2)]));original=copy.deepcopy(data)
        analysis=NetworkAnalysis(data)
        full=analysis.statistics('full',True);lcc=analysis.statistics('mcc',True)
        self.assertEqual((full['nodes'],full['edges'],full['components'],full['isolates']),(5,2,3,2))
        self.assertEqual((lcc['nodes'],lcc['edges'],lcc['components'],lcc['isolates']),(3,2,1,0))
        self.assertEqual(analysis.metrics('full','density')[0]['1'],.5)
        self.assertEqual(analysis.metrics('mcc','density')[0]['1'],1)
        self.assertAlmostEqual(analysis.metrics('full','closeness')[0]['1'],nx.closeness_centrality(analysis.graph)['1'])
        self.assertFalse(analysis.metrics('full','eigenvector')[0])
        self.assertTrue(analysis.metrics('mcc','eigenvector')[0])
        self.assertEqual(data,original)

    def test_communities_and_aggregation_account_for_all_edges(self):
        graph=nx.disjoint_union(nx.complete_graph(6),nx.complete_graph(6));graph.add_edge(0,6)
        analysis=NetworkAnalysis(model(graph));mapping,modularity,_=analysis.communities('full')
        self.assertEqual(len(set(mapping.values())),2);self.assertGreater(modularity,.4)
        aggregate,_=analysis.aggregate('full')
        self.assertEqual(sum(a['size'] for _,a in aggregate.nodes(data=True)),12)
        self.assertEqual(sum(a['internal_edges'] for _,a in aggregate.nodes(data=True))+sum(a['count'] for _,_,a in aggregate.edges(data=True)),31)

    def test_layouts_preserve_topology_and_cache_coordinates(self):
        analysis=NetworkAnalysis(model(nx.disjoint_union(nx.path_graph(6),nx.empty_graph(1))))
        before=list(analysis.graph.edges(data=True))
        for layout in ('circular','force','kamada','hierarchical','community'):
            positions,_=analysis.layout('full',layout)
            self.assertEqual(set(positions),set(analysis.graph))
            self.assertIs(positions,analysis.layout('full',layout)[0])
        self.assertEqual(before,list(analysis.graph.edges(data=True)))

    def test_combined_filters_scores_and_ego_preserve_original(self):
        graph=nx.path_graph(6);nx.set_edge_attributes(graph,.9,'value');graph[2][3]['value']=.1
        data=model(graph);analysis=NetworkAnalysis(data);original=copy.deepcopy(data)
        view=analysis.filtered('full',degree=(2,None),score=(.5,None),ego='2',hops=2)
        self.assertEqual(set(view),{'1','2','3','4'})
        self.assertFalse(view.has_edge('2','3'))
        self.assertEqual(data,original)
        self.assertEqual(analysis.metrics('full','betweenness')[0],nx.betweenness_centrality(analysis.graph,weight=None))

    def test_scores_are_not_assumed_to_be_distance_weights(self):
        graph=nx.path_graph(4);nx.set_edge_attributes(graph,-5,'value')
        analysis=NetworkAnalysis(model(graph))
        self.assertEqual(analysis.metrics('full','closeness')[0],nx.closeness_centrality(analysis.graph))
        self.assertAlmostEqual(sum(analysis.metrics('full','pagerank')[0].values()),1)
        rows=list(csv.DictReader(io.StringIO(analysis.export_csv('full','closeness').decode())))
        self.assertEqual(len(rows),4);self.assertTrue(all(r['universe']=='full' for r in rows))

    def test_cost_limits_and_empty_legacy_results_are_explicit(self):
        analysis=NetworkAnalysis(model(nx.path_graph(1501)))
        self.assertEqual(len(analysis.metrics('full','closeness')[0]),1501)
        self.assertIn('exata',analysis.metrics('full','closeness')[1])
        self.assertIsNone(analysis.layout('full','force')[0])
        self.assertIn('exata',analysis.metrics('full','betweenness')[1])
        empty=NetworkAnalysis(model(nx.Graph()))
        self.assertEqual(empty.communities('full')[0],{})
        self.assertEqual(empty.statistics('full')['nodes'],0)


class AnalyticalViewerTests(unittest.TestCase):
    def test_neighbor_inspector_loads_relationships_progressively(self):
        view=self.viewer(nx.star_graph(80));asyncio.run(view.select('0'))
        relations=lambda:[c for c in view.details.controls if isinstance(c,ft.TextButton) and ' · escore ' in str(c.content)]
        self.assertEqual(len(relations()),50)
        more=next(c for c in view.details.controls if isinstance(c,ft.TextButton) and str(c.content).startswith('Carregar mais'))
        asyncio.run(more.on_click(None));self.assertEqual(len(relations()),80);self.assertTrue(more.disabled)

    def test_dense_network_retains_search_selection_and_neighbor_navigation(self):
        with patch('biomolexplorer.ui.results_viewer.separated_graph_layout',return_value={str(i):(i,0) for i in range(3102)}):
            view=self.viewer(nx.star_graph(3101))
        view.search.value='0';asyncio.run(view.find(None))
        self.assertEqual(len(view.points),3102);self.assertFalse(view.explorer.focus_only.value)
        self.assertTrue(view.raster.visible)
        asyncio.run(view.follow_relation('1'));self.assertEqual(view.selected,'1')
        self.assertEqual(len(view.points),3102)

    def test_search_preserves_full_network_above_vector_render_budget(self):
        with patch('biomolexplorer.ui.results_viewer.separated_graph_layout',return_value={str(i):(i,0) for i in range(3001)}):
            view=self.viewer(nx.path_graph(3001))
        view.search.value='1500';asyncio.run(view.find(None))
        self.assertEqual(view.selected,'1500')
        self.assertEqual(len(view.points),3001)
        self.assertIsNone(view.explorer.ego)
        self.assertTrue(view.raster.visible)
    def test_long_hierarchy_has_bounded_canvas_without_omitting_nodes(self):
        with patch('biomolexplorer.ui.results_viewer.separated_graph_layout',return_value={str(i):(i,0) for i in range(1000)}):
            view=self.viewer(nx.path_graph(1000))
        view.explorer.positions=view.analysis.layout('full','hierarchical')[0]
        view.layout_mode='hierarchical';view.redraw()
        self.assertLessEqual(view.WIDTH,4096);self.assertLessEqual(view.HEIGHT,4096)
        self.assertEqual(len(view.points),1000)
        self.assertIn('4.096',view.explorer.geometry_notice)
    def viewer(self,graph):
        async def call(function,*args,**kwargs):return function(*args,**kwargs)
        async def guard(action):return await action()
        async def noop(*args,**kwargs):pass
        self.ui=SimpleNamespace(token='token',current={'id':'project'},page=SimpleNamespace(width=1440,height=900,update=Mock(),show_dialog=Mock()),
            store=SimpleNamespace(project=lambda *args:None),call=call,guard=guard,notify=Mock(),picker=SimpleNamespace(save_file=noop))
        view=ResultsViewer(self.ui,'project',model(graph));view.build()
        view.viewer.reset=noop;view.viewer.pan=noop
        return view

    def test_metric_scope_mode_filter_and_export(self):
        view=self.viewer(nx.disjoint_union(nx.path_graph(3),nx.empty_graph(1)))
        view.explorer.metric_choice.value='closeness';asyncio.run(view.explorer.change())
        self.assertEqual(set(view.points),set(view.nodes))
        asyncio.run(view.change_mode(SimpleNamespace(control=SimpleNamespace(value='mcc'))))
        self.assertEqual(set(view.points),{'0','1','2'})
        self.assertAlmostEqual(view.analysis.metrics('mcc','closeness')[0]['1'],1)
        view.explorer.bounds[0].value='2';asyncio.run(view.explorer.apply_filters())
        self.assertEqual(set(view.points),{'1'})
        self.assertEqual(len(view.data['nodes']),4)
        self.assertIn(b'<svg',view.explorer.svg())

    def test_communities_ego_aggregation_and_reset(self):
        view=self.viewer(nx.path_graph(6));view.explorer.color_communities.value=True
        asyncio.run(view.explorer.change());self.assertFalse(view.explorer.groups.controls)
        self.assertTrue(view.explorer.analysis.community_cache['full'][0])
        view.selected='2';asyncio.run(view.explorer.ego_selected())
        self.assertEqual(set(view.points),{'1','2','3'})
        view.explorer.aggregate_choice.value=True;asyncio.run(view.explorer.change())
        self.assertFalse(view.points);self.assertIn('agregada',view.summary.value)
        self.assertTrue(view.drawing.shapes)
        view.explorer.aggregate_choice.value=False;asyncio.run(view.explorer.change())
        asyncio.run(view.explorer.reset_filters());self.assertEqual(len(view.points),6)

    def test_large_network_opens_full_molecular_view_without_sampling(self):
        with patch('biomolexplorer.ui.results_viewer.separated_graph_layout',return_value={str(i):(i,0) for i in range(3001)}):
            view=self.viewer(nx.path_graph(3001))
        self.assertEqual(len(view.points),3001)
        self.assertEqual(view.mode,'full')
        self.assertEqual(view.mode_choice.value,'full')
        self.assertFalse(view.component_choice.disabled)
        self.assertFalse(view.explorer.aggregate_mode)
        self.assertTrue(view.raster.visible)
        self.assertTrue(view.raster.src.startswith(b'\x89PNG'))
        self.assertLess(len(view.drawing.shapes),100)
        self.assertEqual(len(view.data['nodes']),3001)
        self.assertIn('clique',view.explorer.notice.value)
        view.selected='1500';asyncio.run(view.explorer.ego_selected())
        self.assertEqual(set(view.points),{'1499','1500','1501'})

    def test_large_disconnected_network_keeps_full_and_mcc_available(self):
        graph=nx.disjoint_union(nx.path_graph(4),nx.empty_graph(3000))
        with patch('biomolexplorer.ui.results_viewer.separated_graph_layout',return_value={str(i):(i,0) for i in range(len(graph))}):
            view=self.viewer(graph)
        self.assertEqual(view.mode,'full')
        self.assertFalse(view.explorer.aggregate_mode)
        self.assertEqual(len(view.points),3004)
        asyncio.run(view.change_mode(SimpleNamespace(control=SimpleNamespace(value='mcc'))))
        self.assertEqual(set(view.points),{'0','1','2','3'})
        self.assertFalse(view.raster.visible)
        self.assertEqual(len(view.analysis.graph),3004)

