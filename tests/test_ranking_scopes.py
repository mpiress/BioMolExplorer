"""Numerical reference checks, sparse workloads and independent ranking scopes."""
import asyncio
import csv
import io
import math
import unittest
from types import SimpleNamespace
from unittest.mock import Mock,patch

import networkx as nx

from biomolexplorer.network_analysis import NetworkAnalysis
from test_network_analysis import model,AnalyticalViewerTests


class CentralityFeasibilityTests(unittest.TestCase):
    def assertValues(self,actual,expected):
        self.assertEqual(set(actual),set(expected))
        for node,value in expected.items():self.assertAlmostEqual(actual[node],value,places=10)

    def test_exact_formulas_and_component_normalization_match_networkx(self):
        graphs=[nx.disjoint_union_all([nx.path_graph(8),nx.star_graph(8),nx.complete_graph(5),nx.empty_graph(2)]),
                nx.disjoint_union(nx.cycle_graph(15),nx.path_graph(5)),nx.balanced_tree(3,3)]
        graphs.extend(nx.gnp_random_graph(20,.15,seed=seed) for seed in range(5))
        for graph in graphs:
            analysis=NetworkAnalysis(model(graph));graph=analysis.graph
            for metric,reference in [('betweenness',nx.betweenness_centrality),('closeness',nx.closeness_centrality),
                                     ('clustering',nx.clustering),('triangles',nx.triangles)]:
                with self.subTest(nodes=len(graph),metric=metric):self.assertValues(analysis.metrics('full',metric)[0],reference(graph))

    def test_isolates_and_small_components_do_not_trigger_false_global_cutoffs(self):
        graph=nx.disjoint_union(nx.path_graph(10),nx.empty_graph(6000))
        analysis=NetworkAnalysis(model(graph))
        for metric in ('betweenness','closeness'):
            values,notice=analysis.metrics('full',metric)
            self.assertEqual(len(values),6010);self.assertGreater(values['4'],0)
            self.assertEqual(values['6009'],0);self.assertIn('exata',notice)
        self.assertAlmostEqual(analysis.metrics('full','closeness')[0]['4'],nx.closeness_centrality(analysis.graph,u='4'))

    def test_tree_above_old_limits_has_exact_path_metrics(self):
        analysis=NetworkAnalysis(model(nx.path_graph(5001)))
        between,notice=analysis.metrics('full','betweenness')
        self.assertIn('exata',notice)
        self.assertAlmostEqual(between['2500'],2500*2500*2/(5000*4999))
        close,_=analysis.metrics('full','closeness')
        self.assertAlmostEqual(close['2500'],5000/(2500*2501))

    def test_betweenness_sampling_is_repeatable_and_covers_nontrivial_components(self):
        data=model(nx.disjoint_union(nx.cycle_graph(501),nx.empty_graph(5001)))
        analysis=NetworkAnalysis(data);values,notice=analysis.metrics('full','betweenness')
        self.assertIn('aproximada',notice);self.assertIn('64',notice)
        raw=nx.betweenness_centrality(analysis.graph.subgraph(map(str,range(501))),k=64,normalized=False,seed=42)
        scale=2/(5501*5500)
        for node,value in raw.items():self.assertAlmostEqual(values[node],value*scale)
        self.assertEqual(values,NetworkAnalysis(data).metrics('full','betweenness')[0])

    def test_sparse_solvers_work_above_old_node_and_edge_caps(self):
        dense=NetworkAnalysis(model(nx.complete_graph(318)))
        self.assertGreater(dense.graph.number_of_edges(),50000)
        vector,notice=dense.metrics('full','eigenvector')
        self.assertEqual(len(vector),318);self.assertIn('esparso',notice)
        for value in vector.values():self.assertAlmostEqual(value,1/math.sqrt(318),places=6)
        page,_=dense.metrics('full','pagerank')
        self.assertAlmostEqual(sum(page.values()),1)
        large=NetworkAnalysis(model(nx.disjoint_union(nx.path_graph(8),nx.empty_graph(20000))))
        page,_=large.metrics('full','pagerank')
        self.assertEqual(len(page),20008);self.assertAlmostEqual(sum(page.values()),1)

    def test_high_degree_sparse_hub_has_feasible_clustering_and_triangles(self):
        analysis=NetworkAnalysis(model(nx.star_graph(5000)))
        for metric in ('clustering','triangles'):
            values,_=analysis.metrics('full',metric)
            self.assertEqual(len(values),5001);self.assertEqual(set(values.values()),{0})

    def test_real_cost_and_convergence_failures_are_reported_without_false_values(self):
        for metric,constant in [('betweenness','BRANDES_WORK'),('closeness','PATH_WORK'),('clustering','TRIANGLE_WORK')]:
            analysis=NetworkAnalysis(model(nx.complete_bipartite_graph(8,8)))
            with patch('biomolexplorer.network_centrality.'+constant,1):
                values,notice=analysis.metrics('full',metric)
            self.assertFalse(values);self.assertIn('orçamento',notice)
        analysis=NetworkAnalysis(model(nx.path_graph(8)))
        with patch('networkx.pagerank',side_effect=nx.PowerIterationFailedConvergence(1000)):
            values,notice=analysis.metrics('full','pagerank')
        self.assertFalse(values);self.assertIn('convergiu',notice)


class RankingScopeTests(unittest.TestCase):
    viewer=AnalyticalViewerTests.viewer

    def calculate(self,view,metric):
        view.explorer.rank_metric.value=metric;asyncio.run(view.explorer.calculate_ranking())
        (analysis,mode),_=view.explorer.ranking_context()
        return analysis.metrics(mode,metric)[0]

    def test_component_ranking_uses_visible_topology_and_keeps_color_metric(self):
        view=self.viewer(nx.disjoint_union(nx.path_graph(8),nx.empty_graph(5001)))
        asyncio.run(view.change_component(SimpleNamespace(control=SimpleNamespace(value='1'))))
        values=self.calculate(view,'betweenness')
        self.assertEqual(set(values),set(map(str,range(8))))
        self.assertEqual(values,nx.betweenness_centrality(view.graph.subgraph(map(str,range(8)))))
        self.assertIn('8 nós / 7 arestas',view.explorer.rank_caption.value)
        self.assertEqual(view.explorer.metric,'degree')
        self.assertNotIn(('full','betweenness'),view.analysis.cache)
        self.assertEqual(set(self.calculate(view,'eigenvector')),set(values))

    def test_universe_scope_preserves_global_normalization(self):
        view=self.viewer(nx.disjoint_union(nx.path_graph(5),nx.path_graph(3)))
        asyncio.run(view.change_component(SimpleNamespace(control=SimpleNamespace(value='1'))))
        visible=self.calculate(view,'closeness')
        view.explorer.rank_scope.value='universe';global_values=self.calculate(view,'closeness')
        self.assertEqual(len(global_values),8);self.assertEqual(len(visible),5)
        self.assertGreater(visible['2'],global_values['2'])
        self.assertIn('8 nós / 6 arestas',view.explorer.rank_caption.value)

    def test_score_filter_changes_shortest_paths_and_invalidates_ranking_cache(self):
        graph=nx.cycle_graph(5);nx.set_edge_attributes(graph,.9,'value');graph[0][1]['value']=.1
        view=self.viewer(graph);before=self.calculate(view,'betweenness')
        view.explorer.bounds[4].value='.5';asyncio.run(view.explorer.apply_filters())
        self.assertIn('Calcular ranking',view.explorer.rank_notice.value)
        after=self.calculate(view,'betweenness')
        expected=nx.Graph(view.explorer.graph())
        self.assertEqual(after,nx.betweenness_centrality(expected));self.assertNotEqual(before,after)

    def test_export_ranking_identifies_selected_scope_and_contains_only_its_nodes(self):
        view=self.viewer(nx.disjoint_union(nx.path_graph(5),nx.empty_graph(3)))
        asyncio.run(view.change_component(SimpleNamespace(control=SimpleNamespace(value='1'))))
        self.calculate(view,'closeness');saved=[]
        async def save(**kwargs):saved.append(kwargs)
        view.ui.picker.save_file=save;asyncio.run(view.explorer.export_ranking())
        rows=list(csv.DictReader(io.StringIO(saved[0]['src_bytes'].decode())))
        self.assertEqual(len(rows),5);self.assertEqual({r['universe'] for r in rows},{'visible:full'})
        self.assertEqual({r['metric'] for r in rows},{'closeness'})

    def test_late_calculation_does_not_overwrite_a_new_scope(self):
        async def scenario():
            view=self.viewer(nx.path_graph(8));explorer=view.explorer
            entered=asyncio.Event();release=asyncio.Event()
            async def call(function,*args,**kwargs):
                if getattr(function,'__name__','')=='metrics':entered.set();await release.wait()
                return function(*args,**kwargs)
            view.ui.call=call;explorer.rank_metric.value='closeness'
            pending=asyncio.create_task(explorer.calculate_ranking());await entered.wait()
            explorer.query.value='0';await explorer.apply_filters()
            release.set();await pending
            self.assertIn('1 nós / 0 arestas',explorer.rank_caption.value)
            self.assertIn('Calcular ranking',explorer.rank_notice.value)
        asyncio.run(scenario())


if __name__=='__main__':unittest.main()
