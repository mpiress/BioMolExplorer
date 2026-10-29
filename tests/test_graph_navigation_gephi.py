"""Round-trip navigation, complete exports and shared appearance semantics."""
import asyncio
import copy
import csv
import io
import unittest
from types import SimpleNamespace
from unittest.mock import AsyncMock,Mock
from xml.etree import ElementTree as ET

import flet as ft
import networkx as nx
from PIL import Image

from test_network_analysis import AnalyticalViewerTests,model
from biomolexplorer.network_analysis import NetworkAnalysis
from biomolexplorer.graph_layouts import forceatlas2
from biomolexplorer.ui.graph_appearance import continuous_color,MISSING


class GraphNavigationTests(unittest.TestCase):
    viewer=AnalyticalViewerTests.viewer

    def test_large_full_community_full_roundtrip_clears_all_restrictions(self):
        view=self.viewer(nx.disjoint_union(nx.path_graph(90),nx.empty_graph(3000)))
        explorer=view.explorer
        asyncio.run(explorer.prepare())
        group=explorer.analysis.community_cache['full'][0]['0']
        asyncio.run(explorer.open_community(group))
        self.assertLess(len(view.points),90);self.assertFalse(view.raster.visible)
        explorer.query.value='0';explorer.bounds[0].value='99'
        explorer.molecular.linked_ids={'0'};explorer.focus_only.value=True;explorer.ego='0'
        explorer.aggregate_mode=True;explorer.aggregate_choice.value=True
        view.component_id=2;view.component_choice.value='2'
        asyncio.run(explorer.restore_full())
        self.assertEqual(view.mode_choice.value,'full')
        self.assertEqual(set(view.points),set(view.nodes))
        self.assertTrue(view.raster.visible)
        self.assertTrue(view.raster.src.startswith(b'\x89PNG'))
        self.assertFalse(explorer.aggregate_mode);self.assertFalse(explorer.aggregate_choice.value)
        self.assertIsNone(view.component_id);self.assertEqual(explorer.query.value,'')
        self.assertEqual(explorer.community_choice.value,'all');self.assertIsNone(explorer.molecular.linked_ids)
        self.assertEqual((view.raster_layer.width,view.raster_layer.height),(view.WIDTH,view.HEIGHT))
        asyncio.run(explorer.restore_full())
        self.assertEqual(set(view.points),set(view.nodes))

    def test_reset_filters_also_clears_component_and_aggregation(self):
        view=self.viewer(nx.disjoint_union(nx.path_graph(5),nx.path_graph(3)))
        view.component_id=2;view.explorer.aggregate_mode=True;view.explorer.aggregate_choice.value=True
        asyncio.run(view.explorer.reset_filters())
        self.assertEqual(len(view.points),8);self.assertIsNone(view.component_id)
        self.assertFalse(view.explorer.aggregate_mode)

    def test_aggregate_with_more_than_500_communities_still_draws_and_drills_down(self):
        view=self.viewer(nx.empty_graph(3001));explorer=view.explorer
        explorer.aggregate_choice.value=True;asyncio.run(explorer.change())
        self.assertTrue(view.raster.visible)
        self.assertEqual(sum(len(cell) for cell in view.index.cells.values()),3001)
        identifier=next(iter(view.analysis.community_cache['full'][0].values()))
        asyncio.run(explorer.open_community(identifier))
        self.assertEqual(len(view.points),1);self.assertFalse(view.raster.visible)
        asyncio.run(explorer.restore_full());self.assertEqual(len(view.points),3001)

    def test_scope_layout_coordinates_do_not_break_full_paint(self):
        view=self.viewer(nx.disjoint_union(nx.path_graph(8),nx.path_graph(3)))
        view.layout_mode='hierarchical'
        view.explorer.positions=view.analysis.layout('mcc','hierarchical')[0]
        view.redraw()
        self.assertEqual(set(view.points),set(view.nodes))
        asyncio.run(view.change_mode(SimpleNamespace(control=SimpleNamespace(value='mcc'))))
        asyncio.run(view.explorer.restore_full())
        self.assertEqual(set(view.points),set(view.nodes))

    def test_restore_flushes_scene_before_fitting_all_nodes(self):
        view=self.viewer(nx.path_graph(80));view.viewport_ready=True
        view.viewport_width=650;view.viewport_height=480
        calls=Mock();calls.attach_mock(view.ui.page.update,'update')
        view.viewer.reset=AsyncMock();view.viewer.zoom=AsyncMock();view.viewer.pan=AsyncMock()
        calls.attach_mock(view.viewer.reset,'reset')
        view.explorer.community_choice.value='99999'
        asyncio.run(view.explorer.apply_filters());self.assertFalse(view.points)
        calls.reset_mock();asyncio.run(view.explorer.restore_full())
        self.assertLess(next(i for i,c in enumerate(calls.mock_calls) if c[0]=='update'),next(i for i,c in enumerate(calls.mock_calls) if c[0]=='reset'))
        scale=view.viewer.zoom.await_args.args[0];dx,dy=view.viewer.pan.await_args.args
        for x,y in view.points.values():
            self.assertTrue(30<=(x+dx)*scale<=620)
            self.assertTrue(30<=(y+dy)*scale<=450)

    def test_svg_and_png_include_dense_base_and_colorbar_without_selection(self):
        view=self.viewer(nx.complete_graph(160));self.assertFalse(view.drawing.shapes)
        svg=ET.fromstring(view.explorer.svg());ns={'s':'http://www.w3.org/2000/svg'}
        self.assertEqual(len(svg.findall('s:line',ns)),12720)
        self.assertGreaterEqual(len(svg.findall('s:circle',ns)),160)
        self.assertIn(b'maior valor',view.explorer.svg())
        image=Image.open(io.BytesIO(view.explorer.image_png()))
        self.assertGreater(image.height,view.HEIGHT)
        self.assertGreater(len(image.getcolors(image.width*image.height)),20)

    def test_gephi_exports_roundtrip_all_properties_scores_and_visual_style(self):
        view=self.viewer(nx.path_graph(6));original=copy.deepcopy(view.data)
        view.nodes['0']['properties'].update(canonical_smiles='CCO',source='PubChem',provenance={'MW':'provider'})
        view.explorer.appearance.palette.value='green_yellow';view.redraw()
        for kind,reader in [('gexf',nx.read_gexf),('graphml',nx.read_graphml)]:
            graph=reader(io.BytesIO(view.explorer.graph_exchange(kind)))
            self.assertEqual(set(graph),set(view.points));self.assertEqual(graph.number_of_edges(),5)
            self.assertEqual(graph.nodes['0']['canonical_smiles'],'CCO')
            self.assertEqual(graph.nodes['0']['source'],'PubChem')
            self.assertEqual(graph.nodes['0']['color'],view.explorer.appearance.colors['0'])
            self.assertNotIn('weight',graph['0']['1'])
        self.assertEqual(view.data['edges'],original['edges'])

    def test_legends_match_numeric_categorical_missing_and_independent_size(self):
        view=self.viewer(nx.star_graph(5));appearance=view.explorer.appearance
        appearance.palette.value='green_yellow';view.redraw()
        self.assertEqual(appearance.colors['0'],continuous_color(5,5,'green_yellow'))
        self.assertEqual(appearance.colors['1'],continuous_color(1,5,'green_yellow'))
        bar=next(c for c in appearance.legend.controls if isinstance(c,ft.Container) and c.gradient)
        self.assertEqual(bar.gradient.colors[-1],'#FDE725')
        view.explorer.metric_choice.value='closeness';asyncio.run(view.explorer.change())
        appearance.attribute.value='degree';view.redraw()
        self.assertIn('closeness',appearance.size_label);self.assertTrue(appearance.title.startswith('Degree'))
        appearance.attribute.value='source'
        view.nodes['0']['properties']['source']='PubChem'
        for n in list(view.nodes)[1:]:view.nodes[n]['properties']['source']='ChEMBL'
        view.redraw()
        self.assertEqual({s[0] for s in appearance.swatches},{'PubChem','ChEMBL'})
        self.assertNotEqual(appearance.colors['0'],appearance.colors['1'])
        appearance.attribute.value='metric';view.explorer.metric='closeness'
        view.analysis.cache[('full','closeness')]=({},'indisponível');view.redraw()
        self.assertEqual(set(appearance.colors.values()),{MISSING})

    def test_data_laboratory_search_pagination_inspection_and_complete_csv(self):
        view=self.viewer(nx.path_graph(120));table=view.explorer.table
        self.assertEqual(len(table.body.controls[0].controls[0].rows),50)
        table.page(1);self.assertEqual(table.offset,50)
        rows=list(csv.DictReader(io.StringIO(table.csv().decode())))
        self.assertEqual(len(rows),120)
        view.nodes['0']['properties']['source']='PubChem'
        table.query.value='PubChem';table.search()
        self.assertEqual(len(table.rows()),1);self.assertEqual(table.offset,0)
        button=table.body.controls[0].controls[0].rows[0].cells[0].content
        asyncio.run(button.on_click(None));self.assertEqual(view.selected,'0')
        table.kind.value='edges';table.query.value='';table.search()
        self.assertEqual(len(table.rows()),119)
        self.assertEqual(len(list(csv.DictReader(io.StringIO(table.csv().decode())))),119)


class GephiAnalysisTests(unittest.TestCase):
    def test_new_layouts_are_deterministic_finite_and_preserve_every_node(self):
        analysis=NetworkAnalysis(model(nx.disjoint_union(nx.path_graph(20),nx.empty_graph(4))))
        original=list(analysis.graph.edges())
        for kind in ('fa2','fa2_linlog','fr_components','random','radial'):
            positions,_=analysis.layout('full',kind)
            self.assertEqual(set(positions),set(analysis.graph))
            self.assertEqual(len(set(positions.values())),len(positions))
            self.assertEqual(positions,NetworkAnalysis(analysis.data).layout('full',kind)[0])
        self.assertEqual(original,list(analysis.graph.edges()))
        self.assertNotEqual(forceatlas2(nx.path_graph(8)),forceatlas2(nx.path_graph(8),linlog=True))

    def test_force_cost_fallback_keeps_large_component_in_view(self):
        analysis=NetworkAnalysis(model(nx.path_graph(1501)))
        positions,notice=analysis.layout('full','fa2')
        self.assertEqual(len(positions),1501);self.assertIn('limite de custo',notice)
        self.assertEqual(len(set(positions.values())),1501)

    def test_new_metrics_and_lcc_path_statistics_match_known_topology(self):
        graph=nx.disjoint_union(nx.path_graph(4),nx.complete_graph(3))
        analysis=NetworkAnalysis(model(graph))
        self.assertEqual(analysis.metrics('full','core')[0],nx.core_number(analysis.graph))
        self.assertEqual(analysis.metrics('full','triangles')[0],nx.triangles(analysis.graph))
        self.assertEqual(analysis.metrics('full','clustering')[0],nx.clustering(analysis.graph))
        statistics=analysis.statistics('full',True)
        self.assertEqual(statistics['diameter_lcc'],3)
        self.assertAlmostEqual(statistics['mean_path_lcc'],nx.average_shortest_path_length(nx.path_graph(4)))

    def test_communities_are_available_above_previous_display_limit(self):
        analysis=NetworkAnalysis(model(nx.empty_graph(10001)))
        mapping,_,_=analysis.communities('full')
        self.assertEqual(len(mapping),10001)


if __name__=='__main__':unittest.main()
