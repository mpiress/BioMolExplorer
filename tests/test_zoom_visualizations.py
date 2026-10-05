"""Regression coverage for shrinking scenes and preserving readable graph nodes."""
import asyncio
import math
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import networkx as nx
import pandas as pd

from biomolexplorer.visualizations import graph_view, egg_view
from caad.graph_results import report_png

try:
    import flet as ft
    import flet.canvas as canvas
    from biomolexplorer.ui.flow_canvas import FlowCanvas
    from biomolexplorer.ui.molecule_3d import Molecule3D
    from biomolexplorer.ui.results_viewer import ResultsViewer
    from biomolexplorer.ui.zoom import MIN_SCALE, MAX_SCALE, zoomable_view
except ImportError:
    ft = None


@unittest.skipIf(ft is None, 'Install the ui extra to test Flet controls')
class ZoomVisualizationTests(unittest.TestCase):
    def setUp(self):
        self.ui = SimpleNamespace(token='test', selected=None,
            current={'pipeline': []}, page=SimpleNamespace(update=lambda: None))
        self.graph = nx.star_graph(50)
        self.graph = nx.relabel_nodes(self.graph, str)
        self.graph.add_node('isolated')
        self.compounds = pd.DataFrame([
            {'molecule_chembl_id': node, 'canonical_smiles': 'CCO',
             'TPSA': 20, 'WLOGP': 1, 'BBB': 'BBB+'} for node in self.graph])
        self.model = graph_view(self.graph, self.compounds, 'Hub and leaves')

    def test_pipeline_graph_egg_and_image_can_shrink_below_the_viewport(self):
        viewers = [FlowCanvas(self.ui, False).viewer,
            ResultsViewer(self.ui, 'project', self.model).viewer,
            ResultsViewer(self.ui, 'project', egg_view(self.compounds, 'EGG')).viewer,
            zoomable_view(ft.Image(src='image.png'), expand=True)]
        for viewer in viewers:
            with self.subTest(viewer=viewer.content.__class__.__name__):
                self.assertLessEqual(viewer.min_scale, .001)
                self.assertGreaterEqual(viewer.max_scale, 100)
                # A small min_scale alone does not fix Flutter's boundary clamp.
                for side in ('left', 'top', 'right', 'bottom'):
                    self.assertEqual(getattr(viewer.boundary_margin, side), math.inf)

    def test_hub_leaf_and_isolate_have_same_small_size_in_full_and_mcc_views(self):
        viewer = ResultsViewer(self.ui, 'project', self.model)
        self.assertEqual(viewer.nodes['0']['properties']['degree'], 50)
        self.assertEqual(viewer.nodes['1']['properties']['degree'], 1)
        self.assertEqual(viewer.nodes['isolated']['properties']['degree'], 0)
        for mode in ('full', 'mcc'):
            viewer.mode = mode
            viewer.redraw()
            circles = [s for s in viewer.base_shapes if isinstance(s, canvas.Circle)
                       and s.paint.style != ft.PaintingStyle.STROKE]
            self.assertEqual(len(circles), 52 if mode == 'full' else 51)
            self.assertEqual(len({circle.radius for circle in circles}), 1)
            self.assertLessEqual(circles[0].radius, 4)
            self.assertEqual(len([s for s in viewer.base_shapes if isinstance(s, canvas.Line)]), 50)
            self.assertNotEqual(circles[0].paint.color, circles[1].paint.color)
            for identifier in ('0', '1'):
                x, y = viewer.points[identifier]
                self.assertEqual(viewer.index.nearest(x, y), identifier)

    def test_exported_mcc_also_uses_small_uniform_markers(self):
        with patch('networkx.draw_networkx_nodes', wraps=nx.draw_networkx_nodes) as draw:
            image = report_png(self.model)
        self.assertTrue(image.startswith(b'\x89PNG'))
        self.assertLessEqual(draw.call_args.kwargs['node_size'], 20)
        self.assertEqual(len(draw.call_args.kwargs['node_color']), 51)

    def test_molecule_3d_shrinks_atoms_bonds_and_coordinates_together(self):
        model = {'atoms': [{'element': 'C', 'x': -1, 'y': 0, 'z': 0},
                           {'element': 'O', 'x': 1, 'y': 0, 'z': 0}],
                 'bonds': [{'a': 0, 'b': 1, 'order': 1}]}
        view = Molecule3D(self.ui.page, model)
        original = [s for s in view.drawing.shapes if isinstance(s, canvas.Circle)]
        distance = math.dist((original[0].x, original[0].y), (original[1].x, original[1].y))
        for scale in (.1, .01, MIN_SCALE, MAX_SCALE):
            with self.subTest(scale=scale):
                view.zoom.value = math.log10(scale)
                view.change_zoom(SimpleNamespace(control=view.zoom))
                circles = [s for s in view.drawing.shapes if isinstance(s, canvas.Circle)]
                line = next(s for s in view.drawing.shapes if isinstance(s, canvas.Line))
                self.assertAlmostEqual(view.scale, scale)
                self.assertAlmostEqual(circles[0].radius, original[0].radius * scale)
                self.assertAlmostEqual(math.dist((circles[0].x, circles[0].y),
                    (circles[1].x, circles[1].y)), distance * scale)
                self.assertAlmostEqual(line.paint.stroke_width, 3 * scale)
        view.reset(None)
        self.assertEqual(view.scale, 1)
        self.assertEqual(view.zoom.value, 0)
        self.assertEqual(view.zoom_label.value, '100%')

    def test_static_image_preview_uses_the_same_zoom_and_reset_controls(self):
        from biomolexplorer.ui.app import WorkspaceUI
        dialogs = []
        async def call(function, *args): return function(*args)
        ui = WorkspaceUI(SimpleNamespace(width=1200, height=900, show_dialog=dialogs.append),
                         SimpleNamespace(read_file=lambda *args: b'image'), None)
        ui.token, ui.current, ui.call = 'test', {'id': 'project'}, call
        asyncio.run(ui.preview_artifact('project', 'plots/mcc.png'))
        content = dialogs[0].content.content
        viewer = content.controls[-1]
        self.assertEqual(viewer.boundary_margin.left, math.inf)
        self.assertEqual(viewer.min_scale, MIN_SCALE)
        toolbar = content.controls[1]
        self.assertEqual([button.label for button in toolbar.controls[:2]], ['Ampliar', 'Reduzir'])
        self.assertEqual(toolbar.controls[-1].content, 'Recentrar')


if __name__ == '__main__':
    unittest.main()
