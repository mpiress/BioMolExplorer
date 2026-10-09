"""Readable complete networks and chemical fragment display across saved runs."""
import asyncio
import json
import math
import unittest
from types import SimpleNamespace

import networkx as nx
import pandas as pd
from rdkit import Chem

from biomolexplorer.visualizations import graph_view, load_view
from caad.graph_results import common_fragment, fragment_molecule

try:
    import flet as ft
    import flet.canvas as canvas
    from biomolexplorer.ui.results_viewer import ResultsViewer
except ImportError:
    ft = None


def fixture():
    graph = nx.disjoint_union_all([nx.star_graph(50), nx.complete_graph(25), nx.path_graph(15)])
    graph = nx.relabel_nodes(graph, lambda value: str(value))
    graph.add_nodes_from(f'isolated-{index}' for index in range(20))
    dataset = pd.DataFrame([{'molecule_chembl_id': node, 'canonical_smiles': 'CCO'} for node in graph])
    return graph, graph_view(graph, dataset, 'Independent molecular components')


class GraphLayoutAndFragmentTests(unittest.TestCase):
    def test_dense_components_and_isolates_have_separated_nodes_and_disjoint_bounds(self):
        graph, model = fixture()
        points = {node['id']: (node['x'], node['y']) for node in model['nodes']}
        self.assertEqual(set(points), set(graph))
        self.assertEqual(len(model['edges']), graph.number_of_edges())
        pairs = list(points.values())
        self.assertGreaterEqual(min(math.dist(a, b) for i, a in enumerate(pairs) for b in pairs[i+1:]), 1-1e-9)
        bounds = []
        for component in nx.connected_components(graph):
            coordinates = [points[node] for node in component]
            bounds.append((min(p[0] for p in coordinates), max(p[0] for p in coordinates),
                           min(p[1] for p in coordinates), max(p[1] for p in coordinates)))
        for i, a in enumerate(bounds):
            for b in bounds[i+1:]:
                self.assertTrue(a[1]+2 <= b[0] or b[1]+2 <= a[0] or a[3]+2 <= b[2] or b[3]+2 <= a[2])

    def test_large_network_keeps_every_vertex_finite_and_distinct(self):
        graph = nx.relabel_nodes(nx.path_graph(1200), str)
        dataset = pd.DataFrame({'molecule_chembl_id': list(graph)})
        model = load_view(json.dumps(graph_view(graph, dataset, 'Large network')).encode())
        points = {(node['x'], node['y']) for node in model['nodes']}
        self.assertEqual(len(points), 1200)
        self.assertTrue(all(math.isfinite(value) for point in points for value in point))
        self.assertEqual(len(model['edges']), 1199)

    def test_fragment_smiles_preserves_atoms_and_bonds_for_aromatic_and_charged_structures(self):
        for structures in (['CCO', 'CCCO'], ['c1ccccc1', 'Cc1ccccc1'],
                           ['c1ccccc1', 'c1ccccn1'], ['[NH3+]CCO', '[NH3+]CCC'], ['[nH]1cccc1']):
            with self.subTest(structures=structures):
                fragment = common_fragment(structures)
                molecule = Chem.MolFromSmiles(fragment['smiles'])
                self.assertIsNotNone(molecule)
                self.assertEqual(molecule.GetNumAtoms(), fragment['atoms'])
                self.assertEqual(molecule.GetNumBonds(), fragment['bonds'])
                self.assertEqual(load_view(json.dumps(dict(fixture()[1], fragment=fragment)).encode())['fragment']['smiles'],
                                 fragment['smiles'])
        self.assertIn('+', common_fragment(['[NH3+]CCO', '[NH3+]CCC'])['smiles'])
        partial = fragment_molecule('[#6]:[#6]:[#6]', 'c1ccccc1')
        self.assertIsNotNone(Chem.MolFromSmiles(Chem.MolToSmiles(partial)))

    def test_complete_fragment_preserves_reference_chirality_and_double_bond_geometry(self):
        for smiles in ('N[C@@H](C)C(=O)O', 'N[C@H](C)C(=O)O', 'F[C@](Cl)(Br)I',
                       'F/C=C/Cl', 'F/C=C\\Cl'):
            with self.subTest(smiles=smiles):
                self.assertEqual(common_fragment([smiles])['smiles'], Chem.MolToSmiles(Chem.MolFromSmiles(smiles)))

    def test_partial_fragment_preserves_stereo_supported_by_retained_subgraph(self):
        for smiles, smarts in (
                ('N[C@@H](C)C(=O)O', '[#7]-[#6](-[#6])-[#6]=[#8]'),
                ('N[C@H](C)C(=O)O', '[#7]-[#6](-[#6])-[#6]=[#8]'),
                ('CC/C=C(/F)Cl', '[#6]-[#6]=[#6](-[#9])-[#17]'),
                ('CC/C=C(\\F)Cl', '[#6]-[#6]=[#6](-[#9])-[#17]')):
            with self.subTest(smiles=smiles):
                reference, query = Chem.MolFromSmiles(smiles), Chem.MolFromSmarts(smarts)
                match = reference.GetSubstructMatch(query)
                bonds = [reference.GetBondBetweenAtoms(match[bond.GetBeginAtomIdx()], match[bond.GetEndAtomIdx()]).GetIdx()
                         for bond in query.GetBonds()]
                expected = Chem.MolFragmentToSmiles(reference, atomsToUse=list(match), bondsToUse=bonds, isomericSmiles=True)
                self.assertEqual(Chem.MolToSmiles(fragment_molecule(smarts, smiles)), expected)
        # The retained alkene has no E/Z geometry after a defining substituent
        # is removed; it must not keep a dangling stereo assignment.
        fragment = fragment_molecule('[#9]-[#6]=[#6]', 'F/C=C/Cl')
        self.assertEqual(Chem.MolToSmiles(fragment), 'C=CF')


@unittest.skipIf(ft is None, 'Install the ui extra to test Flet controls')
class GraphNavigationTests(unittest.TestCase):
    def ui(self):
        async def call(function, *args):
            return function(*args)
        async def guard(action):
            await action()
        return SimpleNamespace(token='owner', call=call, guard=guard,
            store=SimpleNamespace(project=lambda token, project: None),
            page=SimpleNamespace(update=lambda: None))

    def viewer(self, model):
        viewer = ResultsViewer(self.ui(), 'project', model)
        viewer.build()
        self.pans = []
        async def reset():
            return None
        async def pan(x, y):
            self.pans.append((x, y))
        viewer.viewer.reset = reset
        viewer.viewer.pan = pan
        return viewer

    def test_expanding_canvas_preserves_node_spacing_and_component_navigation(self):
        graph, model = fixture()
        viewer = self.viewer(model)
        self.assertTrue(viewer.WIDTH > 720 or viewer.HEIGHT > 460)
        points = list(viewer.points.values())
        self.assertGreaterEqual(min(math.dist(a, b) for i, a in enumerate(points) for b in points[i+1:]), 30-1e-7)
        for component in (2, 3, 4):
            asyncio.run(viewer.change_component(SimpleNamespace(control=SimpleNamespace(value=str(component)))))
            self.assertEqual(set(viewer.points), viewer.components[component])
            lines = [shape for shape in viewer.base_shapes if isinstance(shape, canvas.Line)]
            self.assertEqual(len(lines), graph.subgraph(viewer.components[component]).number_of_edges())
        asyncio.run(viewer.change_component(SimpleNamespace(control=SimpleNamespace(value='all'))))
        self.assertEqual(set(viewer.points), set(graph))

    def test_selecting_node_highlights_relations_and_neighbor_link_centers_target(self):
        _, model = fixture()
        viewer = self.viewer(model)
        asyncio.run(viewer.select('0'))
        highlighted = [shape for shape in viewer.drawing.shapes if isinstance(shape, canvas.Line)
                       and shape.paint.stroke_width == 2]
        self.assertEqual(len(highlighted), 50)
        relations = [control for control in viewer.details.controls if isinstance(control, ft.TextButton) and ' · similaridade ' in str(control.content)]
        self.assertEqual(len(relations), 50)
        asyncio.run(relations[0].on_click(None))
        self.assertEqual(viewer.selected, '1')
        self.assertTrue(self.pans)

    def test_search_locates_compound_outside_current_component(self):
        _, model = fixture()
        viewer = self.viewer(model)
        asyncio.run(viewer.change_component(SimpleNamespace(control=SimpleNamespace(value='2'))))
        self.assertNotIn('0', viewer.points)
        viewer.search.value = '0'
        asyncio.run(viewer.find(None))
        self.assertEqual(viewer.component_id, 1)
        self.assertEqual(viewer.selected, '0')
        self.assertEqual(viewer.component_choice.value, '1')

    def test_saved_smarts_only_result_displays_fragment_smiles_and_regenerates_old_layout(self):
        _, model = fixture()
        model['fragment'] = common_fragment(['CCO', 'CCCO'])
        del model['fragment']['smiles']
        model['layout'] = 'spring'
        for node in model['nodes']:
            node.update(x=0, y=0)
        viewer = self.viewer(model)
        self.assertEqual(len(set(viewer.points.values())), len(model['nodes']))
        fields = [control for control in viewer.fragment_panel.controls if isinstance(control, ft.TextField)]
        self.assertEqual([control.label for control in fields], ['SMILES do fragmento'])
        self.assertEqual(fields[0].value, 'CCO')
        self.assertIsNotNone(Chem.MolFromSmiles(fields[0].value))

    def test_old_fragment_uses_same_reference_order_as_original_image(self):
        graph = nx.Graph()
        graph.add_edge('B', 'A')
        compounds = pd.DataFrame([
            {'molecule_chembl_id': 'B', 'canonical_smiles': 'CCO'},
            {'molecule_chembl_id': 'A', 'canonical_smiles': 'CC[O-]'}])
        model = graph_view(graph, compounds, 'Saved charged compounds')
        model['fragment'] = common_fragment(compounds['canonical_smiles'].tolist())
        original = model['fragment'].pop('smiles')
        self.assertEqual(model['mcc'], ['A', 'B'])
        viewer = self.viewer(model)
        smiles_field=next(control for control in viewer.fragment_panel.controls if isinstance(control,ft.TextField))
        self.assertEqual(smiles_field.value, original)

    def test_singleton_and_flat_component_are_centered_on_unused_axes(self):
        graph = nx.Graph()
        graph.add_node('A')
        compounds = pd.DataFrame([{'molecule_chembl_id': 'A', 'canonical_smiles': 'CCO'}])
        viewer = self.viewer(graph_view(graph, compounds, 'Single compound'))
        self.assertEqual(viewer.points['A'], (viewer.WIDTH/2, viewer.HEIGHT/2))
        graph.add_edge('A', 'B')
        model = graph_view(graph, compounds, 'Horizontal component')
        model['nodes'][0].update(x=0, y=0)
        model['nodes'][1].update(x=2, y=0)
        viewer = self.viewer(model)
        self.assertEqual({point[1] for point in viewer.points.values()}, {viewer.HEIGHT/2})


if __name__ == '__main__':
    unittest.main()
