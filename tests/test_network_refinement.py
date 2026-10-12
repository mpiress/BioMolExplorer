"""Interaction and scientific correctness for the refined network dashboard."""
import asyncio
import threading
import unittest
from types import SimpleNamespace
from unittest.mock import AsyncMock,patch
import networkx as nx
import flet as ft
from biomolexplorer.network_scaffolds import individual,common
from test_network_analysis import AnalyticalViewerTests


class ScaffoldTests(unittest.TestCase):
    def test_individual_ring_scaffold_and_acyclic_case(self):
        result=individual('Oc1ccccc1')
        self.assertEqual(result['atoms'],6);self.assertTrue(result['scaffold']);self.assertEqual(result['smiles'],'c1ccccc1')
        self.assertIsNone(individual('CCCC')['scaffold'])
        self.assertIsNone(individual(None)['image'])

    def test_explicit_common_set_has_full_coverage_and_highlight(self):
        rows=(('a','Oc1ccccc1'),('b','Nc1ccccc1'))
        result=common(rows,2,50,threading.Event())
        self.assertEqual(result['status'],'complete');self.assertEqual(result['coverage'],1)
        self.assertEqual(result['atoms'],6);self.assertEqual(result['ids'],['a','b'])
        self.assertEqual(individual(rows[0][1],result['smarts'])['highlighted_atoms'],6)

    def test_limit_missing_structure_and_cancellation_never_sample(self):
        rows=[(str(i),'c1ccccc1') for i in range(51)]
        self.assertEqual(common(rows,2,50,threading.Event())['status'],'unavailable')
        self.assertEqual(common([('a',None),('b','CC')],2,50,threading.Event())['status'],'unavailable')
        cancelled=threading.Event();cancelled.set()
        self.assertEqual(common(rows[:2],2,50,cancelled)['status'],'cancelled')
        self.assertFalse(common([('a','CCCC'),('b','c1ccccc1')],2,50,threading.Event())['image'])

    def test_partial_mcs_and_progress_cancellation_are_explicit(self):
        cancelled=threading.Event()
        def partial(molecules,parameters):
            self.assertEqual(parameters.Timeout,1);self.assertEqual(parameters.Threshold,1)
            self.assertTrue(parameters.BondCompareParameters.CompleteRingsOnly)
            self.assertTrue(parameters.ProgressCallback(None,parameters))
            return SimpleNamespace(smartsString='c1ccccc1',numAtoms=6,numBonds=6,canceled=True)
        with patch('rdkit.Chem.rdFMCS.FindMCS',side_effect=partial):
            result=common([('a','Oc1ccccc1'),('b','Nc1ccccc1')],1,50,cancelled)
        self.assertEqual(result['status'],'partial');self.assertEqual(result['coverage'],1)
        self.assertIn('não confirmada',result['notice'])
        def interrupted(molecules,parameters):
            cancelled.set();self.assertFalse(parameters.ProgressCallback(None,parameters))
            return SimpleNamespace(canceled=True)
        with patch('rdkit.Chem.rdFMCS.FindMCS',side_effect=interrupted):
            result=common([('a','Oc1ccccc1'),('b','Nc1ccccc1')],1,50,cancelled)
        self.assertEqual(result['status'],'cancelled');self.assertFalse(result['image'])


class RefinedViewerTests(AnalyticalViewerTests):
    def test_empty_canvas_click_clears_only_selection_and_pending_work(self):
        view=self.viewer(nx.path_graph(4));asyncio.run(view.select('1'))
        points=dict(view.points);filters=view.explorer.query.value;version=view.selection_version
        reset=AsyncMock();view.viewer.reset=reset
        asyncio.run(view.click(SimpleNamespace(local_position=SimpleNamespace(x=-100,y=-100))))
        self.assertIsNone(view.selected);self.assertGreater(view.selection_version,version)
        self.assertEqual(view.points,points);self.assertEqual(view.explorer.query.value,filters)
        self.assertFalse(view.explorer.inspector.scaffold.controls);reset.assert_not_awaited()
        self.assertEqual(view.drawing.shapes,view.base_shapes)

    def test_fit_uses_actual_viewport_and_bounds_with_margin(self):
        view=self.viewer(nx.path_graph(4));view.viewport_ready=True;view.viewport_width=310;view.viewport_height=480
        view.ui.page.height=480/.56;view.explorer.graph_box.height=480
        view.viewer.reset=AsyncMock();view.viewer.zoom=AsyncMock();view.viewer.pan=AsyncMock()
        asyncio.run(view.fit_view());scale=view.viewer.zoom.await_args.args[0];dx,dy=view.viewer.pan.await_args.args
        for x,y in view.points.values():
            self.assertGreaterEqual((x+dx)*scale,30);self.assertLessEqual((x+dx)*scale,280)
            self.assertGreaterEqual((y+dy)*scale,30);self.assertLessEqual((y+dy)*scale,450)
        view.viewer.reset.reset_mock()
        asyncio.run(view.resize_graph(SimpleNamespace(width=310,height=480)));view.viewer.reset.assert_not_awaited()
        asyncio.run(view.resize_graph(SimpleNamespace(width=600,height=480)));view.viewer.reset.assert_awaited_once()

    def test_expanded_view_and_escape_preserve_selection_and_controls(self):
        view=self.viewer(nx.path_graph(4));asyncio.run(view.select('1'))
        asyncio.run(view.explorer.toggle_focus());self.assertFalse(view.explorer.inspector_card.visible)
        self.assertFalse(view.explorer.tabs_bar.visible);self.assertEqual(view.explorer.graph_card.col,12)
        asyncio.run(view.explorer.keyboard_handler(SimpleNamespace(key='Escape')))
        self.assertFalse(view.explorer.focus_mode);self.assertEqual(view.selected,'1')
        asyncio.run(view.explorer.keyboard_handler(SimpleNamespace(key='Escape')));self.assertIsNone(view.selected)

    def test_ranking_query_sort_limit_and_visible_selection_keep_camera(self):
        view=self.viewer(nx.path_graph(15));explorer=view.explorer
        explorer.rank_query.value='1';explorer.rank_order.value='asc';explorer.rank_limit.value='10';explorer.refresh_ranking()
        self.assertEqual(len(explorer.ranking.controls),6)
        first=explorer.ranking.controls[0].content.controls[0].content.controls[1]
        view.viewer.reset=AsyncMock();asyncio.run(first.on_click(None))
        self.assertEqual(view.selected,'14');view.viewer.reset.assert_not_awaited()

    def test_source_search_casefold_pagination_and_identity(self):
        view=self.viewer(nx.path_graph(2));molecular=view.explorer.molecular
        molecular.sources={str(i):dict(artifact_id=str(i),name=f'Docking_ache_{i}.csv',stage_name='Producer',stage_id=f'{i:08}',context={}) for i in range(130)}
        molecular.source.value='129';molecular.source_query.value='DOCKING_';molecular.search_sources()
        self.assertEqual(len(molecular.source.options),52);self.assertEqual(molecular.source.value,'129')
        molecular.page_sources(1);self.assertIn('50',[c.key for c in molecular.source.options])
        molecular.source_query.value='not-found';molecular.search_sources(SimpleNamespace())
        self.assertIn('0 correspondências',molecular.source_count.value);self.assertEqual(molecular.source.value,'129')
        self.assertEqual([o.key for o in molecular.source.options],['graph','129'])

    def test_mcs_reuses_snapshot_and_does_not_recompute_on_selection(self):
        view=self.viewer(nx.path_graph(2))
        view.nodes['0']['properties']['canonical_smiles']='Oc1ccccc1'
        view.nodes['1']['properties']['canonical_smiles']='Nc1ccccc1'
        inspector=view.explorer.inspector;inspector.scope.value='full'
        with patch('biomolexplorer.ui.molecular_inspector.common',wraps=common) as calculate:
            asyncio.run(inspector.calculate());asyncio.run(view.select('0'));asyncio.run(view.select('1'));asyncio.run(inspector.calculate())
            self.assertEqual(calculate.call_count,1);self.assertEqual(inspector.result['coverage'],1)
            self.assertEqual(inspector.structure['highlighted_atoms'],6)
        self.assertFalse(view.ui.page.show_dialog.called)
        asyncio.run(inspector.details());self.assertTrue(view.ui.page.show_dialog.called)

    def test_cancelled_pending_mcs_cannot_replace_previous_result(self):
        view=self.viewer(nx.path_graph(2));inspector=view.explorer.inspector
        original_call=view.ui.call
        async def scenario():
            started=asyncio.Event();release=asyncio.Event()
            async def call(function,*args,**kwargs):
                if function is common:
                    started.set();await release.wait()
                    return dict(status='complete',image='',ids=['0','1'],atoms=1,compounds=2)
                return await original_call(function,*args,**kwargs)
            view.ui.call=call
            task=asyncio.create_task(inspector.calculate());await started.wait()
            self.assertTrue(inspector.compute_button.disabled);inspector.cancel()
            await inspector.calculate() # A second search cannot start before the first exits.
            release.set();await task
            self.assertIsNone(inspector.result);self.assertFalse(inspector.compute_button.disabled)
        asyncio.run(scenario())

    def test_molecular_scatter_is_invalidated_on_universe_change(self):
        view=self.viewer(nx.path_graph(4));molecular=view.explorer.molecular
        molecular.scatter_pairs=[('0',1,10)];molecular.draw_scatter()
        asyncio.run(view.change_mode(SimpleNamespace(control=SimpleNamespace(value='mcc'))))
        molecular.scatter_width=500;molecular.draw_scatter()
        self.assertFalse(molecular.scatter_pairs);self.assertFalse(molecular.points)

    def test_modal_relationships_expand_and_navigation_closes_details(self):
        view=self.viewer(nx.star_graph(80));asyncio.run(view.select('0'));asyncio.run(view.explorer.inspector.details())
        inspector=view.explorer.inspector
        asyncio.run(view.explorer.keyboard_handler(SimpleNamespace(key='Escape')));self.assertEqual(view.selected,'0')
        more=next(c for c in inspector.relations_panel.controls if isinstance(c,ft.TextButton) and str(c.content).startswith('Carregar mais'))
        asyncio.run(more.on_click(None))
        buttons=[c for c in inspector.relations_panel.controls if isinstance(c,ft.TextButton) and ' · escore ' in str(c.content)]
        self.assertEqual(len(buttons),80)
        asyncio.run(buttons[0].on_click(None));self.assertIsNone(inspector.dialog);self.assertEqual(view.selected,'1')

    def test_initial_analysis_has_mcc_data_and_communities_without_manual_controls(self):
        view=self.viewer(nx.disjoint_union(nx.path_graph(3),nx.path_graph(2)))
        for i,node in enumerate(view.data['nodes']):node['properties']['MW']=100+10*i
        asyncio.run(view.initialize());molecular=view.explorer.molecular
        self.assertEqual(molecular.scope.value,'mcc');self.assertEqual(set(molecular.points),{'0','1','2'})
        self.assertEqual(molecular.source.value,'graph');self.assertEqual(molecular.property.value,'MW')
        self.assertTrue(view.analysis.community_cache['full'][0]);self.assertTrue(molecular.initialized)
        self.assertEqual(len(view.points),5) # Molecular scope does not silently filter the graph.
        molecular.lower.value='90';molecular.upper.value='200';asyncio.run(molecular.filter_range())
        self.assertEqual(set(view.points),{'0','1','2'})

    def test_mcs_persisted_result_translates_with_dynamic_counts(self):
        from biomolexplorer.ui.localization import LocalizedPage
        view=self.viewer(nx.path_graph(2));inspector=view.explorer.inspector
        inspector.result=dict(scope='LCC · resultado persistido',compounds=71,atoms=17,status='complete',
            notice='Resultado reutilizado do pipeline; os critérios originais estão preservados.',ring_matches_ring_only=True,complete_rings_only=False,timeout=0,image='')
        inspector.show_result();page=LocalizedPage(SimpleNamespace(),'en');page.localize(inspector.result_panel)
        info=inspector.result_panel.controls[0].controls[1].content
        self.assertEqual(info.controls[0].value,'LCC · persisted result')
        self.assertEqual(info.controls[1].value,'71 compounds · 17 atoms')
        self.assertEqual(info.controls[2].value,'Completed')
        self.assertEqual(inspector.result_panel.controls[2].value,'Ring-to-ring matching: Yes')
        page.set_language('pt');page.localize(inspector.result_panel)
        self.assertEqual(info.controls[0].value,'LCC · resultado persistido')
        self.assertEqual(inspector.result_panel.controls[2].value,'Anéis com anéis: Sim')

    def test_close_requests_mcs_cancellation_and_restores_keyboard(self):
        view=self.viewer(nx.path_graph(2));inspector=view.explorer.inspector
        inspector.cancelled=threading.Event();view.close()
        self.assertTrue(inspector.cancelled.is_set());self.assertIsNone(view.ui.page.on_keyboard_event)
