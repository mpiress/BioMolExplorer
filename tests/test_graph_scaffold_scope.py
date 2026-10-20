"""Graph scope, per-compound scaffolds and responsive explicit MCS actions."""
import asyncio
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import test_mcs_ranking_refinements as fixtures
from biomolexplorer.network_scaffolds import common
from biomolexplorer.molecular_network import numeric_properties


class GraphScaffoldScopeTests(unittest.IsolatedAsyncioTestCase):
    def viewer(self):
        return fixtures.MCSConfigurationTests().viewer()

    async def test_graph_mode_recalculates_exact_largest_component_and_full_set(self):
        view=self.viewer();inspector=view.explorer.inspector
        inspector.scope.value='full';await inspector.calculate()
        self.assertEqual(len(inspector.result['ids']),4)
        await view.change_mode(SimpleNamespace(control=SimpleNamespace(value='mcc')))
        self.assertEqual(inspector.scope.value,'mcc')
        self.assertEqual(inspector.result['ids'],['0','1','2'])
        self.assertEqual(inspector.result['atoms'],3)
        self.assertIn('maior componente conectado',inspector.context.value)
        self.assertFalse(view.fragment_panel.visible)
        await view.change_mode(SimpleNamespace(control=SimpleNamespace(value='full')))
        self.assertEqual(inspector.result['ids'],['0','1','2','3'])
        self.assertEqual(inspector.result['atoms'],2)

    async def test_node_click_changes_individual_scaffold_without_changing_group_result(self):
        view=self.viewer();inspector=view.explorer.inspector
        view.nodes['0']['properties']['canonical_smiles']='c1ccccc1O'
        view.nodes['1']['properties']['canonical_smiles']='c1ccncc1'
        await inspector.calculate();result=inspector.result
        await view.select('0',show_popup=False);first=inspector.structure['smiles']
        self.assertIn('0',inspector.scaffold.controls[0].value)
        await view.select('1',show_popup=False)
        self.assertIn('1',inspector.scaffold.controls[0].value)
        self.assertNotEqual(first,inspector.structure['smiles'])
        self.assertIs(inspector.result,result)
        self.assertEqual(inspector.result['ids'],['0','1','2'])
        self.assertIn('apenas deste composto',inspector.scaffold.controls[1].value)

    async def test_filtered_scope_tracks_actual_filters_and_invalidates_empty_groups(self):
        view=self.viewer();explorer=view.explorer;inspector=explorer.inspector
        explorer.molecular.linked_ids={'0','1'}
        inspector.scope.value='filtered';await inspector.change_scope()
        self.assertEqual(inspector.result['ids'],['0','1'])
        explorer.molecular.linked_ids={'2','3'};await explorer.apply_filters()
        self.assertEqual(inspector.result['ids'],['2','3'])
        self.assertEqual(inspector.result['atoms'],2)
        explorer.molecular.linked_ids={'2'};await explorer.apply_filters()
        self.assertIsNone(inspector.result)
        self.assertIn('1 compostos',inspector.notice.value)
        self.assertFalse(inspector.result_panel.controls)

    async def test_calculate_click_shows_progress_during_preparation_and_ignores_duplicates(self):
        view=self.viewer();inspector=view.explorer.inspector;original=view.ui.call
        entered,release=asyncio.Event(),asyncio.Event();calls=[]
        async def call(function,*args,**kwargs):
            if function is numeric_properties:
                entered.set();await release.wait()
            if function is common:calls.append(args)
            return await original(function,*args,**kwargs)
        view.ui.call=call
        first=asyncio.create_task(inspector.compute_button.on_click(None));await entered.wait()
        self.assertTrue(inspector.progress.visible)
        self.assertTrue(inspector.mcs_tile.expanded)
        self.assertTrue(inspector.compute_button.disabled)
        self.assertFalse(inspector.cancel_button.disabled)
        await inspector.compute_button.on_click(None)
        release.set();await first
        self.assertEqual(len(calls),1)
        self.assertEqual(inspector.result['ids'],['0','1','2'])
        self.assertFalse(inspector.progress.visible)
        self.assertFalse(inspector.compute_button.disabled)

    async def test_failed_search_explains_failure_and_can_be_retried(self):
        inspector=self.viewer().explorer.inspector
        with patch('biomolexplorer.ui.molecular_inspector.common',side_effect=RuntimeError('search failed')):
            with self.assertRaisesRegex(RuntimeError,'search failed'):
                await inspector.compute_button.on_click(None)
        self.assertIn('search failed',inspector.notice.value)
        self.assertIsNone(inspector.result)
        self.assertFalse(inspector.progress.visible)
        self.assertFalse(inspector.compute_button.disabled)
        await inspector.compute_button.on_click(None)
        self.assertEqual(inspector.result['status'],'complete')

    async def test_over_limit_click_explains_configured_limit_in_open_panel(self):
        inspector=self.viewer().explorer.inspector
        inspector.scope.value='full';inspector.limit.value='3'
        with patch('biomolexplorer.ui.molecular_inspector.common') as calculate:
            await inspector.compute_button.on_click(None)
        self.assertFalse(calculate.called)
        self.assertTrue(inspector.mcs_tile.expanded)
        self.assertIn('4 compostos; limite configurado: 3',inspector.notice.value)
        self.assertFalse(inspector.progress.visible)

    async def test_changing_molecular_source_removes_old_result_and_highlights(self):
        view=self.viewer();inspector=view.explorer.inspector
        await inspector.calculate();await view.select('0',show_popup=False)
        self.assertGreater(inspector.structure['highlighted_atoms'],0)
        await view.explorer.molecular.load_source()
        self.assertIsNone(inspector.result)
        self.assertEqual(inspector.structure['highlighted_atoms'],0)

    async def test_english_scope_labels_and_feedback_preserve_compound_ids(self):
        from biomolexplorer.ui.localization import Translator
        translate=Translator('en')
        self.assertEqual(translate('MCC / LCC · maior componente conectado\nConjunto do MCS: 3 compostos'),
            'MCC / LCC · largest connected component\nMCS set: 3 compounds')
        self.assertEqual(translate('Bemis–Murcko · CHEMBL220'),'Bemis–Murcko · CHEMBL220')

    async def test_changing_community_recalculates_its_explicit_set(self):
        view=self.viewer();explorer=view.explorer;inspector=explorer.inspector
        view.analysis.community_cache['full']=({'0':1,'1':1,'2':2,'3':2},None,'')
        explorer.community_choice.value='1';inspector.scope.value='community'
        await inspector.change_scope()
        self.assertEqual(inspector.result['ids'],['0','1'])
        explorer.community_choice.value='2';await explorer.apply_filters()
        self.assertEqual(inspector.result['ids'],['2','3'])
        self.assertIn('Comunidade 2',inspector.result['scope'])

    async def test_mode_switch_cancels_pending_result_and_displays_only_latest_set(self):
        view=self.viewer();inspector=view.explorer.inspector;original=view.ui.call
        entered,release=asyncio.Event(),asyncio.Event();calls=[]
        async def call(function,*args,**kwargs):
            if function is common:
                calls.append(args)
                if len(calls)==1:
                    entered.set();await release.wait()
            return await original(function,*args,**kwargs)
        view.ui.call=call;inspector.scope.value='full'
        first=asyncio.create_task(inspector.calculate());await entered.wait()
        await view.change_mode(SimpleNamespace(control=SimpleNamespace(value='mcc')))
        self.assertTrue(calls[0][3].is_set())
        self.assertIsNone(inspector.result)
        release.set();await first
        self.assertEqual(inspector.result['ids'],['0','1','2'])
        self.assertEqual(len(calls),2)
        self.assertFalse(inspector.progress.visible)

    async def test_selection_outside_mcs_set_is_explicitly_identified(self):
        view=self.viewer();inspector=view.explorer.inspector
        await inspector.calculate();await view.select('3',show_popup=False)
        self.assertTrue(any('fora do conjunto' in str(getattr(c,'value','')) for c in inspector.scaffold.controls))
        self.assertEqual(inspector.result['ids'],['0','1','2'])

    async def test_exact_persisted_mcc_fragment_is_reused_above_interactive_limit(self):
        from biomolexplorer.ui.results_viewer import ResultsViewer
        view=self.viewer();await view.explorer.inspector.calculate()
        data=dict(view.data,fragment=dict(view.explorer.inspector.result))
        restored=ResultsViewer(view.ui,'project',data);restored.build()
        inspector=restored.explorer.inspector;inspector.limit.value='2';restored.mode='mcc'
        with patch('biomolexplorer.ui.molecular_inspector.common') as search:
            await inspector.follow_mode()
        self.assertFalse(search.called)
        self.assertEqual(inspector.result['ids'],['0','1','2'])
        self.assertIn('persistido',inspector.result['scope'])
        self.assertIn('reutilizado do pipeline',inspector.notice.value)

    async def test_dropdown_round_trip_restores_saved_scaffold_above_search_limit(self):
        import networkx as nx
        import threading
        import test_network_analysis as analytical
        from biomolexplorer.ui.results_viewer import ResultsViewer
        base=analytical.AnalyticalViewerTests().viewer(nx.path_graph(51))
        for node in base.data['nodes']:
            node['properties']['canonical_smiles']='CCCO'
        saved=common(tuple((n['id'],'CCCO') for n in base.data['nodes']),1,200,threading.Event())
        view=ResultsViewer(base.ui,'project',dict(base.data,fragment=saved));view.build()
        inspector=view.explorer.inspector
        original_image=inspector.result['image']
        self.assertTrue(original_image)
        for scope in ('full','community','filtered','mcc','full','mcc'):
            inspector.scope.value=scope
            await inspector.scope.on_select(None)
            if scope=='mcc':
                self.assertIsNotNone(inspector.result)
                self.assertEqual(inspector.result['image'],original_image)
                self.assertEqual(set(inspector.result['ids']),set(view.analysis.lcc))
                self.assertTrue(inspector.result_panel.controls)
                self.assertFalse(inspector.compute_button.disabled)

    async def test_restoration_rejects_a_persisted_fragment_after_structure_changes(self):
        from biomolexplorer.ui.results_viewer import ResultsViewer
        view=self.viewer();await view.explorer.inspector.calculate()
        restored=ResultsViewer(view.ui,'project',dict(view.data,fragment=dict(view.explorer.inspector.result)));restored.build()
        inspector=restored.explorer.inspector
        restored.nodes['0']['properties']['canonical_smiles']='CCO'
        inspector.scope.value='full';await inspector.change_scope()
        inspector.scope.value='mcc';await inspector.change_scope()
        self.assertEqual(inspector.result['atoms'],2)
        self.assertNotIn('persistido',inspector.result['scope'])

    async def test_full_graph_default_capacity_includes_more_than_fifty_compounds(self):
        import networkx as nx
        import test_network_analysis as analytical
        view=analytical.AnalyticalViewerTests().viewer(nx.path_graph(251))
        for node in view.data['nodes']:node['properties']['canonical_smiles']='CCCO'
        inspector=view.explorer.inspector
        self.assertEqual(int(inspector.limit.value),251)
        inspector.scope.value='full';await inspector.scope.on_select(None)
        self.assertEqual(inspector.result['compounds'],251)
        self.assertEqual(inspector.result['status'],'complete')
        self.assertEqual(set(inspector.result['ids']),set(view.nodes))

    async def test_community_option_offers_local_selection_and_computes_without_external_picker(self):
        view=self.viewer();explorer=view.explorer;inspector=explorer.inspector
        view.analysis.community_cache['full']=({'0':1,'1':1,'2':2,'3':2},None,'')
        inspector.scope.value='community';await inspector.scope.on_select(None)
        self.assertTrue(inspector.community.visible)
        self.assertTrue(inspector.community_row.visible)
        self.assertEqual({o.key for o in inspector.community.options},{'1','2'})
        self.assertEqual(explorer.community_choice.value,'all')
        self.assertEqual(inspector.result['status'],'complete')
        inspector.community.value='1';await inspector.community.on_select(None)
        self.assertEqual(inspector.result['ids'],['0','1'])
        inspector.community.value='2';await inspector.community.on_select(None)
        self.assertEqual(inspector.result['ids'],['2','3'])
        self.assertEqual(explorer.community_choice.value,'all')

    async def test_large_capacity_is_valid_but_an_explicitly_lowered_limit_is_respected(self):
        inspector=self.viewer().explorer.inspector
        inspector.limit.value='201';await inspector.calculate()
        self.assertEqual(inspector.result['status'],'complete')
        inspector.scope.value='full';inspector.limit.value='3';await inspector.change_scope()
        self.assertIsNone(inspector.result)
        self.assertIn('limite configurado: 3',inspector.notice.value)

    async def test_mcs_scaffold_image_is_drawn_directly_without_a_frame(self):
        import flet as ft
        inspector=self.viewer().explorer.inspector
        await inspector.calculate()
        self.assertIsInstance(inspector.result_panel.controls[0],ft.Image)
