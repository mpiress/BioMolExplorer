"""Regression checks for MCS configuration changes and explanatory ranking help."""
import asyncio
import threading
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import networkx as nx
import flet as ft
from biomolexplorer.network_analysis import METRICS,NetworkAnalysis
from biomolexplorer.network_scaffolds import common
from biomolexplorer.molecular_network import numeric_properties
from biomolexplorer.ui.localization import LocalizedPage,Translator
from biomolexplorer.ui.ranking_help import HELP,toy_values
from test_network_analysis import AnalyticalViewerTests


class MCSConfigurationTests(unittest.TestCase):
    def viewer(self):
        fixture=AnalyticalViewerTests()
        view=fixture.viewer(nx.disjoint_union(nx.path_graph(3),nx.empty_graph(1)))
        for node,smiles in zip(view.data['nodes'],('CCCO','CCCN','CCC','CCO')):
            node['properties']['canonical_smiles']=smiles
        return view

    def test_selecting_all_scopes_recalculates_the_explicit_set_and_fragment(self):
        async def scenario():
            view=self.viewer();inspector=view.explorer.inspector
            await inspector.calculate()
            self.assertEqual(inspector.result['atoms'],3)
            self.assertEqual(inspector.result['ids'],['0','1','2'])
            old=inspector.result['image']
            inspector.scope.value='full';await inspector.scope.on_select(None)
            self.assertEqual(inspector.result['ids'],['0','1','2','3'])
            self.assertEqual(inspector.result['atoms'],2)
            self.assertNotEqual(inspector.result['image'],old)
            view.analysis.community_cache['full']=({'0':1,'1':1,'2':2,'3':2},None,'')
            view.explorer.community_choice.value='1'
            inspector.scope.value='community';await inspector.scope.on_select(None)
            self.assertEqual(inspector.result['ids'],['0','1'])
            self.assertEqual(inspector.result['atoms'],3)
            view.explorer.filtered_graph=view.graph.subgraph(['2','3'])
            inspector.scope.value='filtered';await inspector.scope.on_select(None)
            self.assertEqual(inspector.result['ids'],['2','3'])
            self.assertEqual(inspector.result['atoms'],2)
            self.assertEqual(view.data['fragment'] if 'fragment' in view.data else None,None)
        asyncio.run(scenario())

    def test_invalid_set_clears_previous_result_without_sampling(self):
        async def scenario():
            view=self.viewer();inspector=view.explorer.inspector
            await inspector.calculate()
            inspector.limit.value='3';await inspector.limit.on_change(None)
            inspector.scope.value='full'
            with patch('biomolexplorer.ui.molecular_inspector.common',wraps=common) as compute:
                await inspector.scope.on_select(None)
            self.assertIsNone(inspector.result);self.assertFalse(inspector.result_panel.controls)
            self.assertIn('4 compostos',inspector.notice.value);self.assertEqual(compute.call_count,0)
            self.assertFalse(inspector.compute_button.disabled)
        asyncio.run(scenario())

    def test_parameter_changes_remove_stale_highlights_and_have_distinct_cache_entries(self):
        async def scenario():
            view=self.viewer();inspector=view.explorer.inspector
            with patch('biomolexplorer.ui.molecular_inspector.common',wraps=common) as compute:
                await inspector.calculate();await view.select('0',show_popup=False)
                self.assertEqual(inspector.structure['highlighted_atoms'],3)
                inspector.timeout.value='1';await inspector.timeout.on_change(None)
                self.assertIsNone(inspector.result)
                self.assertEqual(inspector.structure['highlighted_atoms'],0)
                await inspector.calculate()
                self.assertEqual(inspector.result['timeout'],1)
                inspector.limit.value='4';await inspector.limit.on_change(None);await inspector.calculate()
                self.assertEqual(inspector.result['limit'],4)
                self.assertEqual(compute.call_count,3)
                self.assertEqual(compute.call_args.args[1:3],(1,4))
                await inspector.calculate();self.assertEqual(compute.call_count,3)
            self.assertFalse(view.ui.page.show_dialog.called)
        asyncio.run(scenario())

    def test_new_scope_cancels_pending_search_and_only_latest_result_is_displayed(self):
        async def scenario():
            view=self.viewer();inspector=view.explorer.inspector;original=view.ui.call
            started,release=asyncio.Event(),asyncio.Event();calls=[]
            async def call(function,*args,**kwargs):
                if function is common:
                    calls.append(args)
                    if len(calls)==1:
                        started.set();await release.wait()
                return await original(function,*args,**kwargs)
            view.ui.call=call
            first=asyncio.create_task(inspector.calculate());await started.wait()
            inspector.scope.value='full';await inspector.scope.on_select(None)
            self.assertTrue(calls[0][3].is_set());self.assertIsNone(inspector.result)
            release.set();await first
            self.assertEqual(len(calls),2)
            self.assertEqual(inspector.result['ids'],['0','1','2','3'])
            self.assertEqual(inspector.result['atoms'],2)
            self.assertFalse(inspector.running);self.assertFalse(inspector.compute_button.disabled)
        asyncio.run(scenario())

    def test_scope_change_during_preparation_cannot_start_a_stale_search(self):
        async def scenario():
            view=self.viewer();inspector=view.explorer.inspector;original=view.ui.call
            started,release=asyncio.Event(),asyncio.Event();blocked=[False];sets=[]
            async def call(function,*args,**kwargs):
                if function is numeric_properties and not blocked[0]:
                    blocked[0]=True;started.set();await release.wait()
                if function is common:sets.append([n for n,_ in args[0]])
                return await original(function,*args,**kwargs)
            view.ui.call=call
            first=asyncio.create_task(inspector.calculate());await started.wait()
            inspector.scope.value='full';await inspector.scope.on_select(None)
            release.set();await first
            self.assertEqual(sets,[['0','1','2','3']])
        asyncio.run(scenario())

    def test_invalid_parameter_values_do_not_start_a_search(self):
        async def scenario():
            inspector=self.viewer().explorer.inspector
            for timeout,limit in [('0','50'),('31','50'),('10','1'),('10','201'),('1.5','50'),('10','')]:
                with self.subTest(timeout=timeout,limit=limit):
                    inspector.timeout.value=timeout;inspector.limit.value=limit
                    await inspector.change_parameters()
                    with patch('biomolexplorer.ui.molecular_inspector.common') as compute:
                        await inspector.calculate()
                    self.assertFalse(compute.called);self.assertIsNone(inspector.result)
        asyncio.run(scenario())

    def test_maximum_time_is_passed_to_rdkit_and_partial_result_retains_parameters(self):
        from rdkit.Chem import rdFMCS
        original=rdFMCS.FindMCS
        with patch('rdkit.Chem.rdFMCS.FindMCS',wraps=original) as search:
            result=common((('A','CCCO'),('B','CCCN')),1,2,threading.Event())
        params=search.call_args.args[1]
        self.assertEqual(params.Timeout,1);self.assertEqual(params.Threshold,1)
        self.assertEqual(result['timeout'],1);self.assertEqual(result['limit'],2)
        with patch('rdkit.Chem.rdFMCS.FindMCS',return_value=SimpleNamespace(smartsString='[#6]-[#6]',numAtoms=2,numBonds=1,canceled=True)):
            partial=common((('A','CCCO'),('B','CCCN')),1,2,threading.Event())
        self.assertEqual(partial['status'],'partial');self.assertEqual(partial['timeout'],1)
        self.assertIn('maximalidade não confirmada',partial['notice'])


class RankingHelpTests(unittest.TestCase):
    def test_every_metric_has_a_translated_explanation_and_matching_toy_values(self):
        self.assertEqual(set(HELP),set(METRICS))
        analysis=NetworkAnalysis({'nodes':[{'id':n} for n in ('A','B','C')],
            'edges':[{'source':'A','target':'B'},{'source':'B','target':'C'}]})
        en,pt=Translator('en'),Translator('pt')
        for metric,paragraphs in HELP.items():
            with self.subTest(metric=metric):
                for text in paragraphs:
                    self.assertNotEqual(en(text),text);self.assertEqual(pt(en(text)),text)
                self.assertEqual(toy_values(metric),analysis.metrics('full',metric)[0])
                self.assertGreater(toy_values(metric)['B'],toy_values(metric)['A'])

    def test_info_button_opens_help_for_the_selected_metric_without_mutating_the_graph(self):
        async def scenario():
            view=AnalyticalViewerTests().viewer(nx.path_graph(3));explorer=view.explorer
            before=(view.mode,set(view.points),explorer.metric,dict(view.analysis.cache))
            for metric in METRICS:
                explorer.rank_metric.value=metric
                await explorer.rank_info.on_tap(None)
                dialog=view.ui.page.show_dialog.call_args.args[0]
                self.assertEqual(dialog.title.value,METRICS[metric].split(' · ')[0])
                texts=[getattr(c,'value',None) for c in dialog.content.content.controls]
                self.assertTrue(all(text in texts for text in HELP[metric]))
                self.assertIn('A — B — C',texts)
                page=LocalizedPage(SimpleNamespace(),'en');page.localize(dialog)
                self.assertIn('What it measures',[getattr(c,'value',None) for c in dialog.content.content.controls])
            self.assertEqual(before,(view.mode,set(view.points),explorer.metric,dict(view.analysis.cache)))
        asyncio.run(scenario())


if __name__=='__main__': unittest.main()
