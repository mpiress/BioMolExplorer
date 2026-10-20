"""Incremental Flet controls, lazy scientific work and dialog request ordering."""
import asyncio
import copy
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import flet as ft
from biomolexplorer.catalog import TITLES, new_stage
from biomolexplorer import flow
from biomolexplorer.ui.flow_canvas import FlowCanvas
from biomolexplorer.ui.compound_catalog import CompoundCatalog
import test_pipeline_editor_regressions as fixtures
from test_frontend_navigation import walk


class EditorResponsivenessTests(unittest.TestCase):
    setUp=fixtures.PipelineEditorRegressionTests.setUp

    def test_add_vina_preserves_existing_controls_and_no_scientific_initialization(self):
        self.ui.schedule_save=lambda:None
        editor=FlowCanvas(self.ui,True)
        original=editor.nodes[self.stage['id']]
        with patch.object(editor,'node',wraps=editor.node) as builds, patch('biomolexplorer.pubchem_retrieval.reference_choices') as parse:
            asyncio.run(editor.add('docking_vina'))
            self.assertIs(original,editor.nodes[self.stage['id']])
            self.assertEqual(builds.call_count,1)
            asyncio.run(self.ui.open_stage_dialog(self.ui.selected))
            parse.assert_not_called()
        self.assertTrue(self.page.dialogs[-1].open)

    def test_replacement_rebinds_gestures_and_removal_releases_controls(self):
        self.ui.schedule_save=lambda:None
        editor=FlowCanvas(self.ui,True)
        old=editor.nodes[self.stage['id']]
        restored=copy.deepcopy(self.ui.current['pipeline'])
        self.ui.current['pipeline']=restored;editor.refresh()
        self.assertIsNot(old,editor.nodes[self.stage['id']])
        header=editor.nodes[self.stage['id']].content.controls[0]
        header.on_pan_update(SimpleNamespace(local_delta=SimpleNamespace(x=8,y=9)))
        self.assertIn('position',restored[0]);self.assertNotIn('position',self.stage)
        self.ui.selected=restored[0]['id'];asyncio.run(editor.delete(None))
        self.assertEqual(editor.nodes,{})
        self.assertEqual(editor._node_states,{})
        self.assertEqual(editor.output_draggables,{})
        asyncio.run(editor.undo(None))
        self.assertEqual(len(editor.nodes),1)

    def test_unbound_vina_does_not_read_history_but_source_selection_loads_ids(self):
        async def scenario():
            source=new_stage('retrieve_compounds');vina=new_stage('docking_vina')
            self.ui.current['pipeline']=[source,vina]
            path=self.store.project_dir(self.pid)/'compounds.csv'
            path.write_text('molecule_chembl_id,canonical_smiles\nCHEMBL1,CCO\n')
            item=dict(source,status='succeeded',artifacts=[str(path)])
            from biomolexplorer.pubchem_retrieval import reference_choices
            with patch.object(self.store,'list_runs',return_value=[dict(id='r',stages=[item])]), patch('biomolexplorer.pubchem_retrieval.reference_choices',wraps=reference_choices) as reader:
                await self.ui.open_stage_dialog(vina['id'])
                reader.assert_not_called()
                dialog=self.page.dialogs[-1]
                source_field=next(c for c in walk(dialog) if isinstance(c,ft.Dropdown) and c.label=='Compostos selecionados')
                source_field.value='stage:'+source['id']
                await source_field.on_select(None)
                self.assertEqual(reader.call_count,1)
                compound=next(c for c in walk(dialog) if isinstance(c,ft.Dropdown) and c.label=='Composto para docking (opcional)')
                self.assertIn('CHEMBL1',[o.key for o in compound.options])
                await source_field.on_select(None)
                self.assertEqual(reader.call_count,1)
                self.page.pop_dialog();self.ui.editing_stage=False
                vina['bindings']['base_selected_mols']={'stage':source['id'],'selector':'auto','compound_id':'CHEMBL1'}
                await self.ui.open_stage_dialog(vina['id'])
                compound=next(c for c in walk(self.page.dialogs[-1]) if isinstance(c,ft.Dropdown) and c.label=='Composto para docking (opcional)')
                self.assertIn('CHEMBL1',[o.key for o in compound.options])
        asyncio.run(scenario())

    def test_catalog_invalidates_changed_files_and_bounds_cache(self):
        async def scenario():
            path=self.store.project_dir(self.pid)/'compounds.csv'
            path.write_text('molecule_chembl_id,canonical_smiles\nCHEMBL1,CCO\n')
            catalog=CompoundCatalog(self.ui,{'source':[path]},lambda:True)
            await catalog.load('source')
            self.assertEqual(set(catalog.references['source'][str(path)]),{'CHEMBL1'})
            path.write_text('molecule_chembl_id,canonical_smiles\nCHEMBL999,CCN\n')
            await catalog.load('source')
            self.assertEqual(set(catalog.references['source'][str(path)]),{'CHEMBL999'})
            self.assertEqual(len(self.ui._compound_reference_cache),1)
            with patch('biomolexplorer.pubchem_retrieval.reference_choices',return_value={'id':'CCO'}):
                for index in range(40):
                    extra=path.with_name(f'{index}.csv');extra.write_text('test')
                    catalog.paths['source']=[extra];await catalog.load('source')
            self.assertLessEqual(len(self.ui._compound_reference_cache),32)
            catalog.valid=lambda:False
            with patch.object(self.ui,'call') as call:
                await catalog.load('source');call.assert_not_called()
        asyncio.run(scenario())

    def test_obsolete_catalog_read_does_not_repopulate_session_cache(self):
        async def scenario():
            path=self.store.project_dir(self.pid)/'compounds.csv';path.write_text('test')
            entered=asyncio.Event();release=asyncio.Event();active=True
            catalog=CompoundCatalog(self.ui,{'source':[path]},lambda:active)
            async def call(function,*args,**kwargs):
                if function==catalog._signature:return function(*args,**kwargs)
                entered.set();await release.wait();return {'CHEMBL1':'CCO'}
            self.ui.call=call
            reading=asyncio.create_task(catalog.load('source'));await entered.wait()
            active=False;self.ui._compound_reference_cache=None
            release.set();await reading
            self.assertIsNone(self.ui._compound_reference_cache)
            self.assertEqual(catalog.references,{})
        asyncio.run(scenario())

    def test_latest_dialog_request_wins_and_does_not_mix_forms(self):
        async def scenario():
            vina=new_stage('docking_vina');self.ui.current['pipeline'].append(vina)
            entered=asyncio.Event();release=asyncio.Event();first=True
            async def call(function,*args,**kwargs):
                nonlocal first
                if first and function==self.store.project:
                    first=False;entered.set();await release.wait()
                return function(*args,**kwargs)
            self.ui.call=call
            old=asyncio.create_task(self.ui.open_stage_dialog(vina['id']))
            await entered.wait()
            latest=asyncio.create_task(self.ui.open_stage_dialog(self.stage['id']))
            await asyncio.sleep(0);release.set();await asyncio.gather(old,latest)
            self.assertEqual(len(self.page.dialogs),1)
            labels=[c.label for c in walk(self.page.dialogs[-1]) if isinstance(c,ft.TextField)]
            self.assertIn('Nome da coleção PDB (opcional)',labels)
        asyncio.run(scenario())

    def test_all_block_forms_open_successively(self):
        async def scenario():
            self.ui.current['pipeline']=[new_stage(operation) for operation in TITLES]
            for _ in range(2):
                for stage in self.ui.current['pipeline']:
                    await self.ui.open_stage_dialog(stage['id'])
                    self.assertTrue(self.page.dialogs[-1].open)
                    self.page.dialogs[-1].on_dismiss(None);self.page.pop_dialog()
                    self.assertFalse(self.ui.editing_stage)
        asyncio.run(scenario())

    def test_bottom_node_drag_keeps_grid_and_uses_targeted_updates(self):
        from unittest.mock import Mock
        self.ui.current['pipeline']=[new_stage('docking_vina') for _ in range(100)]
        self.ui.page=SimpleNamespace(update=Mock())
        editor=FlowCanvas(self.ui,True)
        node=editor.nodes[editor.stages[-1]['id']]
        grid=editor.grid
        header=node.content.controls[0]
        for _ in range(40):
            header.on_pan_update(SimpleNamespace(local_delta=SimpleNamespace(x=1,y=1)))
            self.assertIs(editor.grid,grid)
            self.ui.page.update.assert_called_with(node,editor.edges)
        self.assertGreaterEqual(editor.stack.height,node.top+flow.node_height(editor.stages[-1]))

    def test_producer_rename_invalidates_consumer_labels_only(self):
        source=new_stage('retrieve_compounds');target=new_stage('docking_vina')
        other=new_stage('retrieve_structures')
        self.ui.current['pipeline']=[source,target,other]
        flow.connect(self.ui.current['pipeline'],source['id'],target['id'],'base_selected_mols')
        editor=FlowCanvas(self.ui,True);before=dict(editor.nodes)
        source['name']='Renamed producer';editor.refresh()
        self.assertIsNot(editor.nodes[target['id']],before[target['id']])
        self.assertIs(editor.nodes[other['id']],before[other['id']])
        labels=editor._node_states[target['id']][2]
        self.assertTrue(any('Renamed producer' in label for label in labels))
        self.assertEqual(len(editor.edges.shapes),1)

    def test_repeated_connections_and_disconnections_at_supported_sizes(self):
        async def scenario():
            self.ui.schedule_save=lambda:None
            for count in (10,25,50,100):
                source=new_stage('retrieve_compounds');target=new_stage('docking_vina')
                self.ui.current['pipeline']=[source,target]+[new_stage('retrieve_structures') for _ in range(count-2)]
                self.ui.selected=None;editor=FlowCanvas(self.ui,True)
                untouched=editor.nodes[editor.stages[-1]['id']]
                for _ in range(10):
                    await editor.link(source['id'],target['id'],'base_selected_mols')
                    self.assertEqual(len(editor.edges.shapes),1)
                    flow.disconnect(target,'base_selected_mols');editor.changed()
                    self.assertEqual(len(editor.edges.shapes),0)
                    self.assertIs(editor.nodes[editor.stages[-1]['id']],untouched)
                    self.assertEqual(len(editor.nodes),count)
        asyncio.run(scenario())

    def test_hundred_block_repeated_edits_keep_control_counts_bounded(self):
        async def scenario():
            self.ui.schedule_save=lambda:None
            self.ui.current['pipeline']=[new_stage('docking_vina') for _ in range(100)]
            editor=FlowCanvas(self.ui,True)
            for _ in range(50):
                selected=editor.stages[-1]['id'];self.ui.selected=selected
                await editor.delete(None);await editor.add('docking_vina')
                self.assertEqual(len(editor.nodes),100)
                self.assertEqual(len(editor._node_states),100)
                self.assertEqual(len(editor.output_draggables),100)
                header=editor.nodes[self.ui.selected].content.controls[0]
                for _ in range(3):header.on_pan_update(SimpleNamespace(local_delta=SimpleNamespace(x=1,y=1)))
            self.assertLessEqual(len(editor.history),50)
        asyncio.run(scenario())


if __name__=='__main__':unittest.main()
