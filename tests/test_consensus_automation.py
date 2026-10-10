"""Consensus needs only connected results and sorts the complete score table."""
import asyncio
import csv
import json
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import test_pipeline_execution as execution
import test_compound_tables as compound_tests
from biomolexplorer.catalog import new_stage
from biomolexplorer.pipeline import requires_curation
from biomolexplorer.ui.guided import GuidedForm
from biomolexplorer.ui.docking_results import DockingResultsTable
from biomolexplorer.docking_results import DockingResults
from wrappers.docking import generate_consensus
from biomolexplorer.docking_data import write_csv


class ConsensusAutomationTests(unittest.TestCase):
    setUp=execution.PipelineExecutionTests.setUp
    tearDown=execution.PipelineExecutionTests.tearDown
    upload=execution.PipelineExecutionTests.upload
    imported=execution.PipelineExecutionTests.imported
    save=execution.PipelineExecutionTests.save
    wait=execution.PipelineExecutionTests.wait
    finish_selected=execution.PipelineExecutionTests.finish_selected
    immediate_manager=execution.PipelineExecutionTests.immediate_manager
    # Use existing integration scenario, but never confirm a popup for consensus.
    def test_imported_docking_results_connect_through_the_two_visual_consensus_ports(self):
        with patch.object(self,'finish_selected',side_effect=lambda service,run:self.wait(run)):
            execution.PipelineExecutionTests.test_imported_docking_results_connect_through_the_two_visual_consensus_ports(self)

    def test_configuration_contains_only_the_two_inputs(self):
        vina,dock6,stage=[new_stage(op) for op in ('docking_vina','docking_dock6','consensus')]
        stage['bindings']={'base_vina_path':{'stage':vina['id']},'base_dock6_path':{'stage':dock6['id']}}
        ui=SimpleNamespace(current={'pipeline':[vina,dock6,stage]},service=SimpleNamespace(dock6_path=None),page=SimpleNamespace(update=lambda:None))
        form=GuidedForm(ui,stage,[],True)
        self.assertEqual(set(form.input_editors),{'base_vina_path','base_dock6_path'})
        self.assertNotIn('target',form.field_controls);self.assertNotIn('repulsion_weight',form.field_controls)
        self.assertFalse(hasattr(form,'processing'));self.assertFalse(hasattr(form,'provided'))
        self.assertEqual(form.read()['input_processing'],'merge')
        self.assertFalse(requires_curation(form.read()))


class ConsensusSortingTests(unittest.TestCase):
    tearDown=compound_tests.CompoundTablesTests.tearDown
    page=compound_tests.CompoundTablesTests.page
    def setUp(self):
        compound_tests.CompoundTablesTests.setUp(self)
        self.path.write_text('molecule_chembl_id,canonical_smiles,vina,dock6,min-max,conformer_file\n'+
            ''.join(f'M{i:02},CCO,{-i},{-i*10},{i/40},pose.mol2\n' for i in range(1,41)))
        self.tables=DockingResults(self.store)

    def test_numeric_sorting_precedes_pagination_and_keeps_original_indices(self):
        for column in ('vina','dock6','normalized_score'):
            for ascending in (True,False):
                data=self.page(sort_by=column,ascending=ascending,offset=10,limit=10)
                expected=list(range(1,41))
                expected.sort(key=lambda i:-i if column in ('vina','dock6') else i,reverse=not ascending)
                self.assertEqual([r['id'] for r in data['rows']],[f'M{i:02}' for i in expected[10:20]])
                self.assertEqual([r['index'] for r in data['rows']],[i-1 for i in expected[10:20]])
                self.assertEqual(data['matched'],40)
        self.assertEqual(self.page(sort_by='molecule_chembl_id',ascending=False,limit=1)['rows'][0]['id'],'M40')
        self.assertEqual(self.page(sort_by='vina',query='M0',limit=1)['rows'][0]['id'],'M09')
        with self.assertRaises(ValueError):self.page(sort_by='unknown')

    def test_headers_toggle_order_and_reset_pagination(self):
        async def scenario():
            async def call(fn,*args,**kwargs):return fn(*args,**kwargs)
            async def guard(fn):return await fn()
            ui=SimpleNamespace(store=self.store,token=self.owner,current={'id':self.pid},page=SimpleNamespace(update=lambda:None),call=call,guard=guard)
            viewer=DockingResultsTable(ui,self.pid,self.rid,self.sid,self.tables.tables(self.owner,self.pid,self.rid,self.sid),inline=True)
            await viewer.load()
            self.assertEqual([c.label.value for c in viewer.table.columns[:4]],['Código do composto','Score Vina','Score DOCK6','Score normalizado'])
            viewer.offset=25
            await viewer.table.columns[1].on_sort(None)
            self.assertEqual(viewer.offset,0);self.assertEqual(viewer.data['rows'][0]['id'],'M40')
            await viewer.table.columns[1].on_sort(None)
            self.assertEqual(viewer.data['rows'][0]['id'],'M01')
            self.assertFalse(viewer.table.sort_ascending)
        asyncio.run(scenario())

    def test_normalized_score_is_exported_and_ranked(self):
        root=self.path.parent
        for engine in ('vina','dock6'):
            folder=root/engine;folder.mkdir()
            rows=[]
            for i in (1,2,3):
                pose=folder/f'M{i}.mol2';pose.write_text('pose')
                rows.append(dict(molecule_chembl_id=f'M{i}',canonical_smiles='CCO',receptor_id='1ABC_A',engine=engine,score=-i,conformer_file=pose.name))
            write_csv(folder/'docking_results.csv',rows)
        with patch('wrappers.docking.plot_scatter_comparison'):
            frame=generate_consensus(str(root),str(root/'output'),'consensus',base_vina_path=str(root/'vina'),base_dock6_path=str(root/'dock6'))
        self.assertEqual(frame['molecule_chembl_id'].tolist(),['M3','M2','M1'])
        self.assertEqual(frame['normalized_score'].tolist(),[1.0,0.5,0.0])
        with (root/'output/consensus.csv').open() as stream:
            self.assertIn('normalized_score',next(csv.DictReader(stream)))

if __name__=='__main__':unittest.main()
