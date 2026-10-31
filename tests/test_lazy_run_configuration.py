"""Lazy stage boundaries warm data once without opening or advancing a run."""
import asyncio
import copy
import tempfile
import threading
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

from biomolexplorer.catalog import new_stage
from biomolexplorer.ui.app import WorkspaceUI
from biomolexplorer.ui.file_selection import FileSelection,prepare_selection_data,selection_file_state
from biomolexplorer.ui.redocking_results import RedockingResultsTable


class LazyPreparationTests(unittest.IsolatedAsyncioTestCase):
    async def test_prefetch_is_single_flight_off_loop_and_configuration_scoped(self):
        ui=object.__new__(WorkspaceUI);ui.token='token';ui.store=SimpleNamespace(assets=Mock(return_value=[]))
        stage=new_stage('admet')
        pending=dict(id=stage['id'],status='awaiting_input',configuration=stage)
        run=dict(id='run',project_id='project',status='awaiting_input',stages=[pending])
        original=copy.deepcopy(run);loop_thread=threading.get_ident();threads=[]
        def prepare(*args):threads.append(threading.get_ident());return {'ready':True}
        with patch('biomolexplorer.ui.file_selection.prepare_selection_data',side_effect=prepare) as work:
            first=ui.prepare_pending_inputs(run,pending)
            self.assertIs(first,ui.prepare_pending_inputs(run,pending))
            self.assertEqual(await first,{'ready':True})
            self.assertIs(first,ui.prepare_pending_inputs(run,pending))
            work.assert_called_once()
            self.assertNotEqual(threads[0],loop_thread)
            self.assertEqual(run,original)
            pending['configuration']['parameters']['lipinski']=False
            second=ui.prepare_pending_inputs(run,pending)
            self.assertIsNot(first,second)
            await second
            self.assertEqual(work.call_count,2)

    async def test_prefetch_failure_does_not_change_run_or_open_ui(self):
        ui=object.__new__(WorkspaceUI);ui.token='token';ui.store=SimpleNamespace(assets=Mock(side_effect=OSError('read failed')))
        stage=new_stage('admet');pending=dict(id=stage['id'],status='awaiting_input',configuration=stage)
        run=dict(id='run',project_id='project',status='awaiting_input',stages=[pending]);before=copy.deepcopy(run)
        with patch('biomolexplorer.ui.app.get_logger') as logger:
            task=ui.prepare_pending_inputs(run,pending)
            with self.assertRaises(OSError):await task
            await asyncio.sleep(0)
            logger.return_value.warning.assert_called_once()
        self.assertEqual(run,before)


class PreparedDataTests(unittest.TestCase):
    def test_prepared_redocking_data_preserves_choices_without_parsing_on_open(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);pdb=root/'1ABC.pdb';pdb.write_text('HEADER TEST\n')
            meta=root/'pdb_codes.csv';meta.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN\n1ABC,LIG,2,A\n')
            source=new_stage('retrieve_structures');stage=new_stage('redocking')
            stage['bindings']['base_input_path']={'stage':source['id'],'selector':'1ABC.pdb'}
            entry=dict(id=source['id'],name=source['name'],operation=source['operation'],status='succeeded',artifacts=[str(pdb),str(meta)],configuration=source)
            pending=dict(id=stage['id'],name=stage['name'],operation=stage['operation'],status='awaiting_input',configuration=stage)
            run=dict(id='run',project_id='project',stages=[entry,pending]);before=copy.deepcopy(run)
            with patch('biomolexplorer.result_files.ResultFiles.pdb_review_statuses',return_value={}):
                data=prepare_selection_data(SimpleNamespace(),'token',run,pending,[])
            baseline=FileSelection(run,pending).read()
            with patch('biomolexplorer.redocking_config.metadata_records',side_effect=AssertionError('metadata reparsed')),patch('biomolexplorer.ui.file_selection.artifact_records',side_effect=AssertionError('catalog rebuilt')),patch('biomolexplorer.result_origins.ResultOrigins',side_effect=AssertionError('lineage rebuilt')):
                prepared=FileSelection(run,pending,prepared=data)
                self.assertEqual(prepared.read(),baseline)
            self.assertEqual(run,before)
            meta.write_text(meta.read_text()+'1ABC,OTH,3,B\n')
            self.assertNotEqual(selection_file_state(data['paths']),data['file_state'])


class RedockingPaginationTests(unittest.TestCase):
    def table(self):
        rows=[dict(id=str(i),target='A' if i<30 else 'B',targets=['A' if i<30 else 'B'],
            pdb=f'{i:04d}',ligand='LIG',residue=i,chain='A',rmsd=float(i),status='succeeded') for i in range(61)]
        ui=SimpleNamespace(store=None,language='en',page=SimpleNamespace(update=Mock(),web=False))
        results=SimpleNamespace(ui=ui,valid=lambda:True,writable=False)
        return RedockingResultsTable(results,rows)

    def test_all_targets_every_row_reachable_and_count_matches_page(self):
        table=self.table();seen=[]
        while True:
            seen += [r.cells[1].content.value for r in table.table.rows]
            self.assertIn(f'{table.offset+1}–{table.offset+len(table.table.rows)}',table.count.value)
            if table.next.disabled:break
            table.next.on_click(None)
        self.assertEqual(len(seen),61);self.assertEqual(len(set(seen)),61)
        table.move(1);self.assertEqual(table.offset,50)
        table.target_filter.value='target:A';table.filter_changed()
        self.assertEqual(table.offset,0);self.assertFalse(table.next.disabled)
        table.move(1);self.assertEqual(len(table.table.rows),5);self.assertTrue(table.next.disabled)
        table.search.value='missing';table.filter_changed()
        self.assertEqual(table.table.rows,[]);self.assertIn('0–0 of 0',table.count.value)
        self.assertTrue(table.next.disabled);self.assertTrue(table.previous.disabled)

    def test_navigation_precedes_wide_table_and_is_width_bounded(self):
        import flet as ft
        table=self.table();view=table.build()
        self.assertEqual(view.horizontal_alignment,ft.CrossAxisAlignment.STRETCH)
        navigation=next(i for i,c in enumerate(view.controls) if table.next in (getattr(c,'controls',[]) or []))
        viewport=next(i for i,c in enumerate(view.controls) if table.table in (getattr(c,'controls',[]) or []))
        self.assertLess(navigation,viewport)
