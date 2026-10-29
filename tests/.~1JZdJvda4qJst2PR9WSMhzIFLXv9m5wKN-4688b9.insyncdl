"""Real Flet schema validation, stable typed sorting and read-only inspection."""
import asyncio
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import flet as ft
import test_consensus_automation as consensus
import test_compound_tables as compounds
from biomolexplorer.table_sort import ordered
from biomolexplorer.ui.docking_results import DockingResultsTable
from biomolexplorer.ui.result_table import ResultTableViewer
from biomolexplorer.result_tables import ResultTables


class TypedSortingTests(unittest.TestCase):
    def test_numeric_order_is_stable_and_missing_values_stay_last_both_directions(self):
        values=['-2','-10','2e1',None,'NaN','undefined','',float('inf'),'2']
        self.assertEqual(ordered(values,lambda x:x,True,True)[:4],['-10','-2','2','2e1'])
        self.assertEqual(ordered(values,lambda x:x,False,True)[:4],['2e1','2','-2','-10'])
        self.assertEqual(ordered(values,lambda x:x,False,True)[4:],values[3:8])
        self.assertEqual(ordered(['b','A','a',None],lambda x:x),['A','a','b',None])

    def test_file_sorting_keeps_origin_filter_and_does_not_modify_the_source_list(self):
        from biomolexplorer.ui.file_table import FileTable
        async def scenario():
            files=[dict(name='a.csv',path='a',size=10,origin='A'),dict(name='b.csv',path='b',size=2,origin='A'),dict(name='c.csv',path='c',size=1,origin='B')]
            ui=SimpleNamespace(token='token',current={'id':'project'},page=SimpleNamespace(update=lambda:None))
            viewer=FileTable(ui,'project',files);viewer.origin.value='A';viewer.filter_origin(None)
            await viewer.table.columns[1].on_sort(None)
            self.assertEqual(viewer.table.rows[0].cells[0].content.content.controls[0].value,'b.csv')
            self.assertEqual(viewer.origin.value,'A');self.assertEqual([f['name'] for f in files],['a.csv','b.csv','c.csv'])
        asyncio.run(scenario())


class DockingNavigationTests(unittest.TestCase):
    setUp=consensus.ConsensusSortingTests.setUp
    tearDown=compounds.CompoundTablesTests.tearDown

    def test_every_update_is_a_valid_table_and_sorting_preserves_project_filters_and_cached_data(self):
        async def scenario():
            async def call(fn,*args,**kwargs):return fn(*args,**kwargs)
            async def guard(fn):return await fn()
            project={'id':self.pid};viewer=None;updates=[]
            def update():
                if viewer is not None:
                    viewer.table.before_update();updates.append(len(viewer.table.rows))
            ui=SimpleNamespace(store=self.store,token=self.owner,current=project,tab='Execuções',page=SimpleNamespace(update=update),call=call,guard=guard)
            viewer=DockingResultsTable(ui,self.pid,self.rid,self.sid,self.tables.tables(self.owner,self.pid,self.rid,self.sid),inline=True)
            original=self.path.read_bytes();selector=viewer.selector.value
            with patch.object(viewer.service,'_digest',wraps=viewer.service._digest) as reads:
                await viewer.load();viewer.search.value='M0';await viewer.search_rows(None)
                for index in (1,1,2,3,0):await viewer.table.columns[index].on_sort(None)
                self.assertEqual(reads.call_count,1)
            self.assertEqual(viewer.search.value,'M0');self.assertEqual(viewer.selector.value,selector)
            self.assertIs(ui.current,project);self.assertEqual(ui.tab,'Execuções');self.assertTrue(viewer.active)
            self.assertEqual(self.path.read_bytes(),original);self.assertGreater(len(updates),5)
        asyncio.run(scenario())


class GenericResultSortingTests(unittest.TestCase):
    setUp=compounds.CompoundTablesTests.setUp
    tearDown=compounds.CompoundTablesTests.tearDown

    def test_arbitrary_csv_fields_sort_globally_and_original_indices_remain_valid(self):
        self.path.write_text('code,energy,Alert_PAINS\nB,2,true\nA,-11,false\nC,NaN,true\nD,-2,false\n')
        service=ResultTables(self.store)
        result=service.page(self.owner,self.pid,self.rid,self.sid,str(self.path),sort_by='energy',limit=2)
        self.assertEqual([r['index'] for r in result['rows']],[1,3]);self.assertFalse(result['can_edit'])
        result=service.page(self.owner,self.pid,self.rid,self.sid,str(self.path),sort_by='energy',ascending=False)
        self.assertEqual([r['properties']['code'] for r in result['rows']],['B','D','A','C'])
        original=self.path.read_bytes()
        with self.assertRaisesRegex(ValueError,'somente leitura'):service.remove(self.owner,self.pid,self.rid,self.sid,str(self.path),0,result['version'])
        self.assertEqual(self.path.read_bytes(),original)
        import json
        run=self.store.get_run(self.owner,self.rid);run['stages'][0]['status']='failed'
        with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(run['stages']),self.rid))
        self.assertEqual(service.page(self.owner,self.pid,self.rid,self.sid,str(self.path))['total'],4)
        self.path.write_text('code,energy\nE,-30\n')
        result=service.page(self.owner,self.pid,self.rid,self.sid,str(self.path),sort_by='energy')
        self.assertEqual(result['rows'][0]['id'],'E')

    def test_generic_viewer_commits_columns_and_cells_together_when_schema_changes(self):
        async def scenario():
            async def call(fn,*args,**kwargs):return fn(*args,**kwargs)
            async def guard(fn):return await fn()
            ui=SimpleNamespace(store=self.store,token=self.owner,current={'id':self.pid},page=SimpleNamespace(update=lambda:None),call=call,guard=guard)
            viewer=ResultTableViewer(ui,self.pid,self.rid,self.sid,ResultTables(self.store).tables(self.owner,self.pid,self.rid,self.sid),inline=True)
            ui.page.update=lambda:viewer.table.before_update()
            self.path.write_text('code,energy,Alert_PAINS\nB,2,true\nA,-11,false\n')
            await viewer.load();await viewer.table.columns[1].on_sort(None)
            self.assertEqual(viewer.data['rows'][0]['id'],'A')
            self.path.write_text('code,energy\nX,0\n');await viewer.load()
            self.assertEqual(viewer.table.rows[0].cells[1].content.value,'0')
        asyncio.run(scenario())
