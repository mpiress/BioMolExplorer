"""Regression checks for bounded table views and live access checks."""
import unittest
from unittest.mock import patch
import test_compound_tables as fixtures
from biomolexplorer import compound_tables
from biomolexplorer.workspace import AccessDenied


class TableViewTests(unittest.TestCase):
    setUp=fixtures.CompoundTablesTests.setUp
    tearDown=fixtures.CompoundTablesTests.tearDown
    page=fixtures.CompoundTablesTests.page
    def test_navigation_reuses_sort_and_materializes_only_visible_links(self):
        with patch.object(compound_tables, 'ordered', wraps=compound_tables.ordered) as sorting, \
             patch.object(compound_tables, 'compound_links', wraps=compound_tables.compound_links) as links:
            first=self.page(sort_by='source',limit=1)
            second=self.page(sort_by='source',offset=1,limit=1)
            self.assertEqual(sorting.call_count,1)
            self.assertEqual(links.call_count,2)
            self.assertNotEqual(first['rows'][0]['index'],second['rows'][0]['index'])
            first['rows'][0]['properties']['source']='mutated response'
            self.assertNotEqual(self.page(sort_by='source',limit=1)['rows'][0]['properties']['source'],'mutated response')
            self.path.write_text(self.original.replace('CHEMBL1','CHEMBL99'))
            refreshed=self.page(sort_by='source',limit=1)
            self.assertNotEqual(first['version'],refreshed['version'])
            self.assertEqual(sorting.call_count,2)

    def test_cached_view_does_not_bypass_authorization(self):
        self.page()
        with self.assertRaises(AccessDenied):self.page(self.guest)

    def test_streaming_pages_match_cached_pages(self):
        cached=self.page(offset=1,limit=1)
        from types import SimpleNamespace
        from contextlib import contextmanager
        original=self.tables._reader
        @contextmanager
        def streaming(path):
            with original(path) as (version,reader):
                yield version,SimpleNamespace(fieldnames=reader.fieldnames,rows=iter(reader.rows))
        with patch.object(self.tables,'_reader',streaming):
            streamed=self.page(offset=1,limit=1)
        self.assertEqual(cached,streamed)

    def test_active_run_lookup_uses_partial_index(self):
        with self.store.connect() as db:
            plan=db.execute("EXPLAIN QUERY PLAN SELECT 1 FROM runs WHERE project_id=? AND status IN ('queued','running','awaiting_input')",(self.pid,)).fetchall()
        self.assertIn('runs_active_project',' '.join(str(tuple(row)) for row in plan))

    def test_asset_verification_batches_queries_and_rejects_missing_assets(self):
        from contextlib import contextmanager
        from biomolexplorer import published_results
        ids=[f'asset{i}' for i in range(501)]
        calls=[]
        class Database:
            def execute(self,sql,params):
                calls.append((sql,params))
                return self
            def fetchall(self):
                return [(asset,str(self_path)) for asset in calls[-1][1][1:]]
        self_path=self.path
        project_root=self.store.project_dir(self.pid)
        @contextmanager
        def connect():yield Database()
        entry={'configuration':{},'input_files':{},'source_manifest':{asset:'digest' for asset in ids}}
        with patch.object(published_results,'equivalent',return_value=True), \
             patch.object(published_results,'intact',return_value=True), \
             patch.object(published_results,'asset_references',return_value=ids), \
             patch.object(published_results,'input_digest',return_value='digest'), \
             patch.object(self.store,'connect',connect), \
             patch.object(self.store,'project_dir',return_value=project_root):
            self.assertTrue(published_results.verify(self.store,self.pid,{},entry,[],{}))
            self.assertEqual(len(calls),2)
            self.assertEqual([len(params)-1 for _,params in calls],[500,1])
            with patch.object(Database,'fetchall',return_value=[]):
                self.assertFalse(published_results.verify(self.store,self.pid,{},entry,[],{}))
