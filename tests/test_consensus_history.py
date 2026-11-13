"""Bounded consensus history preserves rollback and portable audit contents."""
import copy
import json
import time
import unittest
from pathlib import Path
from unittest.mock import patch
from uuid import uuid4

from biomolexplorer.catalog import new_stage
from biomolexplorer.project_state import snapshot,record,read_snapshot,_snapshot_runs
from biomolexplorer.workspace import WorkspaceStore
from test_project_portability import ProjectPortabilityTests


class ConsensusHistoryTests(unittest.TestCase):
    setUp=ProjectPortabilityTests.setUp
    tearDown=ProjectPortabilityTests.tearDown

    def completed(self):
        stage=new_stage('consensus')
        output=self.folder/'consensus.csv'
        output.write_text('molecule_chembl_id,vina,dock6,normalized_score\nM1,1,1,1\n')
        item=dict(id=stage['id'],name=stage['name'],operation='consensus',configuration=stage,
                  status='succeeded',artifacts=[str(output)],artifact_records=[],metadata={'source':str(output)})
        rid=uuid4().hex
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,?)',
                (rid,self.pid,self.store.user(self.owner)['id'],'succeeded',json.dumps([item]),time.time(),time.time(),None))
        return rid,item

    def compact_event(self):
        rid,item=self.completed()
        user=self.store.user(self.owner)['id']
        with self.store.connect() as db:
            legacy=record(self.store,db,self.pid,user,'execution','legacy')
            expected=snapshot(self.store,db,self.pid)
            event=record(self.store,db,self.pid,user,'execution','compact',
                         {'_snapshot_reference':legacy,'side':'after'},compact_runs=True)
        return rid,item,expected,event,legacy

    def test_exact_snapshot_and_before_reference_remain_readable(self):
        rid,item,expected,event,legacy=self.compact_event()
        raw=json.loads((self.folder/'.history/snapshots'/f'{event}-after.json').read_text())
        self.assertIsInstance(raw['runs'][0]['stages'],dict)
        actual=read_snapshot(self.store,self.pid,event,'after')
        actual.pop('_compact_runs')
        self.assertEqual(actual,expected)
        self.assertEqual(read_snapshot(self.store,self.pid,event),read_snapshot(self.store,self.pid,legacy,'after'))

    def test_unchanged_runs_reuse_blob_and_changed_runs_preserve_old_snapshot(self):
        rid,item,expected,event,legacy=self.compact_event()
        with self.store.connect() as db:
            with patch('biomolexplorer.project_state.hashlib.sha256',side_effect=AssertionError('Unchanged run rehashed')):
                _snapshot_runs(db,self.pid,self.folder)
            changed=dict(item,metadata={'source':'new'})
            db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps([changed]),rid))
            _snapshot_runs(db,self.pid,self.folder)
        self.assertEqual(read_snapshot(self.store,self.pid,event,'after')['runs'][0]['stages'],expected['runs'][0]['stages'])

    def test_corrupted_blob_is_rejected(self):
        rid,item,expected,event,legacy=self.compact_event()
        raw=json.loads((self.folder/'.history/snapshots'/f'{event}-after.json').read_text())
        key=raw['runs'][0]['stages']['_stage_blob']
        (self.folder/'.history/run-stages'/(key+'.json')).write_text('[]')
        with self.assertRaisesRegex(ValueError,'Corrupt'):
            read_snapshot(self.store,self.pid,event,'after')

    def test_compact_history_exports_in_existing_portable_format(self):
        rid,item,expected,event,legacy=self.compact_event()
        archive=self.store.export_project(self.owner,self.pid)
        other=WorkspaceStore(self.root/'other-machine')
        token=other.register('Other','other@example.org','other-password-2026')
        imported=other.import_project(token,archive,self.root/'imported')
        with other.connect() as db:
            events=db.execute("SELECT id FROM project_history WHERE project_id=? AND summary='compact'",(imported['id'],)).fetchall()
        self.assertEqual(len(events),1)
        restored=read_snapshot(other,imported['id'],events[0][0],'after')
        self.assertIsInstance(restored['runs'][0]['stages'],list)
        self.assertEqual(restored['runs'][0]['stages'][0]['status'],'succeeded')

    def test_preview_matches_authoritative_invalidation_without_full_history(self):
        rid,item=self.completed()
        from biomolexplorer.result_sync import plan
        from biomolexplorer.result_catalog import proof_stage
        stage=copy.deepcopy(item['configuration'])
        stage['parameters']['target']='changed'
        with self.store.connect() as db:
            full=plan(db,self.pid,[stage],self.folder)
            light=plan(db,self.pid,[stage],self.folder,history=[(rid,[proof_stage(item)],'succeeded')])
        for key in ('affected','shared_preserved','blocks','active_conflicts','confirmation'):
            self.assertEqual(full[key],light[key])

    def test_large_execution_payload_is_streamed_with_bounded_python_memory(self):
        import tracemalloc
        rid,item=self.completed()
        item['large_metadata']='x'*(20*1024*1024)
        payload=json.dumps([item])
        with self.store.connect() as db:
            db.execute('UPDATE runs SET stages=? WHERE id=?',(payload,rid))
            tracemalloc.start()
            try:
                _snapshot_runs(db,self.pid,self.folder)
                _,peak=tracemalloc.get_traced_memory()
            finally:tracemalloc.stop()
        self.assertLess(peak,8*1024*1024)


if __name__=='__main__':
    unittest.main()
