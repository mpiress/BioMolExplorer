"""Persistent paging matches CSV semantics and never bypasses live authorization."""
import copy
import csv
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from biomolexplorer.result_index import ResultIndex
from biomolexplorer.result_catalog import ProjectionReader,compact_stage
from biomolexplorer.table_sort import ordered
import test_compound_tables as fixtures


def digest(stream):
    import hashlib
    value=hashlib.sha256(stream.read()).hexdigest();stream.seek(0);return value


class IndexTests(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup)
        self.root=Path(temporary.name);self.path=self.root/'table.csv'
        self.fields=['molecule_chembl_id','canonical_smiles','score','value','description']
        self.rows=[dict(zip(self.fields,values)) for values in [
            ('B','CCO','-2','20','Texto com vírgula, e\nnova linha'),
            ('A','C','-10','2','Å ß'),('C','O','NaN','','missing'),
            ('D','CC','-2','2','é'),('E','N','','-2','other')]]
        self.write(self.rows)

    def write(self,rows):
        with self.path.open('w',encoding='utf-8-sig',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=self.fields);writer.writeheader();writer.writerows(rows)

    def test_typed_sort_missing_values_and_stable_source_ordinals_match_csv(self):
        for column in ('score','value','description'):
            for ascending in (True,False):
                expected=ordered(list(enumerate(self.rows)),lambda item:item[1][column],ascending,column!='description')
                actual=ResultIndex(self.root,self.path,digest).page(0,100,'',column,ascending,True)
                self.assertEqual(actual['visible'],expected)
        with self.path.open('rb') as stream:self.assertEqual(actual['version'],digest(stream))

    def test_literal_search_and_generic_search_have_distinct_views(self):
        index=ResultIndex(self.root,self.path,digest)
        self.assertEqual(index.page(0,100,'vírgula',None,True,True)['matched'],0)
        self.assertEqual(index.page(0,100,'vírgula',None,True,False)['matched'],1)
        self.assertEqual(index.page(0,100,'cco',None,True,True)['visible'][0][0],0)
        self.assertEqual(index.page(0,100,'%',None,True,False)['matched'],0)

    def test_new_service_reuses_index_and_view_without_csv_hash_or_parse(self):
        first=ResultIndex(self.root,self.path,digest).page(0,1,'','score',True,True)
        with patch('biomolexplorer.result_index.csv.DictReader',side_effect=AssertionError('CSV parsed again')):
            index=ResultIndex(self.root,self.path,lambda stream: self.fail('CSV hashed again'))
            second=index.page(1,1,'','score',True,True)
            row=index.row(second['visible'][0][0],first['version'])
        self.assertEqual(row,second['visible'][0][1]);self.assertEqual(row['molecule_chembl_id'],'B')

    def test_same_size_edit_restored_mtime_still_invalidates_index_and_old_preview(self):
        index=ResultIndex(self.root,self.path,digest)
        first=index.page(0,1,'','score',True,True);stamp=self.path.stat()
        self.path.write_bytes(self.path.read_bytes().replace(b'CCO',b'CCC'))
        os.utime(self.path,ns=(stamp.st_atime_ns,stamp.st_mtime_ns))
        second=index.page(0,1,'','score',True,True)
        self.assertNotEqual(first['version'],second['version'])
        with self.assertRaisesRegex(ValueError,'Atualize'):index.row(0,first['version'])

    def test_corrupt_disposable_database_is_rebuilt_from_authoritative_csv(self):
        index=ResultIndex(self.root,self.path,digest);first=index.page(0,25,'',None,True,True)
        index.database.write_bytes(b'corrupt cache')
        self.assertEqual(index.page(0,25,'',None,True,True)['visible'],first['visible'])

    def test_source_change_during_index_build_never_publishes_a_partial_index(self):
        def changed(stream):
            version=digest(stream)
            self.path.write_bytes(self.path.read_bytes().replace(b'CCO',b'CCC'))
            return version
        with self.assertRaisesRegex(ValueError,'durante a indexação'):
            ResultIndex(self.root,self.path,changed).page(0,25,'',None,True,True)
        result=ResultIndex(self.root,self.path,digest).page(0,25,'',None,True,True)
        self.assertEqual(result['visible'][0][1]['canonical_smiles'],'CCC')

    def test_query_cache_is_bounded_and_retains_the_current_view(self):
        import sqlite3
        index=ResultIndex(self.root,self.path,digest)
        for number in range(12):index.page(0,25,str(number),None,True,False)
        result=index.page(0,25,'cco',None,True,True)
        self.assertEqual(result['matched'],1)
        self.assertEqual(result['visible'][0][0],0)
        with sqlite3.connect(index.database) as db:
            self.assertLessEqual(db.execute('SELECT count(*) FROM views').fetchone()[0],8)

    def test_incremental_projection_preserves_unicode_nested_batches_and_shared_origins(self):
        inputs={'source':[dict(path='/source/é.csv',note='\\"\n')]*3}
        stage=dict(id='stage',configuration={'name':'múltiplos alvos'},artifact_records=[
            dict(path=str(i),inputs=inputs,dependency_manifest={'large':'x'*70000},naming={'label':'á'}) for i in range(4)],
            artifact_manifest={'é':'1'},processing_runtime={'large':'x'*90000},execution_manifest={},batches=[{'artifact_records':[]}])
        payload=json.dumps([stage],ensure_ascii=False).encode()
        class Chunks(io.BytesIO):
            def read(self,size=-1):return super().read(min(size,127) if size>=0 else 127)
        reader=ProjectionReader(Chunks(payload))
        self.assertEqual(reader.array(reader.stage),[compact_stage(stage)])
        self.assertEqual(len(compact_stage(stage)['_review_input_groups']),1)


class AuthorizedIndexTests(unittest.TestCase):
    setUp=fixtures.CompoundTablesTests.setUp
    tearDown=fixtures.CompoundTablesTests.tearDown
    page=fixtures.CompoundTablesTests.page

    def test_indexed_response_matches_existing_pager_and_rechecks_permissions(self):
        from biomolexplorer.compound_tables import CompoundTables
        from biomolexplorer.workspace import AccessDenied
        expected=self.page(query='pubchem',sort_by='source')
        with patch('biomolexplorer.result_index.MIN_INDEX_BYTES',0):
            actual=self.page(query='pubchem',sort_by='source')
            self.assertEqual(actual,expected)
            self.tables=CompoundTables(self.store)
            with patch.object(self.tables,'_digest',side_effect=AssertionError('Persistent index was not reused')):
                self.assertEqual(self.page(query='pubchem',sort_by='source'),expected)
                with self.assertRaises(AccessDenied):self.page(self.guest)

    def test_compact_catalog_is_invalidated_even_when_run_timestamp_does_not_change(self):
        first=self.store.get_result_run(self.owner,self.rid)
        with self.store.connect() as db:
            source=json.loads(db.execute('SELECT stages FROM runs WHERE id=?',(self.rid,)).fetchone()[0])
            source[0]['artifacts']=[]
            db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(source),self.rid))
        second=self.store.get_result_run(self.owner,self.rid)
        self.assertTrue(first['stages'][0]['artifacts']);self.assertEqual(second['stages'][0]['artifacts'],[])
        self.assertEqual(first['updated'],second['updated'])

    def test_origins_load_only_the_requested_ancestry(self):
        from biomolexplorer.result_origins import ResultOrigins
        from biomolexplorer.artifacts import records
        source=self.store.get_result_run(self.owner,self.rid)['stages'][0]
        source=dict(source,artifacts=[*source['artifacts'],*[str(self.path.parent/(str(i)+'.csv')) for i in range(1000)]])
        with patch('biomolexplorer.result_origins.records',wraps=records) as reader:
            origins=ResultOrigins([source],self.rid)
            self.assertEqual(reader.call_count,0)
            origins.describe(str(self.path),self.sid,self.rid)
            self.assertEqual(reader.call_count,1)

    def test_cold_review_inside_a_write_transaction_does_not_wait_for_its_own_lock(self):
        with self.store.connect() as db:
            db.execute('BEGIN IMMEDIATE')
            result=self.store.get_result_run(self.owner,self.rid)
            self.assertEqual(result['id'],self.rid)

    def test_review_catalog_preserves_verified_result_selection_and_source_hashes(self):
        self.store.save_pipeline(self.owner,self.pid,[self.stage],0)
        original=self.path.read_bytes()
        full=self.store.available_results(self.owner,self.pid,include_checkpoints=True)
        review=self.store.available_results(self.owner,self.pid,include_checkpoints=True,review=True)
        self.assertTrue(full)
        self.assertEqual([(item['id'],item['artifacts'],item['artifact_manifest']) for item in full],
                         [(item['id'],item['artifacts'],item['artifact_manifest']) for item in review])
        self.assertEqual(self.path.read_bytes(),original)

    def test_progress_only_updates_use_an_overlay_and_remove_obsolete_status_fields(self):
        from biomolexplorer.result_catalog import store_snapshot,read_stage,read_proofs
        with self.store.connect() as db:
            stages=json.loads(db.execute('SELECT stages FROM runs WHERE id=?',(self.rid,)).fetchone()[0])
            stages[0]['progress']={'completed':1}
            stages[0]['artifact_records']=[dict(path=str(self.path),inputs={'source':[{'path':str(self.path)}]})]
            db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),self.rid))
        self.store.get_result_run(self.owner,self.rid)
        stages[0].pop('progress');stages[0]['error']='updated status'
        with patch('biomolexplorer.result_catalog.RecordPacker.add',side_effect=AssertionError('Descriptors repacked')):
            with self.store.connect() as db:
                db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps(stages),self.rid))
                store_snapshot(db,self.rid,stages,unchanged=True)
        versions=[self.store.get_result_run(self.owner,self.rid)['stages'][0],
                  read_stage(self.store,self.rid,self.sid),read_proofs(self.store,self.rid)[0]]
        for stage in versions:
            self.assertNotIn('progress',stage);self.assertEqual(stage['error'],'updated status')


class ExecutionWriteTests(unittest.TestCase):
    """What a progress tick costs, and when the state becomes visible."""
    setUp=fixtures.CompoundTablesTests.setUp
    tearDown=fixtures.CompoundTablesTests.tearDown

    def service(self):
        from biomolexplorer.pipeline import PipelineService
        service=PipelineService(self.store);self.addCleanup(service.close);return service

    def column(self):
        with self.store.connect() as db:
            return db.execute('SELECT stages FROM runs WHERE id=?',(self.rid,)).fetchone()[0]

    def revision(self):
        with self.store.connect() as db:
            return db.execute('SELECT revision FROM result_read_versions WHERE run_id=?',(self.rid,)).fetchone()[0]

    def running(self,service):
        """Put the run in the state a tick arrives in, and index it once."""
        stages=copy.deepcopy(self.store.get_result_run(self.owner,self.rid)['stages'])
        for stage in stages:stage['status']='running'
        service._update(self.rid,'running',stages)
        return stages

    def test_a_progress_tick_never_rewrites_the_execution_column(self):
        service=self.service()
        stages=self.running(service)
        before,revision=self.column(),self.revision()
        stages[0]['progress']={'completed':7,'total':9,'message':'Consultando o ChEMBL'}
        service._run_update_times={}                     # the two-second throttle is not under test
        service._update(self.rid,'running',stages)
        # The row holds a descriptor of every produced file: SQLite would
        # rewrite all of it to advance a counter.
        self.assertEqual(self.column(),before)
        self.assertEqual(self.revision(),revision+1)
        published=self.store.get_result_run(self.owner,self.rid)
        self.assertEqual(published['stages'][0]['progress'],{'completed':7,'total':9,'message':'Consultando o ChEMBL'})
        self.assertEqual(self.store.review_run(self.owner,self.rid)['stages'][0]['progress'],
                         published['stages'][0]['progress'])

    def test_a_state_change_is_always_persisted_in_the_execution_column(self):
        service=self.service()
        stages=self.running(service)
        for change in ({'log_path':'/tmp/run.log'},{'job_id':'job-1'},{'started_at':1.5}):
            with self.subTest(change=change):
                before=self.column()
                stages[0].update(change)
                service._run_update_times={}
                service._update(self.rid,'running',stages)
                self.assertNotEqual(self.column(),before)
                stored=json.loads(self.column())[0]
                for key,value in change.items():self.assertEqual(stored[key],value)
        stages[0].update(status='succeeded',finished_at=2.0)
        service._run_update_times={}
        service._update(self.rid,'succeeded',stages)
        with self.store.connect() as db:
            row=db.execute('SELECT status,stages FROM runs WHERE id=?',(self.rid,)).fetchone()
        self.assertEqual(row['status'],'succeeded')
        self.assertEqual(json.loads(row['stages'])[0]['status'],'succeeded')

    def test_the_final_state_is_committed_before_the_history_is_recorded(self):
        service=self.service()
        stages=self.running(service)
        seen=[]
        def audit(store,db,project_id,user_id,action,summary,before=None,reuse_snapshot=False,**options):
            with self.store.connect() as other:
                seen.append(other.execute('SELECT status FROM runs WHERE id=?',(self.rid,)).fetchone()[0])
        with patch('biomolexplorer.project_state.record',audit):
            service._update(self.rid,'succeeded',stages)
        # A snapshot walks and hashes the whole project; the run must not be
        # reported as running while that happens.
        self.assertEqual(seen,['succeeded'])

    def test_a_failed_audit_never_fails_a_finished_execution(self):
        service=self.service()
        stages=self.running(service)
        with patch('biomolexplorer.project_state.record',side_effect=OSError('disco cheio')):
            service._update(self.rid,'succeeded',stages)
        self.assertEqual(self.store.get_result_run(self.owner,self.rid)['status'],'succeeded')

    def test_following_a_run_reads_its_state_without_the_file_descriptors(self):
        from biomolexplorer.workspace import AccessDenied
        service=self.service()
        # Index the run, so its files really have descriptors to leave out.
        service._update(self.rid,'succeeded',copy.deepcopy(self.store.get_result_run(self.owner,self.rid)['stages']))
        full=self.store.get_result_run(self.owner,self.rid)
        review=self.store.review_run(self.owner,self.rid)
        self.assertEqual((review['id'],review['status'],review['project_id']),
                         (full['id'],full['status'],full['project_id']))
        self.assertEqual([(s['id'],s['status'],s['artifacts']) for s in review['stages']],
                         [(s['id'],s['status'],s['artifacts']) for s in full['stages']])
        self.assertTrue(all(not stage.get('artifact_records') for stage in review['stages']))
        self.assertTrue(any(stage.get('artifact_records') for stage in full['stages']))
        self.assertEqual([run['id'] for run in self.store.list_review_runs(self.owner,self.pid)],[self.rid])
        with self.assertRaises(AccessDenied):self.store.review_run(self.guest,self.rid)
        with self.assertRaises(AccessDenied):self.store.list_review_runs(self.guest,self.pid)
