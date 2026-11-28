"""One content store, the artifact manager and transactional publication.

PHASE-02 (#9): AC-02-01, AC-02-02, AC-02-03 and AC-02-04.
"""
import json
import os
import sqlite3
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from biomolexplorer import artifact_catalog as catalog
from biomolexplorer.artifact_manager import ArtifactManager,STORE_FOLDER
from biomolexplorer.artifact_registry import register_run
from biomolexplorer.artifact_validation import InvalidContent,validate
from biomolexplorer.artifacts import index_run
from biomolexplorer.catalog import new_stage
from biomolexplorer.content_store import ContentStore,OBJECT_MODE
from biomolexplorer.stage_cache import artifact_manifest,file_digest
from biomolexplorer.workspace import WorkspaceStore

COMPOUNDS='molecule_chembl_id,canonical_smiles\nCHEMBL1,CCO\nCHEMBL2,CCC\n'
POSE=('REMARK VINA RESULT: -7.0 0.0 0.0\nROOT\n'
      'ATOM      1  C   UNL A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\n'
      'ENDROOT\nTORSDOF 0\n')


def unlock(path):
    """A stored object is read-only; only its owner can unlock it."""
    os.chmod(path,0o600)


class StoreTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.store=ContentStore(self.root/'store')

    def _source(self,name,content):
        path=self.root/name;path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(content);return path

    def test_identical_content_from_two_producers_is_stored_once(self):
        """AC-02-01: managed storage holds one copy of identical bytes."""
        first=self._source('one/compounds.csv',COMPOUNDS)
        second=self._source('two/compounds.csv',COMPOUNDS)
        digest=self.store.put(first)
        self.assertEqual(self.store.put(second),digest)
        self.assertEqual(self.store.usage(),{'objects':1,'bytes':len(COMPOUNDS.encode())})
        self.assertEqual(self.store.stored(),{digest})

    def test_a_stored_object_is_read_only_and_verifies(self):
        path=self._source('a.csv',COMPOUNDS)
        digest=self.store.put(path)
        stored=self.store.object_path(digest)
        self.assertEqual(stored.stat().st_mode&0o777,OBJECT_MODE)
        with self.assertRaises(PermissionError):stored.open('a')
        self.assertTrue(self.store.verify(digest))
        self.assertEqual(self.store.size(digest),len(COMPOUNDS.encode()))

    def test_publication_leaves_no_partial_object_when_content_changes(self):
        path=self._source('a.csv',COMPOUNDS)
        digest=file_digest(path)
        path.write_text(COMPOUNDS+'CHEMBL3,CCCC\n')
        with self.assertRaises(ValueError):self.store.put(path,digest=digest)
        self.assertFalse(self.store.has(digest))
        self.assertEqual(self.store.stored(),set())

    def test_content_failing_its_contract_is_never_published(self):
        path=self._source('scores.csv','unrelated\n1\n')
        def refuse(_):raise InvalidContent('não satisfaz o contrato')
        with self.assertRaises(InvalidContent):self.store.put(path,validate=refuse)
        self.assertEqual(self.store.stored(),set())

    def test_a_shared_materialization_costs_no_bytes_and_refuses_writes(self):
        path=self._source('receptor.pdbqt',POSE)
        digest=self.store.put(path)
        target=self.root/'job'/'receptor.pdbqt'
        self.store.materialize(digest,target,'shared')
        self.assertEqual(target.stat().st_ino,self.store.object_path(digest).stat().st_ino)
        with self.assertRaises(PermissionError):target.open('a')

    def test_a_writable_materialization_cannot_reach_the_stored_object(self):
        """AC-02-02: execution inputs are isolated."""
        path=self._source('ligand.pdbqt',POSE)
        digest=self.store.put(path)
        target=self.root/'job'/'ligand.pdbqt'
        self.store.materialize(digest,target,'isolated')
        self.assertNotEqual(target.stat().st_ino,self.store.object_path(digest).stat().st_ino)
        target.write_text('rewritten by the engine')
        self.assertTrue(self.store.verify(digest))
        self.assertEqual(self.store.object_path(digest).read_text(),POSE)

    def test_a_tampered_object_is_refused_on_read_and_on_materialization(self):
        path=self._source('a.csv',COMPOUNDS)
        digest=self.store.put(path)
        stored=self.store.object_path(digest);unlock(stored)
        stored.write_text(COMPOUNDS.replace('CCO','CCCC'))
        os.chmod(stored,OBJECT_MODE)
        self.assertFalse(self.store.verify(digest))
        with self.assertRaises(ValueError):self.store.open(digest)
        with self.assertRaises(ValueError):self.store.materialize(digest,self.root/'out.csv')
        self.assertFalse((self.root/'out.csv').exists())

    def test_concurrent_publication_of_the_same_content_is_safe(self):
        sources=[self._source('p%d/a.csv'%index,COMPOUNDS) for index in range(8)]
        digests=[];errors=[]
        def publish(source):
            try:digests.append(self.store.put(source))
            except Exception as error:errors.append(error)
        threads=[threading.Thread(target=publish,args=(source,)) for source in sources]
        for thread in threads:thread.start()
        for thread in threads:thread.join()
        self.assertEqual(errors,[])
        self.assertEqual(len(set(digests)),1)
        self.assertEqual(self.store.usage()['objects'],1)
        self.assertTrue(self.store.verify(digests[0]))

    def test_the_sweep_removes_only_what_nothing_references(self):
        kept=self.store.put(self._source('a.csv',COMPOUNDS))
        orphan=self.store.put(self._source('b.pdbqt',POSE))
        temporary=self.store.object_path(orphan).parent/'.publishing-left-behind'
        temporary.write_text('interrupted')
        self.assertEqual(self.store.sweep({kept}),1)
        self.assertEqual(self.store.stored(),{kept})
        self.assertFalse(temporary.exists())


class ValidationContractTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)

    def _file(self,name,content):
        path=self.root/name;path.write_text(content);return path

    def test_a_table_must_declare_the_columns_of_its_type(self):
        good=self._file('compounds.csv',COMPOUNDS)
        self.assertEqual(validate('compounds',good),good)
        bad=self._file('wrong.csv','molecule_chembl_id\nCHEMBL1\n')
        with self.assertRaises(InvalidContent):validate('compounds',bad)

    def test_a_pose_must_contain_atoms_with_finite_coordinates(self):
        good=self._file('pose.pdbqt',POSE)
        self.assertEqual(validate('molecular_poses',good),good)
        empty=self._file('empty.pdbqt','REMARK only a comment\n')
        with self.assertRaises(InvalidContent):validate('molecular_poses',empty)
        broken=self._file('nan.pdbqt',POSE.replace('0.000   0.000   0.000','  nan     nan     nan'))
        with self.assertRaises(InvalidContent):validate('molecular_poses',broken)

    def test_a_prepared_input_goes_through_the_engine_validator(self):
        good=self._file('ligand.pdbqt',POSE)
        self.assertEqual(validate('prepared_structures',good),good)
        salts=self._file('salts.pdbqt',POSE.replace('ENDROOT','ENDROOT\nROOT'))
        with self.assertRaises(InvalidContent):validate('prepared_structures',salts)

    def test_an_untyped_import_is_checked_for_presence_only(self):
        path=self._file('anything.bin','payload')
        self.assertEqual(validate('other',path),path)
        empty=self._file('empty.bin','')
        with self.assertRaises(InvalidContent):validate('other',empty)

    def test_a_zinc_list_must_parse(self):
        good=self._file('zinc.txt','https://zinc.docking.org/db/bysubset/BA/BAAA.smi\n')
        self.assertEqual(validate('zinc_urls',good),good)
        bad=self._file('bad.txt','not a link at all\n')
        with self.assertRaises(InvalidContent):validate('zinc_urls',bad)


class ManagerFixture(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.workspace=WorkspaceStore(self.tmp.name)
        self.token=self.workspace.register('Owner','store@example.org','owner-password')
        self.pid=self.workspace.create_project(self.token,'Alzheimer store')['id']
        self.root=self.workspace.project_dir(self.pid)
        self.manager=ArtifactManager(self.workspace,self.pid)

    def _run(self,run_id,*,content=COMPOUNDS,name='compounds.csv',operation='retrieve_compounds'):
        stage=new_stage(operation)
        folder=self.root/'runs'/run_id/stage['id'];folder.mkdir(parents=True,exist_ok=True)
        path=folder/name;path.write_text(content)
        stage.update(status='succeeded',artifacts=[str(path)],started_at=1.0,finished_at=2.0,
            artifact_manifest=artifact_manifest([str(path)]),processing_version='digest-'+operation,
            producer_software={'tool':'1.0'},
            configuration={'operation':operation,'parameters':stage['parameters']})
        stages=[stage];index_run(stages,run_id)
        with self.workspace.connect() as db:
            db.execute('INSERT OR REPLACE INTO runs VALUES (?,?,?,?,?,?,?,NULL)',
                (run_id,self.pid,None,'succeeded',json.dumps(stages),time.time(),time.time()))
            register_run(db,self.pid,run_id,stages)
        return stage,path,stage['artifact_records'][0]['occurrence_id']


class ArtifactManagerTests(ManagerFixture):
    def test_registering_fills_the_content_object_left_open_by_phase_one(self):
        stage,path,occurrence=self._run('run-one')
        with self.workspace.connect() as db:
            before=catalog.content_object(db,file_digest(path))
        self.assertIsNone(before['storage_ref'])
        digest=self.manager.register(occurrence)
        with self.workspace.connect() as db:
            after=catalog.content_object(db,digest)
        self.assertEqual(after['storage_ref'],'/'.join((*STORE_FOLDER,'objects',digest[:2],digest)))
        self.assertEqual(after['size'],len(COMPOUNDS.encode()))
        self.assertTrue(self.manager.store.verify(digest))

    def test_two_runs_publishing_the_same_bytes_share_one_managed_object(self):
        """AC-02-01, through the catalog: distinct identities, one object."""
        _,_,first=self._run('run-a')
        _,_,second=self._run('run-b')
        self.assertNotEqual(first,second)
        self.assertEqual(self.manager.register(first),self.manager.register(second))
        self.assertEqual(self.manager.store.usage()['objects'],1)
        with self.workspace.connect() as db:
            shared=catalog.artifacts_by_content(db,self.pid,self.manager.register(first))
        self.assertEqual(len({row['occurrence_id'] for row in shared}),2)

    def test_registering_refuses_content_that_fails_its_contract(self):
        stage,path,occurrence=self._run('run-bad',content='unrelated\n1\n')
        with self.workspace.connect() as db:
            db.execute('''UPDATE artifact_instances SET artifact_type='compounds'
                WHERE project_id=? AND occurrence_id=?''',(self.pid,occurrence))
        with self.assertRaises(InvalidContent):self.manager.register(occurrence)
        self.assertEqual(self.manager.store.stored(),set())

    def test_a_batch_reports_the_refused_artifact_and_publishes_the_rest(self):
        _,good,first=self._run('run-good')
        _,bad,second=self._run('run-bad-batch',content='unrelated\n1\n',name='table.csv')
        with self.workspace.connect() as db:
            db.execute('''UPDATE artifact_instances SET artifact_type='compounds'
                WHERE project_id=? AND occurrence_id=?''',(self.pid,second))
        result=self.manager.publish([first,second])
        self.assertEqual(list(result['registered']),[first])
        self.assertIn(second,result['refused'])
        self.assertIn('colunas',result['refused'][second])
        self.assertEqual(self.manager.store.stored(),{file_digest(good)})

    def test_resolve_returns_the_managed_object_and_detects_tampering(self):
        stage,path,occurrence=self._run('run-resolve')
        digest=self.manager.register(occurrence)
        self.assertEqual(self.manager.resolve(occurrence),self.manager.store.object_path(digest))
        stored=self.manager.store.object_path(digest);unlock(stored)
        stored.write_text('tampered');os.chmod(stored,OBJECT_MODE)
        with self.assertRaises(ValueError):self.manager.resolve(occurrence)

    def test_resolve_falls_back_to_the_published_file_and_checks_its_digest(self):
        stage,path,occurrence=self._run('run-origin')
        self.assertEqual(self.manager.resolve(occurrence),path)
        path.write_text(COMPOUNDS+'CHEMBL9,CC\n')
        with self.assertRaises(ValueError):self.manager.resolve(occurrence)

    def test_materialization_is_shared_by_default_and_isolated_on_request(self):
        """AC-02-02 through the manager."""
        stage,path,occurrence=self._run('run-materialize')
        digest=self.manager.register(occurrence)
        shared=self.manager.materialize(occurrence,self.root/'jobs'/'shared.csv')
        self.assertEqual(shared.stat().st_ino,self.manager.store.object_path(digest).stat().st_ino)
        isolated=self.manager.materialize(occurrence,self.root/'jobs'/'own.csv',writable=True)
        self.assertNotEqual(isolated.stat().st_ino,self.manager.store.object_path(digest).stat().st_ino)
        isolated.write_text('rewritten')
        self.assertTrue(self.manager.store.verify(digest))

    def test_provenance_reads_the_producer_and_its_inputs_from_the_catalog(self):
        """AC-02-03: relationships stay traceable."""
        stage,path,occurrence=self._run('run-provenance')
        record=self.manager.provenance(occurrence)
        self.assertEqual(record['artifact_type'],stage['artifact_records'][0]['artifact_type'])
        self.assertEqual(record['operation']['run_id'],'run-provenance')
        self.assertEqual(record['operation']['operation_type'],'retrieve_compounds')
        self.assertTrue(record['operation']['signature'])
        self.assertEqual({entity['entity_id'] for entity in record['entities']},{'CHEMBL1','CHEMBL2'})
        self.assertFalse(record['managed'])
        self.manager.register(occurrence)
        self.assertTrue(self.manager.provenance(occurrence)['managed'])

    def test_usage_reports_the_store_and_what_the_catalog_claims(self):
        """AC-02-04: the storage impact is measurable."""
        _,_,first=self._run('run-usage-a')
        _,_,second=self._run('run-usage-b')
        self.manager.publish([first,second])
        usage=self.manager.usage()
        self.assertEqual(usage['store']['objects'],1)
        self.assertEqual(usage['store']['bytes'],len(COMPOUNDS.encode()))
        self.assertEqual(usage['registered']['objects'],1)
        self.assertEqual(usage['logical_references'],2)


class TransactionalPublicationTests(ManagerFixture):
    def _journals(self):
        folder=self.root/'.history'/'.result-transactions'
        return sorted(folder.glob('*.json')) if folder.is_dir() else []

    def test_a_complete_publication_leaves_no_journal_and_no_marker(self):
        _,_,occurrence=self._run('run-clean')
        self.manager.publish([occurrence])
        self.assertEqual(self._journals(),[])
        with self.workspace.connect() as db:
            self.assertEqual(db.execute('SELECT COUNT(*) FROM result_file_transactions').fetchone()[0],0)
            self.assertIsNotNone(catalog.content_object(db,file_digest(
                Path(db.execute('SELECT storage_path FROM artifact_instances WHERE occurrence_id=?',
                    (occurrence,)).fetchone()[0])))['storage_ref'])

    def test_an_interruption_before_the_commit_leaves_no_dangling_reference(self):
        stage,path,occurrence=self._run('run-crash-before')
        digest=file_digest(path)
        class Interrupted(Exception):pass
        # The bytes reach the store, and the process dies before the rows commit.
        with patch.object(ArtifactManager,'_record',side_effect=Interrupted):
            with self.assertRaises(Interrupted):self.manager.publish([occurrence])
        self.assertTrue(self.manager.store.has(digest))
        with self.workspace.connect() as db:
            self.assertIsNone(catalog.content_object(db,digest)['storage_ref'])
            self.assertEqual(len(self._journals()),1)
        # Reopening recovers: the journal is resolved and the catalog is intact.
        WorkspaceStore(self.tmp.name)
        self.assertEqual(self._journals(),[])
        with self.workspace.connect() as db:
            rows=[dict(row) for row in db.execute('''SELECT content_sha256 FROM artifact_instances
                WHERE project_id=?''',(self.pid,))]
        self.assertEqual([row['content_sha256'] for row in rows],[digest])
        # The content is still the published artifact's, so recovery keeps it
        # and the next publication deduplicates onto it.
        self.manager.publish([occurrence])
        self.assertEqual(self.manager.store.usage()['objects'],1)
        with self.workspace.connect() as db:
            self.assertIsNotNone(catalog.content_object(db,digest)['storage_ref'])

    def test_an_interruption_before_the_commit_discards_unclaimed_content(self):
        stage,path,occurrence=self._run('run-crash-orphan')
        digest=file_digest(path)
        class Interrupted(Exception):pass
        with patch.object(ArtifactManager,'_record',side_effect=Interrupted):
            with self.assertRaises(Interrupted):self.manager.publish([occurrence])
        self.assertTrue(self.manager.store.has(digest))
        # The run that published those bytes is cleared, so nothing claims them.
        with self.workspace.connect() as db:
            db.execute('DELETE FROM artifact_instances WHERE project_id=?',(self.pid,))
        WorkspaceStore(self.tmp.name)
        self.assertFalse(self.manager.store.has(digest))
        self.assertEqual(self._journals(),[])

    def test_an_interruption_after_the_commit_keeps_the_published_content(self):
        stage,path,occurrence=self._run('run-crash-after')
        digest=file_digest(path)
        from biomolexplorer.result_transactions import ResultTransaction
        with patch.object(ResultTransaction,'finish',side_effect=AssertionError('morreu depois do commit')):
            with self.assertRaises(AssertionError):self.manager.publish([occurrence])
        self.assertEqual(len(self._journals()),1)
        WorkspaceStore(self.tmp.name)
        self.assertTrue(self.manager.store.verify(digest))
        self.assertEqual(self._journals(),[])
        with self.workspace.connect() as db:
            self.assertIsNotNone(catalog.content_object(db,digest)['storage_ref'])

    def test_a_rolled_back_catalog_transaction_leaves_the_store_consistent(self):
        stage,path,occurrence=self._run('run-rollback')
        digest=file_digest(path)
        with patch.object(ArtifactManager,'_record',
                          side_effect=sqlite3.OperationalError('disco cheio')):
            with self.assertRaises(sqlite3.OperationalError):self.manager.publish([occurrence])
        with self.workspace.connect() as db:
            self.assertIsNone(catalog.content_object(db,digest)['storage_ref'])
        self.assertTrue(self.manager.store.verify(digest))


if __name__=='__main__':unittest.main()
