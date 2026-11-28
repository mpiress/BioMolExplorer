"""Content objects, operation identity, dependency edges and catalog migrations.

PHASE-01 (#8): AC-01-01, AC-01-02 and AC-01-03.  The queries must answer from
the catalog alone, so the filesystem is made to fail while they run.
"""
import json
import os
import sqlite3
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from biomolexplorer import artifact_catalog as catalog
from biomolexplorer.artifact_registry import register_run
from biomolexplorer.artifacts import index_run
from biomolexplorer.catalog import new_stage
from biomolexplorer.catalog_schema import MIGRATIONS,apply,pending,schema_version
from biomolexplorer.stage_cache import artifact_manifest
from biomolexplorer.workspace import WorkspaceStore

COMPOUNDS='molecule_chembl_id,canonical_smiles\nCHEMBL1,CCO\nCHEMBL2,CCC\n'


class CatalogFixture(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.store=WorkspaceStore(self.tmp.name)
        self.token=self.store.register('Owner','catalog@example.org','owner-password')
        self.project=self.store.create_project(self.token,'Alzheimer catalog');self.pid=self.project['id']
        self.root=self.store.project_dir(self.pid)

    def tearDown(self):self.tmp.cleanup()

    def _pipeline(self,run_id,*,compounds=COMPOUNDS):
        """A retrieval feeding a fingerprint block, with real files and lineage."""
        producer=new_stage('retrieve_compounds');consumer=new_stage('fingerprints')
        folder=self.root/'runs'/run_id;folder.mkdir(parents=True,exist_ok=True)
        source=folder/producer['id']/'compounds.csv';source.parent.mkdir(parents=True)
        source.write_text(compounds)
        derived=folder/consumer['id']/'fingerprints.csv';derived.parent.mkdir(parents=True)
        derived.write_text('molecule_chembl_id,morgan\nCHEMBL1,0101\n')
        producer.update(status='succeeded',artifacts=[str(source)],started_at=10.0,finished_at=20.0,
            artifact_manifest=artifact_manifest([str(source)]),processing_version='retrieval-digest',
            producer_software={'chembl':'34'},configuration={'operation':'retrieve_compounds','parameters':producer['parameters']})
        consumer.update(status='succeeded',artifacts=[str(derived)],started_at=20.0,finished_at=30.0,
            artifact_manifest=artifact_manifest([str(derived)]),processing_version='fingerprint-digest',
            input_files={'base_input_path':[str(source)]},
            input_manifest={'stage:'+producer['id']:[['compounds.csv',artifact_manifest([str(source)])[str(source.resolve())]]]},
            configuration={'operation':'fingerprints','parameters':consumer['parameters'],
                'bindings':{'base_input_path':{'stage':producer['id'],'selector':'compounds.csv'}}})
        consumer['bindings']={'base_input_path':{'stage':producer['id'],'selector':'compounds.csv'}}
        stages=[producer,consumer]
        index_run(stages,run_id)
        return stages,source,derived

    def _register(self,run_id,stages):
        with self.store.connect() as db:
            db.execute('INSERT OR REPLACE INTO runs VALUES (?,?,?,?,?,?,?,NULL)',
                (run_id,self.pid,None,'succeeded',json.dumps(stages),time.time(),time.time()))
            register_run(db,self.pid,run_id,stages)


class ContentIdentityTests(CatalogFixture):
    def test_two_logical_artifacts_reference_one_content_object(self):
        """AC-01-01: identical bytes published twice are one physical content."""
        stages,source,derived=self._pipeline('run-one')
        twin=self.root/'runs'/'run-one'/stages[1]['id']/'copy.csv';twin.write_text(COMPOUNDS)
        stages[1]['artifacts'].append(str(twin))
        stages[1]['artifact_manifest']=artifact_manifest(stages[1]['artifacts'])
        index_run(stages,'run-one')
        self._register('run-one',stages)
        sha=artifact_manifest([str(source)])[str(source.resolve())]
        with self.store.connect() as db:
            shared=catalog.artifacts_by_content(db,self.pid,sha)
            self.assertEqual({row['storage_path'] for row in shared},{str(source),str(twin)})
            self.assertEqual(len({row['occurrence_id'] for row in shared}),2)
            self.assertEqual(db.execute('SELECT COUNT(*) FROM content_objects WHERE sha256=?',(sha,)).fetchone()[0],1)
            stored=catalog.content_object(db,sha)
        self.assertEqual(stored['size'],len(COMPOUNDS.encode()))
        self.assertEqual(stored['media_type'],'text/csv')
        self.assertIsNone(stored['storage_ref'])

    def test_an_unavailable_file_has_no_content_size(self):
        stages,source,derived=self._pipeline('run-missing')
        self._register('run-missing',stages)
        sha=artifact_manifest([str(source)])[str(source.resolve())]
        source.unlink()
        self._register('run-missing',stages)
        with self.store.connect() as db:
            row=db.execute('SELECT state FROM artifact_instances WHERE project_id=? AND storage_path=?',
                (self.pid,str(source))).fetchone()
            self.assertEqual(row['state'],'unavailable')
            # The object keeps the size measured while the bytes were there.
            self.assertEqual(catalog.content_object(db,sha)['size'],len(COMPOUNDS.encode()))


class OperationIdentityTests(CatalogFixture):
    def test_each_block_execution_becomes_one_operation_with_its_edges(self):
        stages,source,derived=self._pipeline('run-ops')
        self._register('run-ops',stages)
        with self.store.connect() as db:
            operations=catalog.operations_of_run(db,self.pid,'run-ops')
            self.assertEqual([row['operation_type'] for row in operations],['retrieve_compounds','fingerprints'])
            retrieval,fingerprints=operations
            self.assertEqual((retrieval['started'],retrieval['finished']),(10.0,20.0))
            self.assertEqual(fingerprints['state'],'completed')
            self.assertTrue(all(row['signature'] and row['parameters_digest'] and row['software_digest']
                                for row in operations))
            self.assertEqual({row['signature_version'] for row in operations},{catalog.SIGNATURE_VERSION})
            edges=catalog.operation_edges(db,fingerprints['operation_id'])
            self.assertEqual([row['port'] for row in edges['inputs']],['base_input_path'])
            self.assertEqual(edges['inputs'][0]['role'],'artifact')
            produced=stages[0]['artifact_records'][0]['occurrence_id']
            self.assertEqual(edges['inputs'][0]['occurrence_id'],produced)
            self.assertEqual([row['occurrence_id'] for row in edges['outputs']],
                             [stages[1]['artifact_records'][0]['occurrence_id']])

    def test_the_signature_separates_parameters_software_and_input_content(self):
        stages,_,_=self._pipeline('run-signature')
        self._register('run-signature',stages)
        with self.store.connect() as db:
            before=catalog.operations_of_run(db,self.pid,'run-signature')[1]
        stages[1]['configuration']['parameters']=dict(stages[1]['configuration']['parameters'],maccs=True)
        self._register('run-signature',stages)
        with self.store.connect() as db:
            after=catalog.operations_of_run(db,self.pid,'run-signature')[1]
        self.assertEqual(before['operation_id'],after['operation_id'])
        self.assertNotEqual(before['parameters_digest'],after['parameters_digest'])
        self.assertNotEqual(before['signature'],after['signature'])
        self.assertEqual(before['software_digest'],after['software_digest'])

    def test_a_reused_result_points_to_the_operation_that_produced_it(self):
        stages,source,derived=self._pipeline('run-first')
        self._register('run-first',stages)
        produced=stages[0]['artifact_records'][0]['occurrence_id']
        reused=dict(stages[0],reused=True,reused_from_run='run-first',status='succeeded')
        reused['artifact_records']=[dict(record) for record in stages[0]['artifact_records']]
        self._register('run-second',[reused])
        with self.store.connect() as db:
            self.assertEqual(catalog.operations_of_run(db,self.pid,'run-second'),[])
            origin=catalog.producing_operation(db,self.pid,produced)
            self.assertEqual(origin['run_id'],'run-first')
            self.assertEqual(origin['operation_id'],
                catalog.operation_id(self.pid,stages[0]['id'],'run-first'))

    def test_clearing_a_producer_keeps_the_lineage_of_a_run_that_reused_it(self):
        stages,source,_=self._pipeline('run-origin')
        self._register('run-origin',stages)
        produced=stages[0]['artifact_records'][0]['occurrence_id']
        reused=dict(stages[0],reused=True,reused_from_run='run-origin',status='succeeded')
        reused['artifact_records']=[dict(record) for record in stages[0]['artifact_records']]
        self._register('run-reuser',[reused])
        # Clearing the results of a run empties its stage.  The run that reused
        # them is re-indexed in an unspecified order, so the operation must not
        # depend on that order to survive.
        cleared=dict(stages[0],artifacts=[],artifact_records=[],artifact_manifest={},status='succeeded')
        self._register('run-origin',[cleared])
        with self.store.connect() as db:
            origin=catalog.producing_operation(db,self.pid,produced)
            self.assertIsNotNone(origin)
            self.assertEqual(origin['run_id'],'run-origin')
            self.assertEqual([row['occurrence_id'] for row in
                catalog.operation_edges(db,origin['operation_id'])['outputs']],[produced])

    def test_forgetting_a_run_drops_the_operations_nothing_claims(self):
        stages,_,_=self._pipeline('run-alone')
        self._register('run-alone',stages)
        cleared=[dict(stage,artifacts=[],artifact_records=[],artifact_manifest={}) for stage in stages]
        self._register('run-alone',cleared)
        with self.store.connect() as db:
            self.assertEqual(catalog.operations_of_run(db,self.pid,'run-alone'),[])
            self.assertEqual(db.execute('SELECT COUNT(*) FROM operation_outputs').fetchone()[0],0)
            self.assertEqual(db.execute('SELECT COUNT(*) FROM operation_inputs').fetchone()[0],0)

    def test_registering_a_run_again_replaces_its_operations_without_duplicates(self):
        stages,_,_=self._pipeline('run-twice')
        self._register('run-twice',stages)
        self._register('run-twice',stages)
        with self.store.connect() as db:
            self.assertEqual(len(catalog.operations_of_run(db,self.pid,'run-twice')),2)
            self.assertEqual(db.execute('SELECT COUNT(*) FROM operation_inputs').fetchone()[0],1)
            self.assertEqual(db.execute('SELECT COUNT(*) FROM operation_outputs').fetchone()[0],2)

    def test_a_published_file_without_a_descriptor_is_recorded_not_guessed(self):
        producer=new_stage('retrieve_compounds')
        path=self.root/'orphan.csv';path.write_text(COMPOUNDS)
        producer.update(status='succeeded',artifacts=[str(path)],artifact_records=[])
        self._register('run-orphan',[producer])
        with self.store.connect() as db:
            issues=[dict(row) for row in db.execute('SELECT * FROM catalog_backfill_issues')]
            self.assertEqual([row['reason'] for row in issues],['missing_artifact_records'])
            self.assertEqual(issues[0]['stage_id'],producer['id'])
            self.assertEqual(catalog.operations_of_run(db,self.pid,'run-orphan'),[])

    def test_a_run_without_software_identity_has_no_signature(self):
        stages,_,_=self._pipeline('run-legacy')
        for stage in stages:
            stage.pop('processing_version',None);stage.pop('producer_software',None)
            for record in stage['artifact_records']:
                record.update(processing_version=None,producer_software=None,runtime_id=None)
        self._register('run-legacy',stages)
        with self.store.connect() as db:
            operations=catalog.operations_of_run(db,self.pid,'run-legacy')
            self.assertTrue(operations)
            self.assertEqual({row['signature'] for row in operations},{None})
            self.assertTrue(all(row['parameters_digest'] for row in operations))
            reasons={row['reason'] for row in db.execute('SELECT reason FROM catalog_backfill_issues')}
            self.assertEqual(reasons,{'missing_software_identity'})


class CatalogQueryTests(CatalogFixture):
    def test_entities_resolve_to_artifacts_through_the_index(self):
        """AC-01-02: a scientific entity reaches the artifacts that published it."""
        stages,source,derived=self._pipeline('run-entities')
        self._register('run-entities',stages)
        with self.store.connect() as db:
            found=catalog.artifacts_by_entity(db,self.pid,'compound','CHEMBL1')
            # One compound reaches every artifact that published it: the
            # retrieval that found it and the fingerprint table that lists it.
            self.assertEqual({row['storage_path'] for row in found},{str(source),str(derived)})
            self.assertEqual({row['occurrence_id'] for row in found},
                {stage['artifact_records'][0]['occurrence_id'] for stage in stages})
            self.assertEqual([row['storage_path'] for row in
                catalog.artifacts_by_entity(db,self.pid,'compound','CHEMBL2')],[str(source)])
            self.assertEqual(catalog.artifacts_by_entity(db,self.pid,'compound','CHEMBL404'),[])

    def test_producer_and_lineage_queries_answer_from_the_edges(self):
        stages,source,derived=self._pipeline('run-lineage')
        self._register('run-lineage',stages)
        first=stages[0]['artifact_records'][0]['occurrence_id']
        second=stages[1]['artifact_records'][0]['occurrence_id']
        with self.store.connect() as db:
            self.assertEqual([row['storage_path'] for row in
                catalog.artifacts_by_producer(db,self.pid,stages[0]['id'],'run-lineage')],[str(source)])
            self.assertEqual([row['storage_path'] for row in
                catalog.artifacts_by_producer(db,self.pid,stages[0]['id'],port='compounds')],[str(source)])
            self.assertEqual(catalog.artifacts_by_producer(db,self.pid,stages[0]['id'],port='scores'),[])
            self.assertEqual([row['occurrence_id'] for row in catalog.upstream(db,self.pid,second)],[first])
            self.assertEqual([row['occurrence_id'] for row in catalog.downstream(db,self.pid,first)],[second])
            self.assertEqual(catalog.upstream(db,self.pid,first),[])
            self.assertEqual(catalog.downstream(db,self.pid,second),[])
            self.assertEqual(catalog.upstream(db,self.pid,second,depth=0),[])

    def test_indexed_queries_never_read_a_directory(self):
        """AC-01-03: the catalog answers without a recursive scan."""
        stages,source,_=self._pipeline('run-indexed')
        self._register('run-indexed',stages)
        sha=artifact_manifest([str(source)])[str(source.resolve())]
        first=stages[0]['artifact_records'][0]['occurrence_id']
        second=stages[1]['artifact_records'][0]['occurrence_id']
        def refuse(*args,**kwargs):raise AssertionError('A consulta percorreu o sistema de arquivos.')
        with self.store.connect() as db:
            with patch.object(Path,'glob',refuse),patch.object(Path,'rglob',refuse),patch.object(os,'walk',refuse), \
                 patch.object(Path,'iterdir',refuse),patch.object(os,'listdir',refuse),patch.object(os,'scandir',refuse):
                self.assertTrue(catalog.artifacts_by_entity(db,self.pid,'compound','CHEMBL1'))
                self.assertTrue(catalog.artifacts_by_content(db,self.pid,sha))
                self.assertTrue(catalog.artifacts_by_producer(db,self.pid,stages[0]['id']))
                self.assertTrue(catalog.producing_operation(db,self.pid,first))
                self.assertTrue(catalog.operations_of_run(db,self.pid,'run-indexed'))
                self.assertTrue(catalog.upstream(db,self.pid,second))
                self.assertTrue(catalog.downstream(db,self.pid,first))
                self.assertTrue(catalog.content_object(db,sha))

    def test_representative_queries_use_their_indexes(self):
        """AC-01-03: the recorded query plans show no table scan."""
        stages,source,_=self._pipeline('run-plans')
        self._register('run-plans',stages)
        sha=artifact_manifest([str(source)])[str(source.resolve())]
        first=stages[0]['artifact_records'][0]['occurrence_id']
        plans={
            'entity':('''SELECT occurrence_id,qualifiers FROM artifact_entities
                WHERE project_id=? AND entity_type=? AND entity_id=?''',(self.pid,'compound','CHEMBL1'),
                'artifact_entities_entity'),
            'occurrence':('SELECT * FROM artifact_instances WHERE project_id=? AND occurrence_id IN (?)',
                (self.pid,first),'artifacts_occurrence'),
            'content':('SELECT * FROM artifact_instances WHERE project_id=? AND content_sha256=?',
                (self.pid,sha),'artifacts_content'),
            'producer':('SELECT * FROM artifact_instances WHERE project_id=? AND stage_id=? AND indexed_run_id=?',
                (self.pid,stages[0]['id'],'run-plans'),'artifacts_project_stage'),
            'operation':('SELECT * FROM operations WHERE project_id=? AND run_id=?',
                (self.pid,'run-plans'),'operations_project_run'),
            'producing':('''SELECT o.* FROM operation_outputs e
                JOIN operations o ON o.operation_id=e.operation_id
                WHERE e.occurrence_id=? AND o.project_id=?''',(first,self.pid),'operation_outputs_occurrence'),
            'signature':('SELECT * FROM operations WHERE project_id=? AND signature=?',
                (self.pid,'f'*64),'operations_signature'),
            'upstream':(catalog._UPSTREAM,(first,self.pid),'operation_outputs_occurrence'),
            'downstream':(catalog._DOWNSTREAM,(first,self.pid),'operation_inputs_occurrence'),
            # The guard of forget_run runs on every result change, so its cost
            # must follow what the run published, not the size of the project.
            'claimed':(catalog._CLAIMED,(self.pid,'run-plans',self.pid,'run-plans'),'artifacts_occurrence')}
        with self.store.connect() as db:
            for name,(query,values,index) in plans.items():
                with self.subTest(query=name):
                    steps=[row['detail'] for row in db.execute('EXPLAIN QUERY PLAN '+query,values)]
                    report=name+': '+' | '.join(steps)
                    self.assertEqual([step for step in steps if step.startswith('SCAN')],[],report)
                    self.assertTrue(any(index in step for step in steps),report)


class CatalogMigrationTests(CatalogFixture):
    def test_the_schema_records_its_version_and_has_nothing_pending(self):
        with self.store.connect() as db:
            self.assertEqual(schema_version(db),MIGRATIONS[-1][0])
            self.assertEqual(pending(db),[])

    def test_applying_the_migrations_again_changes_nothing(self):
        with self.store.connect() as db:
            before=[tuple(row) for row in db.execute('SELECT version,name,applied FROM schema_migrations')]
            self.assertEqual(apply(db),MIGRATIONS[-1][0])
            self.assertEqual([tuple(row) for row in db.execute('SELECT version,name,applied FROM schema_migrations')],before)

    def test_reopening_the_workspace_keeps_one_row_per_migration(self):
        WorkspaceStore(self.tmp.name)
        with self.store.connect() as db:
            self.assertEqual(db.execute('SELECT COUNT(*) FROM schema_migrations').fetchone()[0],len(MIGRATIONS))

    def test_an_interrupted_migration_completes_on_the_next_open(self):
        database=Path(self.tmp.name)/'interrupted.sqlite3'
        db=sqlite3.connect(database);db.row_factory=sqlite3.Row
        # The base tables are the precondition of a migration; WorkspaceStore
        # creates them before applying any.
        db.executescript('''CREATE TABLE artifact_instances(project_id TEXT,occurrence_id TEXT,
                indexed_run_id TEXT,legacy TEXT);
            CREATE TABLE artifact_entities(project_id TEXT,occurrence_id TEXT,entity_type TEXT,
                entity_id TEXT,qualifiers TEXT);
            CREATE INDEX artifact_entities_lookup ON artifact_entities(project_id,entity_type,entity_id);
            INSERT INTO artifact_instances VALUES ('p','o','r','kept');''')
        statements=MIGRATIONS[0][2]
        for statement in statements[:3]:db.execute(statement)
        db.commit()
        self.assertEqual(schema_version(db),0)
        self.assertEqual([name for _,name in pending(db)],[MIGRATIONS[0][1]])
        self.assertEqual(apply(db),MIGRATIONS[-1][0])
        self.assertEqual(pending(db),[])
        tables={row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        self.assertLessEqual({'content_objects','operations','operation_inputs','operation_outputs'},tables)
        indexes={row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='index'")}
        self.assertIn('artifact_entities_entity',indexes)
        self.assertNotIn('artifact_entities_lookup',indexes)
        self.assertEqual(db.execute('SELECT legacy FROM artifact_instances').fetchone()[0],'kept')
        db.close()

    def test_legacy_rows_and_tables_survive_the_migration(self):
        with self.store.connect() as db:
            rows=[dict(row) for row in db.execute('SELECT * FROM artifact_instances')]
            tables={row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        stages,_,_=self._pipeline('run-survivors')
        self._register('run-survivors',stages)
        WorkspaceStore(self.tmp.name)
        with self.store.connect() as db:
            self.assertLessEqual(tables,{row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='table'")})
            kept=[dict(row) for row in db.execute('SELECT * FROM artifact_instances')]
        self.assertEqual(kept[:len(rows)],rows)
        self.assertEqual(len(kept),len(rows)+2)


class CatalogBackfillTests(CatalogFixture):
    def _legacy_run(self,run_id):
        """A run indexed before the catalog recorded operations."""
        stages,source,derived=self._pipeline(run_id)
        self._register(run_id,stages)
        with self.store.connect() as db:
            db.execute('DELETE FROM catalog_indexed_runs WHERE project_id=? AND run_id=?',(self.pid,run_id))
            catalog.forget_run(db,self.pid,run_id)
            self.assertEqual(catalog.operations_of_run(db,self.pid,run_id),[])
        return stages

    def test_the_backfill_derives_operations_and_runs_only_once(self):
        stages=self._legacy_run('run-backfill')
        with self.store.connect() as db:
            self.assertEqual(catalog.backfill_operations(db),1)
            operations=catalog.operations_of_run(db,self.pid,'run-backfill')
            self.assertEqual(len(operations),2)
            edges=catalog.operation_edges(db,operations[1]['operation_id'])
            self.assertEqual(edges['inputs'][0]['occurrence_id'],stages[0]['artifact_records'][0]['occurrence_id'])
            self.assertEqual(catalog.backfill_operations(db),0)
            self.assertEqual(len(catalog.operations_of_run(db,self.pid,'run-backfill')),2)

    def test_opening_the_workspace_backfills_a_legacy_run(self):
        self._legacy_run('run-opened')
        WorkspaceStore(self.tmp.name)
        with self.store.connect() as db:
            self.assertEqual(len(catalog.operations_of_run(db,self.pid,'run-opened')),2)
            self.assertEqual(db.execute('''SELECT version FROM catalog_indexed_runs
                WHERE project_id=? AND run_id=?''',(self.pid,'run-opened')).fetchone()[0],MIGRATIONS[-1][0])

    def test_an_unreadable_run_is_recorded_and_skipped(self):
        with self.store.connect() as db:
            db.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?,NULL)',
                ('run-broken',self.pid,None,'succeeded','not json',time.time(),time.time()))
            self.assertEqual(catalog.backfill_operations(db),1)
            self.assertEqual([row['reason'] for row in db.execute(
                'SELECT reason FROM catalog_backfill_issues WHERE run_id=?',('run-broken',))],['unreadable_run'])


if __name__=='__main__':unittest.main()
