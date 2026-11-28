"""Ordered, idempotent and non-destructive migrations of the artifact catalog.

The tables ``WorkspaceStore`` creates with ``CREATE TABLE IF NOT EXISTS``
answer whether a table exists, never which version of it the database holds.  A
catalog that gains tables, columns and indexes across releases needs the second
answer, so every change to it is declared here, applied in order, applied once,
and recorded in ``schema_migrations``.

Each migration is a sequence of statements that are safe to run twice, and its
version row is written last.  An interruption therefore leaves the statements
that already ran in place and no version row, and the next open applies the
whole migration again from its start.  No migration drops or rewrites a row:
legacy data is read by the new tables, never replaced by them.
"""
import time

CATALOG_VERSION=1

MIGRATIONS=(
    (1,'content objects, operation identity and dependency edges',(
        # Physical content, identified by the SHA-256 of its bytes (REQ-ID-001).
        # ``storage_ref`` stays empty until PHASE-02 introduces managed
        # storage; a path is a location, never an identity.
        '''CREATE TABLE IF NOT EXISTS content_objects(
            sha256 TEXT PRIMARY KEY,size INTEGER,media_type TEXT NOT NULL,
            storage_ref TEXT,created REAL NOT NULL)''',
        # One execution of one block.  ``run_id`` is the run that produced the
        # artifacts, which is not the run that indexed them when a later run
        # reuses a cached result.
        '''CREATE TABLE IF NOT EXISTS operations(
            operation_id TEXT PRIMARY KEY,project_id TEXT NOT NULL,operation_type TEXT NOT NULL,
            stage_id TEXT NOT NULL,block_id TEXT,run_id TEXT,signature TEXT,signature_version INTEGER NOT NULL,
            parameters_digest TEXT,software_digest TEXT,state TEXT NOT NULL,
            started REAL,finished REAL,updated REAL NOT NULL)''',
        'CREATE INDEX IF NOT EXISTS operations_project_run ON operations(project_id,run_id)',
        'CREATE INDEX IF NOT EXISTS operations_project_stage ON operations(project_id,stage_id,run_id)',
        'CREATE INDEX IF NOT EXISTS operations_signature ON operations(project_id,signature)',
        # Dependency edges as rows.  The JSON column they replace stays in
        # place: it is what the backfill reads, and older releases still read it.
        '''CREATE TABLE IF NOT EXISTS operation_outputs(
            operation_id TEXT NOT NULL,occurrence_id TEXT NOT NULL,port TEXT,
            PRIMARY KEY(operation_id,occurrence_id))''',
        'CREATE INDEX IF NOT EXISTS operation_outputs_occurrence ON operation_outputs(occurrence_id)',
        # An input is an occurrence of another artifact, a project asset or a
        # file outside the project, so the reference columns are nullable and
        # ``role`` says which one was recorded.  ``position`` orders the inputs
        # of one port and keeps the primary key stable across re-derivations.
        '''CREATE TABLE IF NOT EXISTS operation_inputs(
            operation_id TEXT NOT NULL,port TEXT NOT NULL,position INTEGER NOT NULL,
            occurrence_id TEXT,asset_id TEXT,role TEXT NOT NULL,
            PRIMARY KEY(operation_id,port,position))''',
        'CREATE INDEX IF NOT EXISTS operation_inputs_occurrence ON operation_inputs(occurrence_id)',
        # Entity and content lookups reach an occurrence first.  The primary key
        # of ``artifact_instances`` starts with the indexed run, which such a
        # lookup does not know.  ``indexed_run_id`` is part of the index so that
        # asking whether another run still holds an occurrence is answered
        # without reading the row.
        '''CREATE INDEX IF NOT EXISTS artifacts_occurrence
            ON artifact_instances(project_id,occurrence_id,indexed_run_id)''',
        # Looking an entity up returns its occurrences and qualifiers.  With
        # those two columns in the index the planner uses all three equality
        # terms; without them it preferred the primary key on ``project_id``
        # alone, which reads every entity row of the project.
        '''CREATE INDEX IF NOT EXISTS artifact_entities_entity
            ON artifact_entities(project_id,entity_type,entity_id,occurrence_id,qualifiers)''',
        # Redundant once the index above exists: same leading column, fewer
        # terms usable.  Dropping an index removes no row.
        'DROP INDEX IF EXISTS artifact_entities_lookup',
        # What the derivation could not establish is recorded, never guessed,
        # following the policy of ``artifact_migration_issues``.
        '''CREATE TABLE IF NOT EXISTS catalog_backfill_issues(
            project_id TEXT NOT NULL,run_id TEXT NOT NULL,stage_id TEXT NOT NULL,
            reason TEXT NOT NULL,details TEXT NOT NULL,recorded REAL NOT NULL,
            PRIMARY KEY(project_id,run_id,stage_id,reason))''',
        # A run whose operations have been derived.  Reuse gives a run
        # operations that belong to earlier runs, so the presence of rows for
        # the run itself cannot answer whether it has been processed.
        '''CREATE TABLE IF NOT EXISTS catalog_indexed_runs(
            project_id TEXT NOT NULL,run_id TEXT NOT NULL,version INTEGER NOT NULL,updated REAL NOT NULL,
            PRIMARY KEY(project_id,run_id))''',
    )),
)


def apply(db):
    """Bring the catalog schema up to date inside the caller's connection."""
    db.execute('''CREATE TABLE IF NOT EXISTS schema_migrations(
        version INTEGER PRIMARY KEY,name TEXT NOT NULL,applied REAL NOT NULL)''')
    applied=_applied(db)
    for version,name,statements in MIGRATIONS:
        if version in applied:continue
        for statement in statements:db.execute(statement)
        db.execute('INSERT OR REPLACE INTO schema_migrations VALUES (?,?,?)',(version,name,time.time()))
        # The version row is durable only with the statements it describes, so
        # an interruption cannot leave a version claiming absent tables.
        db.commit()
    return schema_version(db)


def _applied(db):
    """Versions recorded as applied; a database predating the mechanism has none."""
    import sqlite3
    try:
        return {row[0] for row in db.execute('SELECT version FROM schema_migrations')}
    except sqlite3.OperationalError:
        return set()


def schema_version(db):
    """Highest migration recorded as applied, or zero for an untouched database."""
    return max(_applied(db),default=0)


def pending(db):
    """Migrations an open database still needs, in order."""
    applied=_applied(db)
    return [(version,name) for version,name,_ in MIGRATIONS if version not in applied]
