# ADR-002 — Persistent metadata catalog technology

- **Status:** Accepted
- **Date:** 2026-11-28
- **Phase:** PHASE-01 (#8)
- **Requirements:** REQ-ARCH-002, REQ-ARCH-003, REQ-CAT-002, REQ-CAT-004, AC-01-03
- **Authors and reviewers:** Michel

## Context

The catalog already exists and is already SQLite: `workspace.py` opens
`workspace.sqlite3`, creates the workspace and artifact-index tables, and
`result_catalog.create_schema` adds the result read model. This phase's decision
is not to choose a new technology but to record why the current one is kept,
what it does not deliver, and how the schema evolves from now on.

What motivates the record:

- **The schema had no version.** Every table was created with
  `CREATE TABLE IF NOT EXISTS`, which answers whether a table exists, never
  which version of it the database holds. Schema repairs lived in
  `artifact_registry.migrate_schema`, triggered by inspecting
  `PRAGMA table_info` — a mechanism that works to rename a column once, but
  neither orders changes nor is idempotent by construction (REQ-CAT-002).
- **The platform is local by requirement.** 01.4 of #6 demands that it work
  without a database server, cloud service or any mandatory infrastructure.
- **Real concurrency is modest.** One host, one interface process and workers in
  subprocesses. Writes are already serialized by `BEGIN IMMEDIATE`, and WAL,
  enabled once per file in `WorkspaceStore._enable_concurrent_reads`, lets the
  interface read a run's progress without waiting for it.

## Alternatives considered

| Alternative | Advantages | Disadvantages | Evidence |
| --- | --- | --- | --- |
| A — SQLite with WAL and ordered migrations *(selected)* | No infrastructure; transactional; verifiable indexes and plans; already in production in this project | One writer at a time; no distributed queue; no planner statistics without `ANALYZE` | Plans in `docs/validation/phase_01_identity_catalog_2026-11-28/README.md` |
| B — PostgreSQL | Real write concurrency; rich types; partitioning | Requires an installed, configured service, against REQ-ARCH-003; a project would stop being a folder the user can copy | 01.4 of #6 |
| C — DuckDB | Excellent for columnar analytical scans of large tables | A concurrency model unsuited to the small, frequent metadata writes of a running execution; it would duplicate the catalog instead of generalizing it (REQ-ARCH-002) | — |
| D — An embedded key-value store (LMDB, SQLite without SQL) | Fast pair writes | No composite-index queries and no verifiable plans; AC-01-03 and REQ-CAT-004 ask for exactly that | — |
| E — JSON files with an in-memory index | No dependency | This is the state F-04 and F-07 measured as the dominant cost; the index is lost on every open | `docs/validation/artifact_architecture_assessment_2026-11-26/README.md` |

For **analytical datasets** — score tables with millions of rows, not metadata —
the choice is a different one. Columnar scan cost is ADR-006's problem
(PHASE-05), and DuckDB or Parquet remain candidates there. This decision covers
the metadata catalog.

## Selected solution

Keep SQLite and give it an explicit mechanism for evolving.

### Ordered, idempotent and non-destructive migrations

`catalog_schema.MIGRATIONS` declares each change as a sequence of statements,
with a permanent number and a name. `catalog_schema.apply` creates
`schema_migrations(version, name, applied)`, applies what is missing in order,
and writes the version row **last**, in the same commit as the statements it
describes.

The three required properties follow from that:

- **Ordered.** Migrations are applied in the order of their numbers, which are
  permanent and never reused.
- **Idempotent.** Every statement is safe to run twice (`IF NOT EXISTS`), and an
  already recorded migration is not reapplied.
- **Resumable.** An interruption leaves the statements that already ran in place
  and no version row. The next open reapplies the whole migration from its
  start, destroying nothing — which is why the statements' idempotence is the
  condition, not an implementation detail.

No migration removes or rewrites a row. The only removal the first migration
performs is of a redundant index (`artifact_entities_lookup`), replaced by one
the planner can use in full.

Migrations run when the workspace opens, right after
`artifact_registry.migrate_schema` and before any indexing, so that the history
backfill writes into the current schema.

### Accepted and recorded limits

| Limit | Practical consequence |
| --- | --- |
| One writer at a time | Two simultaneous runs serialize their index writes. WAL keeps reads free; the computation queue is PHASE-04's subject |
| One host | There is no distributed execution. 01.4 of #6 already establishes this |
| Query variable limit (999) | Queries over a list of occurrences are chunked in blocks of 400 (`artifact_catalog._instances`) |
| No `ANALYZE` | The planner chooses without statistics, so an index is only used with confidence when it covers the returned columns. That is what required ADR-001's entity index |
| `PRAGMA journal_mode=WAL` can be refused | Filesystems without shared memory keep the previous mode, and the workspace opens anyway — existing behavior, preserved |

## Rationale

The current technology meets this phase's verifiable requirements: the query
plans use indexes on every representative access path (AC-01-03), the
transaction is real, and a project remains a folder the user copies. Changing
the database would bring write concurrency the platform does not use, at the
cost of a dependency 01.4 forbids. What was missing was not the database but
schema discipline — and that is what this decision adds.

## Trade-offs

- The catalog does not scale horizontally, and is not meant to.
- Without `ANALYZE`, index choice depends on coverage. Covering indexes cost
  writes and space; the phase report measures that cost.
- `schema_migrations` becomes a point of truth: a published migration is never
  edited, only superseded by another.

## Scientific implications

None. The decision reads, writes or converts no scientific data. Catalog
availability affects what the interface can list, not what the algorithms
compute.

## Migration implications

A database from an earlier version is updated on its first open, additively. An
already updated database opened by an earlier release of the application keeps
working: the new tables are ignored, and the columns that release reads are
still written.

Exported projects (`.bme.zip`) do not change format in this phase. Cache
compatibility is preserved: no identity moved.

## Validation requirements

| Criterion | Verification |
| --- | --- |
| REQ-CAT-002 | `test_artifact_catalog.CatalogMigrationTests` — version recorded, reapplication without effect, an interruption completed on the next open, legacy rows and tables preserved |
| REQ-CAT-004, AC-01-03 | `test_representative_queries_use_their_indexes`, with the plans in the phase report |
| REQ-MIG-001 | `CatalogBackfillTests` — idempotent backfill, an unreadable run recorded and skipped |

## Deviations from the reference specification

None. Sub-issue #8 explicitly asks that SQLite be kept with WAL, short
connections and the existing fallback, which is what this decision records.
