# ADR-001 — Scientific artifact identity model

- **Status:** Accepted
- **Date:** 2026-11-28
- **Phase:** PHASE-01 (#8)
- **Requirements:** REQ-ID-001, REQ-ID-002, REQ-ID-003, REQ-ID-004, REQ-CAT-001, REQ-DEP-001, REQ-PROV-004, AC-01-01, AC-01-02
- **Authors and reviewers:** Michel

## Context

The platform already distinguishes artifact identities, but only three of them
are first-class, and no operation has an identity of its own.

| Identity | Where it lives today | What it answers |
| --- | --- | --- |
| Logical artifact | `artifacts.slot_id` and `artifacts.records` write `artifact_id` and `occurrence_id` into `artifact_instances` | *Which output of which block this is, and from which run* |
| Scientific entity | `artifact_registry.index_entities` and the `artifact_entities` table | *Which compound, structure, receptor or pose this file holds* |
| Content | the `artifact_instances.content_sha256` column, filled from the stage's `artifact_manifest` | *Which bytes were published* |
| Operation | — | *Which execution of which block produced this, and from what* |

The `artifact_id` is the producer plus its qualified output (`slot_id(stage_id,
locator)`); the `occurrence_id` adds the run that produced the file. Neither
depends on a physical path, and the
[PHASE-00 assessment](../validation/artifact_architecture_assessment_2026-11-26/README.md)
confirmed that relocating a project preserves both. That is not the problem.

What is missing was measured in the same report:

- **F-04 — dependencies as JSON.** `artifact_instances.dependencies` holds a
  serialized list of `artifact_id`, and `provenance` holds an object. No index
  can answer *who consumed this occurrence*: the question requires reading and
  deserializing every row of the project.
- **Content is not an object.** As a column, the SHA-256 carries no size, format
  or storage location, and nothing stops the same content from being counted as
  several physical artifacts. F-02, F-03 and F-14 measured exactly that repeated
  copying of identical content across runs.
- **The operation is absent from the catalog.** The stage cache
  (`stage_cache.stage_key`, `operation_digest`) already computes a signature per
  execution, but it lives in a cache key, not in a queryable row. One cannot ask
  *which operations share this signature* or *what did this operation consume*,
  which are the questions of PHASE-03 and PHASE-04.

## Alternatives considered

| Alternative | Advantages | Disadvantages | Evidence |
| --- | --- | --- | --- |
| A — A single content identity: the SHA-256 becomes *the* artifact identity | Trivial deduplication; one identifier | Two distinct experiments that produce identical bytes stop being distinguishable, violating REQ-PROV-004; a pose's provenance would no longer say which redocking produced it | `artifact_instances` already records distinct occurrences with the same `content_sha256` |
| B — Keep dependencies in JSON and query with `json_extract` | No schema migration | No index can serve the JSON content; traversal cost grows with the number of artifacts in the project, which is precisely what F-04 reports | A query plan over `json_each` walks the whole table |
| C — Four separate identities with edges as rows *(selected)* | Each question has an index; reuse is a reference, not a copy; provenance survives deduplication | More tables and a derivation to maintain | Query plans in `docs/validation/phase_01_identity_catalog_2026-11-28/README.md` |
| D — A full provenance model (W3C PROV / RDF) | Standard, interoperable vocabulary | Requires an external dependency and additional storage, against REQ-ARCH-003; the vocabulary does not cover the domain's scientific entities without an extension of its own | 01.4 of #6 |

## Selected solution

Four identities with disjoint roles, none of them derived from a file path.

### 1. Content

`content_objects(sha256, size, media_type, storage_ref, created)`. The identity
is the SHA-256 of the bytes (REQ-ID-001). `media_type` is the format of those
bytes (`chemical/x-pdb`, `text/csv`), which is not the occurrence's scientific
`artifact_type`. `storage_ref` stays empty in this phase: managed storage is
decided in ADR-003 (PHASE-02), and a path is a location, never an identity.

Several logical artifacts reference one content object. That is what
`artifact_catalog.artifacts_by_content` returns, and what AC-01-01 verifies.

### 2. Logical artifact

`artifact_id` and `occurrence_id` stay **unchanged**, with the same computation
and the same columns. No already-persisted identity is recomputed: this phase
rewrites no row of `artifact_instances`.

### 3. Scientific entity

`artifact_entities` stays as it is, with `entity_type` ∈ {`compound`,
`structure`, `receptor`, `pose`, `community`, `molecular_pair`} and the
qualifiers (`chain_id`, `assembly_id`, `ligand_id`) in JSON. The
`artifacts_by_entity` query resolves entity → occurrences → artifacts
(AC-01-02).

The table's index now includes the columns the query returns
(`artifact_entities_entity`), because without them the planner preferred the
primary key using `project_id` alone — that is, it read every entity row of the
project. The previous index, `artifact_entities_lookup`, becomes redundant and
is dropped by the migration; dropping an index removes no row.

### 4. Operation

`operations` records **one execution of one block**:

```
operations(operation_id, project_id, operation_type, stage_id, block_id, run_id,
           signature, signature_version, parameters_digest, software_digest,
           state, started, finished, updated)
```

- `operation_id = SHA-256(['operation', project_id, stage_id, run_id])`, where
  `run_id` is the run that **produced** the artifacts. It is deterministic, so
  the derivation can be repeated without creating new rows.
- `run_id` is not necessarily the run that indexed the artifacts. When a run
  reuses a cached result, its occurrences keep pointing at the operation that
  produced them. That is how reuse appears as a reference rather than a copy
  (REQ-CAT-005).
- `parameters_digest` covers the stage's effective parameters;
  `software_digest` covers `processing_version` (the implementation digest),
  `producer_software` (tool versions) and `runtime_id`.
- `signature` combines operation type, parameters, software and the content
  digests of the inputs, prefixed by `signature_version`. **It is provisional:**
  ADR-004 (PHASE-03) fixes the definitive model, and the version column exists
  so an old signature is never mistaken for a new one.
- `software_digest`, and with it the signature, is `NULL` when the run recorded
  no software identity — the case for every run made before the implementation
  digest existed. A digest of absence would be a signature unable to tell two
  implementations apart, so the gap is recorded in `catalog_backfill_issues`
  instead of filled in.

### Dependency edges

```
operation_outputs(operation_id, occurrence_id, port)
operation_inputs(operation_id, port, position, occurrence_id, asset_id, role)
```

An input is another artifact's occurrence, a project asset or a file outside the
project, so `role` ∈ {`artifact`, `asset`, `external`} says which reference was
recorded and both reference columns are nullable. `position` orders the inputs
of one port and keeps the primary key stable when the derivation is repeated.

Lineage traversal reads those rows: `upstream` goes from an occurrence to the
operation that published it and on to that operation's inputs; `downstream`
does the reverse. The JSON `dependencies` column is still written, because it is
what the derivation reads and what earlier releases read.

### Where the derivation happens

`artifact_registry.register_run` is the single write path of the index, and now
derives the operations of the run it registers. The same derivation serves
history: `artifact_catalog.backfill_operations` reads `runs.stages`, which
already holds the artifact descriptors. A run is processed once per catalog
version, recorded in `catalog_indexed_runs` — the presence of rows in
`operations` cannot serve as the marker, precisely because reuse gives a run
operations that belong to other runs.

Re-indexing a run drops the operations it produced before rebuilding them, under
an explicit invariant: **an operation survives while an occurrence indexed by
another run still declares it as its producer.** Clearing a run's results does
not undo the fact that it ran, and the lineage of a run that reused those
results cannot depend on the order in which the project's runs are re-indexed.

What the derivation cannot establish is recorded, never guessed, following the
policy of `artifact_migration_issues`:

| Reason | Situation |
| --- | --- |
| `missing_artifact_records` | The stage published files with no descriptor; with no occurrence there is no edge |
| `missing_occurrence` | A descriptor without an `occurrence_id` |
| `missing_software_identity` | The run recorded no implementation version; the operation exists, without a signature |
| `unreadable_run` | The run's `stages` column is not valid JSON |

## Rationale

The four questions are different and each needs an index of its own.
Alternative A answers the content question and loses the provenance one: two
docking runs with different parameters that converge on the same pose would
become the same artifact, and the platform would no longer be able to say which
protocol produced each result — a scientific requirement, not a convenience.
Alternative B keeps the current structure and with it the cost F-04 measured.
Alternative D brings a vocabulary that would still need extending for
compounds, chains and crystallographic pairs, at the cost of infrastructure that
01.4 forbids.

Separating the operation, and not only the edges, is what lets the later phases
ask *has this operation already run with this signature* (PHASE-03) and *which
operations depend on this occurrence* (PHASE-04) without reopening the schema.

## Trade-offs

- **Operation granularity is the stage, not the batch.** A stage in individual
  mode processes several batches with their own versions and timings; the
  catalog records one operation per stage and run, consolidating the batches'
  software digests into one. Per-batch granularity is decided in ADR-004, with
  the signature model.
- **A stage that publishes no file produces no operation.** The derivation
  starts from published artifacts. Recording runs without output belongs to
  PHASE-04, which observes execution instead of reconstructing it.
- **This phase's signature is provisional** and must not be used as a cache key
  before ADR-004.
- **Temporary duplication.** Dependencies are stored twice, as JSON and as
  rows, until the older readers are migrated (PHASE-05).
- **One `stat` per published artifact** now feeds the content object's size. It
  is not a new cost: it is the same `stat` that already decided whether the file
  was available.

## Scientific implications

None. No algorithm, protocol, tolerance or rounding is touched. In particular:

- The `artifact_id`, `occurrence_id` and `artifact_type` of already-persisted
  rows stay identical — including the `molecular_poses` of redocking poses,
  deliberately preserved by decision C5 of the
  [test suite triage](../validation/test_suite_port_contract_2026-11-28/README.md).
- No cache identity moves. The 26 identities a source edit could move were
  measured before and after: none changed, so existing preparations, redockings
  and docking runs stay reusable.
- Traceability improves: the operation that produced an occurrence becomes a
  queryable row rather than an inference over filenames.

## Migration implications

The migration is additive: it creates tables and indexes and alters or removes
no row. Projects from earlier versions stay readable, and the history backfill
is idempotent and resumable — an interrupted pass leaves the runs it finished
marked.

An earlier release of the application opening the same database ignores the new
tables and keeps reading the `dependencies` column, which is still written.
`.bme.zip` packages do not change format in this phase.

## Validation requirements

| Criterion | Verification |
| --- | --- |
| AC-01-01 | `test_artifact_catalog.ContentIdentityTests` — two logical artifacts with the same bytes share one content object, with size and media type |
| AC-01-02 | `test_artifact_catalog.CatalogQueryTests.test_entities_resolve_to_artifacts_through_the_index` |
| REQ-DEP-001 | `OperationIdentityTests.test_each_block_execution_becomes_one_operation_with_its_edges` |
| REQ-PROV-004 | `test_a_reused_result_points_to_the_operation_that_produced_it` |
| The no-guessing policy | `test_a_published_file_without_a_descriptor_is_recorded_not_guessed`, `test_a_run_without_software_identity_has_no_signature` |

The full report, with the query plans and the cache identity measurement:
`docs/validation/phase_01_identity_catalog_2026-11-28/README.md`.

## Deviations from the reference specification

Appendix A.2 and sub-issue #8 give reference names for the tables. The
deviations, all documented here, are:

| Deviation | Reason |
| --- | --- |
| `operations` gains `signature_version` and `updated` | Without the version, this phase's provisional signature would be indistinguishable from ADR-004's definitive one |
| `operation_inputs` gains `position`, `asset_id` and `role` instead of `(operation_id, occurrence_id, port, role)` | An input may not be an artifact; `position` gives repeated inputs on the same port a stable primary key |
| The additional `catalog_indexed_runs` and `catalog_backfill_issues` tables | Backfill idempotence, and a record of what cannot be derived |
| `artifact_instances.content_sha256` references `content_objects` with no declared foreign key | Old rows may hold unregistered content; a constraint would make them unreadable, against REQ-MIG-001 |
