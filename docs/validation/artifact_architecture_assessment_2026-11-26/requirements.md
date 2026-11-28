# Requirement register — #6 (PHASE-00 proposal)

Status: **proposed for review** in #6. Once the review is accepted, these
identifiers are permanent and must never be renumbered (02.3 of #6). A
withdrawn requirement keeps its identifier with state `NOT_APPLICABLE` and a
justification.

`REQ-ENG-001 … 010` are defined in section 02.2 of #6 and are not repeated here.
The "Source" column traces every requirement to #6 (sections 01–04 and
Appendix A) or to a PHASE-00 finding (`F-nn`, see `README.md`).

States: `PENDING`, `IN_PROGRESS`, `IMPLEMENTED`, `VALIDATED`, `BLOCKED`,
`NOT_APPLICABLE`. At the end of PHASE-00 every requirement is `PENDING`, except
those listed with another state in the status table at the end.

## REQ-ARCH — Architecture

| ID | Requirement | Source |
| --- | --- | --- |
| REQ-ARCH-001 | Computational modules and presentation components access artifacts through the artifact services and data access layer of 03.2, never through the physical storage layout. | 03.2, 04 |
| REQ-ARCH-002 | There is a single persistent catalog and a single content store. Existing components (`artifact_registry`, `calculation_cache`, `stage_cache`, `preparation_cache`, `result_catalog`) are generalized, not duplicated. | REQ-ENG-001, A.4, F-08 |
| REQ-ARCH-003 | The architecture requires no external database server, cloud service or other mandatory infrastructure. | 01.4 |
| REQ-ARCH-004 | Each significant decision is recorded as an ADR (A.6) in both languages. | REQ-ENG-007, A.6 |

## REQ-ID — Identity

| ID | Requirement | Source |
| --- | --- | --- |
| REQ-ID-001 | Physical content is identified by the SHA-256 of its bytes. | A.2, ADR-001 |
| REQ-ID-002 | Logical artifacts keep independent identities (`artifact_id`, `occurrence_id`) that depend neither on physical paths nor on display names. | 01.2, A.2 |
| REQ-ID-003 | Scientific entities (compound, target, structure/chain, receptor, pose, pair, community, trajectory) have canonical identifiers with explicit qualifiers. | 01.2, F-10 |
| REQ-ID-004 | Each computational operation has a deterministic identity derived from its signature (REQ-CACHE-001). | ADR-001, ADR-004 |

## REQ-CAT — Catalog and indexing

| ID | Requirement | Source |
| --- | --- | --- |
| REQ-CAT-001 | The catalog persists content objects, logical artifacts, occurrences, entities, operations and dependency edges as rows, not as embedded JSON. | 01.2, A.2, F-04 |
| REQ-CAT-002 | Schema changes use ordered, versioned, idempotent and non-destructive migrations. | A.2, REQ-ENG-005 |
| REQ-CAT-003 | Artifacts are retrievable by scientific entity, content hash, producer and dependency without recursive directory traversal. | 01.2, 04 |
| REQ-CAT-004 | Representative queries have documented query plans that use indexes. | A.5 |
| REQ-CAT-005 | Catalog and descriptor overhead is measured (`IO`) and does not grow with the number of times an artifact is reused. A reused artifact is referenced, not copied into each run. | 01.2, F-02, F-03, F-14 |

## REQ-STO — Storage

| ID | Requirement | Source |
| --- | --- | --- |
| REQ-STO-001 | Identical immutable content is stored once in managed storage. | 01.2, 04 |
| REQ-STO-002 | Stored content is immutable and verified by hash before it is reused or restored. | 01.3 |
| REQ-STO-003 | Writable execution inputs are isolated copies. Stored objects are never modified through a materialization. | 01.3 |
| REQ-STO-004 | Publication is atomic (temporary file, fsync, rename) and happens only after validation. | A.2 |
| REQ-STO-005 | Original scientific inputs and essential results are never deleted automatically. | 01.4 |
| REQ-STO-006 | Project history, versioning and execution folders reference content instead of keeping additional full copies of results or of the run state. | 01.2, F-01, F-05, F-06 |

## REQ-CACHE — Computational caching

| ID | Requirement | Source |
| --- | --- | --- |
| REQ-CACHE-001 | Cache signatures are deterministic. They cover operation type and version, effective parameters, input content hashes, tool identity, implementation digest and schema version, and exclude cosmetic fields. | 01.2, A.2 |
| REQ-CACHE-002 | One signature model is shared by all modules. Each existing key format is mapped onto it. | REQ-ENG-001, F-08 |
| REQ-CACHE-003 | Cached results are validated (hashes and artifact contract) before reuse. A failed validation is a miss, never a hit. | 01.3 |
| REQ-CACHE-004 | Users can bypass the cache for one operation and force recomputation of an operation and its descendants. | A.3 |
| REQ-CACHE-005 | Concurrent requests with the same signature compute once (single-flight), and stale locks are recovered. | A.3 |

## REQ-DEP — Dependencies

| ID | Requirement | Source |
| --- | --- | --- |
| REQ-DEP-001 | Every executed operation records its input and output edges. | 01.2 |
| REQ-DEP-002 | A change invalidates only the dependent operations. Invalidation marks artifacts as stale and never deletes content. | 01.2, 01.4 |
| REQ-DEP-003 | Incremental execution plans are deterministic and persisted with the run. | 01.2 |

## REQ-REUSE — Artifact reuse

| ID | Requirement | Source |
| --- | --- | --- |
| REQ-REUSE-001 | Scientifically compatible artifacts are reused across modules, starting with receptors prepared for redocking, Vina and DOCK6. | 01.2, 03.1 |
| REQ-REUSE-002 | Each artifact type has explicit compatibility rules. A mismatch is a rejection with a recorded reason, never a fallback. | 01.3, 01.4 |
| REQ-REUSE-003 | Reuse decisions (reused, rejected with reason, recomputed) are recorded and visible to the user. | 01.3, F-11 |
| REQ-REUSE-004 | Reuse never turns exclusions or failures into successes. | 01.3, `docs/docking_cache.md` |

## REQ-API — Data access

| ID | Requirement | Source |
| --- | --- | --- |
| REQ-API-001 | An artifact repository resolves logical artifacts and entities to validated content. | 04 |
| REQ-API-002 | Scientific query services are independent of the physical storage layout. | 03.2, 04 |
| REQ-API-003 | Queries support stable pagination, batch retrieval and column projection. Large datasets are read selectively. | 01.2, 04 |
| REQ-API-004 | Authorization, project membership and edit permission are checked on every call. | A.2 |
| REQ-API-005 | Query results are equivalent to the existing readers. | 01.3 |

## REQ-VIS — Visualization

| ID | Requirement | Source |
| --- | --- | --- |
| REQ-VIS-001 | Priority dashboards use the standardized interfaces and do not depend on routine directory scans. | 01.2, 04, F-07 |
| REQ-VIS-002 | Molecular structures and poses are retrieved on demand. | 01.2 |
| REQ-VIS-003 | Filtering, sorting and search run in the data layer. | 01.2 |
| REQ-VIS-004 | Derived views are cached by source content, view parameters and renderer version, and are invalidated when their source changes. | 01.2 |

## REQ-INT — Integrity

| ID | Requirement | Source |
| --- | --- | --- |
| REQ-INT-001 | An integrity assessment detects missing content, hash mismatches, orphaned objects, dangling references (including project locations) and catalog/run divergence. | 04, F-15 |
| REQ-INT-002 | The catalog can be rebuilt from authoritative sources, and interrupted publications are recovered. | 04 |
| REQ-INT-003 | Interrupted operations resume from validated outputs only. | 01.2 |

## REQ-PROV — Provenance

| ID | Requirement | Source |
| --- | --- | --- |
| REQ-PROV-001 | Each artifact records its producer operation, inputs, effective parameters, software and run/job. | 01.3 |
| REQ-PROV-002 | Each artifact has a provenance level (`full`, `partial`, `none`). Artifacts with insufficient provenance are never reused automatically. | 01.3, 04 |
| REQ-PROV-003 | Historical provenance survives reuse, migration, export and import. | 01.3 |
| REQ-PROV-004 | Distinct experiments remain distinguishable even when they share content. | 01.3 |

## REQ-MIG — Migration

| ID | Requirement | Source |
| --- | --- | --- |
| REQ-MIG-001 | Legacy projects are indexed in place, idempotently, resumably and non-destructively. | 01.2, 01.4 |
| REQ-MIG-002 | `.bme.zip` export and import stay self-contained and lossless with deduplicated content. | 04 |
| REQ-MIG-003 | Cleanup runs a dry run first, requires explicit confirmation and protects inputs, final results and history references. Retention policies are explicit. | 01.4, F-09 |
| REQ-MIG-004 | Compatibility with projects and caches of earlier versions is stated and tested. | REQ-ENG-005 |

## REQ-TEST — Testing

| ID | Requirement | Source |
| --- | --- | --- |
| REQ-TEST-001 | The existing suite passes after each phase. Until the baseline failures (F-16, F-17) are resolved in separate work, a phase passes when it introduces no failure absent from the PHASE-00 baseline (`tests.json`). | REQ-ENG-008, F-16, F-17 |
| REQ-TEST-002 | Each phase adds unit and integration tests for its components. | A.4 |
| REQ-TEST-003 | Recovery, crash-injection and concurrency scenarios are tested. | A.3 |
| REQ-TEST-004 | Scientific validation uses real engines and is reported separately from mock-based tests. | A.4, REQ-ENG-009 |
| REQ-TEST-005 | Golden or equivalence tests compare new readers and views with the current ones. | 01.3 |

## REQ-PERF — Performance

| ID | Requirement | Source |
| --- | --- | --- |
| REQ-PERF-001 | Baseline (PHASE-00) and final (PHASE-08) measurements follow the methodology of A.5. | A.5 |
| REQ-PERF-002 | `S_comp`, `R_storage`, `R_query`, `CHR`, `ACR` and `IO` are reported for each project size. | A.5 |
| REQ-PERF-003 | Acceptance thresholds are derived from the baseline and reviewed before PHASE-01 (`README.md`, Plan). | A.5 |
| REQ-PERF-004 | Workloads without reuse do not regress beyond the agreed threshold. | A.5 |
| REQ-PERF-005 | Visualization latency improvements are demonstrated on real projects. | 04 |

## Acceptance criteria

`AC-PP-NN` identifies acceptance criterion `NN` of phase `PP`, in the order of the
sub-issue. Each criterion lists the requirements it verifies.

| AC | Criterion | Verifies |
| --- | --- | --- |
| AC-00-01 | Existing behavior is understood and the architecture map is reviewed. | REQ-ARCH-001 |
| AC-00-02 | Scientific invariants are documented. | REQ-ENG-003 |
| AC-00-03 | Reuse and indexing opportunities are identified and ranked by measured impact. | REQ-PERF-001 |
| AC-00-04 | Architectural risks are documented. | REQ-ENG-007 |
| AC-00-05 | Baseline measurements are reproducible from committed scripts. | REQ-PERF-001, REQ-ENG-009 |
| AC-00-06 | No production code was modified. | REQ-ENG-002 |
| AC-01-01 | Multiple logical artifacts can reference one physical content object. | REQ-ID-001, REQ-ID-002, REQ-CAT-001 |
| AC-01-02 | Scientific entities can be associated with artifacts. | REQ-ID-003 |
| AC-01-03 | Indexed queries work without recursive directory scans. | REQ-CAT-003, REQ-CAT-004 |
| AC-01-04 | ADR-001 and ADR-002 are merged in both languages. | REQ-ARCH-004 |
| AC-02-01 | Identical content is stored only once in managed storage. | REQ-STO-001 |
| AC-02-02 | Writable execution inputs are isolated. | REQ-STO-003 |
| AC-02-03 | Artifact relationships remain traceable. | REQ-PROV-001 |
| AC-02-04 | Storage impact, including catalog and index overhead, is measured. | REQ-PERF-002, REQ-CAT-005 |
| AC-03-01 | Compatible prepared receptors can be reused. | REQ-REUSE-001 |
| AC-03-02 | Incompatible receptors are rejected. | REQ-REUSE-002 |
| AC-03-03 | Reuse produces scientifically consistent downstream results. | REQ-TEST-004, REQ-ENG-003 |
| AC-03-04 | Redundant preparation is measurably reduced compared with the baseline. | REQ-PERF-002 |
| AC-03-05 | ADR-004 is merged. | REQ-ARCH-004, REQ-CACHE-001 |
| AC-04-01 | Localized parameter changes affect only dependent operations. | REQ-DEP-002 |
| AC-04-02 | Interrupted runs can resume from validated outputs. | REQ-INT-003 |
| AC-04-03 | Concurrent requests do not duplicate cacheable computations. | REQ-CACHE-005 |
| AC-04-04 | ADR-005 is merged, and the resume behavior change is approved and documented. | REQ-ARCH-004 |
| AC-05-01 | Scientific queries are independent of the physical storage layout. | REQ-API-002 |
| AC-05-02 | Large datasets can be accessed selectively. | REQ-API-003 |
| AC-05-03 | Query results remain equivalent to the existing readers. | REQ-API-005 |
| AC-05-04 | Query plans show index usage, and latency is compared with the baseline. | REQ-CAT-004, REQ-PERF-002 |
| AC-05-05 | ADR-006 and ADR-007 are merged. | REQ-ARCH-004 |
| AC-06-01 | Dashboard operations show measurable improvements against the baseline. | REQ-PERF-005 |
| AC-06-02 | Visualizations preserve scientific correctness. | REQ-TEST-005 |
| AC-06-03 | Routine result retrieval no longer depends on full directory scans. | REQ-VIS-001 |
| AC-06-04 | ADR-011 is merged. | REQ-ARCH-004 |
| AC-07-01 | Existing projects remain accessible. | REQ-MIG-001, REQ-MIG-004 |
| AC-07-02 | Legacy artifacts with insufficient provenance are not silently reused. | REQ-PROV-002 |
| AC-07-03 | Cleanup preserves required scientific content. | REQ-MIG-003, REQ-STO-005 |
| AC-07-04 | Project migration is non-destructive. | REQ-MIG-001 |
| AC-07-05 | ADR-009, ADR-010 and ADR-012 are merged. | REQ-ARCH-004 |
| AC-08-01 | Existing workflows remain operational. | REQ-TEST-001 |
| AC-08-02 | Performance improvements are demonstrated, and regressions are reported. | REQ-PERF-002, REQ-PERF-004 |
| AC-08-03 | Scientific reproducibility and traceability are preserved. | REQ-PROV-003, REQ-TEST-004 |
| AC-08-04 | Every requirement and AC is in a valid state with evidence or justification. | REQ-ENG-010 |

## Status at the end of PHASE-00

| Identifier | State | Evidence |
| --- | --- | --- |
| AC-00-01 … AC-00-06 | `IMPLEMENTED` (awaiting review) | `README.md`, `audit.json`, `baseline.json` |
| REQ-PERF-001 | `IN_PROGRESS` | Baseline done; final measurement in PHASE-08 |
| REQ-PERF-003 | `IN_PROGRESS` | Thresholds proposed in `README.md`, pending review |
| REQ-ARCH-004 | `IN_PROGRESS` | ADR process and template published; no ADR yet |
| REQ-TEST-001 | `IN_PROGRESS` | Baseline suite failed: 45/169 modules, 50 failures, 38 errors, 1 hang (`tests.json`). After the [triage](../test_suite_triage_2026-11-27/README.md), the [decisions](../test_suite_decisions_2026-11-27/README.md), the [code fixes](../test_suite_code_fixes_2026-11-27/README.md) and the [port contract](../test_suite_port_contract_2026-11-28/README.md), the suite is green: 169/169 modules, 0 failures, 3 tracked `expectedFailure` tests. Decisions C1–C8 are closed; the two remaining items are the I1 and I2 investigations, both PHASE-07 (#14). A phase passes when it adds no failure and turns no tracked item into an unexpected success without updating it. |
| All other REQ-* and AC-* | `PENDING` | — |

## Review and status at the end of PHASE-01

The register was approved on 2026-11-28, together with the acceptance thresholds
proposed in `README.md`. The identifiers above are therefore permanent (02.3)
and are not renumbered. PHASE-01 is reported in
[phase_01_identity_catalog_2026-11-28](../phase_01_identity_catalog_2026-11-28/README.md).

| Identifier | State | Evidence |
| --- | --- | --- |
| AC-01-01, AC-01-02, AC-01-03, AC-01-04 | `IMPLEMENTED` | `tests/test_artifact_catalog.py`, query plans and ADR-001/ADR-002 in both languages |
| REQ-ID-001, REQ-ID-002, REQ-ID-003 | `IMPLEMENTED` | `content_objects`; logical identities unchanged; entity queries served by an index |
| REQ-ID-004 | `IN_PROGRESS` | Operation identity recorded; the signature is provisional until ADR-004 (PHASE-03) |
| REQ-CAT-001, REQ-CAT-002, REQ-CAT-003, REQ-CAT-004 | `IMPLEMENTED` | Content, operations and edges as rows; ordered migrations; no plan contains a table scan |
| REQ-CAT-005 | `IN_PROGRESS` | Overhead measured at 670 B per artifact and independent of reuse in the new tables; the `artifact_instances` duplication goes to PHASE-02 |
| REQ-DEP-001 | `IN_PROGRESS` | Edges recorded for operations that publish artifacts; runs without output belong to PHASE-04 |
| REQ-ARCH-002, REQ-ARCH-003 | `IMPLEMENTED` | One catalog generalized, not duplicated; no external service (ADR-002) |
| REQ-ARCH-004 | `IN_PROGRESS` | Two of the twelve ADRs accepted |
| REQ-MIG-001, REQ-MIG-004 | `IN_PROGRESS` | Additive migration and idempotent, resumable backfill, verified against legacy rows; project migration as a whole is PHASE-07 |
| REQ-PERF-003 | `IMPLEMENTED` | Thresholds reviewed and approved on 2026-11-28 |
| REQ-TEST-001, REQ-TEST-002 | `IMPLEMENTED` for this phase | Suite green with 22 new tests, no new failure and no tracked item turned into an unexpected success |

## Status at the end of PHASE-02

Reported in
[phase_02_content_store_2026-11-28](../phase_02_content_store_2026-11-28/README.md).

| Identifier | State | Evidence |
| --- | --- | --- |
| AC-02-01, AC-02-02, AC-02-03, AC-02-04 | `IMPLEMENTED` | `tests/test_content_store.py`; `R_storage` = 0.333 for one receptor published by three blocks, with three hard-linked materializations costing 0 extra blocks |
| REQ-STO-001, REQ-STO-002, REQ-STO-003, REQ-STO-004 | `IMPLEMENTED` | One object per content; verification on read and before reuse; isolated writable materialization; temporary file, fsync, digest check, atomic rename |
| REQ-STO-005 | `IMPLEMENTED` for this phase | Recovery never deletes referenced content |
| REQ-STO-006 | `IN_PROGRESS` | Published results stay in the run folders; ADR-010 (PHASE-07) decides the move, with retention |
| REQ-INT-002 | `IN_PROGRESS` | Interrupted publications recover; rebuilding the catalog belongs to PHASE-07 |
| REQ-CAT-005 | `IN_PROGRESS` | `R_storage` measured; the `artifact_instances` duplication goes to PHASE-05 |
| REQ-ARCH-002 | `IMPLEMENTED` | One store, generalized from `calculation_cache` rather than duplicated |
| REQ-ARCH-004 | `IN_PROGRESS` | Four of the twelve ADRs accepted (ADR-001, ADR-002, ADR-003, ADR-008) |
| REQ-MIG-004 | `IN_PROGRESS` | 16 identities moved by routing `calculation_cache`, registered as an audited equivalence and measured back to their previous values |
| REQ-TEST-003 | `IMPLEMENTED` for this phase | Crash injection at four points and concurrent publication |
| REQ-TEST-001, REQ-TEST-002 | `IMPLEMENTED` for this phase | Suite green with 28 new tests, no new failure and no tracked item turned into an unexpected success |
