# ADR-003 — Content-addressable storage strategy

- **Status:** Accepted
- **Date:** 2026-11-28
- **Phase:** PHASE-02 (#9)
- **Requirements:** REQ-STO-001, REQ-STO-002, REQ-STO-003, REQ-STO-005, REQ-ARCH-002, REQ-ARCH-003, REQ-CAT-005, AC-02-01, AC-02-02, AC-02-04
- **Authors and reviewers:** Michel

## Context

The platform already stored content by its digest, in two places:

| Mechanism | What it does | What was missing |
| --- | --- | --- |
| `calculation_cache` | An `objects` table and directory, `copy_artifact` with copy-on-write, transactional publication after validation, per-key locks | Writable objects; verification of the written bytes left to the caller; no policy for handing content to a consumer |
| `preparation_cache` | Immutable entries published by atomic rename after every required file is validated, with every digest verified on read | A layout of its own, by suffix rather than by content; migrates in PHASE-03 |

Related pieces: `storage.py` (atomic writes and directory fsync),
`result_transactions.py` (file moves recovered around a commit) and the public
names of `artifact_naming.py`.

Three facts measured on this machine decide the policy, and are recorded in the
[phase report](../validation/phase_02_content_store_2026-11-28/README.md):

1. **The filesystem is ext4 and has no reflink support.** `FICLONE` fails with
   `EOPNOTSUPP`, so every copying materialization is a full copy. Appendix A's
   "reflink, then copy" policy degrades to a copy here.
2. **Hard links work.** A link to an object costs no bytes.
3. **An object in mode 0444 refuses writes.** The kernel denies the write open
   instead of letting a consumer corrupt content other runs reference.

Together, (2) and (3) give what (1) denies: a cost-free handover to a read-only
consumer that protects itself.

## Alternatives considered

| Alternative | Advantages | Disadvantages | Evidence |
| --- | --- | --- | --- |
| A — One store, in the layout `calculation_cache` already uses, with read-only objects *(selected)* | Existing caches stay valid without being moved; one implementation of publication, verification and handover | Keeps two indexes (the project's and the catalog's) over one store | The report's measurements |
| B — A third, new store, migrating both current ones | A clean model | Against REQ-ARCH-002 (generalize, do not duplicate) and would invalidate the user's already computed caches | A.4, F-08 |
| C — Writable objects, as today | No permission change | An engine that rewrites its own input corrupts shared content, which is exactly what REQ-STO-003 forbids | — |
| D — Compressing the objects | Fewer bytes | PDB, MOL2 and PDBQT are read by native tools that expect the file; decompressing on every materialization trades space for time on every reuse, and the gain was not measured as necessary | 01.2 |
| E — A hard link for every consumer | Always free | An engine that writes in place corrupts the object; it is only safe for consumers proven read-only | A.2 |

## Selected solution

### Scope and root

The store belongs to the **project**:
`<project>/.biomolexplorer/store/objects`. A project remains a folder the user
copies, which is what supports the portability verified in PHASE-00 and the
`.bme.zip` export. Cross-project sharing stays as it was, through
`BIOMOL_PREPARATION_CACHE` and `BIOMOL_CACHE_DIR`, which set the root of
`calculation_cache`; 01.4 excludes *distributed* caching, not a shared local
directory the user chooses.

### Layout

`objects/<first two characters of the digest>/<digest>` — exactly the layout
`calculation_cache` already used. That is why an existing cache keeps working
without being moved: the implementation changed, the address did not.

### Permissions

An object is written in mode **0444**. Replacing it is governed by the
directory, not by the file, so atomic publication remains possible; writing to
it does not.

### Publication

`ContentStore.put` writes to a temporary file in the destination directory,
fsyncs, **verifies the digest of the written bytes**, applies the read-only
mode and only then renames, syncing the directory. An interruption at any point
leaves either nothing or a complete object — never a partial one. A concurrent
publication of the same content renames an identical object over an identical
one.

`validate` receives the bytes about to be published and refuses them by
raising: **content that fails the contract of its type is never published**.

### Verification

Full verification happens **on read and before every reuse** (REQ-STO-002):
`open`, `materialize` and `ArtifactManager.resolve` recompute the digest and
refuse a divergence rather than treating it as a hit. At publication, an object
already present is trusted — `calculation_cache` passes its own verification
function, with the filesystem stamp cache it already maintained, so its
behavior does not change.

### Materialization

| Mode | How | For whom |
| --- | --- | --- |
| `shared` | A hard link to the read-only object | A consumer that only reads. It costs no bytes and the kernel refuses anyone writing in place |
| `isolated` | An independent writable copy (copy-on-write where available), verified afterwards | An engine that rewrites its own input (REQ-STO-003) |

When a link is not possible (another device, a filesystem without links), the
`shared` mode makes a copy and **keeps it read-only**: the contract of the mode
still holds, the cost no longer does.

### Compression

None by default. The scientific formats are read by native tools that expect
the file; compressing would trade space for decompression on every
materialization. The decision is revisable if a measurement shows the need.

### Why `calculation_cache` is generalized rather than replaced

Because its checkpoints are the user's already computed scientific results.
Moving the storage would invalidate all of them. The change is one of
implementation: `copy_artifact` delegates to the store, `object_path` returns
the same path, and the tables and readers stay identical. The cost of that
decision is in the cache identity section of the report:
`calculation_cache.py` is read directly by
`stage_cache.implementation_digest`, so 16 of the 26 identities moved and were
registered as an audited equivalence.

## Rationale

Alternative B would be cleaner on paper and destructive in practice. C keeps
the risk REQ-STO-003 exists to remove. E is the obvious temptation on a system
without reflink, and it is safe only for a consumer that does not write —
which is why the mode exists, explicitly, rather than being the default.

What the decision adds to what was already there is not a store but the closing
of three gaps: the written bytes are now verified *by the store*, the object
now refuses writes, and the handover to a consumer is now a declared choice
between sharing and isolating.

## Trade-offs

- **Two indexes, one store.** The project index
  (`calculation_cache.objects`) and the workspace catalog (`content_objects`)
  describe content of the same store. Unifying the indexes belongs to PHASE-05.
- **Published results stay in the run folders.** In this phase managed storage
  holds the calculation checkpoints and whatever the manager registers
  explicitly. Bringing every published result into the store would duplicate
  bytes on a system without reflink; doing it by link would require proving
  every consumer read-only, which #9 itself conditions. That decision belongs
  to ADR-010 (PHASE-07), with retention.
- **A content contract is only as strong as its type's declaration.** `scores`
  declares no schema because its columns vary by engine, so its contract checks
  presence and a header, not named columns.
- **Verification on read costs a read.** Reusing fifteen thousand poses reads
  the bytes again to verify them. That is the price of never reusing divergent
  content, and it is what `calculation_cache` already paid.
- **The catalog is not rewritten.** An occurrence indexed by an earlier release,
  with no digest recorded, has its digest **computed** on resolve rather than
  written back.

## Scientific implications

No algorithm, protocol, tolerance or scientific byte changes. In particular:

- The digests, object paths, keys and tables of `calculation_cache` are the
  same, and all 26 cache identities returned to their previous values after the
  audited equivalence was registered: existing preparations, redockings and
  docking runs stay reusable.
- The handover becomes safer, not different: an engine that could previously
  have rewritten a shared object now receives an isolated copy or is refused by
  the kernel.
- Content that fails the contract of its type can no longer be published or
  reused. That is a new refusal, never a different result.

## Migration implications

Additive. Existing caches stay where they are, at the same addresses, and are
read by the new implementation. Objects already written stay writable — only
new ones are born read-only, and nothing is rewritten to change the mode of the
old ones.

An earlier release of the application opening the same project reads the
objects normally, because the layout did not change. `.bme.zip` packages do not
change format in this phase; including the store in the package belongs to
PHASE-07.

## Validation requirements

| Criterion | Verification |
| --- | --- |
| AC-02-01 | `test_content_store.StoreTests.test_identical_content_from_two_producers_is_stored_once`, `ArtifactManagerTests.test_two_runs_publishing_the_same_bytes_share_one_managed_object` |
| AC-02-02 | `test_a_writable_materialization_cannot_reach_the_stored_object`, `test_materialization_is_shared_by_default_and_isolated_on_request` |
| AC-02-04 | `test_usage_reports_the_store_and_what_the_catalog_claims` and the report's `R_storage` measurement |
| REQ-STO-002 | `test_a_tampered_object_is_refused_on_read_and_on_materialization`, `test_resolve_returns_the_managed_object_and_detects_tampering` |
| REQ-STO-004 | `test_publication_leaves_no_partial_object_when_content_changes`, `test_concurrent_publication_of_the_same_content_is_safe` |
| Content contracts | `ValidationContractTests` |

## Deviations from the reference specification

| Deviation | Reason |
| --- | --- |
| The materialization policy is not "reflink, then copy" but "link for a reader, isolated copy for a writer" | Measured: there is no reflink on this filesystem. The reflink attempt stays in the copy, for systems that have it |
| The store does not receive published results in this phase | Without reflink it would duplicate bytes; by link it would require proving every consumer read-only, which #9 itself conditions |
| `preparation_cache` is not routed here | #9 itself defers it to PHASE-03 |
