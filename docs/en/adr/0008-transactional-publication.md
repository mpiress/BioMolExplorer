# ADR-008 — Transactional publication and recovery

- **Status:** Accepted
- **Date:** 2026-11-28
- **Phase:** PHASE-02 (#9)
- **Requirements:** REQ-STO-004, REQ-STO-005, REQ-INT-002, REQ-INT-003, REQ-TEST-003, AC-02-03
- **Authors and reviewers:** Michel

## Context

Publishing an artifact touches three things that have to agree: the **bytes**
in the store, the **rows** of the catalog and the **run state**. Each already
had its own isolated guarantee:

| Part | Existing guarantee |
| --- | --- |
| Bytes | `storage.py` writes to a temporary file with fsync and rename; `content_store.put` adds digest verification before the rename |
| Rows and run state | A single SQLite transaction: `artifact_registry.register_run` runs inside the caller's transaction, which is the same one that writes `runs.stages` |
| File moves | `result_transactions.ResultTransaction` — an on-disk journal plus a marker row in `result_file_transactions`, with `recover` at workspace open |

What was missing was the relation between them. An interruption can land on
either side of the commit, and without a declared order both bad states are
possible: a row pointing at absent bytes, or bytes nothing references.

Of those two, only one is acceptable. A dangling reference breaks the reading
of a scientific result; unreferenced bytes are inert.

## Alternatives considered

| Alternative | Advantages | Disadvantages | Evidence |
| --- | --- | --- | --- |
| A — A declared order (bytes before rows) plus a journal of the content introduced *(selected)* | A dangling reference becomes impossible by construction; recovery only has to handle the inert case | An interruption can leave bytes to clean up | The report's crash-injection tests |
| B — A two-phase commit of its own, with an intention log over all three parts | Full symmetry | Reimplements what SQLite's transaction already gives for two of the three parts, and adds a log of its own to recover | — |
| C — Rows before bytes | A failure leaves the store clean | It leaves exactly the unacceptable state: the catalog promising absent content | 01.3 |
| D — Nothing; rely on re-execution | Simple | A scientific re-execution costs hours, and REQ-INT-003 requires resuming from validated outputs | A.3 |

## Selected solution

### The invariant

> **The bytes reach the store before any row references them.**

It follows that no catalog row can point at absent content. The only possible
residue of an interruption is content nothing references.

### The journal

`ResultTransaction.record_objects` names the content a transaction is about to
introduce, **before** it is written, in the same journal that already recorded
file moves. The journal gains an `objects` key.

The marker row in `result_file_transactions` remains what distinguishes a
finished transaction from an interrupted one: it is inserted in the caller's
transaction, so it exists if and only if that transaction committed.

`ArtifactManager.publish` is the whole sequence: journal, bytes, rows, commit —
and the rows commit with the run state, because they are the same transaction.

### The recovery

`result_transactions.recover`, already called when the workspace opens:

| Interruption | State found | What recovery does |
| --- | --- | --- |
| Before the commit | Journal present, marker absent | Restores the moves and discards the content named in the journal **that nothing references** |
| After the commit, before cleanup | Journal and marker present | Finishes the cleanup and **keeps** the content: the rows reference it |
| During the catalog transaction | SQLite rolls the rows back | The same as the first case; the verified bytes stay in the store |
| After everything | Neither journal nor marker | Nothing to do |

### Recovery is conservative

An object is discarded only when **no** logical artifact of the project claims
that content and **no** registered content object points at it. REQ-STO-005
forbids deleting scientific content automatically, and catalog consistency is
never sufficient reason to delete bytes.

The honest consequence: an interrupted publication whose content belongs to a
still-published artifact leaves managed bytes in the store without the row that
registers them. They are correct content, the next publication deduplicates
onto them, and sweeping orphans at scale belongs to PHASE-07 (REQ-INT-001 and
REQ-INT-002). `ContentStore.sweep` already exists for when it arrives.

## Rationale

The asymmetry between the two bad states is what decides. Since only one of
them is tolerable, ordering the writes is enough to make the other impossible —
and that is stronger than any protocol that tries to repair both afterwards.
The journal does not exist to guarantee consistency; it exists so the cleanup
of the tolerable case is possible instead of permanent.

Alternative B would reimplement, with a new log to recover, what SQLite's
transaction already delivers for the rows and the run state.

## Trade-offs

- An interruption can leave bytes to clean up, and the cleanup is conservative
  by decision.
- The journal is per project and per transaction; it does not replace an
  integrity check of the whole store, which belongs to PHASE-07.
- `publish` opens a transaction of its own. A publication made inside a
  transaction the caller already opened passes `db` and joins it, but then the
  journal is the caller's responsibility.

## Scientific implications

None. The decision computes, converts or discards no result. What it guarantees
is that an already computed result does not become unreadable through an
interruption, which is a traceability requirement (01.3), not a computational
one.

## Migration implications

The `objects` key is additive in the journal:

- a journal written by this release and read by an earlier one is read
  normally, because the earlier release reads `moves` and ignores keys it does
  not know;
- a journal written by an earlier release and read by this one has no
  `objects`, and recovery treats the list as empty.

No new table, no rewritten row.

## Validation requirements

| Criterion | Verification |
| --- | --- |
| Interruption before the commit | `test_content_store.TransactionalPublicationTests.test_an_interruption_before_the_commit_leaves_no_dangling_reference` |
| Cleanup of unclaimed content | `test_an_interruption_before_the_commit_discards_unclaimed_content` |
| Interruption after the commit | `test_an_interruption_after_the_commit_keeps_the_published_content` |
| Catalog transaction rolled back | `test_a_rolled_back_catalog_transaction_leaves_the_store_consistent` |
| A complete publication leaves no residue | `test_a_complete_publication_leaves_no_journal_and_no_marker` |
| REQ-STO-004 | `test_publication_leaves_no_partial_object_when_content_changes` |

## Deviations from the reference specification

None. #9 asks that the content write, the catalog rows and the run-state commit
recover consistently on either side of the commit, and that is what the
invariant plus the journal deliver. Sweeping orphans at scale stays where #6
itself places it, in PHASE-07.
