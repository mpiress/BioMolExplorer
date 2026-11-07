# Artifact refactor status

Updated: 2026-11-06

This file records verified progress for the semantic artifact migration. A
step is marked complete only after its focused regressions pass.

| Step | Status | Evidence |
| --- | --- | --- |
| Preserve future semantic references in the editor | Complete | `InputEditor.read()` now preserves port, selection, current/pinned mode and artifact identity; focused reopen regression passes. |
| Resolve current and pinned references through execution manifests | Complete | Central resolver handles current ports, entity indexes and pinned/shared occurrences; generic, graph, docking and MD adapters consume `_resolved_paths`. |
| Automatically index published entities | Complete | Transactional publication indexes each occurrence/content revision once and records index errors without failing science; focused registry test passes. |
| Complete scientific contracts and semantic selectors | Complete | Required/conditional inputs, batch cardinality, schemas and metadata are declared; new UI links show producer port/mode and hide physical selectors. |
| Integrate specialized resolvers | Complete | Graph, docking, redocking helpers and MD consume manifest-resolved paths before applying their scientific adapters. |
| Fix cache and deferred-configuration regressions | Complete | Cross-block historical cache aliasing and nested async preparation deadlock fixed; both named regressions pass. |
| Migrate legacy/imported projects | Complete | Startup/save/import migration is idempotent, backfills imported runs and records incompatible/ambiguous references without guessing. |
| Alzheimer integration coverage | Complete (architectural) | Small coherent AChE/BChE/BACE1 DAG tests cover branches, convergence, future 4M0E/chain selection, manifests, multiple compounds and reopen persistence. External engines are not scientifically validated by mocks. |
| Portuguese and English documentation | Complete | Architecture guides document contracts, registry, selection, cache/provenance and migration; generated HTML validation passes. |

Consolidated validation: 162 tests passed across semantic artifacts, Alzheimer DAGs, registry/publication, cache, file selection, graphs, project portability/import, pipeline execution, strict ports, docking sources, MD/Amber and redocking. Documentation — 34 generated pages validated. `py_compile` and `git diff --check` pass.

Scientific engine acceptance remains separate from architectural tests. Mocked
or fixture-based docking tests must not be reported as validation of Vina,
DOCK6, AmberTools, or Amber Score.
