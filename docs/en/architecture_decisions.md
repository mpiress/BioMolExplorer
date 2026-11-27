# Architecture decision records

ADRs (Architecture Decision Records) document the decisions of the scientific
artifact architecture described in #6 and in sub-issues PHASE-00 to PHASE-08
(#7 to #15). Each decision records its context, alternatives, solution, rationale,
trade-offs, scientific and migration implications, validation requirements and
deviations from the reference specification.

## Process

- An ADR is created in the phase listed below, with status **Proposed**, from the
  [template](adr/0000-template.md).
- A decision becomes **Accepted** when its AC criteria have verifiable evidence.
  A superseded decision is never deleted: its status becomes
  **Superseded by ADR-XXX**.
- Numbers are permanent and never reused.
- Each ADR has a Portuguese version (`docs/adr/`) and an English one
  (`docs/en/adr/`).

## Index

| ADR | Decision | Phase | Status |
| --- | --- | --- | --- |
| ADR-001 | Scientific artifact identity model | PHASE-01 | Not started |
| ADR-002 | Persistent metadata catalog technology | PHASE-01 | Not started |
| ADR-003 | Content-addressable storage strategy | PHASE-02 | Not started |
| ADR-004 | Computational cache signature model | PHASE-03 | Not started |
| ADR-005 | Dependency tracking and invalidation | PHASE-04 | Not started |
| ADR-006 | Analytical dataset storage | PHASE-05 | Not started |
| ADR-007 | Unified data access interfaces | PHASE-05 | Not started |
| ADR-008 | Transactional publication and recovery | PHASE-02 | Not started |
| ADR-009 | Project migration and backward compatibility | PHASE-07 | Not started |
| ADR-010 | Artifact lifecycle and retention | PHASE-07 | Not started |
| ADR-011 | Visualization data retrieval | PHASE-06 | Not started |
| ADR-012 | Scientific metadata authority and provenance | PHASE-07 | Not started |

## Starting evidence

The PHASE-00 assessment is in
[docs/validation/artifact_architecture_assessment_2026-11-26](../validation/artifact_architecture_assessment_2026-11-26/README.md)
(in Portuguese, following the validation reports convention). It contains the
inventory, the baseline measurements, the requirement register and the
implementation plan that every ADR must take into account.
