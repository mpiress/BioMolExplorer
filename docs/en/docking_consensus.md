# Compound-level docking consensus

Model: `target-balanced-ecr-v1`, implemented in `consensus_model.py`. Each eligible compound produces one output row regardless of the number of targets. This operation reads existing Vina/DOCK6 scores; it does not rerun docking, preparation or footprint calculations.

## Normalization

For each engine e, target t and compound i, retain the lowest finite raw score. Lower raw energy ranks first; neither absolute values nor zero-clipping are applied. Equal energies receive the mean tie rank. With one-based rank r and N observed compounds:

```text
q_iet = (r_iet - 1) / max(N_et - 1, 1)
p_iet = exp(-q_iet / 0.1)
```

A singleton receives q=0 and p=1. If all scores are equal and N>1, every compound receives q=0.5. These factors are relative rankings, not calibrated binding probabilities. There is no universally optimal normalization for every biological target.

## Multiple targets and engine aggregation

T is the union of target identifiers in both inputs. For each engine, coverage c is the fraction of targets with an observation. Missing observations receive worst rank q*=1, without imputed raw energy:

```text
R_ie = c_ie * exp(-sum_t(q*_iet) / (T * 0.1))
C_i = (R_i,Vina + R_i,DOCK6) / 2
```

With complete coverage, R is the geometric mean of target-specific relevance. Higher consensus is better. A compound is eligible only if at least one compound/target pair occurs in both engines. Partial target coverage is retained and penalized. With 1,000 compounds evaluated on the same three targets by both engines, the output has exactly 1,000 rows. No common pair yields the existing skipped-result explanation.

Canonical isomeric SMILES are checked across observations: the same identifier associated with different structures is rejected. Different identifiers are not merged merely because structures resemble each other. Missing stereochemistry is not inferred. RMSD, affinity experiments and toxicity are separate validation criteria, not hidden weights.

## Results and actions

| Display column | CSV field | Meaning |
| --- | --- | --- |
| Compound | `molecule_chembl_id` | One row per eligible identifier |
| Vina relevance | `vina` | Engine relevance across targets |
| DOCK6 relevance | `dock6` | Engine relevance across targets |
| Consensus | `normalized_score`, `consensus` | Mean of the engine relevance factors |

The three relevance columns display three decimal places; computation and CSV data retain precision. `target_count`, `paired_target_count`, coverage fields and the JSON `target_scores` preserve all target-specific raw scores, rank fractions and relevance. `receptor_id`, `vina_raw` and `dock6_raw` describe a representative target only; they are not the global relevance values. `repulsion_weight` remains a legacy no-op argument.

For multiple targets, **3D Vina**, **3D DOCK6**, residue and footprint actions first request a target through a searchable combobox. Only the selected compound/target attachment is resolved and validated. Missing or changed files are not replaced by the representative target. With one target, actions open directly. Removing rows remains an audited curation operation.

## Performance and persistence

The implementation caches SMILES parsing within the calculation, deduplicates observations, materializes only representative structural references and resolves other target attachments lazily from authorized provenance. Published tables are authoritative; removed rows must not reappear from auxiliary tables. Cache compatibility includes the scientific model identity. Engine results and molecular preparation caches are not recalculated by consensus.

For the full derivation and historical performance evidence, see the [Portuguese model guide](../docking_consensus.md). The implemented mathematical policy adapts [Exponential Consensus Ranking](https://doi.org/10.1038/s41598-019-41594-3); geometric target breadth and coverage penalties are explicit application policies, not evidence of universal predictive superiority.
