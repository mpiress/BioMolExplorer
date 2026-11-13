# Docking consensus: mathematical model and interpretation

[Documentation portal](../../index_doc.html) · [Português](../docking_consensus.md) · [API parameters](backend_usage.md#api-consensus)

**Docking consensus** turns existing AutoDock Vina and DOCK6 observations into one priority list **per compound**. The implemented model is `target-balanced-ecr-v1`, in `consensus_model.py`, with fixed `TAU = 0.1`. This guide explains the current calculation; it does not change the algorithm.

**Reading the result:** raw score is an engine's energy-like output; rank is a compound's relative position; relevance transforms rank into a dimensionless factor; consensus combines the two engine factors. **Lower raw score is better; higher relevance/consensus is better.** Consensus is not an energy, binding probability, p-value or calibrated confidence estimate.

## 1. Purpose and scope

Vina and DOCK6 use different scoring functions. Their raw values should not be averaged as if identically calibrated. Target-specific score distributions may also differ. BioMolExplorer ranks compounds **separately for each engine/target**, then aggregates those ranks.

For 1,000 compounds evaluated against the same three targets by both engines, 3,000 Vina and 3,000 DOCK6 observations become **1,000 rows**, provided every compound has at least one target evaluated by both engines. No shared compound/target pair produces the existing no-intersection explanation rather than invented results.

Consensus does not perform docking, structure preparation, minimization or footprint generation. It reads existing scores and structural references.

## 2. Notation and required information

| Symbol | Meaning |
| --- | --- |
| i | Compound identifier, `molecule_chembl_id` |
| e | Engine: Vina or DOCK6 |
| t | Target identifier, `receptor_id` |
| x_iet | Lowest raw score for compound i, engine e, target t |
| N_et | Number of distinct observed compounds in that engine/target population |
| r_iet | One-based mean tie rank of x_iet |
| q_iet | Rank normalized to [0,1] |
| p_iet | Exponential relevance for that engine/target |
| T | Number of targets in the union of both inputs |
| c_ie | Observed fraction of the T targets for that compound/engine |
| R_ie | Engine relevance after target aggregation |
| C_i | Final consensus across the two engines |

Compound/target identifiers and finite scores must be meaningful. Matching uses textual identifiers: the model does not infer that `AChE`, `ACHE` and `4M0E_A` name the same target. Empty/ambiguous target labels prevent sound biological interpretation of that category.

Nonempty supplied isomeric SMILES are canonicalized. Different structures/stereoisomers under the same identifier raise an error. Missing SMILES do not trigger structure reconstruction or establish identity from the identifier alone. Distinct identifiers are not merged automatically.

### Minimal input-format example

To reproduce example 1, tables can reside at `Vina/docking_results.csv` and `Dock6/docking_results.csv` under the input directory. Engine names are `vina` and `dock6`.

Vina table:

```csv
molecule_chembl_id,canonical_smiles,receptor_id,engine,score
A,CCO,T,vina,-9
B,CCC,T,vina,-8
C,CCN,T,vina,-7
```

DOCK6 table:

```csv
molecule_chembl_id,canonical_smiles,receptor_id,engine,score
A,CCO,T,dock6,-30
B,CCC,T,dock6,-40
C,CCN,T,dock6,-20
```

These files demonstrate numeric factors without structural attachments. Real inputs must retain the corresponding files and provenance to open 3D/footprint views. **A missing observation means an absent pair row**, not a blank, NaN or infinite `score` cell; rows with nonnumeric/nonfinite scores are invalid.


## 3. Deduplication and reference populations

For each `(engine, target, compound)` key, retain the lowest raw score. Duplicate poses or repeated imports do not gain extra weight. Equal scores for the same key use structural-path ordering for a deterministic reference, without changing the numeric score.

Ranks use all observed compounds in each engine/target, **including compounds that may later be ineligible for the final output**. Final eligibility requires at least one target shared by both engines for that compound. Different engine populations can therefore produce different relative ranks for a shared compound.

T is the union of target identifiers, not only fully paired targets. A target observed in one engine alone still affects the coverage policy. Check completeness and identifier correspondence before interpreting rankings.

## 4. Ordinal normalization within each engine/target

Sort raw scores in ascending order. The best receives rank 1; ties receive the mean of the positions they occupy:

```text
r_iet = mean occupied position in ascending score order
q_iet = (r_iet - 1) / max(N_et - 1, 1)
```

For N>1, q=0 is best and q=1 worst. Scores `[-10, -8, -8, -6]` have ranks `[1, 2.5, 2.5, 4]` and q values `[0, 0.5, 0.5, 1]`. Compound identifiers do not break energy ties at this stage.

- **Singleton:** N=1, r=1, q=0; this only means no competitor was observed.
- **Constant scores with N>1:** r=(N+1)/2 and q=0.5 for every compound; no discrimination is possible.
- **Positive score:** remains in ascending ordering; no zero-clipping or absolute-value transform is applied.
- **Missing observation:** receives no fabricated raw energy; target aggregation handles it explicitly.
- **NaN/infinity:** invalid nonfinite scores are rejected.

Strictly increasing score transformations preserve ranks. This reduces dependence on energy units/scales and extreme magnitudes but discards energy-gap information. Adding/removing competitors can change N and normalized ranks. Relevance factors from different screening libraries are not automatically comparable.

## 5. Exponential relevance and the ECR connection

The published ECR method combines exponential rank contributions through summation. Its individual term has the form `(1/σ) × exp(-r/σ)`. BioMolExplorer adapts that idea through rescaled ranks and dimensionless relevance, not a literal reproduction of the original density. Its target aggregation is an additional platform policy. [Primary ECR paper](https://doi.org/10.1038/s41598-019-41594-3).

```text
τ = 0.1
p_iet = exp(-q_iet / τ)
```

For N>1, equivalently `p = exp(-(r-1)/(τ × (N-1)))`: the best rank receives 1, and the original ECR density prefactor is not used.

τ is fixed in the current implementation, not an editable form control or a cutoff discarding 90% of molecules. Every eligible compound remains in the output; exponential decay emphasizes better relative ranks.

| q | p = exp(-q/0.1) |
| --- | --- |
| 0.0 | 1.000000 |
| 0.1 | 0.367879 |
| 0.2 | 0.135335 |
| 0.5 | 0.006738 |
| 1.0 | 0.000045 |

Relevance values do not necessarily sum to 1 and are not probabilities. A singleton receives p=1 without demonstrating affinity. Constant scores with N>1 produce `exp(-5) ≈ 0.006738`, not zero. [Pose/rank consensus research](https://pubs.rsc.org/en/content/articlehtml/2021/ra/d1ra05785e) illustrates a 10%-of-library exponential scale; this is not proof of universal superiority for the current adaptation.

## 6. Target aggregation and missing observations

For each compound/engine, count observed targets. A missing target receives **q*=1**, the worst relative rank, without an imputed raw score:

```text
q*_iet = q_iet when observed, otherwise 1
c_ie = observed target count / T
R_ie = c_ie × exp(-sum_t(q*_iet) / (T × τ))
```

Equivalently:

```text
p*_iet = exp(-q*_iet / τ)
R_ie = c_ie × (product_t p*_iet)^(1/T)
```

The implementation sums ranks in the exponent rather than directly multiplying many small values. Targets have equal weight. With complete coverage, R is the geometric mean of target relevance. Missing data both contribute a worst rank and reduce the coverage multiplier.

**This is not an explicit variance or worst-target penalty.** R depends on mean q* and coverage. Different target profiles with equal means and coverage receive equal R; no minimum performance on every target is guaranteed. Example 4 demonstrates this limit.

## 7. Combining Vina and DOCK6

The final factor is an equally weighted arithmetic mean:

```text
C_i = (R_i,Vina + R_i,DOCK6) / 2
```

Higher C means higher priority under this policy. There is no subsequent normalization, hypothesis test or explicit engine-disagreement penalty. A strong engine result can compensate for a weak one: `(0.9 + 0.1)/2` and `(0.5 + 0.5)/2` both equal 0.5. Consensus alone does not establish pose agreement or agreement in physical energies.

Output sorts by descending C, then ascending textual compound identifier for exact ties. The implementation does not estimate confidence intervals, standard errors or statistical significance for C.

## 8. Worked examples

### Example 1 — one target and three compounds

| Compound | Vina raw score | DOCK6 raw score |
| --- | --- | --- |
| A | -9 | -30 |
| B | -8 | -40 |
| C | -7 | -20 |

Both engine populations have N=3:

| Compound | Vina q | DOCK6 q | Vina R | DOCK6 R | Consensus C |
| --- | --- | --- | --- | --- | --- |
| A | 0.0 | 0.5 | 1.000000 | 0.006738 | 0.503369 |
| B | 0.5 | 0.0 | 0.006738 | 1.000000 | 0.503369 |
| C | 1.0 | 1.0 | 0.000045 | 0.000045 | 0.000045 |

A and B tie; textual ordering puts A first. Three compounds produce three rows. The UI displays C's factor as `0.000`, although its internal value is positive. Raw `-9` and `-30` were ranked in different engine populations rather than directly compared.

### Example 2 — AChE, BChE and BACE1

Assume 11 compounds per engine/target, allowing normalized-rank steps of 0.1. Compound X has:

| Target | Vina q | Vina p | DOCK6 q | DOCK6 p |
| --- | --- | --- | --- | --- |
| AChE | 0.0 | 1.000000 | 0.1 | 0.367879 |
| BChE | 0.1 | 0.367879 | 0.2 | 0.135335 |
| BACE1 | 0.2 | 0.135335 | 0.3 | 0.049787 |

Both engines have complete coverage:

```text
R_X,Vina  = exp(-(0.0 + 0.1 + 0.2)/(3 × 0.1)) = exp(-1) = 0.367879
R_X,DOCK6 = exp(-(0.1 + 0.2 + 0.3)/(3 × 0.1)) = exp(-2) = 0.135335
C_X       = (0.3678794412 + 0.1353352832)/2 = 0.2516073622
```

The displayed values are approximately `0.368`, `0.135` and `0.252`. X's three target evaluations become **one row**; target-specific details remain in `target_scores`.

### Example 3 — one missing target

A compound has q=0 for AChE/BChE but no BACE1 observation in one engine:

```text
q* = [0, 0, 1]
c = 2/3
R = (2/3) × exp(-(0 + 0 + 1)/(3 × 0.1)) = 0.0237826622
```

If both engines share this pattern, C is also 0.0237826622. If the other engine instead has complete coverage and q=0 throughout, its R is 1 and C is 0.5118913311. Missing data are penalized but do not automatically exclude a compound with at least one shared target. A compound evaluated by both engines only on disjoint targets is ineligible.

### Example 4 — breadth-policy limitation

| Complete-coverage profile | Target q values | Mean q | R |
| --- | --- | --- | --- |
| Balanced | 0.3, 0.3, 0.3 | 0.3 | 0.049787 |
| Uneven | 0.0, 0.0, 0.9 | 0.3 | 0.049787 |

There is no additional dispersion penalty. To require minimum relevance for every target, inspect target details and apply a separate, justified scientific selection criterion; consensus does not silently impose that rule.

## 9. Final table and traceability

| Field | Meaning |
| --- | --- |
| `molecule_chembl_id` | One reference per eligible compound |
| `vina` | Global Vina R, not raw energy |
| `dock6` | Global DOCK6 R, not raw energy |
| `normalized_score`, `consensus` | The same final C |
| `consensus_model` | `target-balanced-ecr-v1` |
| `target_count` | T, the target union size |
| `paired_target_count` | Targets with both engine observations for that compound |
| `vina_target_coverage`, `dock6_target_coverage` | Engine-specific coverage in [0,1] |
| `target_scores` | JSON target-specific raw scores, rank fractions and relevance; missing values are null |
| `receptor_id`, `vina_raw`, `dock6_raw` | Representative saved context, not global results |
| Structural/footprint fields | Existing attachments for result actions |

The representative target is selected among the compound's paired targets by the greatest sum of engine relevance at that target, with textual target identifiers breaking ties. This supplies saved structural context and does **not** change global factors.

The three relevance columns display three decimal places while calculation, sorting and CSV retain precision. Existing buttons remain: multiple targets first open a searchable target combobox for 3D Vina, 3D DOCK6 and footprint; a single target opens directly. Attachments must match the compound and target. Missing attachments do not justify substituting another target or rerunning docking to create one.

## 10. Alternatives and what is not calculated

The active model uses relative order rather than raw energetic scale. Min-max depends on observed extremes, z-score on population mean/dispersion, and direct energy sums require engine calibration. These are not calculated by the current model. Rank-based exponential normalization is not universally best for every biological setting.

`repulsion_weight` is a legacy no-op in `consensus_rows`; it neither changes engine weights nor corrects the current factor for repulsion. τ and equal engine weights are fixed. No hidden bonus for stereochemistry, lipophilicity, QED, toxicity or similarity is included. Structural identity checking is not an extra stereochemical score.

## 11. Scientific interpretation checklist

1. Check compound/target identifiers, score direction and reference populations.
2. Inspect target count, coverage and shared targets; missing data are not experimental low-affinity measurements.
3. Use unrounded factors for automatic selection; `0.000` may be a positive small value.
4. Inspect per-target details and poses. The global factor does not establish selectivity, contacts or geometric quality.
5. Assess predictive utility using appropriate actives/decoys or experimental references and predefined enrichment criteria. This block does not automatically perform that validation.
6. Do not interpret factors across changed target/library populations as directly calibrated, or translate C into clinical success probabilities or absolute affinity.

## 12. Reproduce aggregation with standard Python

This example receives **already normalized q values**, not raw energies. It reproduces example 2 using only the standard library; it is not a new API and does not run docking. Both engine vectors must represent the same target union in the same order; use `None` for missing observations.

```python
from math import exp, fsum

def relevance(q_values):
    if not q_values:
        raise ValueError("The target union must not be empty")
    if any(q is not None and not 0 <= q <= 1 for q in q_values):
        raise ValueError("Supply normalized ranks between 0 and 1")
    coverage = sum(q is not None for q in q_values) / len(q_values)
    average_rank = fsum(1.0 if q is None else q for q in q_values) / len(q_values)
    return coverage * exp(-average_rank / 0.1)

vina = relevance([0.0, 0.1, 0.2])
dock6 = relevance([0.1, 0.2, 0.3])
consensus = (vina + dock6) / 2
print(vina, dock6, consensus)
```

Approximate output: `0.3678794412 0.1353352832 0.2516073622`. The four examples and their coverage variants were checked against the current `consensus_rows` function without running Vina/DOCK6 or changing application code. This is arithmetic documentation verification, not predictive or experimental validation.

## 13. Sources and implementation checks

- [ECR, Palacio-Rodríguez et al. (2019)](https://doi.org/10.1038/s41598-019-41594-3): the exponential-rank consensus foundation.
- [Pose/rank consensus study (2021)](https://pubs.rsc.org/en/content/articlehtml/2021/ra/d1ra05785e): an example of a library-relative exponential scale.
- `consensus_model.py`, `consensus_io.py` and `wrappers/docking.py`: authoritative sources for the adaptation, deduplication, coverage, export and finalization behavior.

Normalization to [0,1], geometric target aggregation and coverage penalties are explicit BioMolExplorer policies, distinguished here from the original published ECR formula. See the [Portuguese guide](../docking_consensus.md) for dated performance evidence; those timings were not remeasured by this documentation update.
