# Scientific and methodological documentation

This reference documents the scientific behavior of BioMolExplorer 2.2.0 as implemented in the current source. It is the English companion to the fully detailed [Portuguese scientific reference](../scientific_documentation.md); equations, defaults, limitations, implementation paths and bibliography are identical.

## Implemented method map

| Area | Exact active implementation |
| --- | --- |
| ADMET descriptors | RDKit MW, TPSA, Wildman-Crippen logP, HBD/HBA, rotatable bonds, QED, fraction Csp3, rings, charge, fragments and aromatic proportion in `src/caad/admet.py` |
| Screening rules | Exact Lipinski, Veber and extended rule-of-three violation counts; ESOL-style and Potts-Guy-style equations documented in the Portuguese reference |
| BBB/HIA/P-gp | Explicit local ellipse/penalty heuristics; labels are not probabilities |
| ADMET-AI | ADMET-AI 2.x, `include_physchem=False`, no DrugBank percentile, batches of 256, raw provider endpoints prefixed `ML_`, no generic probability threshold |
| Morgan | radius 2, 2048 bits by default, chirality and feature invariants enabled, bond types disabled |
| MACCS | RDKit public 166-key definition in a 167-position vector (bit zero unused) |
| Pharmacophore | RDKit Gobbi Pharm2D topological feature-pair/triplet fingerprint |
| Similarity | Tanimoto, Dice, Cosine, Sokal, Russel, Rogot-Goldberg, AllBit, Kulczynski, McConnaughey, Asymmetric and Braun-Blanquet; values rounded to 3 decimals before thresholding |
| Consensus | target-balanced ECR rank quantiles, `tau=0.1`, explicit target coverage and missing-rank penalty; no z-score or min-max |
| DOCK 6 footprint | Engine-produced residue van der Waals/electrostatic footprint and reference FPS; BioMolExplorer validates and displays rather than recomputes energies |
| 3D contacts | Conservative geometry for H-bond, hydrophobic, pi-parallel, pi-T and van der Waals contacts; no salt-bridge classifier |
| Networks | Unweighted undirected degree, centralities, PageRank 0.85, k-core, clustering, triangles and seeded Louvain |
| Pose validation | PoseBusters-inspired battery `posebusters_inspired_v1`: chemistry, connectivity, coordinates, internal geometry and planarity against an MMFF94/UFF relaxed conformer, internal strain and Bondi van der Waals clash at 0.75 of the radii sum |
| Binding-site surface | Solvent-excluded surface colored by the Kyte-Doolittle residue hydrophobicity scale interpolated between -4.5 and +4.5; a representation, not an energy, buried-area or druggability calculation |

## Core equations

For active-bit sets `A` and `B`, with `a=|A|`, `b=|B|`, `c=|A intersection B|`: `Tanimoto=c/(a+b-c)`, `Dice=2c/(a+b)`, and `Cosine=c/sqrt(ab)`.

For each engine/target group, tied observations occupying zero-based indices `s` through `e-1` among `N` observations receive:

```text
q = ((s + e - 1) / 2) / max(N - 1, 1)
relevance = exp(-q / 0.1)
engine_factor = coverage * exp(-sum_target(q) / (T * 0.1))
consensus = (vina_factor + dock6_factor) / 2
```

Missing target ranks use `q=1` in the penalty. At least one target-compound pair must be present in both engines. `repulsion_weight` is currently a compatibility no-op. This is a ranking policy, not a calibrated binding probability.

The implemented local screening equations are:

```text
ESOL_logS = 0.26121066137801696 - 0.0066138847738667125 MW
            - 0.7416739523408995 WLOGP + 0.003451545565957996 RB
            - 0.42624840441316975 AromaticProportion
logKp = -6.3 + 0.71 WLOGP - 0.0061 MW
```

## 3D geometry and explicit omissions

Hydrogen bonds require `D-A <= 4.1 A`, explicit `D-H` between 0.6 and 1.25 A, and `D-H-A > 100 degrees`. Hydrophobic contacts use 4.0 A. Aromatic centers use 5.5 A, offset at most 2 A, and normal angles at most 30 degrees (parallel) or at least 60 degrees (T-shaped). van der Waals gaps must lie between -0.4 and 0.5 A relative to RDKit radii.

Pose validation measures the delivered geometry: relative bond deviation 0.25, absolute angle deviation 25 degrees, planarity 0.25 A, internal strain 100 kcal/mol and a steric ratio of 0.75 of the van der Waals radii sum, each configurable and recorded with the report. Criteria that depend on bond orders are skipped for proximity-bonded poses, and the steric test is skipped without an associated receptor; both conditions are recorded instead of failing the pose. Internal strain is not a binding or free energy. See [pose validation](pose_validation.md).

The current repository does not implement PCA, t-SNE, UMAP, AUC-ROC, enrichment factor, salt bridges, cation-pi, halogen bonds or metal coordination. These capabilities must not be inferred from generic scientific terminology.

## References and full traceability

The complete bibliography and method-by-method interpretation are maintained in the [Portuguese scientific reference](../scientific_documentation.md#referencias-fundamentais). Source anchors include `src/caad/admet.py`, `src/biomolexplorer/admet_predictors.py`, `src/kernel/descriptors.py`, `src/biomolexplorer/consensus_model.py`, `src/biomolexplorer/footprints.py`, `src/biomolexplorer/docking_interactions.py`, `src/biomolexplorer/pose_geometry.py`, `src/biomolexplorer/molecular_surface.py` and `src/biomolexplorer/network_analysis.py`.
