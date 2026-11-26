# Pose validation: physical and geometric criteria

[Documentation portal](../../index_doc.html) · [Português](../pose_validation.md) · [API parameters](backend_usage.md#api-pose-validation)

**Validate poses** assesses the physical and geometric plausibility of poses produced by AutoDock Vina, DOCK6 or the consensus block before they reach Amber Score, the consensus or a report. It **does not dock, reposition, minimise the delivered pose or change scores**: every pose is measured as received, and the result is a per-pose verdict.

The criteria follow the PoseBusters formulation (Buttenschoen, Morris and Deane, *Chemical Science* **15**, 3130–3139, 2024), implemented in `pose_geometry.py` and identified as `posebusters_inspired_v1`. The implementation is our own, uses RDKit/NumPy/SciPy and does not run the original package; the report therefore records the criteria set identifier rather than a PoseBusters certification.

**Reading the result:** an approved pose satisfied every evaluated criterion. A rejected pose carries the criteria that failed and their measurements. Approval means the delivered geometry is physically plausible — **not** affinity, biological activity or a correct binding mode.

## 1. Inputs and scope

| Input | Accepted content | Use |
| --- | --- | --- |
| `base_input_path` | Receptors prepared by redocking or preparation (`receptor_pdb`, `prepared_structures`, `structures`) | Steric context of the intermolecular test |
| `base_selected_mols` | Poses and tables from Vina, DOCK6 or consensus (`vina`, `dock6`, `molecular_poses`, `scores`) | Poses to validate |

Poses arrive from the connected blocks. The editor lists targets by the `receptor_id` of the prepared receptors and compounds by the molecules that actually produced docking results, read from the artifact registry and the published score tables — the form never loads execution JSON and never walks the pose tree. With no explicit selection, every received pose is validated.

| Selection | Effect |
| --- | --- |
| Selected targets (`selected_receptors`) | Restricts validation to the chosen `receptor_id` values |
| Selected ligands (`selected_ligands`) | Restricts to the chosen compounds, by `compound_id` or `ligand_id` |
| Removed complexes (`excluded_pairs`) | Removes one complex as `compound:receptor`, without removing that compound on other targets |

A standalone CLI run resolves poses from the folder given in `base_selected_mols` or `base_input_path`; the same exclusions are applied over the resolved poses.

## 2. Preparing the measured molecule

Each pose file is read with whatever bond orders are available:

| Format | Reader | Bond orders |
| --- | --- | --- |
| `.sdf` | `SDMolSupplier` without automatic sanitization | Explicit |
| `.mol2` | `MolFromMol2File` without automatic sanitization | Explicit |
| `.pdbqt` | Internal conversion to PDB, selectable model | By proximity |

The distinction matters: without explicit bond orders every bond would be read as single, and the force field and the planarity test would measure a different molecule. The internal-geometry, planarity and energy criteria are therefore **evaluated only when bond orders do not come from proximity**. The report's `bond_orders` field records which reading was used, and a criterion without a measurement is reported as unavailable instead of silently failing the pose.

## 3. Evaluated criteria

| Criterion | Measurement | Rejection |
| --- | --- | --- |
| Chemistry (`valid_chemistry`) | RDKit sanitization of the molecule as read | Invalid valence or aromaticity |
| Connectivity (`valid_connectivity`) | Number of fragments | More than one disconnected fragment |
| Coordinates | Finite and distinct positions | Missing, repeated or non-finite coordinates |
| Internal geometry (`valid_geometry`) | Largest relative bond-length deviation and largest absolute angle deviation against the relaxed conformer | Deviation above `bond_tolerance` or `angle_tolerance` |
| Planarity (`valid_planarity`) | Largest distance to the fitted plane of aromatic rings and non-ring double bonds | Deviation above `planarity_tolerance` |
| Internal energy (`valid_energy`) | ΔE between the docked pose and the same ligand minimised | ΔE above `delta_energy_limit` |
| Steric clash (`valid_clash`) | Smallest ratio between the observed distance and the sum of van der Waals radii, over heavy atoms of pose and receptor | Ratio below `vdw_clash_scale` |

### 3.1 Relaxed internal reference

The geometric and energetic reference is a **copy of the same connectivity**, minimised with MMFF94 and, when the field does not parameterise the molecule, with UFF. The effective method is recorded in `energy_method`. The original pose is never modified: minimisation runs on the copy.

ΔE = E(pose) − E(relaxed conformer), in kcal/mol, with the same force field on both evaluations. ΔE measures the **internal strain of the delivered conformation**; it is not a binding energy, not a free energy and not an interaction term with the receptor.

### 3.2 Internal deviations

The bond deviation is relative: |d_pose − d_relaxed| / d_relaxed, taking the largest value in the molecule. The 0.25 default accepts deviations of up to 25% per bond. The angle deviation is absolute, in degrees, with a 25° default. Both compare the pose with the relaxed conformer of the same molecule, not with tabulated values.

### 3.3 Planarity

For every aromatic ring and every non-ring double bond, the group plane is fitted and the largest atom distance to that plane is measured. The number of evaluated groups appears in `planar_groups`.

### 3.4 Intermolecular steric clash

The test uses Bondi van der Waals radii (unlisted elements use the carbon radius, 1.70 Å). Solvent is discarded, and the site's own crystallographic ligand is excluded through `exclude_residues`, so the crystallographic reference is not counted as a clash. The decision uses the smallest distance/radii-sum ratio of the set; the report also carries the absolute overlap of that worst pair (`min_vdw_overlap`), the closest heavy-atom distance (`min_heavy_distance`) and how many pairs fell below the sum of their radii (`clashing_contacts`).

The 0.75 default keeps hydrogen bonds and tight polar contacts valid, as in the PoseBusters intermolecular distance criterion. **Without an associated receptor the criterion is not evaluated**: the pose is not rejected for that, and the report records the absence in `resolved_by` and in the summary (`unresolved_inputs`).

## 4. Tolerances and defaults

| Parameter | Default | Accepted domain |
| --- | --- | --- |
| `strict_filter` | `true` | boolean |
| `vdw_clash_scale` | `0.75` | 0.3 to 1.0 |
| `delta_energy_limit` | `100.0` | 0 to 100000 kcal/mol |
| `bond_tolerance` | `0.25` | 0 to 2 (fraction) |
| `angle_tolerance` | `25.0` | 0 to 180 degrees |
| `planarity_tolerance` | `0.25` | 0 to 5 Å |
| `minimization_steps` | `400` | 1 to 20000 iterations |

With `strict_filter` enabled, only approved poses are handed to the next blocks. Disabled, every pose remains available with the validation metadata and the decision stays with whoever reads the report. The applied tolerances are written into the run summary, so an older report stays interpretable with the limits that produced it.

## 5. Published results

| File | Content |
| --- | --- |
| `pose_validation.csv` | One row per evaluated pose, with identity, engine, source score, per-criterion verdicts, measurements, rejection reasons and source files |
| `pose_validation.json` | Summary: target, total poses, approved, rejected, published, strict filter, counts per reason, applied tolerances, unresolved inputs and criteria identifier |
| `poses/` | Copies of the poses handed to the next blocks |
| `receptors/` | Receptors used in the steric test, as PDB |

Output ports are `validated_poses` (`molecular_poses`, `vina`, `dock6`), `validation_report` (`scores`) and `receptors` (`receptor_pdb`). Every published pose keeps its own file: when the engine names a pose after the compound alone — as DOCK6 does — the copy is qualified by the receptor, so three targets of the same compound never collapse into a single file.

Report ordering: rejected first, then by target, compound and pose. The order exists to put the poses that need a decision first.

## 6. Report columns and their origin in the reference work

### 6.1 Tests defined by the work

PoseBusters groups validity into three categories. The thresholds below are the work's, not this implementation's.

| Category | Test | Definition and threshold in the work |
| --- | --- | --- |
| Chemical | File loads | RDKit can read the molecule |
| Chemical | Sanitisation | The molecule passes RDKit sanitisation: valence, aromaticity, radicals and hybridisation |
| Chemical | Molecular formula | Matches the true ligand, compared through standard InChI, without isotopes and with charges neutralised where possible |
| Chemical | Bonds | Match the true ligand through the InChI connectivity layer |
| Chemical | Tetrahedral chirality | Matches the true ligand through the InChI stereo layer |
| Chemical | Double bond stereochemistry | Matches the true ligand through the InChI double-bond layer |
| Intramolecular | Bond lengths | Between 0.75 of the lower bound and 1.25 of the upper bound of the RDKit distance-geometry bounds (25% tolerance) |
| Intramolecular | Bond angles | The same 25% tolerance over the distance-geometry bounds |
| Intramolecular | Internal steric clash | Non-covalently bound pairs above 0.8 of the lower distance bound |
| Intramolecular | Aromatic ring flatness | 5- and 6-membered aromatic rings with every atom within 0.25 Å of the closest shared plane |
| Intramolecular | Double bond flatness | Non-ring, non-aromatic C=C bonds: the two carbons and their four neighbours within 0.25 Å of the plane |
| Intramolecular | Internal energy | An ensemble of 50 conformers generated with ETKDGv3 and relaxed with UFF (up to 200 iterations); the ratio between the pose's UFF energy and the ensemble mean may not exceed 100 — less strict than the ratio of 7 used by Wills et al. |
| Intermolecular | Minimum protein–ligand distance | Heavy-atom distance above 0.75 times the sum of the van der Waals radii |
| Intermolecular | Minimum distance to organic cofactors | Above 0.75 times the sum of the van der Waals radii |
| Intermolecular | Minimum distance to inorganic cofactors | Above 0.75 times the sum of the covalent radii |
| Intermolecular | Volume overlap | Below 7.5% with the protein and with organic cofactors (van der Waals radii scaled by 0.8) and below 7.5% with inorganic cofactors (scale 0.5), through RDKit's ShapeTverskyIndex |

A pose that passes every test is **PB-valid**. Accuracy is measured separately: a minimum symmetry-aware heavy-atom RMSD (`GetBestRMS`) of at most 2 Å to the closest crystallographic ligand. Validity and accuracy are distinct criteria — a pose can be valid and wrong, or close to the reference and physically invalid.

### 6.2 Columns of `pose_validation.csv`

| Column | What it records here | Correspondence in the work |
| --- | --- | --- |
| `pose_id` | Pose identifier, from the producer or the file name | Identification; not a test |
| `receptor_id` | Target of the pose, by prepared receptor | Identification |
| `ligand_id` | Ligand reported by the source | Identification |
| `compound_id` | Evaluated compound | Identification |
| `engine` | Engine that produced the pose: `vina`, `dock6` or consensus | Identification |
| `score` | Raw score of the source engine, preserved untransformed | None; the work does not use the score to judge validity |
| `pose_valid` | Approval on every criterion actually evaluated | Equivalent to **PB-valid** restricted to the implemented subset |
| `valid_chemistry` | RDKit sanitization of the molecule as read | Sanitisation test |
| `valid_connectivity` | Molecule in a single fragment | Connectivity: the work's bond consistency requires the true ligand and is not reproduced |
| `valid_geometry` | Lengths and angles within the configured tolerances | Bond length and bond angle tests, with a different reference (see 6.4) |
| `valid_planarity` | Flatness of aromatic rings and double bonds | Aromatic ring and double bond flatness |
| `valid_energy` | Internal strain within the configured limit | Internal energy, in a different formulation (see 6.4) |
| `valid_clash` | Intermolecular contact above the configured fraction of the radii | Minimum distance to protein and to organic cofactors |
| `min_vdw_ratio` | Smallest distance/van der Waals radii-sum ratio; the decision metric | The same quantity as the minimum distance test, whose threshold in the work is 0.75 |
| `min_vdw_overlap` | Absolute overlap of the worst pair, in Å; negative means closer than the radii sum | Auxiliary measurement; not a test of the work |
| `min_heavy_distance` | Closest heavy-atom distance between pose and receptor, in Å | Raw quantity of the minimum distance test |
| `clashing_contacts` | How many pairs fell below the radii sum | Auxiliary measurement; the work decides by the worst pair and by volume overlap |
| `delta_energy` | ΔE between the pose and the relaxed conformer, in kcal/mol | Internal energy test, in a different formulation (see 6.4) |
| `energy_method` | Effective force field: `MMFF94` or `UFF` | The work uses UFF |
| `max_bond_deviation` | Largest relative bond-length deviation | Bond lengths test |
| `max_angle_deviation` | Largest absolute angle deviation, in degrees | Bond angles test |
| `max_planarity_deviation` | Largest distance to the fitted plane, in Å | Flatness tests, with the same 0.25 Å tolerance |
| `planar_groups` | How many planar groups were evaluated | Auxiliary coverage measurement |
| `atoms` | Atoms read in the pose | Context of the measurement |
| `bond_orders` | Origin of the bond orders: `explicit` or `proximity` | Applicability condition of the intramolecular tests |
| `resolved_by` | How the receptor was associated: `binding`, `identifier` or `unresolved` | Applicability condition of the intermolecular test |
| `rejection_reasons` | Keys of the failed criteria, separated by `;` | Set of failed tests |
| `rejection_details` | Readable description of the same reasons | Presentation |
| `error` | Reading or evaluation diagnosis, when present | Equivalent to a file-loading failure |
| `pose_file` | Pose file evaluated | Provenance |
| `receptor_file` | Receptor used in the intermolecular test | Provenance |
| `source_table` | Score table the pose came from | Provenance |

### 6.3 Columns of the displayed table

| Column | Content | Origin |
| --- | --- | --- |
| Pose | `pose_id` | Identification |
| Target | `receptor_id` | Identification |
| Compound | `compound_id`, with `ligand_id` as fallback | Identification |
| Engine | `engine` in upper case | Identification |
| Score | `score` of the source engine | Does not take part in validity |
| State | `pose_valid` as approved or rejected | PB-valid of the implemented subset |
| Failed criteria | Labels of the criteria with a false verdict | Failed tests |
| Min. contact · vdW fraction | `min_vdw_ratio` | Intermolecular minimum distance; threshold 0.75 in the work |
| Internal ΔE | `delta_energy` | Internal energy, in a different formulation |
| Actions | Per-pose verdict and receptor–pose complex in 3D | Presentation |

Each row's verdict repeats the six criteria with their outcomes, every measurement with its unit and the tolerances applied in that run.

### 6.4 Declared divergences from the work

| Point | The work | This implementation |
| --- | --- | --- |
| Reference for lengths and angles | RDKit distance-geometry bounds | A conformer of the same molecule relaxed with MMFF94 or UFF; the 25% relative tolerance coincides, the reference does not |
| Internal energy | Ratio between the pose's UFF energy and the mean of 50 ETKDGv3 conformers, threshold 100 | ΔE difference against a single relaxed conformer, default limit of 100 kcal/mol; an energy, not a dimensionless ratio |
| Internal steric clash | Non-bonded pairs above 0.8 of the lower bound | Not evaluated |
| Volume overlap | Below 7.5%, with 0.8 and 0.5 scaling | Not evaluated; the intermolecular decision uses distance only |
| Identity consistency | Formula, bonds, chirality and stereochemistry against the true ligand, through InChI | Not evaluated: a prospective screen has no true ligand. The compound SMILES, when available, only serves as a template to assign bond orders |
| Cofactors and waters | Separate tests for organic cofactors, inorganic cofactors and waters, with covalent radii for ions | A single set of receptor atoms with van der Waals radii, discarding solvent and the site's crystallographic ligand |
| Ring flatness | 5- and 6-membered aromatic rings | Every aromatic ring, with no size restriction |
| Accuracy | RMSD ≤ 2 Å to the crystallographic ligand | Out of scope for this block; redocking computes RMSD against the crystallographic reference |

Because of these differences the report identifies the criteria as `posebusters_inspired_v1` and claims no numerical equivalence with the work's PB-valid. A result from this implementation must not be published as a PB-valid rate.

## 7. Presenting the results

On the **Runs** tab the block presents the validation table as its result:

- Summary cards with evaluated, approved, rejected and delivered poses, plus the count per rejection reason and the criteria identifier.
- A paged table with pose, target, compound, engine, score, state, failed criteria, minimum contact (fraction of the van der Waals radii) and internal ΔE, with a filter across all, approved and rejected poses.
- Per row, the complete verdict — the six criteria, each measurement with its unit and the applied tolerances — and the receptor–pose complex in 3D, using the copies the block itself published.
- Below the table, collapsed, the block's files for download, preview or removal.

## 8. Interpretation limits

- Approval indicates physical plausibility of the delivered geometry. It does not indicate affinity, activity, selectivity or a correct binding mode.
- ΔE is conformational internal strain in the force field used. It is not a binding or free energy, and values from different molecules are not directly comparable.
- The steric clash depends on the associated receptor and the excluded residues. Without a receptor that criterion has no measurement.
- Poses read by proximity (PDBQT) do not receive the criteria that depend on bond orders. The report records that condition in `bond_orders`.
- Tolerances are configurable: one result set is comparable with another only when the recorded tolerances and criteria match.
- The block does not repair poses. Rejecting a pose neither repositions nor minimises it for the next block.

## 9. Verification

```bash
PYTHONPATH=src:tests python -m unittest test_pose_validation test_pose_validation_flow
```

`tests/test_pose_validation.py` covers reading per format, each criterion with approved and rejected cases, tolerances, the strict filter, published pose identity and the parameters the block declares. `tests/test_pose_validation_flow.py` covers the target and compound catalog, the form, canvas connections, pose resolution during execution, the identified refusal of an unconfigured block and the results table with its files.

## 10. References

- Buttenschoen, M., Morris, G. M., Deane, C. M. (2024). PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences. *Chemical Science*, 15, 3130–3139. DOI 10.1039/D3SC04185A. The tests, categories and thresholds quoted in section 6 come from this work.
- Bondi, A. (1964). van der Waals Volumes and Radii. *The Journal of Physical Chemistry*, 68(3), 441–451.
- Halgren, T. A. (1996). Merck molecular force field. I. Basis, form, scope, parameterization, and performance of MMFF94. *Journal of Computational Chemistry*, 17(5–6), 490–519.
- Rappé, A. K. et al. (1992). UFF, a full periodic table force field for molecular mechanics and molecular dynamics simulations. *Journal of the American Chemical Society*, 114(25), 10024–10035.
- Landrum, G. RDKit: Open-source cheminformatics. https://www.rdkit.org
