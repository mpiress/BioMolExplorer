# Structural support for native DOCK6 6.11 Amber Score

This is target-independent preprocessing implemented in `structure_sanitation.py`, `md_residue_support.py`, `md_cofactor_support.py`, `md_torsion_support.py` and integrated by `md_engine.py`. The 2WSL target is a validation case, not a branch in production code.

## Execution contract

Immutable inputs → execution-local sanitation → explicit residue variants → coherent protein profile → retained-component parameters → native topology preparation → topology/torsion validation → native or local-compatible executable → implicit-solvent refinement/analysis → verified per-pair cache.

Target PDB and a target-matched coordinate pose are required. Bare SMILES are insufficient. The protocol is native Amber Score, not conventional trajectory-producing molecular dynamics; no RMSF, equilibrium trajectory or universal pKa/tautomer prediction is supplied.

## Sanitation and biological roles

`HETATM` is not a deletion rule. Peptide geometry, LINK/CONECT records, recognized cofactors and explicit retention choices distinguish covalent residues from removable additives. Recognized unbound additives/water and isolated amino-acid HETATMs may be removed; covalently connected PTMs, recognized cofactors/metals and ambiguous components are retained for parameter diagnosis. In Redocking, this also runs immediately before Chimera because the complex PDB may not exist when scripts are generated. The selected reference ligand is excluded from the receptor copy; covalent reference ligands require an appropriate bonded protocol.

Alternate conformations are chosen consistently per residue by occupancy with deterministic ties. Shared blank-altLoc atoms remain. Nonfinite coordinates, unresolved duplicates, duplicate serials and unselected multiple models fail before engines. Standard-protein aliases OT1/O, OT2/OXT, HN/H and ILE CD/CD1 are recorded. Existing explicit HIS hydrogens, disulfide annotations and backbone gaps control variant/TER normalization without moving observed atoms or inventing loops.

Classic Chimera DockPrep does not reconstruct missing side chains in the current helper (`incompleteSideChains=False`) and mutation of modified residues is disabled. AMBER templates can complete supported missing side-chain atoms when the required backbone exists; this is a separate stage.

## Protein profiles and unknown chemistry

`residue_parameter_policy=auto_curated` retains native ff94 for standard targets. Supported internal PTMs switch the entire protein to ff03 plus published Forcefield_PTM parameters, avoiding a mixture of incompatible backbone-charge profiles. `native_ff94` retains the native protein profile; `ff03_ffptm` forces the uniform ff03 profile.

There are 32 bundled models, including CSS, SEP, TPO and PTR. Each model defines its own atom graph/types and chemical charge state; CSS is neutral while SEP/TPO/PTR have charge −2. These are fixed published states, not pH predictions. Hashes, graph, charge, origin and license are verified. Two duplicate hydrogen-name errors in source data (DA1 and SEN) were renamed without changing charges, types or bonded parameters; the exact adaptations are in the manifest.

Post-preparation checks require complete model atoms, expected bonds/types/charges and preservation of observed heavy coordinates within 0.002 Å. Incompatible termini and special crosslinks are rejected. Additional validated ff03 profiles can be provided through `BIOMOL_AMBER_PTM_PROFILES`, with `manifest.json`, PREP/frcmod files, graphs, charges, hashes, source URL and license; built-in models cannot be overridden.

Isolated organic components without native templates can use CCD-backed RDKit chemistry plus Chimera-bundled `antechamber`/`parmchk2`, GAFF and AM1-BCC. Graph/element inventories, stereochemistry, explicit charge/protonation evidence, finite charges and unresolved frcmod terms are checked. Parameters and source metadata are cached. A detached small-molecule GAFF model is never substituted for an unknown covalently embedded PTM. MSE without a validated bonded profile, unsupported metals/spin states and special crosslinks still require appropriate parameters. Selenium is not silently replaced by sulfur.

## Local torsion compatibility

Topologies are inspected before execution. When the original kernel cannot represent their periodicities/phases, sources from the licensed DOCK6 6.11 installation are copied and compiled locally with the full expression `E = k * (1 + cos(n*phi - phase))` and forward-mode coordinate derivatives. Terms are not truncated or fitted. The original installation remains unchanged. Native pi rounding (3.141594, tolerance 1e-5 rad) retains the native fast path for standard topologies.

This requires licensed DOCK6 sources, C/C++ compilers, make, flex and zlib. Build results have content-addressed cache identities. `requirements.yml` declares the compiler dependencies; the Python package does not install DOCK6 or Chimera. Unsupported source layouts fail explicitly.

## Cache, results and validation limits

Completed-pair reuse and independent receptor/ligand/complex preparation checkpoints are separate. Empty diagnostic logs from reused preparations are excluded from pair-cache objects; empty scientific topologies remain errors. Execution/protocol/support JSON files record decisions, model sources and executable hashes.

A dated 2WSL + CHEMBL92663 run completed automatically and its repeat restored cache with `run_pair` blocked in the test. The 32 libraries were checked through native topology generation; numerical tests compared 100 energies and 1,200 derivatives. These are focused checks, not a guarantee for every chemistry or complete desktop acceptance. A positive score does not demonstrate favorable binding. See the [detailed report and measured evidence](../amber_score_structural_support.md).

## License and reproduction

Code retains MIT licensing; bundled [Forcefield_PTM](https://doi.org/10.1021/ct400556v) data have **CC BY-NC 4.0** terms documented in LICENSE and the parameter NOTICE/manifest. CCD data and explicit local chemical definitions retain their respective provenance/license metadata.

```bash
PYTHONPATH=src python scripts/validate_amber_support.py --dock6-root /path/to/dock6 --output /tmp/amber-library-validation
```

This invokes native tools and writes validation products under `--output`. Optional `--pairs-json` accepts an exported nonempty `md_inputs` array and validates full execution plus completed-pair cache reuse. The example installation path must be replaced with the actual licensed installation.
