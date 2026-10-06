# UCSF Chimera replacement assessment

[Documentation](../README.md) · English · [Português](../chimera_migration.md)

## Decision

The October 5, 2026 review did not establish a complete replacement preserving scientific protocols, accepted inputs and editable templates. Chimera remains required. No scientific backend or functionality was changed. A Python implementation is possible in principle, but it needs implementation and comparative validation before replacing the existing workflow. See the [detailed assessment in Portuguese](../chimera_migration.md) for the full inventory and reasoning.

Code, workflows, resources, installation, UI, tests and documentation were inspected, including searches for notebooks. Official documentation was consulted. `chimera` was unavailable on this session's PATH; OpenMM, PDBFixer, OpenFF and ParmEd were absent from the inspected scientific environment. The project file search found no versioned PDB/MOL2/PDBQT reference files. No experimental comparison between Chimera and proposed replacements was performed.

In redocking, effective selections and options are generated for each validated pair. The base template still removes solvent and hydrogens when used directly by standalone preparation; redocking moves these operations to receptor and ligand stages to respect their respective settings.

## Inventory

Actual resources live in `src/biomolexplorer/resources/chimera/`; `paths.resolve_path` redirects legacy `src/scripts/chimera/` paths to packaged or worker-specific resources.

| Resource or entry point | Function | Outputs or constraint |
| --- | --- | --- |
| `prepare_complex.template` | Keep the selected chain and explicitly declared cofactors, delete the inverse selection | `{PDB}_{CHAIN}.complex.pdb`; no cofactor is included by default; receptor/ligand stages handle solvent and hydrogens |
| `prepare_receptor.template` | Remove ligands, keep protein, add H, assign ff14SB charges with `method gas`, minimize, export MOL2, reopen and remove H | `.dockprep.mol2` and `.noH.pdb` |
| `prepare_ligand.template` | Isolate residue/chain, add H, assign selected charges, minimize, export PDB and reopen for MOL2 export | `.lig.pdb` and `.lig.mol2`; attribute preservation through reopening needs comparison |
| `prepare_better_conform.template` | Prepare, charge and minimize the extracted first Vina pose | `.lig.mol2` used by DOCK6 |
| `prepare_md.template` | Remove solvent/H, add H, write PDB | Published in the editor; no execution caller identified; does not run MD |
| `prepare_on_chimera` | Run `chimera --nogui --silent`, propagate command failures and remove scripts only after success | Implemented in both `caad/docking.py` and `caad/redocking.py` |
| `perform_consensus` | Instruct manual Chimera loop refinement when `pdb_code` is absent | Human step followed by workflow interruption |

The active redocking wrapper imports classes from `caad.docking`. The duplicate `caad.redocking` module also needs coverage for its consumers. Preparation executes complex, receptor and ligand groups in sequence, parallelizing each group; redocking uses the same selected chain for receptor and ligand. Ligand preparation and conformation share options, and explicitly selected cofactors are preserved.

The catalog publishes all five templates for preparation, redocking, Vina and DOCK6. Template validation accepts the scientific commands `open`, `delete`, `select`, `write`, `close`, `addh`, `addcharge`, and `minimize`, without restricting every argument to its original setting. The UI can disable H addition, minimization and removals, and change charge methods. Replacing only default templates would leave existing customizations unsupported.

Open Babel subsequently generates PDBQT and applies the requested pH. Chimera templates do not receive that pH value. PyMOL calculates centers of mass and performs structural analyses. Vina runs docking/redocking, the native DMS port computes molecular surfaces, and DOCK6/accessory programs perform spheres, boxes, grids, minimization, docking and footprint. Removing Chimera would leave these dependencies; Python libraries also remain third-party software.

## Candidates and limits

| Candidate | Intended role | Compatibility concern |
| --- | --- | --- |
| Biopython, already declared | PDB/mmCIF reading and structural selection | Must define models, altlocs, insertion codes and residue/solvent classification; no charges/minimization |
| RDKit, already declared | Chemical identity, stereochemistry, Gasteiger, conformers | MMFF/UFF do not reproduce the current Amber protocol |
| Open Babel, already declared | Python conversion/perception/charge APIs | Exporting MOL2 does not establish chemical equivalence |
| PDBFixer | Missing atom/residue repair | Repair changes structures and should be an explicit validated option |
| OpenMM | Amber protein parameterization, H placement, minimization | Requires suitable templates; different minimizer and protonation decisions |
| OpenFF Toolkit + openmmforcefields | Ligand parameterization and OpenMM integration | Requires known bond orders/formal charges; alternative force fields change the protocol |
| AmberTools | AM1-BCC/GAFF through Antechamber | Retains compiled executables even if orchestrated from Python |
| ParmEd | Parameter transport and MOL2 I/O | Output atom types still need DOCK6 compatibility validation |
| Meeko | AutoDock preparation and chemical pose recovery | Does not independently replace Amber receptor preparation, minimization and DOCK6 MOL2 |

These are future candidates, not activated dependencies. No untested version set was selected. [OpenMM model editing](https://docs.openmm.org/latest/userguide/application/03_model_building_editing.html), [PDBFixer](https://github.com/openmm/pdbfixer), [RDKit force-field helpers](https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html), [Open Babel charges](https://open-babel.readthedocs.io/en/latest/Charges/charges.html), [openmmforcefields](https://github.com/openmm/openmmforcefields), [ParmEd MOL2](https://parmed.github.io/ParmEd/html/api/parmed/parmed.formats.mol2.html) and [Meeko](https://github.com/forlilab/Meeko) document the respective capabilities. My assessment is that they could form a future backend, subject to the following constraints.

## Unresolved requirements

1. **Chemical identity:** `MolConverter.extract_pdb_to_pdbqt` extracts the first pose and truncates coordinate records at column 66. The recovery path does not pass a chemically complete reference molecule. Requiring SDF/SMILES/CCD or additional metadata would restrict current inputs unless a compatible alternative were validated. [OpenFF's FAQ](https://docs.openforcefield.org/en/latest/faq.html) explains why PDB coordinates alone cannot reliably define the necessary chemistry.
2. **Both charge methods:** the public API exposes `gas` and `am1`; providing Gasteiger alone is incomplete. `am1` selects the Chimera AM1-BCC route. Protein standard residues use ff14SB charges; `method gas` does not mean assigning Gasteiger to the entire protein. Compare perception, protonation and per-atom charges across implementations. See [Chimera minimization](https://www.rbvi.ucsf.edu/chimera/docs/UsersGuide/midas/minimize.html) and [Chimera team's Gasteiger guidance](https://mail.cgl.ucsf.edu/mailman/archives/list/chimera-users%40cgl.ucsf.edu/thread/MRTULUWKFOL6TJWDAQVU2WFPWI3UQSH2/). The specific addcharge page was blocked by the source website in this session.
3. **Minimization:** Chimera uses MMTK with Amber/Antechamber parameters, structural preparation and steepest descent followed by conjugate gradient. [OpenMM's minimizer](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.LocalEnergyMinimizer.html) uses L-BFGS. Differences require geometric/downstream validation; they do not alone prove damage.
4. **Downstream contracts:** MOL2 charges/types affect DOCK6 grids, scoring and footprint; noH PDB feeds the native surface generator, and minimized ligand centers define Vina boxes. Matching names/extensions alone is insufficient.
5. **Customization and loops:** preserve editable selection/minimization arguments and saved configurations, and define how the manual loop step is replaced. No equivalent Chimera scientific-command interpreter was implemented.
6. **Reference evidence:** mocked-engine tests cannot establish equivalence for both charge methods and all supported cases.

## Environment and validation

`requirements.yml` and `requirements.yml` were reviewed and already have identical dependency lists, including RDKit, Biopython, Open Babel, PyMOL and Vina. Their dependencies were preserved and review comments added. No unused candidate libraries were added, no general package upgrade was performed, and no Conda solve/install was attempted.

Existing `test_docking_handoffs`, `test_workspace` and `test_stage_default_isolation` check artifact handoffs, template constraints, isolation and workspace behavior. All three modules passed during this review: 31 tests in 18.151 seconds, using the BioMolExplorer environment Python and `PYTHONPATH=src:tests`. Relevant scientific engines are mocked; passing these tests is software-contract evidence, not chemical equivalence. HTML pages were regenerated and `git diff --check` reported no whitespace errors.

A future migration needs paired references covering real project inputs, chains/FAD, termini/histidines/disulfides, modified residues, aromatic/charged ligands and imported poses; chemical reference/atom mapping; both charge methods and editable options; declared scientific tolerances for atom/bond identity, charges, MOL2 types, geometry, centers, RMSD, docking scores, footprint and ranking; real-engine comparisons and clean Python 3.12 installation. Only then should Chimera installation, code, templates and instructions be removed.

Potential additions include preparation provenance, optional structural repair, explicit protonation states and chemically faithful pose recovery. These are proposals and were not implemented in this review.
