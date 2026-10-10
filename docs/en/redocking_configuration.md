# Redocking configuration

[Documentation](../README.md) · English · [Português](../redocking_configuration.md)

## Select and prepare pairs

After retrieving PDB structures, open the redocking stage configuration. In **Input Data**, explicitly select each PDB / ligand / residue pair. No pair is selected automatically; at least one is required before applying configuration.

Each selection has a single **Chain**, used for both receptor and ligand, and its own preparation options. Resolution follows the selected metadata and has no editable control. For user-supplied PDBs without metadata, use the manual pair option.

Enable the cofactor option and enter residue codes separated by commas, such as `FAD, MG`. Without cofactors, complex selection is `select #0:.{chain}`. Declared cofactors are preserved as part of the receptor, including those outside the selected chain. The redocking ligand cannot also be a cofactor.

Configure solvent removal, removal of existing hydrogens, hydrogen addition, minimization and charge method for receptor and ligand. Ligand preparation and ligand conformation share the same settings. pH is configured once for the stage.

Pairs using the same PDB and chain share their prepared receptor file and must use identical receptor options and cofactors. Their ligand settings may differ. Individual processing passes only pairs and settings for the corresponding structure.

Configured pairs are reused without another selection prompt. Before calculations, the pipeline validates receptor, ligand, residue, chain and cofactors against the structure. Prepared complexes require valid PDBQT files and ligand centers. Chimera, Open Babel and Vina must be available in the scientific environment; prepared inputs require only Vina. Verbosity is always zero.

## Diagnosing a failed run

Review the reported executable and process message in the execution log. Failed Chimera command scripts are retained. A missing or invalid ligand center and incomplete prepared files prevent docking. Correct the configuration or installation and retry; see [pipeline validation](pipeline_validation.md) for test coverage and native execution limits.

## Backend configuration

Pass `pdb_codes=[["4M0E", "1YL", 604, "A", 2.0]]` for the selected records. `preparation_pairs` is a dictionary keyed by `4M0E|1YL|604|A`, with `cofactors`, `receptor` and `ligand` settings. Allowed preparation options are `remove_solvent`, `remove_hydrogens`, `add_hydrogens`, `minimize` and `charge_type` (`gas` or `am1`). Omitted options use the defaults. The fifth record element is resolution metadata. Existing projects with a separate `ligand_chain` use the shared `CHAIN` value.

Selecting a subset does not delete other PDB structures or prepared files from the input collection. When a pair does not specify a charge method, ligand preparation retains the stage-wide charge method.

## Preparation failures with classic Chimera

Chimera can exit with status zero even when a `.com` command fails. The pipeline checks those messages and expected files before continuing. Each output must exist, be nonempty and contain atoms; Open Babel conversions follow the same checks. Invalid outputs stop the stage and retain the Chimera script for diagnosis.

Paths in `open` and `write` commands are passed without shell quotes: classic Chimera interprets those quotes as part of the filename. Paths with spaces were verified with the real executable. Tool stdout and stderr are recorded in the scientific logs.

## Inspect redocking results

Expanding a completed stage in **Runs** shows a table of PDB, ligand, residue, chain and **RMSD (Å)**. **View simulation** opens a dialog listing files for that selection. The list distinguishes the prepared receptor, reference ligand, poses and metadata. Files shared by the receptor or collection accompany their corresponding simulations.

Use each row's download button for an individual file, or **Download all (ZIP)** for the simulation bundle. The ZIP preserves directories, including when the reference ligand and Vina output share a filename. Files from other simulations are excluded.

**View 3D structure** opens the browser viewer, following the same flow as PDB inspection. It accepts PDB, PDBQT and MOL2 for inspecting receptors, ligands and poses before downloading. Use **Model** when opening a file containing multiple poses. The scene opened by the table’s **3D** button automatically selects the lowest-scoring pose, or the first unscored pose; that scene does not offer model switching. Labels and actions follow the selected language. Readers can also inspect and download results.

The table summarizes finite, nonnegative values from the `RMSD` column of `pdb_codes.csv`, displayed to three decimal places; the file retains its original precision. RMSD compares the reference ligand with the redocking result. The interface does not impose an automatic scientific acceptance threshold. PDB / ligand / residue / chain identifies each simulation. Failed or running stages retain the general file list for diagnosis.

[Backend usage](backend_usage.md) includes a complete parameter dictionary and the output file layout.

## Reuse the receptor for candidate docking

Connect redocking to **PDB receptor (retrieval or redocking)** in **Prepare for docking** and select the `.dockprep.pdbqt` receptor under **Receptor to use**. Ligand-only files are hidden. Receptor preparation controls are disabled while companion formats and binding centers are reused. Candidates from one or more ChEMBL, PubChem, ZINC or uploaded sources still undergo preparation. Choose Vina, DOCK6 or both, then connect receptor and compounds to each selected engine. See the [user manual](user_manual.md#preparation-and-redocking) for individual and merge processing.

## Styles and 3D interactions

The **3D** scene compares pose, receptor and reference in their original
coordinates. Choose independent ligand/receptor styles and use **Ligand
hydrogens** to show only H present in the file. Classified interactions have
colored dashed traces, filters and hover identification. PDBQT requires an
associated SMILES to recover topology; redocking results lacking it retain
distance contacts and explain chemical classification unavailability. See the
[viewer guide](molecular_viewer.md) for colors, requirements and classification
limits.
