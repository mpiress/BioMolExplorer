# Redocking configuration

[Documentation](../README.md) · English · [Português](../redocking_configuration.md)

## Select and prepare pairs

After retrieving PDB structures, open the redocking stage configuration. In **Input Data**, explicitly select each PDB / ligand / residue pair. Initial retrieval does not select execution pairs without a user choice; at least one is required before applying configuration.

An explicit choice in **Results → Ligands → Save ligands** also updates pairs in redocking stages connected to that retrieval. This includes the first ligand selected for a PDB without a configured reference. If redocking inputs were restricted to other files, the reviewed structure is added to the input selection. Other PDB pairs and their settings are retained; a new ligand on the same chain inherits shared receptor preparation settings. Saving the same list again does not duplicate pairs.

Then choose **Run pipeline** and **Yes, reuse** to calculate the new pairs and incorporate them into accumulated results. Curation itself does not start calculations. To repair a ligand choice saved by an older version without reaching redocking, reopen **Ligands** and save the list again.

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

Chimera can exit with status zero even when a `.com` command fails. The pipeline checks those messages and expected files before continuing. Each output must exist, be nonempty and contain atoms; Open Babel conversions follow the same checks. Redocking records per-pair preparation failures and continues with the other pairs, retaining the Chimera script for diagnosis. Global installation or filesystem failures still stop the stage.

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

## Reuse and duration of subsequent batches

Valid redocking preparation is reused for candidate docking rather than repeated per compound. Vina retains maps per configuration within a batch when its Python API supports the template. DOCK6 shares surface, spheres and grids per receptor/site. Candidate preparation remains necessary without a compatible prepared input.

There is no implicit 300-second scientific timeout. Calculations continue until completion, error or cancellation; HTTP timeouts are separate. Individual candidate problems appear in **Excluded compounds**. See [Vina reuse](user_manual.md#reusing-receptors-in-vina), [DOCK6 validation](user_manual.md#dock6-validation-and-exclusions) and [backend duration](backend_usage.md#execution-duration).

## Partial preparation failures and source deletion

Redocking isolates preparation errors per receptor/ligand pair, reusing successfully prepared receptors. A receptor failure skips its dependent pairs; a ligand failure leaves the other ligands eligible. Partial pair outputs are cleaned. If all pairs fail preparation, the stage completes with a warning without starting Vina. Other preparation functions retain their strict behavior, and global installation/configuration/access/storage errors still fail the operation.

`pdb_codes.csv` keeps all selected pairs with `STATUS`, `FAILED_STEP` and `ERROR`. Failed pairs have `STATUS=preparation_failed` and no RMSD. `redocking_failures.csv` records the errors; failed pairs are excluded from valid prepared-pair choices. The interface reports the failure count and highlights error rows in red with **Preparation failed**, RMSD **—**, and **View error**. Pose/contact actions are disabled for those rows.

**Delete failure and source complex** requires confirmation and removes the failed row, the source PDB and execution snapshot, associated source metadata and verified input-cache copies. Source identity follows exact input/batch provenance; unrelated collections with identical filenames are preserved. A source changed since execution cannot be deleted through this stale result.

Deleting one PDB makes all its pairs unavailable for the next execution of that redocking block, while preserving completed historical results and their prepared files. The confirmation explains this shared-file effect. Retrieval sources add the ID to `excluded_pdb_ids`, visible under **PDB IDs excluded from future retrievals**, to prevent automatic download on rerun. Clear that exclusion for an intentional retry. Imported source assets and import references are removed instead.

Deletion requires owner/editor access and a finished execution. File/configuration changes and cache invalidation are audited and rolled back together if the operation fails. Historical restoration copies are not active redocking inputs. Preparation failures never trigger automatic source deletion.

## Incremental redocking and accumulated results

Choose **Yes, reuse** when starting the pipeline to retain compatible completed calculations. A new ligand for an existing PDB/chain reuses the prepared receptor and runs only ligand preparation and its redocking. A new receptor is fully prepared. Previously completed compatible pairs are not calculated again.

The new result table combines new pairs with compatible previous results from the same input source, including previously evaluated receptors and ligands that were not selected again. The original PDB must remain available: deleted inputs are never automatically restored from history. Each pair appears once. Reused poses and prepared structures are copied into the new run so downloads and 3D views remain available.

Compatibility uses PDB contents, chain, cofactors, preparation options, pH, charges, templates, Vina parameters and scientific implementation. Changing the structure or receptor preparation invalidates receptor reuse; changing docking parameters requires a new pose calculation. SHA-256 hashes verify reused files. Older native results can also be adopted when their configuration, prepared files, pose, RMSD and integrity manifest are available.

Only successfully completed pipeline runs supply reusable data. Cancelled, interrupted, failed or awaiting-input runs never supply incremental results, even if an individual stage finished. Partial files may remain in history for diagnostics but are discarded as sources for subsequent calculations. Preparation failures and removed outputs are also excluded. **No, run again** bypasses reuse.

## Redocking engine and Chimera preparation

Choose **Redocking engine**: AutoDock Vina (default) or DOCK6. Search controls
follow the selected engine. Older projects default to Vina when
`redocking_engine` is absent. Both engines use Chimera preparation and the
administrator supplies the DOCK6 installation.

DOCK6 redocks the pair's reference ligand with a shared site/grid and flexible
or rigid search. Poses follow the existing redocking result contract. Its RMSD
uses heavy atoms and molecular symmetry in the receptor coordinate frame;
it does not align the pose onto the reference.

Structural cleanup precedes hydrogens, charges and minimization, including
ligand alternate locations. Local chemical definitions and repaired export
bond orders participate in traceability and caching. See
[molecular preparation](molecular_preparation.md). Complete Chimera output is
retained in logs. Obsolete failure reports are removed after a successful retry.
