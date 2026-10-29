# Redocking configuration

[Documentation](../README.md) · English · [Português](../redocking_configuration.md)

## Select and prepare pairs

After retrieving PDB structures, open the redocking stage configuration. In **Input Data**, explicitly select each PDB / ligand / residue pair. Initial retrieval does not select execution pairs without a user choice; at least one is required before applying configuration.

An explicit choice in **Results → Ligands → Save ligands** also updates pairs in redocking stages connected to that retrieval. This includes the first ligand selected for a PDB without a configured reference. If redocking inputs were restricted to other files, the reviewed structure is added to the input selection. Other PDB pairs and their settings are retained; a new ligand on the same chain inherits shared receptor preparation settings. Saving the same list again does not duplicate pairs.

Then choose **Run pipeline** and **Yes, reuse** to calculate the new pairs and incorporate them into accumulated results. Curation itself does not start calculations. To repair a ligand choice saved by an older version without reaching redocking, reopen **Ligands** and save the list again.

The **Selected pairs** table contains PDB, Ligand, Residue number, Chain and Actions. **Edit pair** opens a separate form; Cancel discards changes. **Remove pair** removes the row from the selection without deleting its source PDB. Resolution follows the metadata and is not editable. Use the manual pair option for PDBs without metadata.

Expand **Receptor cofactors (optional)** and enable the cofactor option and enter residue codes separated by commas, such as `FAD, MG`. Without cofactors, complex selection is `select #0:.{chain}`. Declared cofactors are preserved as part of the receptor, including those outside the selected chain. The redocking ligand cannot also be a cofactor.

Configure solvent removal, removal of existing hydrogens, hydrogen addition, minimization and charge method for receptor and ligand. **Shared preparation for the collection** contains one panel for all receptors and another for all ligands. Ligand preparation and conformation share these options. PDBQT conversion pH appears under Vina settings. Panels align side by side on wide screens and stack on narrow screens; the cofactor input is hidden until enabled.

Pairs using the same PDB and chain share their prepared receptor. The interface applies common receptor and ligand settings to the entire set. Legacy per-pair projects remain executable through the backend; the interface warns before replacing heterogeneous settings with shared preparation when saving. Individual processing passes only pairs and settings for the corresponding structure.

Configured pairs are reused without another selection prompt. Before calculations, the pipeline validates receptor, ligand, residue, chain and cofactors against the structure. Prepared complexes require valid binding centers and engine-specific formats: PDBQT for Vina; receptor and ligand MOL2 plus receptor `.noH.pdb` for DOCK6. Preparation uses classic Chimera and Open Babel; prepared inputs skip these steps. The selected engine must be installed. DOCK6 also requires sphgen, sphere_selector, showbox, grid and parameter files from its 6.11 installation. Verbosity is always zero.

## Diagnosing a failed run

Review the reported executable and process message in the execution log. Failed Chimera command scripts are retained. A missing or invalid ligand center and incomplete prepared files prevent docking. Correct the configuration or installation and retry; see [pipeline validation](pipeline_validation.md) for test coverage and native execution limits.

## Choose the engine

Under **Redocking method**, select **AutoDock Vina** or **DOCK6 6.11**. Only active engine fields are shown and validated; switching retains the other engine's values for later editing.

- **Vina:** box X/Y/Z in Å, conversion pH, search effort and pose count. The reference ligand supplies the box center. **Vina · engine options** contains seed and scoring function.
- **DOCK6:** administrator-managed installation, flexible/rigid search, SES density, probe radius, sphere selection distance, maximum orientations, saved poses and seed. Expandable sections contain engine, grid, box and sphere options. The workflow uses the prepared ligand directly, without additionally minimizing the crystallographic reference.

```text
Select and review pairs in the table.
Set one preparation for receptors and one for ligands.
Choose Vina or DOCK6 6.11 and adjust its fields.
Save, run the pipeline and review RMSD and diagnostics.
```

New runs with MOL2 available calculate symmetry-corrected heavy-atom RMSD without geometric alignment. The reference restores ligand coordinates from the original PDB; if no PDB exists, results explicitly record `RMSD_REFERENCE=prepared`. Legacy Vina inputs containing only PDBQT retain the legacy calculation, identified by `RMSD_METHOD`. These methods are not interchangeable. CSV export includes `Engine`, `RMSD_method` and `RMSD_reference`. A finite RMSD is not automatic protocol acceptance.

Preparation uses the component CCD graph before Chimera adds hydrogens and charges. IDATM types are fixed from that graph to prevent incorrect geometric reinterpretation. MOL2 is written without reopening PDB. After PDBQT conversion, names are restored by matching coordinates before docking. The `.export.json` or `.input.json` records the definition and missing atoms; observed fragments are not completed with invented heavy atoms. Recoverable Vina/RMSD pair failures are recorded without aborting the remaining pairs.

## Backend configuration

Pass `pdb_codes=[["4M0E", "1YL", 604, "A", 2.0]]` for the selected records. `preparation_options` specifies common `receptor`, `ligand` and `cofactors` settings and takes precedence. For compatibility, `preparation_pairs` is a dictionary keyed by `4M0E|1YL|604|A`, with `cofactors`, `receptor` and `ligand` settings. Allowed preparation options are `remove_solvent`, `remove_hydrogens`, `add_hydrogens`, `minimize` and `charge_type` (`gas` or `am1`). Omitted options use the defaults. The fifth record element is resolution metadata. Existing projects with a separate `ligand_chain` use the shared `CHAIN` value.

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

Connect redocking to **PDB receptor (retrieval or redocking)** in **Prepare for docking** and select the `.dockprep.pdbqt` receptor (or `.dockprep.mol2` for a DOCK6-only preparation) under **Receptor to use**. Ligand-only files are hidden. Receptor preparation controls are disabled while companion formats and binding centers are reused. Candidates from one or more ChEMBL, PubChem, ZINC or uploaded sources still undergo preparation. Choose Vina, DOCK6 or both, then connect receptor and compounds to each selected engine. See the [user manual](user_manual.md#preparation-and-redocking) for individual and merge processing.

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


## Target and best RMSDs per target

The redocking table includes **Target**, obtained from the retrieval block that
supplied each PDB through saved file provenance. Merged inputs match the PDB to
its specific source. Imported metadata may also provide `Target` or `TARGET`.
When the source cannot be identified, the cell shows **—**; the redocking block
name is not used to invent an association.

Use **Therapeutic target** to filter and search for target, PDB or ligand. The
initial ordering groups targets and prioritizes smaller RMSDs. Click **Target**
or **RMSD (Å)** to sort. **Highlight best RMSDs per target** marks each minimum and
its ties in green (numerical tolerance 1e-9 Å). Failures remain red and are excluded
from comparison. Minima use all results for the target, including other pages and
rows outside the text search.

**Download RMSD table** exports filtered pairs with `Target` and `Best_for_target`,
revalidating access. Labels do not modify poses, RMSDs, original scientific
metadata or caches. Earlier results can display the column when opened if their
provenance is available.

```text
Open results and select ACHE in the target filter.
Enable Highlight best RMSDs per target.
Compare green rows and open the desired simulation.
Download the RMSD table to retain the selection and its origin.
```

## Native preparation integrity

Ligands receive the observed component's CCD topology before Chimera adds hydrogens.
After protonation, the molecular graph is validated and its formal charge is supplied
explicitly to native AM1-BCC or Gasteiger charging. Minimization uses this already
prepared graph without repeating DockPrep or replacing the chosen charges. The
`.lig.preparation.json` report records charge, method and validation.

For Chimera's known `accPhiPsi` FindHBond error, only hydrogen placement is retried
using native geometric AddH; charging and minimization remain enabled. Free
crystallization amino acids recorded as HETATM are excluded from the protein unless
explicitly selected as cofactors. Peptide-bonded amino acids remain in the receptor.

RMSD compares heavy atoms in the crystal reference and best pose without alignment
that would hide displacement. A high RMSD is an evaluation result, not an execution
failure. Topology, tool and configuration errors remain visible with their failing
step and diagnostics, while other pairs continue independently.
