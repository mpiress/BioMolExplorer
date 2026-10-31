# User manual: from first sign-in to results

[Documentation](../README.md) · English · [Português](../user_manual.md)

This manual follows a study from choosing a language to exporting the project. Button names below refer to the English interface. For technical parameters and APIs, see the [Flet workspace](frontend.md), [backend usage](backend_usage.md), [projects and versions](projects.md) and [pipeline validation](pipeline_validation.md) guides.

[Quick examples for each feature](#quick-examples-for-each-feature)

## 1. Prepare the environment and start the application

Download the [project from GitHub](https://github.com/mpiress/BioMolExplorer) and enter its source code folder:

```bash
git clone https://github.com/mpiress/BioMolExplorer.git
cd BioMolExplorer
```

Without Git, use **Code → Download ZIP**, extract the archive and open a terminal inside the extracted folder. The [installation and configuration guide](installation.md) explains downloading, external tools and paths.

**Configure UCSF Chimera, Open Babel and Vina in the main scientific environment.** DOCK6 is required for its stages. Amber Score uses only the native DOCK6 6.11 bundle; no external AmberTools environment is required.

1. Create the scientific environment from the downloaded source code root.
2. Activate the environment and install the interface.
3. Configure external tool directories on `PATH` and replace `/path/to/dock6-6.11` with the actual DOCK6 root containing `bin/` and `parameters/`.
4. Start browser or desktop mode after checking the executables.

```bash
conda env create -f requirements.yml
conda activate BioMolExplorer
python -m pip install -e '.[ui]'
biomolexplorer-ui --web --language en --dock6-path /path/to/dock6-6.11
```

In your browser, open `http://127.0.0.1:8550`. For desktop mode, omit `--web`. Without `--language`, the application starts in English. Use `--language pt` to start in Portuguese; this option also works in desktop mode.

On a server, use `--web --no-browser --port 8550`. `--data-dir /path/workspace` selects storage for accounts and project indexes; use the same path on subsequent launches. `--worker-python /path/python` selects another Python interpreter for calculations. `--dock6-path /path/dock6` specifies the DOCK6 installation.

Fingerprints, similarity, graphs and ADMET use the local scientific environment. Retrieval queries external providers. Preparation and docking need the corresponding engines: Chimera, Open Babel, Vina and, for DOCK6, that protocol's executables and resources. Installing the interface does not install these engines. The backend guide explains the dependencies; validate a docking protocol against a reference case before using it in your study.

## 2. Choose a language, create an account and sign in

1. In the top bar, click the flag in the right corner to open **Select language**.
2. Choose the United States flag for English or the Brazilian flag for Portuguese. The screen updates immediately; fields already filled in are preserved.
3. For your first visit, click **I do not have an account yet**.
4. Enter your name, email, password and password confirmation. The password must contain at least ten characters; both password fields must match.
5. Click **Create account**. Subsequently, use the same email and password under **Sign in to workspace**.

The language applies to that interface session: menus, forms, selections, progress and recognized messages use the chosen language. User-provided names, files, molecular identifiers, SMILES, scientific parameters and JSON keys are preserved. External engine logs retain the language produced by the engine. To change language during use, open the same flag menu in the top bar without signing out. The next launch starts in the language specified by the startup command.

Accounts belong to the workspace selected at startup. This version has no email password recovery. **Sign out** ends the session; an expired session returns to authentication.

## 3. Create a project and choose its folder

1. In the workspace, click **New project**.
2. Enter a name identifying the study; add a description and choose the card color from visual swatches. The Tags field has been removed.
3. Under **Project folder**, click the folder button inside the field and select the parent folder. Enter the project name first: the project will live in a subfolder with that name. An existing destination requires replacement confirmation; removal happens only when saving. The field cannot be typed into. Desktop and browser mode open visible folder navigation on the computer running BioMolExplorer. Create folder asks for a name, creates it inside the current location and selects it automatically.
4. Confirm creation and open the project.

Linux example: `/home/researcher/studies/enzyme-a`. In browser mode, this path belongs to the backend computer. A path on your laptop works only if the backend runs there or can access the folder. The application needs write permission. Folders belonging to other projects, overlapping folders and the accounts directory cannot be reused for a new project.

The selected folder stores `project.json`, files in `assets/`, runs in `runs/`, versions in `.history/` and project working copies. The central database retains accounts, sessions, permissions and indexes; general infrastructure logs use the configured logging directory. Results and run-associated logs are available through the project. Restarting with the same workspace preserves the folder association.

Search by name or description to find a project. Archiving removes it from the active list; restoring allows you to use it again. Deletion requires confirmation and permanently removes the project folder, including inputs, results and history. Export a backup first; history cannot undo project deletion.

## 4. Understand the project areas

| Area | When to use it |
| --- | --- |
| Pipeline | Add, connect, configure and execute blocks |
| Files | Upload your own data, inspect inputs and download files |
| Runs | Follow statuses, open results and inspect errors and logs |
| Share | Invite collaborators and review permissions |
| History | Inspect changes and, as owner, restore versions |

The pipeline defines dependencies. Block positions help visual organization; they do not determine calculation order. Independent branches can run even when another branch fails. Descendants of a failed stage cannot proceed until their input becomes available.

## 5. Upload and check your own files

1. Add **Import my files** to the pipeline and open its configuration.
2. Choose **Type of files to add** and click **Select files from disk**. Select all required files; each file is validated before entering the project.
3. Reuse existing uploads through **File already uploaded to the project** or **Add all project files**.
4. Review the **Type / File / Actions** table. You can change each row's type (with fresh validation) or use **Remove** in the **Actions** column. Removing a row preserves the file in the project library.
5. Include metadata and companion formats for prepared receptors. Set the target folder for structures and apply the configuration. One block can publish several types, connected to the corresponding consumer inputs.

A minimal compound CSV is:

```csv
name,smiles
ETHANOL,CCO
BENZENE,c1ccccc1
PYRIDINE,c1ccncc1
```

`molecule_chembl_id,canonical_smiles` is also accepted. Missing codes are generated from the structure. Use stable codes suitable for filenames, without slashes or paths. A code must not represent different molecules. Use commas to separate fields and a decimal point for numbers; save as UTF-8. Values containing commas, such as fingerprint lists, must be quoted in the CSV.

| Type | Required content |
| --- | --- |
| Compounds | CSV with `canonical_smiles` or `smiles`; an identifier is recommended |
| Fingerprints | CSV with `molecule_chembl_id,fingerprint`; equally sized lists of 0/1 bits |
| Similarities | CSV with `source,target,value`; identifiers and finite values between 0 and 1 |
| PDB complexes | `.pdb` with ATOM/HETATM and valid coordinates; ligand/chain records in the form or `pdb_codes.csv` |
| Prepared receptors | Complete receptor, ligand and metadata files compatible with the protocol |
| Vina | `.pdbqt` with coordinates and `REMARK VINA RESULT` |
| DOCK6 | `*_scored.mol2` with MOLECULE/ATOM sections and a numeric `Grid_Score` |
| Scores | CSV with an identifier (`molecule`, `molecule_chembl_id` or `id`) and numeric scores |
| ZINC download list | TXT/URI file or exported tranche-browser script with SMI/MOL2 links, including compressed files |
| Other | General data files; choose a specific type to validate scientific formats |

For prepared receptors, preserve `<PDB>_<CHAIN>.dockprep.pdbqt`, `pdb_codes.csv` with `PDB_CODE,LIGAND,RESNUM,CHAIN` and `centers.csv` containing three coordinate rows per complex. DOCK6 also uses `.dockprep.mol2` and `.noH.pdb`; redocking uses reference ligands. Auxiliary files produced by preparation automatically accompany receptors selected through pipeline connections. For an external import, upload the complete dataset.

## 6. Build and configure a pipeline

1. Open **Pipeline**. Use a preset or click a category title in the library to expand its blocks. Categories start collapsed; search opens categories containing matches.
2. Drag a block into the grid or click **+**. Drag its header to reposition it.
3. Connect the source output to the consumer input. You can also click the output first and then the input.
4. Use input names to distinguish receptors, compounds, fingerprints and poses.
5. Open configuration by double-clicking the header or using **Configure selection**.
6. In **Visual** mode, review the name, enabled state, parameters, inputs and filters. **Add input** allows several sources in one field.
7. Click **Apply configuration**. **Cancel** keeps the previous configuration.
8. Use **Save** when you need to confirm a draft; ordinary edits are also saved after a short pause.

Incompatible inputs, cycles and dependencies on disabled blocks are rejected. **Advanced** allows JSON and template editing; use it for options outside the form. Preserve technical parameter names and required template placeholders. Correct any validation error before applying.

The **ⓘ** icon on each block, in the library and canvas, opens a short explanation of accepted inputs and outputs. Incompatible connections show only “Incompatible data. Choose an output of the type expected by this input.” and preserve the pipeline. For similarity, use **ChEMBL / PubChem / ZINC → Generate fingerprints → Calculate similarity**.

**Ctrl+Z** undoes, **Ctrl+Y** redoes and **Delete** removes the selection when no configuration form is open. **Escape** cancels a connection in progress. Duplicate a block to compare parameters and rename the alternatives. Automatic layout, navigation, zoom and recentering help explore large workflows.

## 7. Understand the available connections

| Consumer | Permitted canvas input |
| --- | --- |
| Expand similar compounds | ChEMBL retrieval with original downloads; a generic compound CSV does not replace this source |
| Retrieve ZINC | Download list of type `zinc_urls` or `other`, with tranche links |
| Prepare for docking | Raw PDB or redocking receptor, plus one or more compound sources |
| Evaluate ADMET / Generate fingerprints | Compounds from retrieval, import, ADMET or graphs |
| Calculate similarity | Compatible fingerprints; do not mix types or widths |
| Filter by graphs | Output from **Calculate similarity**; an external CSV can be selected in the form |
| Run redocking | Raw PDBs with preparation enabled; prepared receptors with preparation disabled |
| Docking with Vina | Prepared receptors and selected compounds |
| Docking with DOCK6 | Prepared receptors and selected compounds or poses; optional Vina poses |
| Docking consensus | Vina poses and DOCK6 results with matching identifiers |

Compound and structure retrieval can begin independent branches. **Import my files** publishes the types of the selected files. A block's **ready-made results** mode validates existing outputs and skips its calculation; it does not turn any CSV into that block's result. See the [validation matrix and tests](pipeline_validation.md) for additional conditions.

## 8. Execute, select files and choose individual or merge

1. Execute the full pipeline or the selected stage with its dependencies.
2. If completed results are saved, answer the reuse question described in the next section.
3. Follow the execution window. When a connected input still needs automatic file selection, execution pauses and opens **Select files**. Explicit selections saved in the block skip this question.
4. Check only the files you want for each input. The list shows compatible files from connected sources, without automatically selecting all newly produced outputs.
5. Choose **Process individually** or **Merge files**.
6. To adjust parameters, click **Configure stage**. Choices already made are included in the form, avoiding another selection of the same dataset.
7. Click **Continue with selected files**. Repeat this decision for subsequent stages that need processing.

| Choice | Effect | Example |
| --- | --- | --- |
| Individual | Separate processing and outputs for each file | `alpha.csv` and `beta.csv` produce their own fingerprints and remain distinguishable |
| Merge | Combines the selected data for a joint execution | Compounds from `alpha.csv` and `beta.csv` participate in the same analysis |

For operations with distinct inputs, such as receptors and compound tables, individual processing forms combinations. DOCK6 restricts these combinations to receptors/compounds with corresponding poses. Consensus pairs Vina and DOCK6 by pose identifier, preventing combinations between different molecules or references. Merge still requires compatibility: fingerprint types, identifiers and metadata must agree. Identical filenames with different contents must be disambiguated or renamed before combining.

Structure metadata accompanies selections in both modes. **Select later** closes the window while keeping the run paused; **Configure inputs and continue** reopens it. Closing and reopening during the session preserves checked files and processing mode. Confirmed choices remain in the block configuration when you open it later. When filenames repeat, select the path distinguishing each result.

Selection occurs before consuming blocks that need execution; blocks without input do not need this decision. Reused stages do not repeat selection. The final stage makes results available for inspection without requiring another consumer. Readers follow progress; editors and owners confirm choices.

In Vina or DOCK6 docking, choosing a compound file reveals **Compound for docking (optional)**. Select an identifier to use only that molecule; **All compounds** or an empty selection uses the entire table. Each file has its own choice. Changing the file clears the previous compound. Options become available after the source produces its data; uploaded or imported files can be configured before execution.

When receptor and compound inputs already identify explicit files in the Vina/DOCK6 configuration, execution proceeds without asking again. Connections using **Identify automatically**, without a chosen file, still pause for input configuration. Files and identifiers are revalidated at execution; a choice that no longer exists produces a validation message.

## 9. Reuse data after restarting

1. Close and reopen the application with the same `--data-dir`.
2. Sign in and open the existing project without creating another project over its folder.
3. Click execute. If intact results exist, **Reuse data and results?** appears before submission.
4. Choose **Yes, reuse** to immediately complete compatible stages. This includes previous selections and avoids repeating their popups.
5. Choose **No, run again** to recalculate the requested stages and confirm inputs again where necessary.

Reuse checks parameters, connections, templates, inputs and result integrity. Moving or renaming a block does not require recalculation. Changing data, individual/merge mode or scientific options invalidates affected stages; others can still be reused. Adding a new block allows you to use completed preceding stages. Reused artifacts remain in the run that produced them.

A calculation interrupted during engine execution does not resume the engine's internal instruction. Start a new run and consider reusing previous completed stages. A run paused for selection can be reopened and continued. Cancelling the reuse question does not start a new execution.

## 10. Your first complete study with local files

This walkthrough exercises import, fingerprints, similarity, graphs and ADMET without provider queries or docking.

1. Create a project in a new folder.
2. Save the example CSV from section 5 as `molecules.csv` and upload it as **Compounds**.
3. Add **Import my files** and select that CSV.
4. Add **Generate fingerprints**; connect the import to its compound input.
5. Choose Morgan, radius 2 and 2,048 bits. These are example settings, not a validation of your study.
6. Add **Calculate similarity** and connect the fingerprints. Choose Tanimoto and a 0% threshold to observe the available relationships in this small example.
7. Add **Filter by graphs** and connect the similarity output.
8. Add **Evaluate ADMET** and connect the graph-selected compounds.
9. Execute. At each selection, check the previous block's output and confirm individual mode.
10. Under **Runs**, open fingerprints, similarities, graphs and ADMET. Check identifiers and files to verify that the dataset traversed the intended stages.
11. Execute again, choose reuse and check the **Reused** indication.
12. To compare datasets, upload another CSV, add it to the import and repeat with individual processing. Then compare with merge, observing results and graph composition.

A small dataset may produce isolated components or be entirely removed by ADMET filters. Check reports; an absence of results is not experimental evidence of efficacy or toxicity.

## 11. Configure each scientific operation

### Retrieval and expansion

Under **Retrieve compounds**, choose a target search by name/IDs/UniProt/text or a direct compound search by name/IDs/similarity/substructure. Review optional filters and limits, and choose ChEMBL/PubChem expansion. Threshold and maximum similar compounds control expansion. After execution, inspect `compounds.csv` and available ChEMBL downloads; use the appropriate table for the next analysis.

**Expand similar compounds** uses ChEMBL downloads from existing retrieval. Connect that source and review threshold/maximum; do not substitute an ordinary CSV for the ChEMBL context. **Retrieve ZINC** accepts the TXT/URI file or download script exported by ZINC tranches. Upload it directly under **ZINC tranche download list** or connect an import of type **ZINC download list**. SMI/MOL2 tranches, including `.gz` and `.bz2`, generate `compounds.csv`; MOL2 conformations travel with the table. Connect it to ADMET/fingerprints or preparation and docking. See the [retrieval guide](retrieval.md#zinc-2d-and-3d-tranches).

**Retrieve PDB structures** accepts text, PDB IDs, UniProt accessions and ligand codes; EC and collection names are optional. Organism, resolution, polymer and method refine the search. See the [retrieval guide](retrieval.md) and retrieval_report.json for criteria and outcomes. Use each PDB’s **Ligands** button to review, remove or add records before preparation or redocking. The browser-based 3D explorer offers mouse rotation, zoom and representation options; the RCSB button opens its official page. `pdb_codes.csv` and `retrieval_report.json` remain on disk but are hidden in this listing. If retrieval finds no compatible structures, revise the filters rather than proceeding with an empty dataset.

### Preparation and redocking

Under **Prepare for docking**, connect **PDB receptor (retrieval or redocking)** and **Selected compounds** separately. For a raw retrieval PDB, enter `[PDB, reference ligand, residue, chain]` records or provide `pdb_codes.csv`; receptor preparation options remain available. For redocking receptors, select the `.dockprep.pdbqt` file under **Receptor to use**: ligand-only files are hidden and receptor preparation options are disabled. Companion receptor files and binding centers follow the selection, and the receptor is reused without preparation.

Add one or more **ChEMBL**, **PubChem**, **ZINC** or uploaded compound sources under **Selected compounds**. **Ligand preparation and conformation** configures these external candidates and remains available. The PDB reference ligand defines the docking site and is not a candidate. Choose **Prepare output for → Vina, DOCK6 or Vina and DOCK6**. Each candidate is prepared once and selected formats are exported from the same conformation. Connect this block to both the receptor and compound inputs of each selected engine. Use **Merge files (merge)** to combine sources in one execution; individual mode prepares each selected combination separately. Inspect the generated receptors, candidates and centers.

Use a separate block for each receptor mode: do not combine raw PDBs and prepared receptors in the same block. pH remains available for candidate preparation when reusing a receptor. Before preparing compounds, the block checks binding centers (three finite coordinates) and receptor files required by the selected output. Vina requires `.dockprep.pdbqt`; DOCK6 requires `.dockprep.mol2` and `.noH.pdb`. Selection accepts receptor MOL2 if PDBQT is unavailable, and lists only PDBQT when both formats exist. A Vina consumer can export its rigid format in a working copy without repeating Chimera preparation. Missing files stop the process with the required filename, without repeating receptor preparation.

Under **Run redocking**, keep **Prepare complexes before redocking** enabled for raw PDBs. Disable it only for a prepared dataset. Choose Vina or DOCK6 6.11 and check its search, site and pose settings. Examine RMSDs and logs; establish the study's protocol acceptance criteria before using its receptors to dock candidates.

In **Input Data**, select and validate at least one PDB / ligand / residue / **Chain** pair; both structures use that chain. Review the table using edit/remove actions. Define common receptor and ligand preparation, with optional cofactors in a collapsed section. Resolution comes from metadata and ligand options also apply to conformation. Configured pairs are reused without another prompt. See [Redocking configuration](redocking_configuration.md) for options, validation and diagnosis.

### Inspect redocking results

Expanding a completed stage in **Runs** shows a table of PDB, ligand, residue, chain and **RMSD (Å)**. **View simulation** opens a dialog listing files for that selection. The list distinguishes the prepared receptor, reference ligand, poses and metadata. Files shared by the receptor or collection accompany their corresponding simulations.

Use each row's download button for an individual file, or **Download all (ZIP)** for the simulation bundle. The ZIP preserves directories, including when the reference ligand and Vina output share a filename. Files from other simulations are excluded.

**View 3D structure** opens the browser viewer, following the same flow as PDB inspection. It accepts PDB, PDBQT and MOL2 for inspecting receptors, ligands and poses before downloading. Use **Model** when opening a file containing multiple poses. The scene opened by the table’s **3D** button automatically selects the lowest-scoring pose, or the first unscored pose; that scene does not offer model switching. Labels and actions follow the selected language. Readers can also inspect and download results.


The **3D** button on each result opens an overlay of the best Vina pose in cyan, the crystallographic ligand in gold and the docking receptor. Enable the full crystallographic complex in the layers panel. Coordinates are preserved without realignment that could conceal displacement. If the original PDB is unavailable, an existing prepared reference is explicitly labeled as such.

**Residues** lists minimum heavy-atom distances between receptor and ligand for the pose and reference. The default cutoff is 4 Å, adjustable between 2 and 8 Å. Download the table as CSV or select a residue in the viewer to highlight it. This table represents geometric proximity. Type classification appears separately in viewer traces when chemical data allows it, as described in the next section.

### Styles, hydrogens and pose interactions

**Receptor style** and **Ligand style** are independent: for example, show receptor
ribbons with ligand spheres. **Ligand hydrogens** shows or hides only H present in
the file; it does not generate H absent from Vina output or change preparation or
scores.

**Show 3D interactions** draws each classified relation as a dashed trace: green
for hydrogen bonds, pink/coral for parallel/T-shaped π–π, lilac for hydrophobic
contacts and lime green for geometric van der Waals contacts. Hover over a trace
for type, residue, chain and distance. The legend and filters use matching
colors; click the list to highlight a residue.

Classification requires valid chemical topology and residue identities;
hydrogen bonds require explicit H. If necessary data is missing, the interface
explains the limitation and retains distance contacts. The 3D analysis is
geometric and differs from **Footprint** per-residue energies. See the
[viewer guide](molecular_viewer.md) for layers, pose selection, controls,
interpretation and updates.

### Fingerprints, similarity and graphs

**Generate fingerprints** chooses one type per block: Morgan, MACCS or pharmacophore. Radius and bits apply to Morgan. Duplicate the block to compare types and keep results identifiable.

**Calculate similarity** reads fingerprints and identifies the type of native outputs. For your own files, check the type selector. Metric and threshold belong to this stage. Approximate mode changes the search strategy; compare it with exact analysis according to your protocol.

**Filter by graphs** preserves the incoming relationships and selects the MCC, the largest connected component. Configure the MCS timeout and ring options to inspect the common fragment. SMILES come from the lineage producing the similarities; for external data, supply a corresponding molecular table when you want structures displayed.

### ADMET

**Evaluate ADMET** starts with broad exploration. Choose oral, CNS, peripheral, fragment, natural-product or custom hypotheses; annotate all compounds or filter them. BBB/HIA classes, rule violation tolerances, alerts and inclusive descriptor/model ranges are editable. RDKit provides local descriptors and estimates; optional ADMET-AI adds distinct `ML_` outputs after installation in the worker environment.

Inspect the responsive EGG, click points for molecular structures, and review `_evaluated.csv`, selection/rejection reasons and `_admet_report.json`. Local BBB/HIA and the displayed ellipses are not a validated reproduction of published BOILED-Egg. See the [ADMET guide](admet_research.md) for methods, setup, units, profiles, exports and limits.

### Vina, DOCK6 and consensus

**Docking with Vina** and **Docking with DOCK6** accept compounds from any block publishing molecules: ChEMBL, PubChem, ZINC, ADMET, fingerprints, graphs, consensus or user imports. Select a table at the compound input and connect the prepared receptor separately. Graphs, ADMET and the other engine are optional.

Reuse conformations by connecting **Vina → DOCK6** or **DOCK6 → Vina** at the compound input and selecting `docking_results.csv` or a specific pose. Identity and SMILES are retained; conversion uses the selected conformation instead of generating one from SMILES. Without an input pose, preparation generates a 3D conformer. DOCK6 uses the receptor site center to select spheres and position newly generated ligands; its optional Vina pose input remains available for refining corresponding candidates.

Review Vina box, pH and search effort and DOCK6 charges, surface, distance, radius, flexible/rigid search and footprint. Both require prepared receptors, metadata and centers; DOCK6 also needs receptor MOL2 and hydrogen-free PDB. Results show compound ID, receptor, score, SMILES, **3D** of the calculated pose and **Remove**. **3D** overlays the best pose with its docking receptor; an available prepared reference appears in gold. **Residues** provides the same contact table as redocking. New runs store receptor and reference copies alongside results to preserve this association.

The DOCK6 **Footprint** button opens an integrated popup showing per-residue Van der Waals and electrostatic energies comparing the reference with the final docked pose. Zoom and drag to explore the chart; **Recenter** restores the initial view and **Download PDF** saves the original document. The same popup is available in consensus results when a footprint is attached. Generation validates energy data and stops the stage if the engine produces an empty file. Older results may have evaluated the minimized ligand before docking; the interface flags that origin. Rerun those experiments to obtain a final-pose graph.

**Docking consensus** gathers selected results from both engines, including multiple batches and branches. It computes only the intersection by compound ID and receptor, choosing the lowest score for repeated poses. Matching identifiers with different structures are rejected. With no intersection, the block is skipped and explains that the engines did not evaluate common compounds against the same receptor.

After connecting Vina and DOCK6 results, consensus runs automatically in the pipeline without file confirmation or additional parameters. The table shows compound code, Vina score, DOCK6 score and normalized score (0 to 1). Click a column header to toggle ascending and descending order across the complete table before pagination. The normalized score averages the two engines’ min-max normalizations; higher values rank better. The CSV also retains z-score and min-max. Each row has **3D Vina**, **3D DOCK6** and **Remove** buttons per row. Consensus uses the original scores from both engines without requiring receptor, reference ligand or pose files. Available structural files are optional visualization attachments; missing files do not block computation. Singleton or constant scores normalize to zero. Removal requires editing permission, is audited and affects subsequent inputs; completed analyses are not recalculated. Each batch summary is authoritative, preventing deleted rows from returning from auxiliary tables.

### Native Amber Score · DOCK6 6.11

Under **Molecular Dynamics**, **Amber Score · DOCK6** runs native Amber Score only. It does not run conventional molecular dynamics or invoke an external AmberTools environment. UCSF Chimera continues to prepare docking/redocking inputs. AMBER components distributed within Chimera and DOCK6 remain required parts of those tools.

1. Connect receptor PDB files and ligand 3D poses, then explicitly select pairs.
2. Configure pre/post minimization, internal refinement steps, temperature, implicit GB solvent (1/2/5), cutoff and movable region.
3. Enable preservation of input charges to use Chimera-prepared MOL2 charges. DOCK6 supplies GAFF typing, topologies and complex preparation; receptor force field is ff94. No silent substitution for ff14SB/GAFF2 occurs.
4. Inspect the four energy terms, structures, logs and `protocol.json`. Valid pairs reuse checkpoints; individual failures preserve other pairs.

Native refinement uses a fixed 0.001 ps step, Langevin (gamma 2), without RATTLE constraints. The `nothing` mode evaluates energy without moving atoms. Explicit solvent boxes, pressure, conventional heating/equilibration and trajectory controls are unavailable because they are not Amber Score 6.11 parameters. Implicit solvent is not a water box.

When DOCK6 writes a final complex, results include backbone/ligand RMSD against the prepared initial complex after backbone alignment, initial/final mass-weighted complex radius of gyration, and hydrogen atom count. **Hydrogen atom count is not hydrogen-bond count.** Plots are structural comparisons, not time series. RMSF, hydrogen-bond occupancy and equilibrium claims are unavailable without a trajectory. Energy-only mode reports unavailable final-structure metrics explicitly.

Historical `dock6_sander` protocols require explicit reconfiguration; earlier results remain readable and are never relabeled as native calculations. Missing parameters for metals, cofactors or unusual residues require explicit resolution, without deleting components or approximating scores.

```bash
PYTHONPATH=src python scripts/validate_native_structures.py \
  --dock6 /path/to/dock6 --steps 2 \
  --output /tmp/amber_score_validation.json
```

This short check uses the installed 1LGU example and is not scientific acceptance for every system.


## 12. Explore results and navigate graphs

1. Open **Runs** and expand the stage, or use **Results** on the block.
2. If several tables or analyses exist, select the one you want.
3. Use search and pagination to find compounds; select 10, 25, 50 or 100 items per page where available.
4. **2D** shows the structure; **3D** generates a local SMILES conformer with rotation, zoom and recentering. This is not a docking pose.
5. Download the file to preserve the complete analysis data.

In the graph viewer, choose **Full graph** to see all components. Components occupy separate areas and vertices are spaced apart. Use the component selector to focus a region; click a vertex to inspect its details, neighbors and relationship weights. Searching a code locates the node and navigates to its component. Use dragging and zoom; **Fit to view** and **Recenter** restore an overview.

Each analysis writes one indexed full-graph CSV with `graph_analysis_id`, `graph_component` and `graph_is_mcc`. Components are filters over that file and the interactive artifact, rather than physical copies. Export CSV, SVG, PNG or metrics on demand. Older results with one CSV per component remain readable; new runs avoid that file multiplication.

Choosing **MCC** restricts exploration to the largest connected component. The panel shows **Fragment SMILES**, its structure and the common-fragment search status when available. The SMILES is extracted from a reference molecule for the MCS pattern; shared connectivity does not establish stereochemical identity across all compounds. Interpret time-limited searches according to their reported status. Without corresponding SMILES, topology remains explorable but no molecular fragment structure can be displayed.

Older results may contain only PNG. To generate interactive artifacts, execute again and choose **No, run again**. `*.biomol-view.json` visualizations larger than the 32 MB opening limit can be downloaded. In BOILED-Egg, overlapping points offer compound selection; search also distinguishes them.

## 13. Check exclusions, errors and cancellation

During processing, invalid rows are removed from working copies while originals are preserved. A stage reports the exclusion count and offers `molecule_exclusions.json`, including source, row, identifier and reason. Inspect this report before interpreting the final dataset size.

| Situation | Recommended action |
| --- | --- |
| CSV rejected | Check UTF-8, header, delimiter, field count and expected format shown in the form |
| No compatible file in the popup | Check the source's published type and outputs; select the correct source |
| Ambiguous filename | Select the full path distinguishing the run or batch |
| Incompatible fingerprints | Use the same type and width; process alternatives separately |
| Metadata/centers do not match | Check PDB, ligand, residue and chain; upload the complete prepared dataset |
| Poses do not match compounds/receptors | Check identifiers and select Vina/candidates for the same reference |
| External query failure | Read query context in the log; review filters and retry after checking the provider |
| Engine missing or command failed | Check the worker environment, executables, input files and command log |
| Descendant stage prevented from running | Correct the failed source and execute again |
| Results were not reused | Check whether inputs, configuration, processing mode or output integrity changed |

**Minimize** keeps progress in a project strip; **Track run** reopens the window. To interrupt, use **Cancel run** and wait for the final state. Cancellation may take time while the backend stops the worker. Completed results from previous stages can be reused in a new run; a partial stage is not treated as completed.

## 14. Curation, collaboration and history

**Remove** in a compound table asks for confirmation and removes the row only from that CSV. It does not automatically synchronize other tables. Explicitly choose the table feeding the next stage. Removing a file requires disconnecting inputs using it. Editing/removal requires an editor or owner and no active pipeline; originals and authorship remain in history. Changed data invalidates dependent analyses.

To share, the owner opens **Share**, enters the email of an existing account and chooses **Editor** or **Reader**. The recipient must accept in the workspace. Invitations are internal, without email delivery. Readers inspect and download; editors also configure, select files and execute; owners manage collaborators and project deletion.

Collaborator changes update in the open project. Forms and drafts being edited are preserved; a conflict on the same field requires review before reapplying. Use block names and descriptions to make alternatives understandable to the team.

In **History**, inspect the date, user and change. As owner, **Restore before** returns configurations, files, results and permissions to the state before the selected record; read the confirmation and cancel any active run first. Restoration is also recorded. See [projects and versions](projects.md) for effects on collaborators and imported versions.

## 15. Export, back up and migrate

1. Complete or cancel active runs.
2. Use **Export** on the project and preserve the `.bme.zip` package.
3. For a complete installation backup, close the application and also copy the accounts workspace and project folders, including hidden files.
4. On another computer, install the environment, create/sign in to an account and use **Import project**.
5. Choose a new or empty folder on the backend computer and import the package.
6. Open the project, inspect inputs/results and engine installations. On execution, decide whether to reuse compatible data.

`project.json` alone does not contain results or inputs. The package includes data, configurations and history; accounts and passwords are not exported. The importing user becomes owner. Collaborators need new acceptance at the destination. Large packages generated in browser mode remain in `.exports/` for direct copying; see the limits in the projects guide.

## 16. Glossary and reference workflow

| Term | Meaning in this application |
| --- | --- |
| Block / stage | A configured pipeline operation |
| Input / binding | A file or result connected to an operation field |
| Artifact | A file produced by a stage |
| Fingerprint | Molecular representation used for comparison |
| Similarity | Numeric relationship between two molecular identifiers |
| Connected component | Group of vertices joined by paths in a graph |
| MCC | Largest connected component |
| MCS | Common substructure search, which may have a timeout |
| SMILES | Textual representation of a molecular structure |
| Redocking | Repositioning a reference ligand to evaluate a protocol |
| Reuse | Using complete, compatible, previously persisted results |

A common workflow is retrieve/import compounds → fingerprints → similarity → graphs → ADMET, alongside retrieve/import PDBs → preparation/redocking. The branches converge at Vina → DOCK6 → consensus. Choose files and processing mode at each handoff; record parameters, exclusions and scientific criteria with your study. The [validation report](pipeline_validation.md) explains tests and verification limits for this chain.

## Flexible retrieval

Use the [retrieval guide](retrieval.md) to search PDB by text, IDs, UniProt or ligand codes without requiring EC. ChEMBL supports target and direct compound searches by name, IDs, similarity and substructure, with optional filters and explicit limits. Review retrieval_report.json before passing selected files downstream.

## Checks before execution

Pair selection is mandatory and uses one **Chain** for both receptor and ligand. A configured redocking stage reuses its selected pairs; the pipeline requests selection only when it is missing. Resolution comes from metadata and has no editable field. The table lists pairs with edit/remove actions. Cofactors and solvent, hydrogen, minimization and charge options are defined once for all receptors and once for all ligands. Ligand preparation and conformation share those options.

Before calculations, the system checks residues, chain and cofactors in the PDB files. Choose **AutoDock Vina** or **DOCK6 6.11** and fill in its displayed settings. Prepared inputs require PDBQT for Vina or MOL2 and hydrogen-free receptor PDB for DOCK6, plus valid centers. Chimera and Open Babel are needed for preparation; the selected engine and its tools must be installed. The scientific interpreter's directory is also included in subprocess PATH.

Tool failures report the executable, exit code and process message. Failed Chimera scripts remain available for diagnosis. Correct the reported input or installation and retry. There is no verbose control: Vina verbosity stays at zero.

See [Redocking configuration](redocking_configuration.md) for the complete workflow.

### Compounds excluded during preparation

Candidate preparation for Vina and DOCK6 handles compounds independently. If conformation generation, conversion, charge assignment or minimization fails, that compound is excluded from docking inputs while the others continue. Partial compound outputs are removed. If no candidate can be prepared, the engine is not started and the stage explains why.

In **Runs → Results by stage**, use **Excluded compounds** to inspect the paginated list and **Download excluded compound list** to download `preparation_failures.csv`. The report records identifiers, submitted SMILES, input file/row, preparation context and failure reasons. These failures also appear in the general exclusion report. Original inputs are preserved. Infrastructure failures, such as a missing engine or full disk, remain stage errors.

## Reusing receptors in Vina

Connecting prepared redocking outputs to docking reuses the receptor and binding center across batches. The engine's required files and center coordinates are checked before docking. Missing prepared files cause an input error rather than silently triggering preparation. Raw structure preparation converts each PDB/chain receptor once, even when several reference ligands share it.

The standard workflow uses Vina's Python API to retain the receptor and affinity maps in memory during a batch. Maps are computed before loading any ligand to cover every atom type, following the [official Vina documentation](https://autodock-vina.readthedocs.io/en/latest/docking_python.html). Receptor, center, box dimensions or initialization parameter changes create a new session. Template settings remain respected; unsupported special options or a missing Python API use the CLI executable and log this fallback.

Events `receptor.reused`, `ligand.preparation.started`, `vina.maps.started` and `vina.maps.reused` distinguish the operations. `prepare.com` prepares an external candidate, which can be slow with AM1 charges or minimization; it does not prepare the receptor. Already prepared candidates are reused. The configured seed and search parameters are preserved, with compounds processed in stable order. Reuse applies within each batch; in-memory maps are not shared between workers.


## DOCK6 validation and exclusions

The workflow checks the installation, global parameters, prepared receptor, site center, spheres, box and grid files. Each ligand MOL2 is checked for atom/bond counts, identifiers, finite coordinates and charges, chemical topology and connectivity. Disconnected or chemically invalid ligands are excluded without guessing replacement bonds. This matters when converting Vina PDBQT poses, which lack the full chemical information of MOL2. A requested Vina pose must match the selected compound and receptor; another input is not silently substituted.

The surface, spheres and grid are shared by compounds using the same receptor/site. Minimization, docking and footprint use independent processes per compound. Native failures (including SIGSEGV, exit code -11), missing/empty outputs and invalid scores exclude that compound from subsequent stages and final results while others continue. Footprint plot failures are also recorded individually. Generated configurations and tool logs are preserved for diagnosis. Identifiers containing dots remain intact.

The **Excluded compounds** list now covers preparation and execution failures, including step, identifier, SMILES, input and reason. It retains the `preparation_failures.csv` filename for compatibility with the Vina report. If every compound is excluded, the stage reports that no usable result was produced. Global errors — missing installation, incomplete parameters, unavailable receptor/site/grids or filesystem failures — still stop the stage.

The review follows the [official DOCK6 manual](https://dock.compbio.ucsf.edu/DOCK_6/dock6_manual.htm), [rigid/flexible docking tutorial](https://dock.compbio.ucsf.edu/DOCK_6/tutorials/ligand_sampling_dock/ligand_sampling_dock.html) and [grid generation tutorial](https://dock.docking.org/DOCK_6/tutorials/grid_generation/generating_grid.html). Structural checks cannot guarantee the native engine is free of defects; compound isolation allows the batch to continue and identifies affected inputs.

## Navigation and result maintenance

Workspace actions are in the sidebar. Inside a project, **Project navigation** and **Pipeline and execution** are collapsible, with contextual help. Use the expanded canvas, zoom and **Arrange blocks** for large pipelines. The [interface guide](frontend.md#navigation-canvas-and-completed-runs) explains controls, saved selections, completed-run cards and result removal.

In 3D views, **Show residue names** allows cleaner exports while preserving interactions. See the [viewer guide](molecular_viewer.md#hide-residue-names).

[View a DOCK6 footprint](../../imgs/examples/dock6-footprint.png) and its [original PDF](../../imgs/examples/dock6-footprint.pdf): the 1ABE tutorial example compares reference/pose energy contributions after Vina → DOCK6. These energy profiles are distinct from geometric 3D contacts.


## Exporting the pipeline as PNG

The camera in the canvas controls saves the **entire pipeline as PNG**, including blocks outside the visible area. The image preserves positions, names, colors, ports and connections independently of zoom. The file includes 300 dpi metadata and increased resolution, capped at 32 million pixels to control memory use. Viewers can also export the diagram for articles and presentations.

## Recognizing dataset origins and following new inputs

Results display the original stage and target identifier, such as **ACHE · CHEMBL220**, **BCHE · CHEMBL1914** and **BACE1 · CHEMBL4822**. Exact file provenance follows ADMET, PubChem, fingerprints, similarity and graphs, including cached results and repeated names such as `compounds.csv`.

Search input selectors by target or ChEMBL identifier and filter result tables by origin. ADMET and graph dashboards show **Data source** and **File and pipeline path**. Graph molecular sources, compound inspection and exported reports retain this identity. Companion files belong to the selected batch rather than every analysis with the same filename.

Adding the three direct ADMET BBB+ datasets alongside three PubChem datasets produces six separate fingerprint datasets. An older similarity or graph selection that covered every previous output follows this growth; the effective selection is audited and saved. The result is six graphs, one direct ADMET dataset and one PubChem dataset per target.

**Follow new datasets from the source** includes compatible outputs from connected producers. Disable it to use only explicitly selected files. Older partial selections are not automatically expanded. Use **Process individually** to keep the six datasets separate; merge deliberately combines selected inputs.

Source labels are reconstructed from actual file inputs and execution history, including earlier runs that produced reused files. They are not guessed from filename patterns. Existing results receive the improved labels when opened.


## Quick examples for each feature

The following four-line examples describe the current interface. IDs in parentheses match the catalog; import is a workspace block. Use data and criteria appropriate to your study. External searches require network access; preparation and docking require the engines listed in the installation guide.

### Import my files (`import_results`)

```text
Add the import block and choose the data type.
Upload compounds.csv with molecule_chembl_id and canonical_smiles.
Validate the file and connect its output to Evaluate ADMET.
Run the consumer block and inspect its result table.
```

### Retrieve ChEMBL (`retrieve_compounds`)

```text
Add Retrieve ChEMBL and search for target CHEMBL220.
Review bioactivity filters and similar-compound expansion.
Save the configuration and run the block.
Select molecules or similars at the next block input.
```

### Retrieve PubChem (`retrieve_pubchem`)

```text
Add Retrieve PubChem and choose a manual reference.
Select CID and enter 2244; set cutoff 90 and maximum 10.
Save and run the search, which requires network access.
Connect the compound table to ADMET or fingerprints.
```

### Expand similar compounds (`expand_similar_compounds`)

```text
Connect a ChEMBL retrieval to Expand similar compounds.
Choose the source collection and the desired threshold.
Save and run the expansion.
Inspect new compounds and relationships before selecting its output.
```

### Retrieve PDB structures (`retrieve_structures`)

```text
Add Retrieve PDB structures and enter ID 1P0I.
Set the target name and review filters and ligand screening.
Run the block and open the retrieved records.
Inspect chain, ligand and 3D structure before redocking.
```

### Retrieve ZINC (`retrieve_zinc`)

```text
Obtain a list of ZINC SMI or MOL2 tranche URLs.
Upload the list and select it in Retrieve ZINC.
Set concurrent downloads and run the block.
Inspect the report and connect the compound table to the next block.
```

### Prepare for docking (`prepare_structures`)

```text
Connect the receptor and compound table at separate inputs.
Choose Vina, DOCK6 or both as output formats.
Review cofactors, hydrogens, gas/am1 charges and minimization.
Run and connect prepared receptor and candidates to the selected engine.
```

### Evaluate ADMET (`admet`)

```text
Connect a compound table to Evaluate ADMET.
Choose an oral profile, annotation mode and local RDKit method.
Review BBB/HIA, rules and alerts; save and run.
Open the EGG and tables; select a compatible output to continue.
```

### Generate fingerprints (`fingerprints`)

```text
Connect a table containing identifiers and valid SMILES.
Choose Morgan with radius 2 and 2048 bits.
Save and run; inspect the generated collection.
Connect the output to Calculate similarity.
```

### Calculate similarity (`similarity`)

```text
Connect previously calculated Morgan fingerprints.
Choose Tanimoto and a 70% threshold.
Review approximate mode before running.
Connect relationship CSVs to Filter by graphs.
```

### Filter by graphs (`graphs`)

```text
Connect calculated similarity relationships to the graph block.
Review common-fragment options and run.
Open Full Graph or MCC and inspect components and molecules.
Select the analysis compound table for ADMET or docking.
```

### Run redocking (`redocking`)

```text
Connect PDBs and select receptor/ligand/residue/chain pairs.
Validate curation and review cofactors and Chimera preparation.
Choose Vina or DOCK6 6.11, review its options and run.
Inspect RMSD, logs, Residues and the 3D button for completed pairs.
```

### Docking with Vina (`docking_vina`)

```text
Connect the receptor and compounds at their respective inputs.
Select the receptor, molecular table and required preparation.
Review box, pH, search effort and pose count.
Run and inspect scores and the pose/receptor 3D comparison.
```

### Docking with DOCK6 (`docking_dock6`)

```text
Configure the DOCK6 installation containing bin/ and parameters/.
Connect a compatible receptor and compounds; Vina poses are optional.
Review charges, surface, spheres and flexible or rigid search.
Run and inspect scores, 3D poses and available footprints.
```

### Docking consensus (`consensus`)

```text
Connect Vina and DOCK6 results at their respective inputs.
Use matching compounds and receptors in both result sets.
Review repulsion weight and run consensus.
Sort normalized scores and inspect the available poses.
```

### Projects, language and collaboration

```text
Choose English at sign-in and access your account.
Create a project and select its destination folder.
Open Share and invite an existing account as viewer or editor.
Inspect history to identify changes and versions.
```

### Execution, cancellation and reuse

```text
Configure and save blocks; connect their dependencies.
Run the pipeline and monitor stages and logs in Runs.
Use Cancel run to stop an active run.
Reopen the project and reuse compatible results when continuing.
```

### Input selection and merge

```text
Connect two compound sources to a consumer block.
Select only the produced outputs you want to process.
Use individual processing to keep collections separate.
Choose merge when compatible inputs should be combined.
```

### 2D, 3D and contact inspection

```text
Open a compound table and inspect its 2D structure.
In docking or redocking results, click 3D.
In the viewer tab, rotate the scene and change styles and layers.
Open Residues to inspect geometric contacts.
```

### Results, cleanup and history

Click data-column headers to toggle ascending and descending order. The arrow
shows the direction; missing values remain last. Sorting applies to the filtered
dataset before pagination and preserves filters and the open stage. Use **View
table** to inspect every CSV column, including ADMET properties and alerts,
fingerprints, similarities and reports. This inspection is read-only; sorting
neither modifies files nor executes the pipeline.

Two actions with red icons appear beside **Refresh Results**:

| Action | Removes | Preserves |
|---|---|---|
| **Remove All** | Generated results and computational caches in the project | Blocks, connections and every setting |
| **Clear All** | Those results, input selections, scientific settings and custom templates | Blocks and connections; settings return to defaults |

Both preserve original library files, logs and history. Review the dialog before
confirming; canceling leaves data intact. Finish or cancel active runs before
cleaning. Configure the blocks again after Clear All.

With reuse enabled, valid batches remain available after restarting. Configuration,
input and result equivalence is checked when saving and before execution; no manual
synchronization is needed. Adding independent inputs or actions preserves earlier
results and processes only the complement. Changing parameters or removing inputs
or actions requires confirmation before saving. The dialog lists affected blocks,
including downstream dependents. Confirm saves the change and deletes obsolete
results; Cancel restores the previous configuration and preserves files. Other
inputs and shared files survive. A run still using affected data must finish or
be cancelled before that change.

For example, adding MACCS to completed Morgan work executes only MACCS. Changing
Morgan's radius afterwards preserves MACCS. Removing one file from an individual
set removes only its batch and dependent batches. A merged result represents the
entire set: changing its members requires rebuilding that combination, without
approximating the calculation using incompatible individual batches.

```text
Finish or cancel the run before removing its data.
Open Results by stage and choose Remove all results.
Review the confirmation and apply the intended removal.
Use history to restore result files if needed.
```

### Export and migrate

```text
Use the canvas camera to export the pipeline as PNG.
Download tables, figures or the simulation ZIP from Results.
Export the project from its card to preserve configuration and data.
Import the package in the destination workspace and check its tools.
```


## Identify ADMET output collections

New ADMET outputs include their collection origin before the result type:

| Tag | Origin |
| --- | --- |
| `MOL` | Reference molecules from ChEMBL selection or a file identified as molecules |
| `SIM` | Similar molecules from ChEMBL selection or a file identified as similars |
| `PUBCHEM` | Compounds retrieved from PubChem |
| `ZINC` | Compounds retrieved from ZINC |
| `MIX` | Merged categories or a consolidated collection |
| `INPUT` | A file whose category cannot be determined safely |

Examples: `ACHE_molecules_MOL_BBB+.csv`, `ACHE_similars_SIM_HIA-.csv` and
`ACHE_molecules_MOL_egg.png`. The same identification accompanies evaluated,
selected/rejected tables, interactive plots, method reports and exclusion reports
when produced. Pipeline lineage supplies the target name; standalone files use
conservative filename-based identification. Unsafe characters and long names are
normalized with a short identifier when needed.

Readers and downstream blocks still recognize `_evaluated`, `_BBB+`, `_HIA-`,
`_egg` and `_admet_report` endings. Naming applies to new outputs; historical files
are neither renamed nor required to be recalculated to remain readable.
