# User manual: from first sign-in to results

[Documentation](../README.md) · English · [Português](../user_manual.md)

This manual follows a study from choosing a language to exporting the project. Button names below refer to the English interface. For technical parameters and APIs, see the [Flet workspace](frontend.md), [backend usage](backend_usage.md), [projects and versions](projects.md) and [pipeline validation](pipeline_validation.md) guides.

## 1. Prepare the environment and start the application

Download the [project from GitHub](https://github.com/mpiress/BioMolExplorer) and enter its source code folder:

```bash
git clone https://github.com/mpiress/BioMolExplorer.git
cd BioMolExplorer
```

Without Git, use **Code → Download ZIP**, extract the archive and open a terminal inside the extracted folder. The [installation and configuration guide](installation.md) explains downloading, external tools and paths.

**Install UCSF Chimera 1.17 and DOCK6 6.11 before starting the application.** They must be available on the computer running calculations, including in web mode. Installing the interface and Conda environment does not install these tools; stages requiring them will fail without them. Complete the installation guide's checks before running the startup commands.

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

Search by name or description to find a project. Archiving removes it from the active list; restoring allows you to use it again. Deleting removes normal project access while preserving its files on disk.

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
4. Review the **Type / File / Actions** table. You can change each row's type (with fresh validation) or use **Remove from list**. Removing a row preserves the file in the project library.
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
3. Follow the execution window. When a stage needs inputs from its sources, execution pauses and opens **Select files**.
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

Use a separate block for each receptor mode: do not combine raw PDBs and prepared receptors in the same block. pH remains available for candidate preparation when reusing a receptor. Before preparing compounds, the block checks binding centers (three finite coordinates) and receptor files required by the selected output. Vina requires `.dockprep.pdbqt`; DOCK6 additionally requires `.dockprep.mol2` and `.noH.pdb`. PDBQT remains the receptor selection file for either engine. Missing files stop the process with the required filename, without repeating receptor preparation.

Under **Run redocking**, keep **Prepare complexes before redocking** enabled for raw PDBs. Disable it only for a prepared dataset. Check records, box dimensions, search effort and pose count. Examine RMSDs and logs; establish the study's protocol acceptance criteria before using its receptors to dock candidates.

In **Input Data**, select and validate at least one PDB / ligand / residue / **Chain** pair; both structures use that chain. Configure cofactors and receptor/ligand preparation for that selection. Resolution comes from metadata and ligand options also apply to conformation. Configured pairs are reused without another prompt. See [Redocking configuration](redocking_configuration.md) for options, validation and diagnosis.

### Inspect redocking results

Expanding a completed stage in **Runs** shows a table of PDB, ligand, residue, chain and **RMSD (Å)**. **View simulation** opens a dialog listing files for that selection. The list distinguishes the prepared receptor, reference ligand, poses and metadata. Files shared by the receptor or collection accompany their corresponding simulations.

Use each row's download button for an individual file, or **Download all (ZIP)** for the simulation bundle. The ZIP preserves directories, including when the reference ligand and Vina output share a filename. Files from other simulations are excluded.

**View 3D structure** opens the browser viewer, following the same flow as PDB inspection. It accepts PDB, PDBQT and MOL2 for inspecting receptors, ligands and poses before downloading. Use **Model** in the viewer for files containing multiple poses. Labels and actions follow the selected language. Readers can also inspect and download results.

### Fingerprints, similarity and graphs

**Generate fingerprints** chooses one type per block: Morgan, MACCS or pharmacophore. Radius and bits apply to Morgan. Duplicate the block to compare types and keep results identifiable.

**Calculate similarity** reads fingerprints and identifies the type of native outputs. For your own files, check the type selector. Metric and threshold belong to this stage. Approximate mode changes the search strategy; compare it with exact analysis according to your protocol.

**Filter by graphs** preserves the incoming relationships and selects the MCC, the largest connected component. Configure the MCS timeout and ring options to inspect the common fragment. SMILES come from the lineage producing the similarities; for external data, supply a corresponding molecular table when you want structures displayed.

### ADMET

**Evaluate ADMET** calculates properties and applies the available molecular filters. Inspect tables, exclusions and the BOILED-Egg plot. BBB/HIA and other indicators are computational estimates and rules; they do not replace experimental validation. Consult the backend guide for the calculations actually implemented.

### Vina, DOCK6 and consensus

**Docking with Vina** and **Docking with DOCK6** accept compounds from any block publishing molecules: ChEMBL, PubChem, ZINC, ADMET, fingerprints, graphs, consensus or user imports. Select a table at the compound input and connect the prepared receptor separately. Graphs, ADMET and the other engine are optional. PubMed provides literature references; molecules obtained from those references must be imported with identifiers and SMILES or valid structures.

Reuse conformations by connecting **Vina → DOCK6** or **DOCK6 → Vina** at the compound input and selecting `docking_results.csv` or a specific pose. Identity and SMILES are retained; conversion uses the selected conformation instead of generating one from SMILES. Without an input pose, preparation generates a 3D conformer. DOCK6 uses the receptor site center to select spheres and position newly generated ligands; its optional Vina pose input remains available for refining corresponding candidates.

Review Vina box, pH and search effort and DOCK6 charges, surface, distance, radius, flexible/rigid search and footprint. Both require prepared receptors, metadata and centers; DOCK6 also needs receptor MOL2 and hydrogen-free PDB. Results show compound ID, receptor, score, SMILES, **3D** of the calculated pose and **Remove**. This viewer opens the structure produced by docking.

**Docking consensus** gathers selected results from both engines, including multiple batches and branches. It computes only the intersection by compound ID and receptor, choosing the lowest score for repeated poses. Matching identifiers with different structures are rejected. With no intersection, the block is skipped and explains that the engines did not evaluate common compounds against the same receptor.

The consensus table shows Vina and DOCK6 scores, SMILES, receptor and z-score/min-max normalizations, with **3D Vina**, **3D DOCK6** and **Remove** buttons per row. The DOCK6 consensus score is `min(0, Grid_Score + repulsion_weight × Internal_energy_repulsive)`; missing repulsion is zero. Singleton or constant scores normalize to zero. Removal requires editing permission, is audited and affects subsequent inputs; completed analyses are not recalculated. Each batch summary is authoritative, preventing deleted rows from returning from auxiliary tables.

## 12. Explore results and navigate graphs

1. Open **Runs** and expand the stage, or use **Results** on the block.
2. If several tables or analyses exist, select the one you want.
3. Use search and pagination to find compounds; select 10, 25, 50 or 100 items per page where available.
4. **2D** shows the structure; **3D** generates a local SMILES conformer with rotation, zoom and recentering. This is not a docking pose.
5. Download the file to preserve the complete analysis data.

In the graph viewer, choose **Full graph** to see all components. Components occupy separate areas and vertices are spaced apart. Use the component selector to focus a region; click a vertex to inspect its details, neighbors and relationship weights. Searching a code locates the node and navigates to its component. Use dragging and zoom; **Fit to view** and **Recenter** restore an overview.

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

Pair selection is mandatory and uses one **Chain** for both receptor and ligand. A configured redocking stage reuses its selected pairs; the pipeline requests selection only when it is missing. Resolution comes from metadata and has no editable field. Cofactors and solvent, hydrogen, minimization and charge options are configured per pair. Ligand preparation and conformation share those options.

Before calculations, the system checks residues, chain and cofactors in the PDB files. Prepared inputs require PDBQT files and three finite coordinates in `Prepared/centers.csv`. The scientific environment must provide Chimera, Open Babel (`obabel`) and Vina; without preparation, only Vina is required. The scientific interpreter's directory is also included in subprocess PATH.

Tool failures report the executable, exit code and process message. Failed Chimera scripts remain available for diagnosis. Correct the reported input or installation and retry. There is no verbose control: Vina verbosity stays at zero.

See [Redocking configuration](redocking_configuration.md) for the complete workflow.
