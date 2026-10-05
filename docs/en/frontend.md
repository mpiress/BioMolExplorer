# Flet workspace

[Documentation](../README.md) · English · [Português](../frontend.md)

The interface brings together authentication, private projects, invitations,
a file library, a pipeline editor and execution monitoring. It uses the backend
Python services; scientific calculations continue to run in separate processes.

## Portable collaborative projects

Creating a project requires selecting a folder through the button inside **Project folder**; the field cannot be typed into. On desktop, use the system picker; on the web, browse and create folders on the computer running BioMolExplorer. Choose a new or empty folder; its
inputs, configuration, history and results are stored there. In web mode, this
path refers to the computer running the backend. The workspace offers **Import project**,
**Export** and **History**, including complete version restoration. Configuration
is saved automatically and collaborator updates appear in open projects. Blocks
accept multiple inputs and validated completed results. See
[the projects and versions guide](projects.md) for instructions and formats.

## Getting started

Before starting, download the source code from [GitHub](https://github.com/mpiress/BioMolExplorer) and install **UCSF Chimera 1.17, DOCK6 6.11 and DMS** on the computer running calculations. Follow the [installation and configuration guide](installation.md) to prepare the environment, configure `PATH` and check executables. Installing the interface does not install these tools; stages requiring them will fail if they are absent.

From the downloaded source code root, with the scientific environment prepared:

```bash
conda activate BioMolExplorer
python -m pip install -e '.[ui]'
biomolexplorer-ui --web --language en --dock6-path /path/to/dock6-6.11
```

Open `http://127.0.0.1:8550`. To serve the interface without opening a browser
automatically:

```bash
biomolexplorer-ui --web --no-browser --data-dir ./workspace-data --port 8550 --dock6-path /path/to/dock6-6.11
```

`--worker-python /path/to/python` selects the calculation environment.
`--dock6-path /path/to/dock6` sets the DOCK6 installation available to stages.
Chimera, Open Babel, Vina, DOCK6 and their dependencies remain necessary for the
respective protocols. Installing the interface does not install those engines.

The default data directory is `~/.local/share/biomolexplorer`. It contains the
SQLite account/project database and temporary uploads. Inputs, configuration,
versions and results reside in each selected project folder. Use a persistent directory
and back it up while the application is stopped. One instance supervises each
directory; two processes cannot run the same workspace simultaneously.

## Choose a language

Start with `--language en` or `--language pt`; without this option the interface starts in English. In the top bar, the flag menu in the right corner selects English (United States) or Portuguese (Brazil) for the session and preserves filled fields. User names, filenames, molecular data and scientific parameter keys retain their original values. To switch during use, open the same flag menu without signing out. See the [detailed user manual](user_manual.md).

## Accounts and collaboration

Create an account on the initial screen. Passwords have at least 10 characters
and are stored using an individual salt and scrypt. Sessions last eight hours
and are invalidated on logout.

Each project starts private. In the **Compartilhar** (Share) tab, the owner enters
the email of an already registered account and chooses a permission. An invitation
appears in the recipient's workspace and must be accepted. The owner can revoke it.
Invitations show whether the role is **Editor** or **Leitor** (Viewer). Declined
or revoked invitations cannot be accepted later. Changing a collaborator's role
creates a new invitation and suspends previous access until acceptance. Inviting
an active collaborator again with the same role does not remove access.
If the role changes while an invitation is displayed, the application refreshes
it before allowing acceptance. Authorization is checked on requests: when access
loss is detected, the interface closes project viewers and returns to the
workspace. Expired sessions return to authentication. Responses started in a
previous session cannot reopen projects or charts after logout.
The interface supports English and Portuguese. This reference guide retains some
Portuguese labels with English explanations; the user manual follows the English controls.

| Permission | Capabilities |
| --- | --- |
| Viewer | View pipelines, files, results and logs; download artifacts |
| Editor | In addition to viewing, configure stages, upload files and run/cancel analyses |
| Owner | In addition to editing, invite/revoke collaborators and delete the project |

Invitations are delivered inside the application. This version does not send
email, verify ownership of email addresses, recover passwords or integrate with
an institutional identity provider. Before publishing an institutional service
on the internet, configure HTTPS, verified identity and infrastructure controls.
Local workers share the operating system user: isolation is implemented through
application authorization and path restrictions, not containers.

## Projects and pipelines

Create projects with a name, description, tags and color. Search by name or tags.
Projects can be archived, restored or deleted. Deletion removes access and keeps
the files in storage; it does not physically purge data.

In **Pipeline**, click a category title in the stage library to expand its blocks; categories start collapsed. Search expands categories containing matches. Drag blocks onto the grid, or use
**+**. Drag a block header to move it. Drag an output to another block's input,
or click the output followed by the input. Colors identify retrieval, analysis,
docking and imports. Incompatible connections and cycles are rejected without
changing the graph.

Double-click a header or use **Configurar seleção** (Configure selection). The
**Visual** mode provides forms, lists, selectors and controls for ChEMBL filters,
PDB records, docking boxes and engine options. Scripts are generated from these
choices. **Avançado** (Advanced) allows JSON parameter/connection editing and full
scientific templates. **Aplicar configuração** validates and commits the draft;
**Cancelar** keeps the previous block configuration.

After a stage completes, execution pauses before the next connected block and
automatically opens a file selection popup. Check the compatible files to use for
each input, then click **Continuar com os arquivos selecionados**. Only outputs
from this run are listed. Choose **Processar individualmente** (the popup default)
to keep each file and its results separate, or **Mesclar arquivos (merge)** to
combine selected inputs. Multiple input ports form file combinations in individual
mode. DOCK6 retains only matching receptor/compound/pose combinations; consensus pairs matching Vina/DOCK6 identifiers. Required structure metadata accompanies selected receptors in both modes.
The choice is recorded for each stage. Stages that need processing require
confirmation; reused stages complete without repeating file selection.
Closing the popup keeps the run paused and preserves the selected files and mode.
Reopen it from the progress window or the runs tab; **Configurar etapa** opens
the full configuration with those choices already filled in. Confirmed selections
are also saved in the block. Resuming keeps
completed outputs and follows the dependencies defined by the connections.

Use zoom, canvas navigation, automatic layout, duplication, undo and redo.
The pipeline, charts and images support zoom from 0.1% to 10,000%, including
when the content is smaller than its viewport. Reset the view to restore its
initial position and scale. The 3D molecule zoom control uses the same range,
with a logarithmic slider for fine adjustments at small and large scales.
**Ctrl+Z**, **Ctrl+Y** and **Delete** operate while configuration is closed;
**Escape** cancels a pending connection. Undo history belongs to the editing
session; positions and connections autosave after a short editing pause.
**Salvar** (Save) also commits changes. The selection
menu and buttons provide alternatives to dragging. Visual positions do not define
execution order: connections define dependencies. An enabled stage cannot depend
on a disabled stage.

Run the entire pipeline or the selected stage with its dependencies. Each run
records a configuration snapshot. When saved results exist, **Reaproveitar dados
e resultados?** asks whether to reuse them. **Sim, reaproveitar** completes
compatible stages from those results, including after restarting the application
or adding new blocks. The system compares parameters, templates,
inputs and output integrity. Moving or renaming blocks does not repeat work;
changed parameters or inputs recalculate affected stages. Select **Não, executar
novamente** to recalculate all stages. Reused
outputs remain in the directory of the run that produced them.

Starting a run opens a dialog showing the current stage, elapsed time, every
block's state and the number of completed stages. ChEMBL retrieval also reports
the query phase and the number of processed molecular records. The animated
indicator shows ongoing work without estimating a percentage for calculations
whose duration is unknown. **Minimizar** (Minimize) keeps a project status banner;
**Acompanhar execução** (Track run) reopens the dialog. Completion, cancellation
or failure is shown with access to stage logs and results. Editors and owners can
request **Cancelar execução** (Cancel run); the request remains visible while
the backend stops the process. Readers can track runs.

**Execuções** (Runs) shows states, errors, logs and downloadable files. Cancellation stops the worker;
after an application interruption, the previous run is recorded as interrupted.
A new run reuses compatible completed stages and executes the remaining ones.
The progress dialog and history label these stages **Reaproveitado** (Reused).
Every two seconds, open projects load collaborator changes when there is no
open configuration form or local draft. If permission changes to viewer, editing drafts and private dialogs are
discarded when the change is detected.

Templates include retrieval → ADMET, graph selection → ADMET, user-supplied
compounds → ADMET, user-supplied PDBs → preparation and a complete route through
consensus. They are starting points: review targets, inputs and scientific choices
before running. DOCK6 requires a configured installation and a reference complex.

## Exploring charts and compounds

In **Execuções** (Runs), expand **Recuperar compostos** (Retrieve compounds) or
**Expandir similares** (Expand similar compounds). The compound table appears
directly within the stage. The table
selector lists available `<target>_FULL`, `<target>_MOLS`, `<target>_SIMS` CSVs and
`compounds.csv`, the integrated ChEMBL/PubChem dataset. Columns show compound codes,
SMILES and actions, with search and a selected-table download. **Itens por página**
(Items per page) offers 10, 25, 50 or 100 rows; the default is 25.
Only generated files appear; runs without similar compounds may omit some CSVs.

**2D** displays the molecular structure. **3D** generates a local conformer from
SMILES with drag rotation, zoom and reset; it is not a docking pose. Readers can
also use these previews.

**Remover** (Remove) asks for confirmation and deletes only the selected CSV row,
preserving its other columns. Editors and owners can curate tables while no
pipeline is active in the project. The original file and authorship record remain
in private curation history (`.curation/` and SQLite table `compound_edits`).
Retrieval remains reusable after curation; stages consuming changed data are
recalculated on the next run. Other CSVs are not synchronized: select the table
to use in the next block's input connection.

After running a stage, click the block's **Resultados** (Results) icon, or open
**Execuções** (Runs) and expand the stage. Fingerprints, similarity and other
file lists use paginated tables with download, preview when available, and remove
actions on the right. **Arquivos** (Files) uses the same layout. Removing files
requires editor access and no active pipeline; unlink input files from blocks
before deleting them. Project history can restore removed files. Stages with
removed outputs cannot be reused as complete results on the next execution.

In **Filtrar por grafos** (Graph filtering), **Similaridades** is the only
canvas input. Connect one or more **Calcular similaridade** (Compute similarity)
stages and select their output files. Each CSV produces an independent analysis.
Fingerprint, compound retrieval and generic import stages cannot connect directly.
Metric, fingerprint type and threshold belong to the similarity stage; graph
filtering preserves the relationships and weights it receives.

Optionally, **Enviar meus arquivos** uploads external UTF-8 CSVs with
`source,target,value` columns, valid identifiers and finite numeric weights
between 0 and 1. Example: `MOL1,MOL2,0.85`. Invalid headers, missing identifiers
and invalid weights produce messages explaining the expected format.

Pipeline compound SMILES are identified automatically. For external inputs,
**Arquivo externo de compostos e SMILES · opcional** accepts a companion CSV
with `molecule_chembl_id,canonical_smiles` and matching edge identifiers. This
provides molecular 2D structures, common-fragment search and isolated compounds.
Without the companion table, topology remains available and missing structures
are explained; these CSVs cannot feed stages requiring SMILES until the compound
table is supplied.

Older projects preserve similarity connections and saved experiment results.
Direct fingerprint connections and graph calculation parameters are removed
when opening the project. If no similarity connection remains, connect a
**Calcular similaridade** stage before execution. Older archives receive the
same update on import.

Expand the stage in **Execuções** and select **Análise de grafos · entrada e
origem**. The graph appears inline. Switch between the full graph, including
isolates, and the MCC (largest connected component). Blue outlines highlight
MCC nodes in the full graph. Nodes have a small, uniform size, including in the
exported MCC presentation. Viridis colors indicate degree
(connection count). The color selector also offers normalized connectivity,
defined as degree / (n − 1) for the displayed dataset. The legend gives the
range. This per-node measure differs from global edge density, recorded in the
result model. Components occupy separate areas and vertices retain spacing as
the canvas expands. **Navegar entre componentes** selects a component;
**Ajustar à tela** fits the visible graph. Pan, zoom, recenter or choose circular
organization. Selecting a node highlights its neighbors and relations; the
neighbor list follows connections and displays similarity values. Code searches
can locate compounds in other components. Toggle node
labels as needed. Hover to identify a compound; click to open its molecular 2D
structure and properties.

In MCC mode, the right panel shows the common molecular fragment, its image,
**SMILES do fragmento**, atom and bond counts, and search status. The image and
SMILES depict the matched fragment of a reference molecule. SMARTS remains in
the result model for matching. Older saved graphs receive this display and
layout update when opened. Configure the time limit
(default: 30 seconds) and ring-matching options in the block. The search requires
a match in every MCC compound. When the limit is reached, the interface and
exported figure label the best fragment found as partial; its maximum size is
unconfirmed. This distinction follows the [RDKit FindMCS documentation](https://www.rdkit.org/docs/source/rdkit.Chem.rdFMCS.html).

**Baixar apresentação MCC** exports a PNG with the degree-colored MCC, fragment
image, degree rank, histogram and degree distributions. **Baixar compostos do
MCC** exports the selected analysis CSV with original identifiers and metadata.
**Arquivos desta análise** lists files for the selected result in a paginated table (25 items by
default). Every analysis has its own CSV, figure, interactive model and edges.
`Molecules/molecules.csv` retains the union of MCCs for existing pipelines. If
independent analyses reuse an identifier for different structures, this union
uses suffixed identifiers and preserves `original_code`; individual analyses
keep their original identifiers. Select a specific output in the next block
to work with one analysis.

Equal-sized components are resolved deterministically by compound code. With
no edges, one isolated compound forms an MCC of size one. Networks larger than
1,000 nodes use concentric rings within each component to reduce layout cost, preserving all nodes and
relationships.

For **Avaliar ADMET** (Evaluate ADMET), choose a result CSV and one of four
datasets: all evaluated compounds, BBB+, BBB− or HIA+. The interactive EGG appears
within the stage and shows TPSA × WLOGP and HIA/BBB regions. Axis bounds include
outliers. Clicking a point opens a popup with its molecular 2D structure and
properties. Two buttons on the right download the selected EGG chart or matching
CSV. Filtered CSV downloads preserve compound properties. JSON files remain
internal visualization data and are hidden from the ADMET file list.

**Gerar fingerprints** (Generate fingerprints) offers one algorithm per block:
Morgan, MACCS or pharmacophore. Radius and bit-count fields appear only for Morgan.
**Calcular similaridade** (Compute similarity) automatically identifies and locks
the algorithm for selected application-generated inputs. The selector remains
editable for user-supplied inputs. Different fingerprint types cannot be combined.
For older projects generating several types per block, select a specific output
file before configuring similarity.

In both interactive viewers:

1. Hover over a point to see the compound code.
2. Click to open its 2D structure, SMILES and available properties.
   If several compounds occupy the same point, a list lets you choose one;
   the `(+N)` indicator shows how many other compounds overlap there.
3. Search by code to select a compound with the keyboard, including overlapping
   points. In MCC mode, search is restricted to that component.
4. Drag the background to pan, use zoom and **Recentrar** (Reset view).

Viewers can explore and download charts; editing and execution remain restricted
to editors and owners. Authorization is checked when opening files and compound
details. Versioned `*.biomol-view.json` artifacts are produced alongside PNGs,
without external services or HTML/script execution. Opening is limited to 32 MB;
larger artifacts can be downloaded. Older results retain their PNGs: select
**No, run again** in the reuse prompt, then rerun the stage
to generate interactive artifacts.

## Using your own files

Each block supports **Enviar meus arquivos** (Upload my files), **Adicionar entrada**
(Add input) and a separate completed-results mode, with expected-format validation.
The import block remains available:

1. Open **Arquivos** (Files), choose the type and upload files.
2. Add **Importar meus arquivos** (Import my files) to the pipeline and select the uploads.
3. Configure the import type and, for structures, the target directory.
4. Connect the import output to the stage that will consume the data.

Each file can be up to 200 MB. In the browser, uploads use a signed, temporary URL;
files enter the library only after the session, permission and size are checked
again. The server limits transfers to the configured size, with another check
before incorporating a file into the project. Files are not served from a public
assets directory; downloads go through authorization.

Compound CSVs accept `canonical_smiles` and `molecule_chembl_id`, or `smiles` and
`name`. Import normalizes column names and generates missing identifiers.
Identifiers must be safe for filenames. The uploaded original is preserved.
To select a particular result, use **Resultado usado nesta entrada** (Input result).
The selector lists imported files and outputs from previous runs; automatic
bindings can start automatic, but execution asks you to confirm files. Advanced mode accepts other `selector` values.

For your own PDBs, use the **Complexos PDB** (PDB complexes) import type, followed
by **Preparar meus complexos** (Prepare my complexes), and add records through the
form (PDB, ligand, residue and chain). You can also upload `pdb_codes.csv` with
`PDB_CODE,LIGAND,RESNUM,CHAIN`. PDB retrieval or redocking is then unnecessary.

For already prepared receptors, choose **Receptores preparados** (Prepared
receptors) and upload the bundle required by the existing protocol:
`<PDB>_<CHAIN>.dockprep.pdbqt` files, `centers.csv` and `pdb_codes.csv`.
Import creates `<target>/Prepared` and puts metadata in `<target>/pdb_codes.csv`.
Connect this import directly to Vina, and connect compounds to the other input.
A raw PDB alone still needs molecular preparation.

Fingerprints, similarity relationships and Vina/DOCK6 results can also be
imported. They must follow backend formats: fingerprints contain
`molecule_chembl_id,fingerprint`, relationships contain `source,target,value`,
Vina results include `.lig.pdbqt`, and DOCK6 results include `_scored.mol2`.
Consensus accepts separate inputs from the two engines, including imports.

## Diagnostics and ChEMBL retrieval

The login screen uses original project logos from `imgs`, bundled with the package,
and adapts to smaller screens. Block configuration uses identification, input and
parameter cards. Filters and engine options expand into spaced, responsive forms.

Errors are written under the repository's `logs` directory. Outside a source
checkout, the default is `logs` in the working directory. Use UI option
`--log-dir /path/logs` or environment variable `BIOMOL_LOG_DIR` to override it.

| File | Contents |
| --- | --- |
| `logs/frontend.log` | Interface actions, authentication and event failures |
| `logs/backend.log` | Supervisor/pipeline failures with execution IDs |
| `logs/errors.log` | Aggregated main-process errors with tracebacks |
| `logs/jobs/<job-id>/` | Scientific logs, errors and an `execution.log` copy for each worker |
| `logs/chembl-check.log` | Live retrieval probe diagnostics |
| `logs/chembl-check-report.json` | Latest probe status, target and error/counts |

Diagnostic files rotate at 5 MB, retaining three previous copies. `execution.log`
is copied when the worker exits; the original private log remains available in
**Execuções** while running. Do not serve `logs` as public assets: diagnostics
include paths and project identifiers.

ChEMBL retrieval uses paginated REST endpoints with bounded retries and timeouts,
without initializing `/spore`. An unavailable API remains an explicit error, rather
than an empty dataset. Activity thresholds use `standard_value`/`standard_units`;
natural-product filtering applies to both Yes and No.

The **Recuperar compostos** block calls `wrappers.crawlers.retrieve_compounds` in a
separate Python worker. `workflow/1-InformationRetrieval/retrieve_compounds.py` is
a manual example calling the same function; its options and output directory may
differ from the project's configuration. Worker startup records its interpreter
in `logs/jobs/<job-id>/backend.log`. The private job's `request.json` records the
operation and supplied parameters. Block filters can also come from saved stage
templates.

ChEMBL endpoint HTTP 500 responses or retrieval timeouts are reported to the
frontend with query context after automatic retries. Downloading some files before
a failure does not complete retrieval: partial files remain in that run's private
directory and dependent stages are not executed. The dialog and **Ver log da etapa**
(View stage log) help distinguish this situation from a worker startup failure.

To probe CHEMBL220 using only IC50, without similarity expansion:

```bash
PYTHONPATH=src python scripts/check_chembl.py --target CHEMBL220 --limit 5
```

This probe downloads at most five IC50 activity records in nM, up to 5000 nM and
with pChEMBL, plus their compounds. Successful results are exported under
`datasets/chembl-check/CHEMBL220`. HTTP 500, timeouts and empty results are recorded
in the report and produce exit code 1. Offline tests use controlled responses to
verify filtering and pagination independently of live service availability.

## Advanced configuration and extension

Each stage exposes its wrapper's arguments in guided forms. JSON and template
editors are available in Advanced mode. Inputs and outputs
are managed by the project to preserve its scope. ChEMBL filters and Chimera,
Vina and DOCK6 templates can be edited per stage; the worker receives its own copy.
Keep formatting placeholders intact. Editing accepts scientific configuration
and rejects external commands and unrestricted paths.

The catalog is in `biomolexplorer/catalog.py`; scientific contracts and execution
are in `operations.py`; authentication/files are in `workspace.py`; ordering and
execution are in `pipeline.py`; the interface is in `ui/app.py`, canvas in
`ui/flow_canvas.py`, forms in `ui/guided.py` and visual contracts in `flow.py`. A newly registered
operation can use the same parameter editor, with a title and help added to the catalog.

The current executor uses SQLite and a bounded local queue, with two concurrent
runs and at most 20 queued/active runs. For multiple servers, the interface can
retain these services as a contract while storage and supervision are replaced
with distributed services. This version does not expose a public HTTP API for
scientific services.

For contracts and test evidence, see [pipeline validation](pipeline_validation.md).
