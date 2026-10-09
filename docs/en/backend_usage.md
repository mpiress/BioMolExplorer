# Application layer usage

[Documentation](../README.md) · English · [Português](../backend_usage.md)

For the ready-to-use interface with authentication and projects, see the
[Flet workspace guide](frontend.md). This document covers the CLI and the Python
services used by the interface.

Before running the CLI or starting the interface, follow the [installation and configuration guide](installation.md): download the [GitHub repository](https://github.com/mpiress/BioMolExplorer) and install **UCSF Chimera 1.17 and DOCK6 6.11** on the computer running calculations. Configure their executables on `PATH`. Stages requiring these tools will fail without them.

Create or activate the Conda environment and install the project from the downloaded source code root:

```bash
conda env create -f requirements.yml
conda activate BioMolExplorer
python -m pip install -e . --no-deps --no-build-isolation
```

Scientific dependencies are managed by the Conda files. The Python package
provides the application layer, scientific modules and resources; installing
the package alone does not install Chimera 1.17 or DOCK6 6.11.
Open Babel and Vina are declared in the Conda environment. Update an existing
environment as needed and check those executables too.

## CLI

```bash
biomolexplorer retrieve_compounds \
  --parameters examples/retrieve_compounds.json \
  --output ./datasets
```

The equivalent command is `python -m biomolexplorer ...`. The CLI runs the
operation to completion and returns JSON containing the artifacts. Without
installing the package, run `PYTHONPATH=src python -m biomolexplorer ...` from the
project root. Scripts in `workflow/` remain available and locate `src` relative
to their own files. `/datasets` retains its meaning as a workspace-relative
directory; `/tmp/study`, for example, is a normal absolute path.

ChEMBL filters are stored in `src/biomolexplorer/resources/crawlers`.
You can also supply a `chembl_filters` dictionary in the parameters, with
`target`, `bioactivity`, `molecules` and `similars` keys. Each supplied key replaces
the filter set for that stage; omitted stages use their defaults. This supports
configuring UI runs without changing global files.

## Jobs for Flet

```python
from pathlib import Path
from biomolexplorer.config import AppConfig
from biomolexplorer.jobs import JobManager

config = AppConfig(
    workspace=Path('/tmp/my-study'),
    max_jobs=1,
    cpu_workers=2,
    max_queued_jobs=100,
    job_timeout=86400,
)
manager = JobManager(config)
job = manager.submit('retrieve_compounds', {
    'search_term': 'CHEMBL220',
    'include_pubchem': True,
    'pubchem_threshold': 75,
    'pubchem_max_records': 1000,
})
# This call returns without waiting for scientific retrieval to finish.
status = manager.get(job['id'])
# manager.cancel(job['id'])
# When shutting down the application: manager.close()
```

Use one long-lived `JobManager` instance per workspace. A `with` block closes
the manager on exit and cancels pending jobs; do not create a new manager for
every UI click.

`examples/flet_controller.py` provides a controller independent of widgets.
Its methods use `asyncio.to_thread` for management calls, and `updates` yields
snapshots for updating controls in the UI loop:

```python
job = await controller.submit('admet', {
    'base_input_path': '/tmp/my-study/datasets/compounds/CHEMBL220',
    'input_file': 'compounds.csv',
})
async for snapshot in controller.updates(job['id']):
    # Update Flet controls here; snapshots are serializable dictionaries.
    print(snapshot['status'], snapshot['error'])
```

Terminal states are `succeeded`, `failed`, `cancelled` and `interrupted`.
`result.artifacts` provides file paths; `log_path` points to the worker log.
Execution states are available, but stage progress percentages are not yet
provided. The interface must handle validation/submission errors and the `failed`
state of an already accepted job.

If the interface uses another Python environment, configure
`AppConfig(worker_python=Path('/path/conda/envs/BioMolExplorer/bin/python'), ...)`.
Scientific processes use that interpreter. Workers do not require Flet
dependencies. The controller does not build a visual interface.

## Available operations

See `biomolexplorer.operations.OPERATIONS` for the complete parameters.

| Operation | Main inputs | Outputs |
| --- | --- | --- |
| `retrieve_compounds` | `search_term`, PubChem options and ChEMBL filters | Source-specific data and `compounds/<target>/compounds.csv` |
| `expand_similar_compounds` | `search_term`, `base_input_path` containing `ChEMBL/molecules` and `ChEMBL/similars` | New compounds, relationships and consolidated dataset |
| `retrieve_structures` | Text/IDs/UniProt/ligands or attribute filters; optional `target` | PDB files, `pdb_codes.csv`, retrieval report |
| `retrieve_zinc` | `base_input_path`, tranche-list `filename` | `compounds.csv`, individual MOL2 files and report |
| `admet` | `base_input_path`, optional `input_file` | Assessment CSVs, subsets and plot |
| `fingerprints` | `base_input_path`, algorithms, `chunk_size` | Fingerprint CSVs |
| `similarity` | Fingerprint `base_input_path`, metric, threshold and `approximate` | Edge CSVs in the `Similarity` subdirectory |
| `graphs` | `similarity_path`; optional compound `base_input_path`; MCS options | Independent graph/MCC models, fragments and CSVs |
| `prepare_structures` | Raw or prepared PDB receptor, compounds, `target`, PDB records and role options | Receptors, centers, `pdb_codes.csv` and candidates in selected formats |
| `redocking` | PDB directory, `target`, preparation parameters | Working copies, Vina results and RMSD |
| `docking_vina` | Prepared complexes, selected compounds, `mol_filename` | PDBQT files and Vina results |
| `docking_dock6` | Prepared complexes, compounds, `base_vina_path`, `pdb_code`, DOCK6 installation | Conformations, scores and footprints |
| `consensus` | Directory containing `Vina` and `Dock6`, or separate `base_vina_path` and `base_dock6_path` inputs, `target` | Consensus CSV and plot |

Enums accept names or values, such as `TanimotoSimilarity`/`Tanimoto` and
`Morgan`/`morgan`. PDB filter lists also use serializable strings.
The Vina `pdb_code` parameter follows the existing function: a list of
`[PDB_CODE, LIGAND, RESNUM, CHAIN]` records; DOCK6 receives one record with those
four fields. Redocking requires explicitly selected `pdb_codes`: a list of
`[PDB_CODE, LIGAND, RESNUM, CHAIN]`, with optional resolution as its fifth value,
 supplemented by `pdb_codes.csv` metadata. `preparation_pairs` stores options
under `PDB|LIGAND|RESNUM|CHAIN` keys; receptor and ligand use the same chain.
Legacy `verbose` arguments are ignored; verbosity stays at zero.
See [redocking configuration](redocking_configuration.md).

Each job receives a new output directory. Chain jobs by passing their artifacts
to the next operation. `graphs` accepts only ready similarity CSVs or directories
through `similarity_path`. Each file produces an independent analysis.
`base_input_path` supplies an optional companion compound CSV or directory.
Edge files use `source,target,value`, with weights between 0 and 1, and need no
metric-based filename prefix. Compound tables with
`molecule_chembl_id,canonical_smiles` preserve isolated nodes and enable structures
and common-fragment search. Without them, only graph topology is available.

Metric, fingerprint and threshold belong to `similarity`. `graphs` does not
recalculate or filter weights. `mcs_timeout` defaults to 30 seconds,
`mcs_ring_matches_ring_only` to true and `mcs_complete_rings_only` to false.
See the [workspace guide](frontend.md) for selection, interaction and downloads.

Advanced callers may pass `graph_inputs`: dictionaries with `kind: "similarity"`,
`file`, `compound_files` (an empty list when structures are absent), and optional
`label`, `metric` and `fingerprint` provenance metadata. These fields describe
inputs rather than configure computation. `fingerprints` entries are rejected.
Workspace paths must belong to the project. Compound SMILES are resolved through
upstream stages, and each analysis remains separate under `plots/`, `Molecules/`,
`data/maxcomp/` and `centroids/`. MCS queries remain internal SMARTS; the model and interface also provide
fragment SMILES extracted from the reference molecule and search status.

## Redocking example

The example requires an existing `/tmp/study/PDB/Estruturas/4M0E.pdb` containing residue 1YL / 604 in chain A. Adapt paths and the pair to your collection; this call does not retrieve the PDB.

```python
from biomolexplorer.operations import execute_operation

parameters = {
    "base_input_path": "/tmp/study/PDB",
    "target": "Estruturas",
    "pdb_codes": [["4M0E", "1YL", 604, "A", 2.0]],
    "prepare_complex": True,
    "pH": 7.4,
    "sizeof_box": [24, 24, 24],
    "exhaustiveness": 20,
    "num_modes": 10,
    "charge_type": "gas",
    "preparation_pairs": {
        "4M0E|1YL|604|A": {
            "cofactors": [],
            "receptor": {
                "remove_solvent": True,
                "remove_hydrogens": True,
                "add_hydrogens": True,
                "minimize": True,
                "charge_type": "gas",
            },
            "ligand": {
                "remove_solvent": True,
                "remove_hydrogens": True,
                "add_hydrogens": True,
                "minimize": True,
                "charge_type": "gas",
            },
        },
    },
}
result = execute_operation("redocking", parameters, "/tmp/study/redocking-output")
```

When pair options are omitted, solvent removal, existing hydrogen removal, hydrogen addition and minimization default to `True`; receptor charges use `gas`, while omitted ligand charges inherit stage-wide `charge_type`. No cofactor is selected by default. With `prepare_complex=False`, retain receptors `<PDB>_<CHAIN>.dockprep.pdbqt`, reference ligands `<PDB>_<LIGAND>_<RESNUM><CHAIN>.lig.pdbqt` and three finite coordinates per complex in `Prepared/centers.csv`.

`execute_operation` prepares a copy in `structures/<target>/`, preserving the original input. Metadata with `RMSD` is saved in `structures/<target>/pdb_codes.csv`, prepared files in `structures/<target>/Prepared/`, and poses in `<target>/` within the output directory. `result.artifacts` lists exportable files. Direct calls to `wrappers.redocking.perform_redocking` do not provide the application layer's isolation. The CLI accepts the same dictionary through `--parameters`, without the interface's inspection dialogs.

## Tests

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

Run tests in the scientific environment. See the [architecture review](architecture.md)
for corrected behavior, verification and current limitations.

## Portable projects and multiple inputs

`WorkspaceStore.create_project(..., directory="/new/folder")` selects the project
file directory. `project_dir(id)` resolves both existing projects and selected
folders. `export_project(token,id)` and `import_project(token,archive,directory)`
transfer configuration, results and versions with integrity verification; accounts
and sessions are excluded. `history(token,id)` lists changes, while
`rollback(token,id,event_id)` restores the state preceding an event, requiring
ownership and no active execution.

Bindings accept a single reference or `{"sources":[...references...]}`. Each
reference uses `stage` or `asset`, with an optional `selector`. Selected CSVs are
normalized into working copies. `input_processing="individual"` separates files;
`input_processing="merge"` combines selected inputs. `provided_results` contains `kind`
and `asset_ids` to bypass calculation with validated results; completed ADMET
requires identifiers, SMILES, TPSA and WLOGP. See [projects and versions](projects.md).


## Selection, resume and interface language

`PipelineService.submit(token, project_id, reuse_results=True)` reuses compatible results. The interface checks `existing_results` and asks for the user's decision; `reuse_results=False` forces recalculation. Stages needing data pause as `awaiting_input`; `resume(token, run_id, configuration)` confirms explicit references and mode. Completed results persist across restarts. The CLI runs one operation and does not offer workspace popups.

`biomolexplorer-ui --language en` (default) starts in English; `--language pt` starts in Portuguese. Login allows session-specific switching. `ui/localization.py` applies catalogs in `resources/i18n/` only to presentation, preserving parameters, identifiers and editable data. See the [manual](user_manual.md) and [validation](pipeline_validation.md) guides.

## Flexible retrieval parameters

See the [retrieval guide](retrieval.md) for all ChEMBL `search_mode` values, direct compound search examples, optional PDB collection names, limits, filters and report files. Existing operation names and downstream CSV contracts remain available.

See [Logs and diagnostics](logging.md) for the common format, execution context, failure codes, job summary and `python -m biomolexplorer.log_report` command.


## Interoperable docking contract

**Docking with Vina** and **Docking with DOCK6** accept compounds from any block publishing molecules: ChEMBL, PubChem, ZINC, ADMET, fingerprints, graphs, consensus or user imports. Select a table at the compound input and connect the prepared receptor separately. Graphs, ADMET and the other engine are optional. PubMed provides literature references; molecules obtained from those references must be imported with identifiers and SMILES or valid structures.

Reuse conformations by connecting **Vina → DOCK6** or **DOCK6 → Vina** at the compound input and selecting `docking_results.csv` or a specific pose. Identity and SMILES are retained; conversion uses the selected conformation instead of generating one from SMILES. Without an input pose, preparation generates a 3D conformer. DOCK6 uses the receptor site center to select spheres and position newly generated ligands; its optional Vina pose input remains available for refining corresponding candidates.

Review Vina box, pH and search effort and DOCK6 charges, surface, distance, radius, flexible/rigid search and footprint. Both require prepared receptors, metadata and centers; DOCK6 also needs receptor MOL2 and hydrogen-free PDB. Results show compound ID, receptor, score, SMILES, **3D** of the calculated pose and **Remove**. This viewer opens the structure produced by docking.

**Docking consensus** gathers selected results from both engines, including multiple batches and branches. It computes only the intersection by compound ID and receptor, choosing the lowest score for repeated poses. Matching identifiers with different structures are rejected. With no intersection, the block is skipped and explains that the engines did not evaluate common compounds against the same receptor.

The consensus table shows Vina and DOCK6 scores, SMILES, receptor and z-score/min-max normalizations, with **3D Vina**, **3D DOCK6** and **Remove** buttons per row. The DOCK6 consensus score is `min(0, Grid_Score + repulsion_weight × Internal_energy_repulsive)`; missing repulsion is zero. Singleton or constant scores normalize to zero. Removal requires editing permission, is audited and affects subsequent inputs; completed analyses are not recalculated. Each batch summary is authoritative, preventing deleted rows from returning from auxiliary tables.

Native outputs use `docking_results.csv`: `molecule_chembl_id,canonical_smiles,receptor_id,engine,score,conformer_file`. Pose paths are relative to the table. Import the conformation files too. For direct backend calls, `base_selected_mols` is the table directory and `mol_filename` its name without `.csv`; DOCK6 `base_vina_path` is optional. Consensus accepts directories through `base_vina_path` and `base_dock6_path`, exports both poses into `poses/` and returns `skipped_reason` for an empty intersection.


`prepare_structures` accepts `preparation_options` containing `receptor`, `ligand` and `cofactors`. Each role accepts boolean `remove_solvent`, `remove_hydrogens`, `add_hydrogens`, `minimize` and a `charge_type` of `gas` or `am1`. The wrapper applies these settings to explicit records or records inferred from `pdb_codes.csv`. pH remains shared across roles. Legacy calls without this parameter retain previous preparation behavior; the UI uses role settings and migrates earlier configurations while retaining custom templates.


`prepare_structures` also accepts `base_selected_mols`, `mol_filename` (default `compounds`), `receptor_prepared` and `docking_engines` (`vina`, `dock6` or `both`). The pipeline detects prepared redocking receptors automatically. Prepared receptors are copied with companion files and binding centers without repeated preparation; raw PDBs use the redocking preparation procedure. External candidates are prepared under `Target/Compounds/compounds.csv`, preserving identifiers and SMILES with `prepared_pdbqt` and/or `prepared_mol2` columns. Vina and DOCK6 reuse these files without further minimization. The manifest records available formats; an engine cannot use an output that excludes its format. The CSV and all referenced files must travel together. Legacy calls without `base_selected_mols` still prepare structures only.

Use a separate block for each receptor mode: do not combine raw PDBs and prepared receptors in the same block. pH remains available for candidate preparation when reusing a receptor. Before preparing compounds, the block checks binding centers (three finite coordinates) and receptor files required by the selected output. Vina requires `.dockprep.pdbqt`; DOCK6 additionally requires `.dockprep.mol2` and `.noH.pdb`. PDBQT remains the receptor selection file for either engine. Missing files stop the process with the required filename, without repeating receptor preparation.

## Per-file imports and ZINC tranches

`import_results` accepts `asset_types`, a mapping `{asset_id: type}` covering exactly the selected `asset_ids`. A block can combine compound CSVs, structures and ZINC lists, publishing their respective types. Without `asset_types`, legacy `kind` still applies to every file. The pipeline revalidates files and bundles before execution; prepared receptors still require metadata. Removing a file from the form table changes the block's selection only.

`download_workers` configures parallel downloads in `retrieve_zinc` and `wrappers.crawlers.load_zinc`: integer from 1 to 16, default 4; 1 runs sequentially. Each task owns its HTTP session. Prefetch is bounded by the thread count; molecular normalization follows list order. The setting is validated before downloads and the report records the effective `download_workers`, limited to the number of links.

`retrieve_zinc` calls `wrappers.crawlers.load_zinc`. `base_input_path` identifies the list folder and `filename` its name (default `zinc_urls.txt`); the interface resolves both from selected files. The output has the stable name `compounds.csv` regardless of the list name. Individual MOL2 files reside in `Conformers/`, with authorized paths in the manifest. `conformer_origin=library` identifies conformations still requiring binding-site placement; the preparation block preserves this origin as `prepared_origin`.

```python
from biomolexplorer.operations import execute_operation

execute_operation('retrieve_zinc', {
    'base_input_path': '/path/to/lists',
    'filename': 'zinc-download.uri',
}, '/path/to/output')
```

See [2D and 3D tranches](retrieval.md#zinc-2d-and-3d-tranches) for formats and artifact organization.

## Specific compounds and preconfigured docking inputs

Input references under `base_selected_mols` in `docking_vina` and `docking_dock6` accept optional `compound_id` alongside `stage`/`asset` and `selector`. Without this key, every compound in the file participates. Filtering occurs per reference before deduplication and preserves conformations and prepared files. A missing code stops the stage with a validation message. The selection is part of configuration, caching and materialized input provenance.

```json
{"stage": "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", "selector": "compounds.csv", "compound_id": "CHEMBL1"}
```

`requires_curation` skips repeat confirmation for Vina/DOCK6 when receptor and compounds are configured and every reference identifies files (`asset` or an explicit `selector`). Stage references with `selector=auto` still require selection. Configuration is revalidated before execution. Incompatible connections use the same standard message in the canvas and `validate_pipeline`.
