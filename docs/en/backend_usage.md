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
| `retrieve_zinc` | `filename`, `base_input_path` containing the URI file | ZINC CSV |
| `admet` | `base_input_path`, optional `input_file` | Assessment CSVs, subsets and plot |
| `fingerprints` | `base_input_path`, algorithms, `chunk_size` | Fingerprint CSVs |
| `similarity` | Fingerprint `base_input_path`, metric, threshold and `approximate` | Edge CSVs in the `Similarity` subdirectory |
| `graphs` | `similarity_path`; optional compound `base_input_path`; MCS options | Independent graph/MCC models, fragments and CSVs |
| `prepare_structures` | User-supplied PDBs, `target`, PDB records, pH and charges | Prepared structures, centers and `pdb_codes.csv` |
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
