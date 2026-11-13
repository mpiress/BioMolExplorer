# Current source reference

2026-11-13 · BioMolExplorer 2.2.0

Source of truth: the files currently present in this checkout, including uncommitted implementations. This is a documentation audit, not a release certification or a claim that changes have been installed. Historical validation reports describe their dated state.

## Entry points

```bash
biomolexplorer --help
biomolexplorer-ui --help
python -m biomolexplorer.log_report --help
```

`biomolexplorer OPERATION --parameters FILE.json --output DIRECTORY` runs one registered operation. `--parameters` and `--output` are required. JSON stdout describes operation, artifacts and details; scientific progress is redirected to stderr. Failures exit nonzero. `base_output_path` is supplied by `--output`, not an extra JSON key.

### Workspace CLI

| Flag | Default / meaning |
| --- | --- |
| `--data-dir` | `~/.local/share/biomolexplorer` |
| `--web` | False; enable web hosting |
| `--no-browser` | False; suppress browser launch in web mode |
| `--host` | `127.0.0.1` |
| `--port` | `8550` |
| `--log-dir` | `BIOMOL_LOG_DIR` or context-dependent `logs/` |
| `--worker-python` | Current Python when omitted |
| `--dock6-path` | Explicit installation root or resolver discovery |
| `--language` | `en`; choices from packaged `resources/i18n/locales.json` |

The scientific CLI is not a REST server. Web routes expose authenticated report/viewer resources plus the Flet host; internal worker requests are not user-facing CLI operations.

## Active operations and defaults
 
There are 15 registered scientific operations and the workspace-only `import_results` block. Unknown JSON keys are rejected. A callable may retain legacy arguments that the public operation registry does not accept. Required values below are API requirements; typed pipeline bindings can supply inputs at execution. `md_inputs`, `md_seed_path` and `md_cache_epoch` are runtime/cache metadata, not ordinary scientific controls.

### dynamics_dock6

`wrappers.dynamics.perform_dynamics_dock6`

| Parameter | Required | Callable default |
| --- | --- | --- |
| `base_input_path` | yes | `—` |
| `base_selected_mols` | yes | `—` |
| `dock6_app_path` | no | `""` |
| `md_inputs` | no | `null` |
| `md_seed_path` | no | `null` |
| `selected_pairs` | no | `null` |
| `excluded_pairs` | no | `null` |
| `amber_score_type` | no | `"dock6"` |
| `amber_score_movable_region` | no | `"ligand"` |
| `amber_score_minimization_rmsgrad` | no | `0.01` |
| `amber_score_before_md_minimization_cycles` | no | `100` |
| `amber_score_md_steps` | no | `3000` |
| `amber_score_after_md_minimization_cycles` | no | `100` |
| `amber_score_gb_model` | no | `5` |
| `amber_score_nonbonded_cutoff` | no | `18.0` |
| `amber_score_temperature` | no | `300.0` |
| `amber_score_movable_distance_cutoff` | no | `10.0` |
| `amber_score_receptor_movable_atom_expr` | no | `"::"` |
| `amber_score_complex_movable_atom_expr` | no | `"::"` |
| `amber_score_abort_on_unprepped_ligand` | no | `true` |
| `use_existing_ligand_charges` | no | `false` |
| `residue_parameter_policy` | no | `"auto_curated"` |
| `md_cache_epoch` | no | `""` |

### retrieve_pubchem

`biomolexplorer.pubchem_retrieval.retrieve_pubchem`

| Parameter | Required | Callable default |
| --- | --- | --- |
| `target` | no | `"PubChem"` |
| `reference_source` | no | `"compounds"` |
| `reference_type` | no | `"smiles"` |
| `reference` | no | `""` |
| `selection_mode` | no | `"all"` |
| `compound_id` | no | `""` |
| `threshold` | no | `75` |
| `max_records` | no | `1000` |
| `base_input_path` | no | `null` |
| `input_file` | no | `null` |

### retrieve_compounds

`wrappers.crawlers.retrieve_compounds`

| Parameter | Required | Callable default |
| --- | --- | --- |
| `search_term` | yes | `—` |
| `include_pubchem` | no | `false` |
| `pubchem_threshold` | no | `75` |
| `pubchem_max_records` | no | `1000` |
| `chembl_filters` | no | `null` |
| `search_mode` | no | `"target"` |
| `max_targets` | no | `25` |
| `max_records` | no | `1000` |
| `expand_chembl` | no | `true` |
| `similarity_threshold` | no | `70` |

### expand_similar_compounds

`wrappers.crawlers.expand_similar_compounds`

| Parameter | Required | Callable default |
| --- | --- | --- |
| `search_term` | yes | `—` |
| `threshold` | no | `75` |
| `max_records` | no | `1000` |
| `base_input_path` | no | `null` |
| `input_file` | no | `null` |

### retrieve_structures

`wrappers.crawlers.load_pdb`

| Parameter | Required | Callable default |
| --- | --- | --- |
| `target` | no | `"Estruturas"` |
| `pdb_query` | no | `null` |
| `pdb_ids` | no | `null` |
| `uniprot_ids` | no | `null` |
| `ligand_ids` | no | `null` |
| `max_records` | no | `100` |
| `pdb_ec` | no | `null` |
| `organism` | no | `null` |
| `PolymerEntityTypeID` | no | `null` |
| `ExperimentalMethodID` | no | `null` |
| `max_resolution` | no | `null` |
| `must_have_ligand` | no | `true` |
| `ligand_selection` | no | `"recommended"` |
| `excluded_pdb_ids` | no | `null` |

### retrieve_zinc

`wrappers.crawlers.load_zinc`

| Parameter | Required | Callable default |
| --- | --- | --- |
| `base_input_path` | yes | `null` |
| `filename` | no | `"zinc_urls.txt"` |
| `verbose` | no | `0` |
| `download_workers` | no | `4` |

### prepare_structures

`wrappers.redocking.prepare_structures`

| Parameter | Required | Callable default |
| --- | --- | --- |
| `base_input_path` | yes | `—` |
| `target` | yes | `—` |
| `pdb_codes` | no | `null` |
| `pH` | no | `7.4` |
| `charge_type` | no | `"gas"` |
| `preparation_options` | no | `null` |
| `base_selected_mols` | no | `null` |
| `mol_filename` | no | `"compounds"` |
| `receptor_prepared` | no | `false` |
| `docking_engines` | no | `"both"` |

### admet

`wrappers.admet.ADMETWrapper`

| Parameter | Required | Callable default |
| --- | --- | --- |
| `base_input_path` | yes | `null` |
| `input_file` | no | `null` |
| `verbose` | no | `0` |
| `admet_profile` | no | `"legacy"` |
| `admet_mode` | no | `null` |
| `prediction_backend` | no | `"rdkit"` |
| `bbb_target` | no | `null` |
| `hia_target` | no | `null` |
| `rule_set` | no | `null` |
| `max_rule_violations` | no | `null` |
| `alert_policy` | no | `null` |
| `alert_catalogs` | no | `null` |
| `property_limits` | no | `null` |
| `dataset_origin` | no | `null` |

### fingerprints

`wrappers.molecular_analyzer.generate_fingerprints`

| Parameter | Required | Callable default |
| --- | --- | --- |
| `base_input_path` | yes | `—` |
| `morgan_n_bits` | no | `2048` |
| `radius` | no | `2` |
| `files` | no | `null` |
| `morgan` | no | `true` |
| `maccs` | no | `true` |
| `pharmacophore` | no | `true` |
| `chunk_size` | no | `100` |

### similarity

`wrappers.molecular_analyzer.compute_similarity`

| Parameter | Required | Callable default |
| --- | --- | --- |
| `base_input_path` | yes | `—` |
| `metric` | no | `"similarityFunctions.TanimotoSimilarity"` |
| `fingerprint` | no | `"fingerprints.Morgan"` |
| `filename` | no | `null` |
| `threshold` | no | `null` |
| `approximate` | no | `true` |

### graphs

`wrappers.molecular_analyzer.analyze_graphs`

| Parameter | Required | Callable default |
| --- | --- | --- |
| `similarity_path` | no | `null` |
| `base_input_path` | no | `null` |
| `mcs_timeout` | no | `0` |
| `mcs_ring_matches_ring_only` | no | `true` |
| `mcs_complete_rings_only` | no | `false` |
| `graph_inputs` | no | `null` |

### redocking

`wrappers.redocking.perform_redocking`

| Parameter | Required | Callable default |
| --- | --- | --- |
| `base_input_path` | yes | `—` |
| `target` | yes | `—` |
| `pdb_codes` | no | `null` |
| `pH` | no | `7.4` |
| `sizeof_box` | no | `[24, 24, 24]` |
| `exhaustiveness` | no | `20` |
| `num_modes` | no | `10` |
| `prepare_complex` | no | `true` |
| `charge_type` | no | `"gas"` |
| `preparation_pairs` | no | `null` |
| `preparation_options` | no | `null` |
| `redocking_engine` | no | `"vina"` |
| `dock6_app_path` | no | `null` |
| `density` | no | `0.5` |
| `radius` | no | `1.4` |
| `distance` | no | `10.0` |
| `conformer_search_type` | no | `"flex"` |
| `max_orientations` | no | `1000` |
| `num_scored_conformers` | no | `1` |
| `dock6_seed` | no | `73` |
| `dock6_mode` | no | `"redock"` |

### docking_vina

`wrappers.docking.perform_docking_vina`

| Parameter | Required | Callable default |
| --- | --- | --- |
| `base_input_path` | yes | `—` |
| `target` | yes | `—` |
| `base_selected_mols` | yes | `—` |
| `mol_filename` | yes | `—` |
| `pdb_code` | no | `null` |
| `pH` | no | `7.4` |
| `sizeof_box` | no | `[24, 24, 24]` |
| `exhaustiveness` | no | `20` |
| `num_modes` | no | `10` |
| `preparation_options` | no | `null` |
| `receptor_prepared` | no | `true` |

### docking_dock6

`wrappers.docking.perform_docking_dock6`

| Parameter | Required | Callable default |
| --- | --- | --- |
| `base_input_path` | yes | `—` |
| `target` | yes | `—` |
| `base_selected_mols` | yes | `—` |
| `dock6_app_path` | yes | `—` |
| `charge_type` | yes | `—` |
| `mol_filename` | yes | `—` |
| `pdb_code` | no | `null` |
| `base_vina_path` | no | `null` |
| `density` | no | `0.5` |
| `radius` | no | `1.4` |
| `distance` | no | `10.0` |
| `conformer_search_type` | no | `"flex"` |
| `plot_max_residues` | no | `50` |
| `pH` | no | `7.4` |
| `preparation_options` | no | `null` |
| `receptor_prepared` | no | `true` |
| `dock6_mode` | no | `"auto"` |

### consensus

`wrappers.docking.generate_consensus`

| Parameter | Required | Callable default |
| --- | --- | --- |
| `base_input_path` | yes | `—` |
| `target` | yes | `—` |
| `repulsion_weight` | no | `1.0` |
| `base_vina_path` | no | `null` |
| `base_dock6_path` | no | `null` |

## Runtime configuration and tools
 
`AppConfig` controls `workspace`, `max_jobs=1`, `cpu_workers=2`, `max_queued_jobs=100`, `job_timeout=None`, `worker_python=None`, `resource_dir=None`, `cache_dir=None`, `reuse_results=True`, `preparation_cache=None`, `reuse_preparations=True`. These are Python API options, not extra workspace CLI flags.

Python >=3.12 is required by packaging. `requirements.yml` pins Python 3.12.3 and the scientific/UI stack. `pyproject.toml` deliberately declares no base scientific dependencies; the `ui` extra adds Flet desktop/web. External Chimera and DOCK6 are not installed by pip. Compiler, make, flex and zlib dependencies support execution-local AMBER compatibility builds. The optional ADMET-AI setup is separate.

### Environment names used by source

| Name | Implementation references |
| --- | --- |
| `BIOMOL_ADMET_PYTHON` | `biomolexplorer/admet_predictors.py` |
| `BIOMOL_AMBER_PTM_PROFILES` | `biomolexplorer/md_residue_support.py` |
| `BIOMOL_CACHE_DIR` | `biomolexplorer/jobs.py`, `biomolexplorer/ligand_mol2.py`, `biomolexplorer/molecular_cache.py`, `biomolexplorer/preparation_cache.py`, `crawlers/pubchem.py` |
| `BIOMOL_CCD_PATH` | `biomolexplorer/chimera_charge_repair.py`, `biomolexplorer/structure_selection.py` |
| `BIOMOL_CPU_WORKERS` | `biomolexplorer/jobs.py`, `biomolexplorer/paths.py` |
| `BIOMOL_DOCK6_ROOT` | `biomolexplorer/docking_tools.py` |
| `BIOMOL_LOG_CONTEXT` | `biomolexplorer/diagnostics.py`, `biomolexplorer/jobs.py` |
| `BIOMOL_LOG_DIR` | `biomolexplorer/diagnostics.py`, `biomolexplorer/jobs.py`, `biomolexplorer/ui/app.py` |
| `BIOMOL_MD_DEBUG` | `biomolexplorer/md_engine.py` |
| `BIOMOL_MD_DOCK_BINARY` | `biomolexplorer/md_engine.py` |
| `BIOMOL_PREPARATION_CACHE` | `biomolexplorer/calculation_cache.py`, `biomolexplorer/chimera_charge_repair.py`, `biomolexplorer/jobs.py`, `biomolexplorer/ligand_cache.py` |
| `BIOMOL_PREPARATION_DEBUG` | `biomolexplorer/molecular_core.py` |
| `BIOMOL_PREPARATION_WORKERS` | `biomolexplorer/preparation_resources.py` |
| `BIOMOL_PROGRESS_FILE` | `biomolexplorer/jobs.py`, `biomolexplorer/progress.py` |
| `BIOMOL_RESOURCE_DIR` | `biomolexplorer/jobs.py`, `biomolexplorer/paths.py` |
| `BIOMOL_REUSE_PREPARATIONS` | `biomolexplorer/calculation_cache.py`, `biomolexplorer/jobs.py`, `biomolexplorer/ligand_cache.py` |
| `BIOMOL_REUSE_RESULTS` | `biomolexplorer/jobs.py`, `biomolexplorer/molecular_cache.py`, `biomolexplorer/preparation_cache.py` |
| `BIOMOL_WORKSPACE` | `biomolexplorer/jobs.py`, `biomolexplorer/paths.py` |

Variables in this inventory include internal worker-propagated values. Their presence does not make each one a supported user control; consult the indicated implementation and the preparation/installation guides.

## Supported formats and persistence
 
Inputs depend on semantic ports: compound CSV/SMILES and compatible molecular structures, receptor PDB and prepared PDBQT/MOL2 bundles, fingerprint binaries, similarity/graph CSV and authorized result manifests. Native Amber Score requires target PDB and a target-matched coordinate pose, not bare SMILES. Engine-native files, compatibility CSVs and images/reports remain durable; SQLite indexes support lookup, lineage, table pagination and scientific caches and are not a replacement for all structural artifacts.

## Compatibility and removed behaviors
 
- No standalone AmberTools backend is selectable for docking preparation; `auto`/`python` preparation identifiers normalize to the Chimera path.
- Native Amber Score is not sander/cpptraj MD; conventional trajectories, RMSF, equilibration controls and arbitrary metal parameterization are not provided.
- Consensus z-score/min-max generation and the old clipped-energy consensus formula are retired; legacy readers remain compatible.
- Automatic updates install a validated package; manual ZIP extraction is not the update button workflow. A restart remains manual.
- There is no public relational-output analysis operation in `OPERATIONS`; SQLite implementation files/indexes must not be presented as a new scientific block.
- Deprecated timing/trajectory arguments retained in the Python callable are not accepted JSON parameters unless listed in the registry above.

## Automation scripts

Flags below are extracted from each script parser. Scripts with no parser are Python utilities with implementation-defined entry points, not extra scientific operations. Some scripts start engines, access providers, build derived indexes or drive browser sessions; inspect their options and use a separate output directory before running them. They were not executed merely to validate this documentation.

| Script | Flags declared in source |
| --- | --- |
| `scripts/audit_canvas_ui.py` | `--port` |
| `scripts/audit_dependencies.py` | `--output` |
| `scripts/audit_docking_reuse.py` | `--database`, `--run-id`, `--output`, `--consolidate-caches` |
| `scripts/audit_localization.py` | `--output` |
| `scripts/audit_settings_ui.py` | `--port` |
| `scripts/benchmark_canvas_scaling.py` | `--baseline`, `--frames`, `--output` |
| `scripts/benchmark_editor_interactions.py` | `--baseline`, `--repeats` |
| `scripts/benchmark_molecular_preparation.py` | `--sizes`, `--legacy-max`, `--output`, `--engines` |
| `scripts/benchmark_pipeline_editor.py` | `--baseline`, `--frames` |
| `scripts/benchmark_pipeline_save.py` | `--repeats`, `--files` |
| `scripts/benchmark_report.py` | No argparse flags |
| `scripts/benchmark_result_index.py` | `--rows`, `--database`, `--run-id`, `--output` |
| `scripts/build_docs.py` | No argparse flags |
| `scripts/check_chembl.py` | `--target`, `--limit`, `--output` |
| `scripts/compare_antechamber_migration.py` | `--dock6`, `--output`, `--charge-only` |
| `scripts/prepare_result_indexes.py` | `--database`, `--project-id`, `--run-id`, `--limit`, `--report` |
| `scripts/setup_admet_ai.py` | `--directory` |
| `scripts/validate_amber_support.py` | `--dock6-root`, `--output`, `--pairs-json` |
| `scripts/validate_canvas_browser.py` | `--port`, `--output` |
| `scripts/validate_dock6_modes.py` | `--prepared-input`, `--target`, `--poses-table`, `--receptor-id`, `--compound-id`, `--dock6-root`, `--output`, `--cache-dir`, `--require-existing-grid` |
| `scripts/validate_docking.py` | `--dock6-root`, `--output`, `--prepare-inputs` |
| `scripts/validate_docking_viewer.py` | `--chrome`, `--pdb`, `--structure`, `--redocking-root`, `--docking-root`, `--smiles`, `--screenshot` |
| `scripts/validate_docking_warm_cache.py` | `--cache-root`, `--output` |
| `scripts/validate_docs.py` | No argparse flags |
| `scripts/validate_native_structures.py` | `--dock6`, `--output`, `--steps` |
| `scripts/validate_native_surface.py` | `--dms-executable` |
| `scripts/validate_pdb_viewer.py` | `--chrome`, `--pdb`, `--structure`, `--smiles`, `--screenshot` |
| `scripts/validate_preparation_migration.py` | `--dock6`, `--output` |
| `scripts/validate_settings_browser.py` | `--output` |
