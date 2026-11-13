# Referência do código atual

2026-11-13 · BioMolExplorer 2.2.0

Fonte da verdade: os arquivos presentes neste checkout, incluindo implementações ainda não commitadas. Esta é uma auditoria documental, não uma certificação de release nem uma afirmação de instalação das alterações. Relatórios históricos descrevem o estado na data indicada.

## Pontos de entrada

```bash
biomolexplorer --help
biomolexplorer-ui --help
python -m biomolexplorer.log_report --help
```

`biomolexplorer OPERAÇÃO --parameters ARQUIVO.json --output DIRETÓRIO` executa uma operação registrada. `--parameters` e `--output` são obrigatórios. O JSON em stdout descreve operação, artefatos e detalhes; o progresso científico vai para stderr. Falhas retornam código diferente de zero. `base_output_path` vem de `--output`, não de uma chave adicional no JSON.

### CLI do workspace

| Flag | Padrão / significado |
| --- | --- |
| `--data-dir` | `~/.local/share/biomolexplorer` |
| `--web` | False; habilita hospedagem web |
| `--no-browser` | False; não abre navegador no modo web |
| `--host` | `127.0.0.1` |
| `--port` | `8550` |
| `--log-dir` | `BIOMOL_LOG_DIR` ou `logs/` conforme o contexto |
| `--worker-python` | Python atual quando omitido |
| `--dock6-path` | Raiz explícita da instalação ou descoberta pelo resolvedor |
| `--language` | `en`; escolhas em `resources/i18n/locales.json` |

A CLI científica não é um servidor REST. As rotas web oferecem recursos autorizados de relatório/visualização e o host Flet; requisições internas dos workers não são operações públicas da CLI.

## Operações ativas e valores padrão
 
Há 15 operações científicas registradas e o bloco `import_results`, exclusivo do workspace. Chaves JSON desconhecidas são recusadas. Uma função pode conservar argumentos legados que o registro público não aceita. Campos obrigatórios abaixo são exigências da API; vínculos tipados podem fornecê-los na execução. `md_inputs`, `md_seed_path` e `md_cache_epoch` são metadados de execução/cache, não controles científicos comuns.

### dynamics_dock6

`wrappers.dynamics.perform_dynamics_dock6`

| Parâmetro | Obrigatório | Padrão da função |
| --- | --- | --- |
| `base_input_path` | sim | `—` |
| `base_selected_mols` | sim | `—` |
| `dock6_app_path` | não | `""` |
| `md_inputs` | não | `null` |
| `md_seed_path` | não | `null` |
| `selected_pairs` | não | `null` |
| `excluded_pairs` | não | `null` |
| `amber_score_type` | não | `"dock6"` |
| `amber_score_movable_region` | não | `"ligand"` |
| `amber_score_minimization_rmsgrad` | não | `0.01` |
| `amber_score_before_md_minimization_cycles` | não | `100` |
| `amber_score_md_steps` | não | `3000` |
| `amber_score_after_md_minimization_cycles` | não | `100` |
| `amber_score_gb_model` | não | `5` |
| `amber_score_nonbonded_cutoff` | não | `18.0` |
| `amber_score_temperature` | não | `300.0` |
| `amber_score_movable_distance_cutoff` | não | `10.0` |
| `amber_score_receptor_movable_atom_expr` | não | `"::"` |
| `amber_score_complex_movable_atom_expr` | não | `"::"` |
| `amber_score_abort_on_unprepped_ligand` | não | `true` |
| `use_existing_ligand_charges` | não | `false` |
| `residue_parameter_policy` | não | `"auto_curated"` |
| `md_cache_epoch` | não | `""` |

### retrieve_pubchem

`biomolexplorer.pubchem_retrieval.retrieve_pubchem`

| Parâmetro | Obrigatório | Padrão da função |
| --- | --- | --- |
| `target` | não | `"PubChem"` |
| `reference_source` | não | `"compounds"` |
| `reference_type` | não | `"smiles"` |
| `reference` | não | `""` |
| `selection_mode` | não | `"all"` |
| `compound_id` | não | `""` |
| `threshold` | não | `75` |
| `max_records` | não | `1000` |
| `base_input_path` | não | `null` |
| `input_file` | não | `null` |

### retrieve_compounds

`wrappers.crawlers.retrieve_compounds`

| Parâmetro | Obrigatório | Padrão da função |
| --- | --- | --- |
| `search_term` | sim | `—` |
| `include_pubchem` | não | `false` |
| `pubchem_threshold` | não | `75` |
| `pubchem_max_records` | não | `1000` |
| `chembl_filters` | não | `null` |
| `search_mode` | não | `"target"` |
| `max_targets` | não | `25` |
| `max_records` | não | `1000` |
| `expand_chembl` | não | `true` |
| `similarity_threshold` | não | `70` |

### expand_similar_compounds

`wrappers.crawlers.expand_similar_compounds`

| Parâmetro | Obrigatório | Padrão da função |
| --- | --- | --- |
| `search_term` | sim | `—` |
| `threshold` | não | `75` |
| `max_records` | não | `1000` |
| `base_input_path` | não | `null` |
| `input_file` | não | `null` |

### retrieve_structures

`wrappers.crawlers.load_pdb`

| Parâmetro | Obrigatório | Padrão da função |
| --- | --- | --- |
| `target` | não | `"Estruturas"` |
| `pdb_query` | não | `null` |
| `pdb_ids` | não | `null` |
| `uniprot_ids` | não | `null` |
| `ligand_ids` | não | `null` |
| `max_records` | não | `100` |
| `pdb_ec` | não | `null` |
| `organism` | não | `null` |
| `PolymerEntityTypeID` | não | `null` |
| `ExperimentalMethodID` | não | `null` |
| `max_resolution` | não | `null` |
| `must_have_ligand` | não | `true` |
| `ligand_selection` | não | `"recommended"` |
| `excluded_pdb_ids` | não | `null` |

### retrieve_zinc

`wrappers.crawlers.load_zinc`

| Parâmetro | Obrigatório | Padrão da função |
| --- | --- | --- |
| `base_input_path` | sim | `null` |
| `filename` | não | `"zinc_urls.txt"` |
| `verbose` | não | `0` |
| `download_workers` | não | `4` |

### prepare_structures

`wrappers.redocking.prepare_structures`

| Parâmetro | Obrigatório | Padrão da função |
| --- | --- | --- |
| `base_input_path` | sim | `—` |
| `target` | sim | `—` |
| `pdb_codes` | não | `null` |
| `pH` | não | `7.4` |
| `charge_type` | não | `"gas"` |
| `preparation_options` | não | `null` |
| `base_selected_mols` | não | `null` |
| `mol_filename` | não | `"compounds"` |
| `receptor_prepared` | não | `false` |
| `docking_engines` | não | `"both"` |

### admet

`wrappers.admet.ADMETWrapper`

| Parâmetro | Obrigatório | Padrão da função |
| --- | --- | --- |
| `base_input_path` | sim | `null` |
| `input_file` | não | `null` |
| `verbose` | não | `0` |
| `admet_profile` | não | `"legacy"` |
| `admet_mode` | não | `null` |
| `prediction_backend` | não | `"rdkit"` |
| `bbb_target` | não | `null` |
| `hia_target` | não | `null` |
| `rule_set` | não | `null` |
| `max_rule_violations` | não | `null` |
| `alert_policy` | não | `null` |
| `alert_catalogs` | não | `null` |
| `property_limits` | não | `null` |
| `dataset_origin` | não | `null` |

### fingerprints

`wrappers.molecular_analyzer.generate_fingerprints`

| Parâmetro | Obrigatório | Padrão da função |
| --- | --- | --- |
| `base_input_path` | sim | `—` |
| `morgan_n_bits` | não | `2048` |
| `radius` | não | `2` |
| `files` | não | `null` |
| `morgan` | não | `true` |
| `maccs` | não | `true` |
| `pharmacophore` | não | `true` |
| `chunk_size` | não | `100` |

### similarity

`wrappers.molecular_analyzer.compute_similarity`

| Parâmetro | Obrigatório | Padrão da função |
| --- | --- | --- |
| `base_input_path` | sim | `—` |
| `metric` | não | `"similarityFunctions.TanimotoSimilarity"` |
| `fingerprint` | não | `"fingerprints.Morgan"` |
| `filename` | não | `null` |
| `threshold` | não | `null` |
| `approximate` | não | `true` |

### graphs

`wrappers.molecular_analyzer.analyze_graphs`

| Parâmetro | Obrigatório | Padrão da função |
| --- | --- | --- |
| `similarity_path` | não | `null` |
| `base_input_path` | não | `null` |
| `mcs_timeout` | não | `0` |
| `mcs_ring_matches_ring_only` | não | `true` |
| `mcs_complete_rings_only` | não | `false` |
| `graph_inputs` | não | `null` |

### redocking

`wrappers.redocking.perform_redocking`

| Parâmetro | Obrigatório | Padrão da função |
| --- | --- | --- |
| `base_input_path` | sim | `—` |
| `target` | sim | `—` |
| `pdb_codes` | não | `null` |
| `pH` | não | `7.4` |
| `sizeof_box` | não | `[24, 24, 24]` |
| `exhaustiveness` | não | `20` |
| `num_modes` | não | `10` |
| `prepare_complex` | não | `true` |
| `charge_type` | não | `"gas"` |
| `preparation_pairs` | não | `null` |
| `preparation_options` | não | `null` |
| `redocking_engine` | não | `"vina"` |
| `dock6_app_path` | não | `null` |
| `density` | não | `0.5` |
| `radius` | não | `1.4` |
| `distance` | não | `10.0` |
| `conformer_search_type` | não | `"flex"` |
| `max_orientations` | não | `1000` |
| `num_scored_conformers` | não | `1` |
| `dock6_seed` | não | `73` |
| `dock6_mode` | não | `"redock"` |

### docking_vina

`wrappers.docking.perform_docking_vina`

| Parâmetro | Obrigatório | Padrão da função |
| --- | --- | --- |
| `base_input_path` | sim | `—` |
| `target` | sim | `—` |
| `base_selected_mols` | sim | `—` |
| `mol_filename` | sim | `—` |
| `pdb_code` | não | `null` |
| `pH` | não | `7.4` |
| `sizeof_box` | não | `[24, 24, 24]` |
| `exhaustiveness` | não | `20` |
| `num_modes` | não | `10` |
| `preparation_options` | não | `null` |
| `receptor_prepared` | não | `true` |

### docking_dock6

`wrappers.docking.perform_docking_dock6`

| Parâmetro | Obrigatório | Padrão da função |
| --- | --- | --- |
| `base_input_path` | sim | `—` |
| `target` | sim | `—` |
| `base_selected_mols` | sim | `—` |
| `dock6_app_path` | sim | `—` |
| `charge_type` | sim | `—` |
| `mol_filename` | sim | `—` |
| `pdb_code` | não | `null` |
| `base_vina_path` | não | `null` |
| `density` | não | `0.5` |
| `radius` | não | `1.4` |
| `distance` | não | `10.0` |
| `conformer_search_type` | não | `"flex"` |
| `plot_max_residues` | não | `50` |
| `pH` | não | `7.4` |
| `preparation_options` | não | `null` |
| `receptor_prepared` | não | `true` |
| `dock6_mode` | não | `"auto"` |

### consensus

`wrappers.docking.generate_consensus`

| Parâmetro | Obrigatório | Padrão da função |
| --- | --- | --- |
| `base_input_path` | sim | `—` |
| `target` | sim | `—` |
| `repulsion_weight` | não | `1.0` |
| `base_vina_path` | não | `null` |
| `base_dock6_path` | não | `null` |

## Configuração e ferramentas
 
`AppConfig` controla `workspace`, `max_jobs=1`, `cpu_workers=2`, `max_queued_jobs=100`, `job_timeout=None`, `worker_python=None`, `resource_dir=None`, `cache_dir=None`, `reuse_results=True`, `preparation_cache=None`, `reuse_preparations=True`. São opções da API Python, não flags adicionais da CLI do workspace.

O empacotamento exige Python >=3.12. `requirements.yml` fixa Python 3.12.3 e a pilha científica/UI. `pyproject.toml` não declara dependências científicas de base; o extra `ui` adiciona Flet desktop/web. Chimera e DOCK6 externos não são instalados pelo pip. Compiladores, make, flex e zlib dão suporte à compilação AMBER em cópia local. O preparo opcional do ADMET-AI é separado.

### Variáveis de ambiente utilizadas pelo código

| Nome | Referências de implementação |
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

O inventário inclui valores internos propagados pelos workers. A presença de uma variável não a torna automaticamente um controle de usuário; consulte a implementação indicada e os guias de instalação/preparação.

## Formatos e persistência
 
As entradas dependem das portas semânticas: CSV/SMILES e estruturas compatíveis de compostos, receptor PDB e conjuntos preparados PDBQT/MOL2, binários de fingerprints, CSVs de similaridade/grafos e manifests de resultados autorizados. Amber Score exige target PDB e pose com coordenadas correspondente ao alvo, não apenas SMILES. Arquivos nativos, CSVs de compatibilidade e imagens/relatórios permanecem persistidos; índices SQLite atendem consulta, linhagem, paginação e caches científicos e não substituem todos os artefatos estruturais.

## Compatibilidade e comportamentos descontinuados
 
- Não há backend AmberTools independente selecionável para preparar docking; identificadores `auto`/`python` de preparo são normalizados para Chimera.
- Amber Score nativo não é dinâmica sander/cpptraj; não oferece trajetória convencional, RMSF, controles de equilíbrio ou parametrização arbitrária de metais.
- A geração de z-score/min-max no consenso e a antiga fórmula com energia truncada foram descontinuadas; leitores legados continuam compatíveis.
- Atualizações automáticas instalam um pacote validado; a extração manual de ZIP não é o fluxo do botão. O reinício continua manual.
- Não existe operação pública de análise relacional em `OPERATIONS`; arquivos/índices SQLite não devem ser apresentados como novo bloco científico.
- Argumentos legados de tempo/trajetória mantidos na função Python não são parâmetros JSON aceitos se não constarem no registro acima.

## Scripts de automação

As flags abaixo foram extraídas dos parsers. Scripts sem parser são utilitários Python com pontos de entrada próprios, não operações científicas adicionais. Alguns iniciam motores, acessam provedores, constroem índices derivados ou conduzem sessões de navegador; confira as opções e use um diretório de saída separado antes de executá-los. Não foram executados apenas para validar esta documentação.

| Script | Flags declaradas no código |
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
| `scripts/benchmark_report.py` | Sem flags argparse |
| `scripts/benchmark_result_index.py` | `--rows`, `--database`, `--run-id`, `--output` |
| `scripts/build_docs.py` | Sem flags argparse |
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
| `scripts/validate_docs.py` | Sem flags argparse |
| `scripts/validate_native_structures.py` | `--dock6`, `--output`, `--steps` |
| `scripts/validate_native_surface.py` | `--dms-executable` |
| `scripts/validate_pdb_viewer.py` | `--chrome`, `--pdb`, `--structure`, `--smiles`, `--screenshot` |
| `scripts/validate_preparation_migration.py` | `--dock6`, `--output` |
| `scripts/validate_settings_browser.py` | `--output` |
