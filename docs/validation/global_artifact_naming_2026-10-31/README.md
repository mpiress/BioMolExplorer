# Convenção global — validação em 31/10/2026

## Implementação

Padrão: `TG[ALVO]-SRC[ORIGEM]-OP[OPERAÇÃO]-TIPO-ID[hash].ext`.
A implementação compartilhada cobre os 16 tipos de bloco, mantém apenas a origem imediata no nome e conserva a linhagem completa nos descritores. O sufixo determinístico distingue ocorrências, inclusive entre execuções. Não há novas bibliotecas.

A convenção é aplicada como **nome público persistente**, usado na exibição, seleção e download. Os caminhos físicos operacionais dos programas e os nomes privados da publicação analítica anterior permanecem compatíveis. Não foi realizada uma renomeação global de arquivos físicos: isso romperia referências de poses e pacotes científicos. Novas saídas registram o mapeamento atomicamente em `.artifact-publication.json`; resultados legados são indexados sem escrever em seus arquivos. ZIPs de simulação contêm um índice adicional e mantêm os caminhos internos nativos; backups preservam descritores também de execuções legadas.

Docking/Amber Score propagam o alvo do receptor, e não o alvo da coleção do ligante. Coleções sem alvo biológico conhecido usam `TG[NA]`. Categorias intermediárias incluem `MOL`, `SIM`, `BBB+`, `Fingerprints-Morgan`, `REC` e `PDB`; misturas e conjuntos extensos têm resumo curto com os detalhes nos metadados. Tipos representam o conteúdo efetivo; nenhuma trajetória/métrica de dinâmica é criada pelo mecanismo de nomes.

## Exemplos reais do projeto existente

Estes nomes foram obtidos da execução `d5250afecb2c4a04aae1426ddc37699b`, em modo somente leitura. São aliases dos arquivos existentes, não resultados recalculados.

| Bloco | Nome público real |
| --- | --- |
| retrieve_compounds | `TG[ACHE]-SRC[ChEMBL]-OP[ChEMBL]-MOL-ID[c195c00a].csv` |
| retrieve_structures | `TG[ACHE]-SRC[PDB]-OP[PDB]-Complex-4M0E-ID[4ccd13ef].pdb` |
| retrieve_pubchem | `TG[ACHE]-SRC[BBB+]-OP[PubChem]-SIM-ID[410a0b13].csv` |
| admet | `TG[ACHE]-SRC[MOL]-OP[ADMET]-BBB+-ID[24113d53].csv` |
| fingerprints | `TG[ACHE]-SRC[BBB+]-OP[Fingerprints]-Fingerprints-Morgan-ID[d8c7d19d].csv` |
| similarity | `TG[ACHE]-SRC[Fingerprints-Morgan]-OP[Similarity]-Similarity-Tanimoto-ID[7b97ff4d].csv` |
| graphs | `TG[ACHE]-SRC[Similarity-Tanimoto]-OP[Graphs]-Molecules-LCC-ID[ab7509e3].csv` |
| redocking | `TG[BCHE]-SRC[PDB]-OP[Redocking]-Pose-1P0I_BUA_606A-ID[a56de9e9].pdbqt` |

Os testes de contrato também cobrem Docking-Vina, Docking-DOCK6, Consensus e AmberScore-DOCK6 com artefatos controlados. Não houve nova execução nativa de Vina, DOCK6 ou Amber Score nesta tarefa. A validação de nomenclatura desses produtores não equivale à validação científica de um novo experimento.

## Integridade, desempenho e compatibilidade

- 1.620 artefatos existentes, 1.620 ocorrências e 1.620 nomes distintos; todos os 1.620 arquivos encontrados tiveram seus bytes comparados por SHA-256 antes/depois do catálogo. Nenhuma escrita foi feita no projeto original.
- Catálogo completo: aproximadamente 0,46 s em cada uma de duas medições locais após eliminar decoração duplicada e leituras repetidas de cabeçalhos. Esse tempo mede a construção do catálogo, não todo o fluxo de abertura da interface.
- Os 16 identificadores efetivos de versão científica foram comparados com a captura anterior à alteração e preservados por equivalências exatas auditadas. O identificador da normalização de compostos também foi preservado. Alterações futuras não recebem equivalência automaticamente.
- O metadado de publicação não entra no manifesto científico de diretórios. A alteração de bytes científicos continua alterando o manifesto; há teste explícito para ambas as condições.
- Reinicialização, reúso sem novos workers, exportação/importação com remapeamento de execução e projeto legado sem descritores/sidecars foram exercitados em projetos temporários.
- Navegador Chrome/Flet real: 360, 1024 e 1440 px, nomes canônicos, ordenação sem fechar o painel, botão de fechar acessível. Zero exceções/erros JavaScript; a mensagem informativa “Flutter app loaded” aparece no console. As tabelas usam sua rolagem interna nos espaços estreitos.

## Testes executados

A bateria principal de aceitação inicial executou 48 testes com sucesso; uma execução posterior de 12 testes específicos inclui a exportação de projeto legado sem descritores, também aprovada. O subconjunto de regressão de cache executou **49 testes com sucesso** em 72,49 s (`cache-regression.log`), com as sete exclusões explícitas reproduzidas na linha de base. O script `cache_regression_runner.py` permite repetir essa seleção. A bateria integrada final executou **145 testes com sucesso em 173,11 s** (`integrated-suite.log`). Total das duas baterias finais: **194 testes aprovados**. O popup de footprint e seus contratos de acesso/conteúdo também passaram em uma bateria adicional de 11 testes (`footprint-download.log`).

Após os últimos ajustes de classificação dos componentes e alvo explicitamente informado, os **15 testes de metadados/subconjuntos** foram repetidos com sucesso (`final-metadata.log`); o caso adicional de alvo explícito foi aprovado (`explicit-target.log`). A investigação de visualizações executou 41 casos: 40 passaram e um ainda esperava `Graph-mcc` para uma tabela que contém arestas; a expectativa correta `Edges-LCC` foi aprovada na repetição isolada (`graph-download-recheck.log`). Essa divergência de expectativa e a ausência de arquivos redundantes de componente na fixture anterior também estão documentadas nos logs de investigação.

Os testes atualizados mantêm as verificações de conteúdo, IDs, permissões e integridade, substituindo expectativas de nomes antigos pela interpretação da nova convenção. Os ZIPs agora têm um membro adicional (`artifact-names.json`). Fixtures que selecionavam uma tabela ADMET agregada obsoleta passaram a selecionar a tabela efetivamente produzida `_evaluated.csv`; nenhum algoritmo foi alterado para satisfazer testes.

Para repetir as baterias a partir da raiz do repositório, no ambiente BioMolExplorer existente:

```bash
PYTHONPATH=src:tests python -m unittest test_global_artifact_names test_artifact_publication test_artifact_registry test_stable_artifact_selectors test_dataset_identity test_project_portability test_result_persistence test_published_results test_results_display test_redocking_results test_result_sorting test_graph_input_contract
PYTHONPATH=src:tests python docs/validation/global_artifact_naming_2026-10-31/cache_regression_runner.py
PYTHONPATH=src:tests python -m unittest test_footprint_popup test_docking_scene test_graph_subset_exports test_graph_explorer
```

As execuções registradas usaram o Python do ambiente existente, `PYTHONDONTWRITEBYTECODE=1` e `MPLCONFIGDIR=/tmp/biomol-mpl`. As evidências de browser foram obtidas com servidor local temporário e Chrome já instalado, sem novas dependências.

### Falhas preexistentes identificadas separadamente

Uma cópia intacta de `HEAD`, extraída para `/tmp`, executou a mesma bateria de 56 testes de cache/redocking: 4 falhas e 3 erros. Os sete casos abaixo não são declarados aprovados; estão registrados em `baseline-regression.log` e foram excluídos somente da bateria comparativa final dos outros 49 casos:

- `test_adding_unconnected_block_and_moving_canvas_preserves_existing_results`: expectativa de novo worker para bloco equivalente, enquanto o sistema compartilha o resultado.
- `test_changed_processing_mode_after_restart_requires_new_selection`: fixture salva alteração destrutiva sem confirmação exigida pelo sistema.
- `test_changed_uploaded_input_invalidates_consumer_after_upgrade`: fixture submete invalidação sem confirmação.
- `test_declining_reuse_runs_collection_and_waits_for_new_file_selection`: expectativa de seleção novamente, enquanto a escolha confirmada foi preservada.
- `test_modified_consumer_output_cannot_be_reused_after_upgrade`: expectativa de espera difere do fluxo atual de retomada.
- `test_tampered_original_cache_key_is_not_adopted_after_upgrade`: expectativa antiga sobre adoção de chave difere do comportamento atual.
- `test_three_old_pubchem_selections_expand_to_six_real_graphs_with_exact_origins`: fixture histórica contém seleção não disponível/ambígua no contrato atual. Seu mock de publicação foi adaptado para aceitar os novos argumentos contextuais, mas a incompatibilidade de seleção já ocorre em `HEAD`.

A bateria completa do repositório e uma campanha científica nativa de todos os blocos não foram executadas. Esses limites não afetam a evidência de preservação dos arquivos atuais e não são apresentados como aprovação global de release.

## Arquivos modificados

- `docs/en/user_manual.md`
- `docs/user_manual.md`
- `src/biomolexplorer/artifact_choices.py`
- `src/biomolexplorer/artifact_naming.py`
- `src/biomolexplorer/artifacts.py`
- `src/biomolexplorer/compound_tables.py`
- `src/biomolexplorer/docking_inputs.py`
- `src/biomolexplorer/fingerprint_selection.py`
- `src/biomolexplorer/input_validation.py`
- `src/biomolexplorer/pipeline.py`
- `src/biomolexplorer/project_archive.py`
- `src/biomolexplorer/redocking_results.py`
- `src/biomolexplorer/resources/processing_compatibility.json`
- `src/biomolexplorer/result_files.py`
- `src/biomolexplorer/result_origins.py`
- `src/biomolexplorer/stage_cache.py`
- `src/biomolexplorer/ui/app.py`
- `src/biomolexplorer/ui/compound_table.py`
- `src/biomolexplorer/ui/docking_results.py`
- `src/biomolexplorer/ui/docking_scene.py`
- `src/biomolexplorer/ui/file_selection.py`
- `src/biomolexplorer/ui/file_table.py`
- `src/biomolexplorer/ui/footprint.py`
- `src/biomolexplorer/ui/graph_explorer.py`
- `src/biomolexplorer/ui/graph_table.py`
- `src/biomolexplorer/ui/input_editor.py`
- `src/biomolexplorer/ui/redocking_results.py`
- `src/biomolexplorer/ui/result_table.py`
- `src/biomolexplorer/ui/stage_results.py`
- `tests/test_analysis_origins.py`
- `tests/test_artifact_publication.py`
- `tests/test_dataset_identity.py`
- `tests/test_footprint_popup.py`
- `tests/test_graph_explorer.py`
- `tests/test_graph_subset_exports.py`
- `tests/test_project_portability.py`
- `tests/test_redocking_results.py`
- `tests/test_results_display.py`
- `tests/test_stable_artifact_selectors.py`
- `tests/test_global_artifact_names.py` (novo)
- `docs/artifact_naming.md` (novo)
- Este relatório e evidências em `docs/validation/global_artifact_naming_2026-10-31/`.

O mecanismo e as regras para usuários estão descritos em [Convenção global](../../artifact_naming.md). Logs finais contêm comandos/contextos de suíte; as capturas não usam dados de credenciais nem modificam o experimento original.
