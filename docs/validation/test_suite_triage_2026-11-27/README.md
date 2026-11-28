# Triagem da suíte de testes de referência — 27/11/2026

Tarefa separada, proposta na PHASE-00 de #6 (achados F-16, F-17 e F-13 de
[`artifact_architecture_assessment_2026-11-26`](../artifact_architecture_assessment_2026-11-26/README.md)).
Branch `fix/test-suite-baseline`, criado a partir de `e88f751`.

## Resultado

| Medida | Linha de base (PHASE-00) | Depois desta triagem |
| --- | --- | --- |
| Módulos sem falhas | 124 de 169 | **156 de 169** (13 restantes = C e I abaixo) |
| Testes executados | 1.663 | 1.671 (2 pulados por ambiente) |
| Falhas / erros | 50 / 38 | 22 / 6 |
| Execução completa (`unittest discover -s tests`) | **trava** (F-17) | **termina**: 1.671 testes em 756 s, mesmos resultados da execução por módulo |
| Módulos que passaram a falhar | — | **0** |

Cada módulo foi executado separadamente, com 6 em paralelo e limite de 600 s. Os módulos
que falhavam também foram repetidos em série, com o mesmo resultado.

Cada falha foi classificada em uma de quatro categorias:

- **A** — teste desatualizado em relação a uma mudança intencional e documentada. O teste
  foi atualizado sem enfraquecer o que ele verifica.
- **B** — bug no código. Foi corrigido com uma mudança mínima.
- **C** — divergência que exige decisão do mantenedor. Nenhuma mudança foi feita, e o
  teste continua falhando pelo motivo correto.
- **I** — exige investigação mais profunda. Ficou como está.

## B — Bugs corrigidos no código

| # | Arquivo | Problema | Efeito |
| --- | --- | --- | --- |
| B1 | `workspace.py` (`save_pipeline`) | Um salvamento sem mudança de pipeline retornava cedo e ignorava `prune_results=True`. | Depois de uma curadoria ou edição em disco, a confirmação de invalidação feita na interface (`ui/app.py:1237`) não tinha efeito e a execução ficava bloqueada. |
| B2 | `result_sync.py` (`plan`) | Mudar o modo de `individual` para `both` (ou o inverso) invalidava os lotes individuais. | Recálculo desnecessário de lotes cientificamente idênticos. Lotes combinados continuam sendo invalidados. |
| B3 | `pipeline.py` | "Não, executar novamente" passava `reuse_results=False`, mas não `reuse_preparations=False`, regressão de `9bdb316`. | Os checkpoints de `CalculationCache` e `LigandCache` (ligantes, Vina, DOCK6, pares de redocking) eram reaproveitados mesmo quando o usuário pedia recálculo, contrariando `docs/architecture.md`. |
| B4 | `network_centrality.py` (`betweenness`) | O sorteio das fontes usava a ordem de um `set` interno, que depende do `PYTHONHASHSEED`. | **A betweenness aproximada mudava entre execuções** (por exemplo 0,001998 / 0,001931 / 0,001858 para o mesmo nó), violando a reprodutibilidade (01.3). Agora o valor é idêntico com qualquer semente de hash. As comunidades Louvain foram verificadas e estavam estáveis. |
| B5 | `ui/app.py` (`prepare_pending_inputs`) | O pré-carregamento de entradas usava o adaptador de processamento em primeiro plano. | Podia abrir a tela de processamento ou falhar sem página. Agora roda como tarefa de segundo plano (`background`), como o polling. |
| B6 | `resources/i18n/*.json` | 35 mensagens sem tradução para o inglês e 44 pares inversos ausentes no catálogo português. | As mensagens apareciam em português na interface em inglês. As entradas foram acrescentadas sem alterar as existentes. |
| — | `docs/artifact_naming.md` | A regra de rota semântica de SRC (`95f1b10`) não estava documentada. | Documentação alinhada ao código. |

## A — Testes atualizados (mudanças intencionais)

| Mudança no produto | Testes |
| --- | --- |
| Sincronização automática de resultados (`9bdb316`): mudança que afeta resultados exige confirmação, e a interface salva antes de `resume` | `test_analysis_handoff`, `test_pipeline_guidance`, `test_pubchem_workflow`, `test_curated_redocking` (só a confirmação) |
| Referências semânticas tipadas (`output_port`) e migração ao reabrir | `test_canvas_drag_regression`, `test_guided_input_preservation`, `test_automatic_result_sync`, `test_result_persistence` (recuperação após falha) |
| Convenção global de nomes (nomes nativos e etiqueta de origem) | `test_input_processing`, `test_pubchem_workflow`, `test_run_progress`, `test_chembl_rest`, `test_result_visualizations`, `test_global_artifact_names` (rota SRC) |
| Modelo de consenso `target-balanced-ecr-v1` (scores brutos em `*_raw`) | `test_workspace`, `test_docking_handoffs` |
| Preparação exclusiva pelo Chimera (26/10) | `test_native_structures` (minimização via `sander` removida), `test_chimera_charge_replacement` (rota direta pelo Antechamber removida), `test_zinc_tranches`, `test_ligand_preparation_performance` (o motor não refaz a preparação compartilhada) |
| MOL2 nativo sempre guardado para reuso, e preparação em paralelo | `test_docking_preparation`, `test_docking_interoperability` |
| Patches de `subprocess` globais demais, ou de uma referência que o código não usa mais | `test_pubchem_workflow`, `test_ligand_preparation_performance`, `test_docking_interoperability`, `test_zinc_tranches`, `test_pipeline_failure_recovery` (Vina via `vina_batch`) |
| Fluxo de execução com tela de carregamento e confirmação da ordem de execução | `test_run_reuse_ui` (**era o teste que travava a suíte**), `test_configuration_loading`, `test_lazy_run_configuration`, `test_guided_ui`, `test_mcs_ranking_refinements` |
| Revisão consolidada de resultados (sem histórico recolhível) e métrica de arquivos ainda existentes | `test_workspace_actions`, `test_frontend_navigation` |
| Tipos de importação (`receptor_pdb`), parâmetros internos ocultos e saídas ADMET particionadas | `test_import_files`, `test_workspace`, `test_molecule_quality` |
| Texto de validação | `test_redocking_chimera_repair` (mensagem) |

> **Continuação:** as decisões C1–C8 e I1–I2 abaixo foram tomadas e aplicadas a
> testes e documentação em
> [`test_suite_decisions_2026-11-27`](../test_suite_decisions_2026-11-27/README.md),
> que também traz o censo de exclusões de preparação pedido por C8 e dois
> achados novos (C1a, C1b). As correções de código de C1, C2, C3 e C7c estão em
> [`test_suite_code_fixes_2026-11-27`](../test_suite_code_fixes_2026-11-27/README.md).
> A suíte está verde, e restam 6 falhas esperadas rastreadas: C5, C8 e I1–I2.

## C — Decisões necessárias

| # | Divergência | Testes | Opções |
| --- | --- | --- | --- |
| C1 | **Contrato de saídas ADMET.** `docs/admet_output_selection.md` proíbe `evaluated` e diz que `alerts` só é gerado com catálogos; o código grava os dois. Os testes se contradizem: um exige `evaluated`, outro proíbe; um espera o `compounds_INPUT.csv` primário; outro espera que o perfil `legacy` filtre compostos. Projetos reais consomem `*_evaluated.csv`. | `test_admet_standard_outputs` (7), `test_admet_research`, `test_application` | (a) a documentação vale: remover `evaluated` e o `alerts` vazio, com migração; (b) o código vale: atualizar a documentação e os testes |
| C2 | **Invalidação por par no redocking.** Editar os ligantes de um único PDB muda o `pdb_codes.csv` compartilhado e invalida **todos** os pares, que são recalculados. | `test_curated_redocking` | Tornar a invalidação parcial (`pair_removals`) sensível à linha do PDB alterado; relacionado a PHASE-03/04 |
| C3 | **Modo `both` com arquivos enviados** executa como `merge` (um lote só). B2 corrige o plano de invalidação, não a execução. | `test_input_processing` (`both`) | Corrigir a execução ou retirar o suporte |
| C4 | **Reparo de metadados a partir de pares curados.** O teste nunca passou, e a documentação promete o caminho inverso (pares a partir dos metadados). | `test_input_processing` (redocking) | Implementar ou retirar o teste |
| C5 | **Matriz de conexões.** Poses de redocking entrando no Amber Score e na validação de poses; tipos `other` de importação. | `test_connection_matrix`, `test_redocking_receptor_outputs` | Definir o contrato científico das portas |
| C6 | **Cancelar na configuração de uma execução pausada.** Antes cancelava a execução (`Cancelar execução`); desde 27/10 só fecha a caixa. | `test_execution_cancellation` (2) | Confirmar o comportamento atual e atualizar os testes, ou restaurar o cancelamento |
| C7 | **Pausa para seleção após reinício ou recusa de reuso.** A seleção curada agora é persistida, então não há nova pausa. Uma `cache_key` adulterada é aceita, porque o reuso confere a requisição e os hashes, não a chave gravada. | `test_restart_reuse` (6), `test_file_selection` (2) | Confirmar a semântica atual (e atualizar os testes) ou restaurar as pausas e a verificação da chave |
| C8 | **Reparo automático de cargas no Chimera.** O reparo com tipos atômicos explícitos após falha do Antechamber foi desligado em `9bdb316`; `chimera_charge_repair.charge_commands` ficou sem uso. | `test_redocking_chimera_repair` (1) | Religar o reparo (pode reduzir exclusões de preparação) ou remover o código morto e o teste |

## I — Investigação pendente

| # | Problema | Testes | Relação com #6 |
| --- | --- | --- | --- |
| I1 | Cache v1 legado não é adotado ao acrescentar uma entrada: o lote antigo é recalculado (3 chamadas em vez de 2). | `test_result_persistence` (2) | PHASE-07 (migração) e REQ-CACHE |
| I2 | Seletores antigos de saídas PubChem não são resolvidos ao expandir a seleção. | `test_analysis_origins` | PHASE-07 (projetos legados) |

## Higiene (F-13)

Os `.pyc` versionados em `src/**/__pycache__` continuam no controle de versão. Nesta
triagem os testes rodaram com `PYTHONDONTWRITEBYTECODE=1`, para não alterá-los. A remoção
desses arquivos do repositório fica como recomendação, fora desta mudança.

## Reprodução

```sh
PY=~/progs/miniconda3/envs/BioMolExplorer/bin/python
export PYTHONPATH=src:tests PYTHONDONTWRITEBYTECODE=1
$PY -m unittest discover -s tests
```
