# Decisões C1–C8 e I1–I2: aplicação em testes e documentação — 27/11/2026

Continuação da [triagem da suíte de testes](../test_suite_triage_2026-11-27/README.md).
Branch `fix/test-decisions-pr1`, a partir de `281abe4` (PR #16).

Este é o primeiro de três passos combinados. **Aqui só mudam testes e
documentação**: nenhuma linha de `src/` foi alterada. O objetivo é deixar a
suíte verde e cada divergência registrada de forma rastreável.

| Medida | Triagem (PR #16) | Agora |
| --- | --- | --- |
| Módulos sem falhas | 156 de 169 | **169 de 169** |
| Falhas / erros | 22 / 6 | **0 / 0** |
| `unittest discover -s tests` | 22 falhas e 6 erros | **OK**, 1.673 testes em 738 s |
| Falhas esperadas (`expectedFailure`) | 0 | **10**, cada uma com a decisão ou investigação que a resolve |
| Mudanças em `src/` | 10 arquivos | **nenhuma** |

Critério aplicado: **o que o código faz hoje é o contrato, desde que não viole
os invariantes científicos (01.3 de #6), não descarte dados e não recalcule sem
necessidade.** Quando o código viola, o código é corrigido (passo 2). Quando
não viola, testes e documentação passam a descrever o comportamento real.

## Decisões aplicadas (testes e documentação)

### C1 — Contrato das saídas ADMET

**Decisão: o código é o contrato, com duas limpezas pendentes.** A tabela
`evaluated` é mantida e passa a ser documentada: é a única visão com todos os
compostos válidos e o motivo de cada exclusão (`ADMET_selected`,
`ADMET_reasons`), e é consumida por projetos reais.

[`docs/admet_output_selection.md`](../../admet_output_selection.md) foi
atualizado: a tabela de avaliados entra na lista publicada, a proibição de
`evaluated` saiu, e a seção **Limites conhecidos** registra o que será corrigido
no passo 2.

Dois achados novos, encontrados ao escrever os testes:

| # | Achado | Evidência |
| --- | --- | --- |
| C1a | A tabela de regras é gravada **duas vezes**: `_rules-<regra>.csv` (canônica, só a regra) e `_rules_<regra>.csv` (legada, também filtrada pela hipótese BBB/HIA e com a coluna `Rule_set`). Não são cópias idênticas. | 8 CSVs por entrada, em vez dos 7 do contrato |
| C1b | **Defeito.** Na publicação, a tabela legada é renomeada para `compounds.csv` e o relatório passa a apontar o papel da regra para esse arquivo, em vez da tabela canônica. | `test_published_report_points_the_rule_role_at_the_rule_table` (falha esperada) |

Também ficou registrado que a tabela de alertas é gravada mesmo sem catálogos
selecionados, o que a documentação já proibia.

Correções no passo 2: remover a tabela legada (o que resolve C1a e C1b) e não
gravar alertas sem catálogo.

### C4 — Reparo de metadados a partir de pares curados

**Decisão: o recurso não existe e não será criado.** Selecionar uma estrutura
cujo PDB não tem registro correspondente nos metadados continua sendo recusado
explicitamente, em vez de inventar o registro. O teste passou a verificar a
recusa e que os metadados em disco permanecem intactos.

### C6 — Cancelar na configuração

**Decisão: fechar uma caixa de configuração nunca é destrutivo.** Cancelar uma
execução é uma ação explícita do painel de execução, por `WorkspaceUI.cancel_run`.

Ao escrever os testes, confirmou-se que **nenhuma** das duas caixas cancela mais
a execução: a de configuração tem `Cancelar` (só fecha) e a de seleção tem
`Selecionar depois` (mantém a execução pausada). Três testes cobrem agora os
três caminhos, incluindo o cancelamento real, que fecha as caixas, descarta o
rascunho, atualiza o painel e avisa o usuário.

### C7 — Reuso após reinício

**Decisão: a semântica atual está correta**, com uma distinção que o teste
antigo não fazia:

- **Não há nova pausa** quando a execução é apenas recalculada: a seleção curada
  é configuração do projeto e é reaproveitada.
- **Há pausa** quando a mudança altera o que cada lote consome, como trocar o
  modo para `individual`. O teste correspondente estava certo; só faltava
  confirmar a invalidação antes de salvar.
- **A `cache_key` gravada é metadado derivado.** O reuso recalcula a assinatura a
  partir da requisição salva e confere os SHA-256 dos artefatos, então uma chave
  adulterada não decide nada por si. O teste passou a verificar essa política
  (`verified_completed_request_v1`), que é mais forte que confiar na chave.

Um defeito real foi separado e continua registrado: **C7c**, em que o consumidor
é recalculado após o reinício, sem nada ter mudado, quando o produtor publicou
mais de um arquivo e um deles foi selecionado. É o recálculo redundante medido na
PHASE-00.

Em `test_file_selection`, duas falhas tinham outra causa: o fixture não definia
`ui.base_pipeline`, que a interface real sempre define ao abrir um projeto. Vale
registrar que, sem ele, `save_draft` falha e a execução **não** é retomada, sem
mensagem — o que merece um aviso explícito na interface (candidato ao passo 2).

### C8 — Reparo automático de cargas do Chimera

**Decisão: adiada, agora com dados.** Censo somente de leitura das 77 tabelas
`preparation_failures.csv` do projeto Alzheimer, 12.254 linhas (uma por
composto e receptor):

| Causa | Linhas | Compostos distintos |
| --- | --- | --- |
| Topologia ou tipos atômicos MOL2 inválidos | 10.521 | **694** |
| PDBQT com múltiplos fragmentos ROOT (contraíons) | 576 | 48 |
| Falhas de ferramenta (ANTECHAMBER e Chimera) | 108 | **6** |
| PDBQT com coordenadas, cargas ou tipos inválidos | ~90 | 5 |

O reparo de cargas atacaria **6 dos 809 compostos excluídos (menos de 1%)**, e
não é prioridade. A decisão entre religar o reparo e remover o código morto
continua aberta, com custo baixo nos dois casos.

**Achado novo e mais relevante:** a causa dominante das exclusões é a topologia
ou os tipos atômicos MOL2, com **694 compostos perdidos** no projeto Alzheimer.
Isso merece investigação própria, fora de #6, e não tem relação com C8.

## Pendências registradas como falhas esperadas

Cada uma usa `@unittest.expectedFailure` com a decisão ou a investigação no
docstring. A suíte fica verde e, quando a causa for corrigida, o próprio
unittest acusa "unexpected success".

| Marcador | Teste | Resolve em |
| --- | --- | --- |
| C1b | `test_admet_standard_outputs.test_published_report_points_the_rule_role_at_the_rule_table` | passo 2 |
| C2 | `test_curated_redocking.test_1p0i_bua_saved_after_first_run_is_docked_and_accumulated` | passo 2 |
| C3 | `test_input_processing.test_both_reuses_individual_jobs_and_merges_vectors_without_a_worker` | passo 2 |
| C5 | `test_connection_matrix.test_every_native_imported_and_supplied_result_against_every_port` | passo 3 (contrato das portas) |
| C5 | `test_redocking_receptor_outputs.test_redocking_cannot_connect_to_any_compound_or_score_input` | passo 3 |
| C7c | `test_restart_reuse.test_restart_restores_only_the_previously_selected_file` | passo 2, após investigação |
| C8 | `test_redocking_chimera_repair.test_chimera_charge_failure_retries_once_with_explicit_types` | passo 3 |
| I1 | `test_result_persistence`, dois testes de migração do cache v1 | PHASE-07 (#14) |
| I2 | `test_analysis_origins.test_three_old_pubchem_selections_expand_to_six_real_graphs_with_exact_origins` | PHASE-07 (#14) |

## Próximos passos combinados

- **Passo 2, correções de código: concluído** em
  [`test_suite_code_fixes_2026-11-27`](../test_suite_code_fixes_2026-11-27/README.md).
  C1a e C1b (tabela legada e alertas sem catálogo removidos), C2 (invalidação
  por par, com execução real de redocking em dois PDBs), C3 (execução do modo
  `both`) e C7c, cujo diagnóstico apontou uma causa diferente da suposta: a
  migração para referências tipadas mudava a chave de cache de todo consumidor
  ao reabrir o projeto. Inclui o aviso explícito quando a execução não é
  retomada.
- **Passo 3: concluído** em
  [`test_suite_port_contract_2026-11-28`](../test_suite_port_contract_2026-11-28/README.md).
  C5 divergia em uma única célula por teste: uma importação sem tipo
  alimentando a recuperação ZINC (o código estava certo, o teste e a
  documentação errados) e as poses do redocking alimentando dinâmica e
  validação de poses. Decidiu-se que o redocking entrega apenas os receptores
  que preparou e validou, declarados por `INTERNAL_PORTS` — alinhando o
  contrato tipado a três documentos que já diziam isso. C8 foi removido: o
  redesenho em duas fases já **detecta** o grafo de camada aberta que o reparo
  contornava, e a remoção foi registrada como equivalência auditada para não
  invalidar preparações já calculadas.

## Reprodução

```sh
PY=~/progs/miniconda3/envs/BioMolExplorer/bin/python
export PYTHONPATH=src:tests PYTHONDONTWRITEBYTECODE=1
$PY -m unittest discover -s tests
```

O censo de exclusões usou os scripts somente de leitura desta avaliação,
registrados em `preparation_failures.json`.
