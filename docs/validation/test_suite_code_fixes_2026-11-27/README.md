# Correções de código das decisões C1–C3, C7 — 27/11/2026

Passo 2 dos três combinados após a
[triagem](../test_suite_triage_2026-11-27/README.md). O passo 1 está nas
[decisões aplicadas a testes e documentação](../test_suite_decisions_2026-11-27/README.md).
Branch `fix/test-decisions-pr2`, a partir de `431e137`.

Aqui o código muda. Cada correção tem um teste que falhava antes e passa agora,
e quatro das dez falhas esperadas registradas no passo 1 foram resolvidas.

| Medida | Passo 1 | Agora |
| --- | --- | --- |
| `unittest discover -s tests` | OK, 1.673 testes | **OK**, 1.674 testes em 734 s |
| Falhas esperadas rastreadas | 10 | **6** (C5 ×2, C8, I1 ×2, I2) |
| Arquivos de `src/` alterados | 0 | 11, dos quais 2 são catálogos de tradução |

## C1a e C1b — uma tabela por papel no ADMET

**Causa.** `ADMETWrapper.export_results` gravava, logo antes de
`admet_outputs.export_standard`, as mesmas quatro partições BBB/HIA e a mesma
tabela de alertas, que eram sobrescritas em seguida. O único arquivo que
sobrevivia era a tabela de regras legada `_rules_<regra>`, filtrada também pela
hipótese BBB/HIA e com a coluna `Rule_set` — uma tabela diferente da canônica,
com o mesmo papel.

O C1b vinha daí: `artifact_naming.output_role` só reconhecia o sufixo legado
`_rules_` (sublinhado). A tabela canônica `_rules-<regra>` (hífen) caía no papel
genérico `compounds` e era publicada como `B01__admet__compounds.csv`, de modo
que o relatório passava a apontar o papel da regra para esse nome, enquanto
`B01__admet__rules-<regra>.csv` ficava com a tabela legada.

**Correção.**

| Arquivo | Mudança |
| --- | --- |
| `src/wrappers/admet.py` | `export_results` removido com sua chamada. Toda saída vem de `export_standard`, que já gravava os mesmos arquivos depois. Elimina também a gravação duplicada de cada partição. |
| `src/biomolexplorer/admet_outputs.py` | A tabela de alertas só é gravada com catálogos selecionados. Sem catálogo ela não tem nenhuma coluna `Alert_*` e apenas duplicaria a tabela de avaliados — o que a documentação já proibia. |
| `src/biomolexplorer/artifact_naming.py` | `output_role` reconhece `_rules-<regra>` e `_rules_<regra>`. O nome legado continua lido, para projetos antigos em disco. |
| `src/biomolexplorer/result_files.py` | Os dois pontos que removem o sufixo para achar o relatório da entrada aceitam as duas formas. |

Por entrada ficam seis arquivos, não oito: quatro partições, a tabela da regra
configurada e a de avaliados, mais a de alertas quando há catálogo.
[`docs/admet_output_selection.md`](../../admet_output_selection.md) foi
atualizado.

## C3 — execução do modo `both`

**Causa.** `_variants` já montava os três lotes do modo `both` (um por entrada
mais o conjunto completo) e `_cache_key` já enumerava os três. Só o despacho da
execução testava `input_processing=='individual'`, em dois pontos de
`pipeline.py`, então `both` caía no caminho de lote único e executava como
`merge`, descartando os lotes individuais.

**Correção.** Os dois testes passaram a aceitar `('individual','both')`. O modo
agora faz o que seu rótulo na interface diz, *Individual + conjunto completo
(ALL)*: os lotes individuais são reaproveitados e só o combinado é calculado.

Ao escrever o teste apareceu um problema de identidade: os três lotes produziam
nomes públicos idênticos a menos do hash, porque compartilham a mesma categoria
de origem. `name_records` passa a qualificar o lote combinado com `ALL`
(`SRC[ALL]`, ou `SRC[<origem>-ALL]`), de modo que o conjunto completo se
distingue das entradas individuais na lista de resultados.

## C7c — recálculo redundante após reabrir o projeto

**Causa, que exigia investigação.** Nada a ver com a seleção de arquivos. Ao
reabrir um projeto criado antes das referências semânticas tipadas,
`semantic_references.migrate_pipeline_references` acrescenta `output_port` a
cada ligação legada — uma migração que a própria docstring descreve como não
destrutiva. Mas `stage_cache.stage_key` aplicava hash às ligações cruas, então o
acréscimo mudava a chave de **todo** consumidor e recalculava um pipeline inteiro
que não havia mudado. Medido no diagnóstico: a entrada resolvida era
byte-idêntica, e só o campo `bindings` diferia entre as duas chaves.

**Correção.** `stage_cache.request_bindings` remove `output_port` das ligações
antes do hash. A porta é um rótulo tipado calculado a partir da operação
produtora e do contrato da entrada consumidora (`artifact_contracts.output_port`),
nunca escolhido pelo usuário; a identidade do dado continua fixada por
`stage`/`asset`, pelo seletor e pelo manifesto de entradas resolvidas. Uma porta
que venha a ser selecionável pelo usuário terá de voltar à chave.

Este é o recálculo redundante medido na PHASE-00, e vale para qualquer projeto
reaberto, não só para o caso de vários arquivos publicados.

## C2 — invalidação por par no redocking

**Causa.** Curar os ligantes de um PDB reescreve o `pdb_codes.csv` compartilhado
por todos os pares daquela estrutura. Esse arquivo é uma das entradas
materializadas de cada lote, então `result_state.inputs_changed` passava a
acusar mudança em todos eles. A via parcial de `result_sync.plan`, que preserva
os arquivos dos pares sobreviventes com `pair_removals`, exigia justamente
`not inputs_changed(batch)` e nunca era alcançada: todo par era descartado e
redocado.

**Correção.** A comparação passou a dizer *quais* entradas mudaram, não apenas
se alguma mudou:

| Arquivo | Mudança |
| --- | --- |
| `src/biomolexplorer/stage_cache.py` | `selected_identity` virou função de módulo, para que o nome gravado no manifesto possa ser comparado fora do cálculo da chave. |
| `src/biomolexplorer/result_state.py` | `changed_input_names` devolve os nomes das entradas cujo conteúdo não corresponde mais, ou `None` quando a comparação é impossível (manifestos antigos sem versão, entrada ilegível), caso que continua tratado como mudança total. |
| `src/biomolexplorer/result_sync.py` | `pair_index_only` aceita a via parcial quando o índice de pares é a única entrada alterada. `pair_removals` então decide par a par, a partir do índice gravado com os resultados. |

A decisão por par continua sendo a de `pair_removals`, que já verificava o
digest do `redocking_cache.json` e preserva apenas os arquivos de pares ainda
selecionados e com a mesma preparação. A chave do par
(`PDB|LIGANTE|RESNUM|CADEIA`) carrega toda a requisição científica do par, então
uma edição nas linhas de outro PDB não pode alterá-la. Nenhuma chave de cache
mudou: resultados de redocking já existentes continuam reaproveitáveis.

Verificado com uma execução real de redocking em dois PDBs
(`test_curated_redocking.CuratedRedockingExecutionTests`): 1ABC/LIG é docado uma
única vez, 1P0I/BUA entra depois da curadoria, e a recuperação de estruturas não
se repete.

## Aviso quando a execução não é retomada

Em dois pontos de `ui/app.py`, um `save_draft` recusado fechava a caixa e
retornava em silêncio, deixando a execução pausada — exatamente com a aparência
de uma execução que seguiu. Agora o usuário é avisado nos dois casos, na
configuração do bloco e na seleção de arquivos, com as mensagens acrescentadas
aos catálogos `en`/`pt` sem alterar as existentes. Coberto por
`test_file_selection.test_unsaved_selection_leaves_the_run_paused_and_says_so`.

## Pendências que seguem rastreadas

| Marcador | Teste | Resolve em |
| --- | --- | --- |
| C5 | `test_connection_matrix`, `test_redocking_receptor_outputs` | passo 3 (contrato das portas) |
| C8 | `test_redocking_chimera_repair` | passo 3 (religar ou remover) |
| I1 | `test_result_persistence`, dois testes de migração do cache v1 | PHASE-07 (#14) |
| I2 | `test_analysis_origins` | PHASE-07 (#14) |

Fora de #6, a investigação aberta no passo 1 continua: a topologia ou os tipos
atômicos MOL2 respondem por 694 compostos excluídos no projeto Alzheimer.

Uma otimização identificada e **não** feita aqui: no modo `both`, o lote
combinado de fingerprints poderia ser montado a partir dos lotes individuais,
sem worker, já que os fingerprints são por molécula. Isso é execução
incremental, pertence a REQ-REUSE e às fases PHASE-03/04 de #6, e não a esta
correção.

## Reprodução

```sh
conda activate BioMolExplorer
export PYTHONPATH=src:tests PYTHONDONTWRITEBYTECODE=1
python -m unittest discover -s tests
```

O ambiente precisa estar ativado, e não apenas o seu interpretador invocado por
caminho: `test_molecular_core` e `test_ligand_preparation_performance` chamam
`obabel` por nome, então sem o `bin` do ambiente no `PATH` os dois acusam
`FileNotFoundError: 'obabel'` — uma falha de ambiente, não do código.
