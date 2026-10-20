# Ligantes revisados no PDB entram no redocking — 2026-10-19

## Causa reproduzida

O diálogo **Ligantes** gravava a escolha apenas em `pdb_codes.csv` da recuperação. A execução do redocking usa uma seleção separada, `parameters.pdb_codes`, e seus arquivos de entrada explícitos. Esses dois componentes do plano permaneciam anteriores à edição. Assim, uma nova combinação podia ser omitida antes mesmo de chegar ao cache incremental; o cache reaproveitava corretamente o pedido antigo, mas o pedido novo nunca era enviado.

O teste de regressão falhou antes da correção: após salvar um ligante de uma estrutura que não estava selecionada, a configuração persistida continha somente o par anterior.

## Correção

`ResultFiles.set_pdb_ligands` reconcilia a escolha explícita com os blocos de redocking conectados àquela recuperação. Acrescenta os pares e, quando necessário, o seletor do novo arquivo PDB. Retém pares de outras estruturas e suas configurações. Para um ligante novo do mesmo receptor/cadeia, conserva as opções compartilhadas de receptor e cofatores, sem copiar indevidamente as opções específicas de outro ligante.

A remoção de um ligante da curadoria retira somente a escolha correspondente do plano. A repetição da mesma edição não duplica pares. Uma escolha já salva pela versão anterior, mas ausente do plano, é reconciliada ao salvar a lista novamente. Estruturas sem nenhum ligante, inclusive metadados com uma linha vazia, permitem a inclusão do primeiro par.

Metadados, configuração do pipeline, revisão, manifestos de integridade e histórico são atualizados na mesma operação auditada. Uma falha restaura os metadados e desfaz as alterações do banco. Blocos sem ligação com a recuperação e resultados fornecidos pelo usuário são preservados.

O front end salva o rascunho antes da curadoria e atualiza o projeto após a confirmação. A atualização faz merge das alterações independentes no canvas, preservando a identidade do projeto e sua revisão. Isso evita exibir ou salvar novamente a seleção antiga depois de a nova já ter sido persistida. A mensagem orienta executar o pipeline para processar as combinações atualizadas.

## Integração e limites

O teste usa uma recuperação chamada **BChE PDB**, com um par anterior e uma estrutura identificada como **1P0I** sem ligante selecionado. Salva **BUA / cadeia A**, executa novamente e confirma:

- o novo par chega ao plano persistido e à seleção de arquivos;
- o redocking anterior não é recalculado;
- somente 1P0I/BUA é acrescentado aos cálculos;
- a tabela contém os resultados anterior e novo, uma vez cada;
- a recuperação PDB não é repetida.

São dados sintéticos com esses identificadores e resíduo 3. O teste não representa um download ou uma nova validação molecular do PDB real 1P0I. Armazenamento, curadoria, agendador, adaptador de operações, seleção de lotes, cache e serviço de resultados são reais; Chimera, Open Babel e Vina são simulados.

## Verificações executadas

Ambiente Python 3.12 do Conda BioMolExplorer, `PYTHONPATH=src:tests`, `MPLCONFIGDIR=/tmp/biomol-mpl`.

```text
python -m unittest test_curated_redocking test_retrieval_usability
  test_redocking_incremental test_redocking_configuration
  test_redocking_execution_review test_redocking_partial_failures
  test_stage_cache test_restart_reuse test_run_reuse_ui
  test_pipeline_editor_regressions test_localization -q
```

Resultado: **107 testes aprovados em 52,514 segundos**. Cobertura adicional inclui idempotência, reparação de escolha já salva, primeira escolha com metadados vazios, preparação compartilhada, remoção localizada, rollback e atualização do diálogo sem perder edições do canvas.

```text
python -m unittest test_results_display test_artifact_registry
  test_redocking_results test_canvas_drag_regression test_pipeline_png -q
```

Resultado adicional: **52 testes aprovados em 23,471 segundos**, totalizando **159 testes aprovados** nas duas baterias. `git diff --check` sem erros de espaçamento.

A documentação de configuração do redocking foi atualizada em português e inglês, com HTML regenerado e instrução para reconciliar escolhas feitas pela versão anterior.
