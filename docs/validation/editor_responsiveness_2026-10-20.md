# Auditoria de responsividade do editor — 20/10/2026

## Resultado e alcance

Foram corrigidos três gargalos reproduzíveis no editor Flet, sem alterar algoritmos científicos, APIs, identificadores ou o formato dos pipelines:

1. Reconstrução de todos os controles ao incluir ou alterar um bloco.
2. Leitura e validação molecular de CSVs históricos sem relação com a entrada selecionada ao abrir Vina/DOCK6.
3. Reconstrução da grade e atualização da página inteira a cada pixel de arraste na extremidade inferior do grafo.

Também foi protegida a abertura concorrente de configurações. As evidências e medições abaixo se referem ao código Python e aos callbacks reais dos controles, com persistência temporária real. O congelamento persistente da sessão gráfica relatada pelo usuário não foi reproduzido em um navegador conectado. Portanto, os resultados comprovam a remoção desses gargalos, mas não certificam isoladamente a eliminação de toda possível causa de travamento no cliente Flutter.

Revisão original: `4ed01455a4d7b2e7c0109443183b764e973a38dd`.
Ambiente: Linux x86_64, Python 3.12.3, Flet 1.0.3, ambiente Conda `BioMolExplorer` já instalado. Nenhuma dependência nova.

## Fluxo investigado e diagnóstico

### Canvas, eventos e estado

`FlowCanvas.add()` cria metadados com `new_stage()`, encontra uma posição livre, registra undo, inclui o bloco e chama `changed()`. Não foi encontrado disparo de docking ou criação de worker científico nesse caminho. `catalog.operation_fields()` já usa AST e cache limitado; não foi necessário substituir esse mecanismo.

Antes, `FlowCanvas.refresh()` recriava todos os nós e as opções de seleção. A instrumentação de `node()` confirmou 100 construções ao incluir o centésimo bloco. Essa substituição também renovava callbacks e componentes gráficos não afetados. O código já preservava os controles de saída de arraste por UUID; essa proteção foi mantida.

`node.move()` já enviava apenas nó e arestas quando a grade não mudava. Contudo, a altura da grade era calculada diretamente a partir da coordenada do bloco mais abaixo. Em 100 blocos, um deslocamento de um pixel mudava a altura, reconstruía as linhas da grade e caía em `page.update()` global. O teste da extremidade inferior confirma agora a mesma instância de grade e atualizações restritas durante 40 deslocamentos consecutivos.

Seleção e configurações seguem caminhos distintos: a seleção muda bordas; duplo clique, botão de configuração e configuração da seleção abrem o formulário pelo UUID. As validações transacionais de conexões em `flow.connect()` e `pipeline.validate_pipeline()` foram preservadas.

### Abertura das configurações e trabalho molecular

`WorkspaceUI.open_stage_dialog()` consulta acesso ao projeto, arquivos e execuções, monta os catálogos e constrói `GuidedForm`. Antes, tanto Vina quanto DOCK6 percorriam todos os CSVs dos resultados mais recentes de cada produtor e dos arquivos enviados, antes de mostrar o formulário, inclusive quando o novo bloco não tinha conexão.

`pubchem_retrieval.reference_choices()` lê tabelas com pandas e chama `molecule_quality.clean_dataframe()`, que valida/canonicaliza SMILES com RDKit. Esse trabalho já era encaminhado por `WorkspaceUI.call()` para `asyncio.to_thread()`: não foi diagnosticado como processamento síncrono na thread de eventos. O problema confirmado era sua execução antecipada e repetida, com espera de cada arquivo antes de abrir o painel.

O cenário controlado com cinco CSVs de mil linhas confirma cinco leituras científicas por abertura na versão original e zero na versão atual para Vina sem conexão. Escolher uma origem continua carregando e validando seus compostos. Entradas previamente configuradas continuam disponíveis na primeira renderização do formulário.

A versão anterior permitia sobreposição das aberturas, com `selected`, `artifact_items` e catálogos compartilhados modificados por corrotinas distintas. Um teste com consulta atrasada e nova solicitação comprova agora que apenas a solicitação mais recente publica um formulário.

### Persistência, execução e tarefas adjacentes

Foram examinados o salvamento com revisão e merge (`_save_draft()`), debounce de 600 ms, polling de projeto a cada dois segundos e polling de execução a cada segundo. O salvamento já possui lock e proteção contra cancelamento de uma gravação iniciada; as respostas antigas de polling já são rejeitadas por revisão. Essas proteções foram conservadas e seus testes executados.

`PipelineService` usa executor limitado, estado protegido por lock, workers de jobs e espera por evento de 300 ms no monitoramento. A inclusão de blocos não passa por `submit()`. Não houve evidência que justificasse eliminar threads ou substituir o agendamento científico. Nenhuma alteração foi feita nesses algoritmos ou parâmetros.

Histórico undo/redo continua limitado a 50 snapshots e pipelines a 100 etapas. Não se removeu validação de dependências ou descoberta persistida de resultados para acelerar a interface.

## Modificações e decisões

| Arquivo | Alteração e finalidade |
|---|---|
| `src/biomolexplorer/ui/flow_canvas.py` | Reutiliza controles não afetados; invalida por identidade do objeto, conteúdo, posição efetiva, rótulos de entrada e seleção. Atualiza somente as bordas anterior/nova e seletor ao selecionar. Reutiliza opções sem alterações de IDs/nomes. |
| mesmo arquivo | Dimensiona grade em faixas de 400 unidades e evita reduzi-la durante o arraste, eliminando reconstruções por pixel. Expansão necessária ainda atualiza o canvas. |
| `src/biomolexplorer/ui/compound_catalog.py` | Separa descoberta de origens de leitura molecular; carrega somente origens solicitadas. Cache de sessão limitado a 32 arquivos e 20.000 IDs; arquivos maiores continuam funcionais, sem entrar nesse cache. |
| mesmo arquivo | Invalida por caminho, tamanho, `mtime_ns` e `ctime_ns`; verifica mudanças durante a leitura. Mutações de cache ocorrem no event loop; resultados obsoletos não repovoam a sessão. |
| `src/biomolexplorer/ui/app.py` | Indexa caminhos autorizados, pré-carrega somente as entradas de docking já configuradas e protege publicação por geração, token e identidade do projeto. Aberturas são serializadas; solicitações ultrapassadas não publicam formulários. Limpa o cache ao sair do projeto. |
| `src/biomolexplorer/ui/input_editor.py` | Escolha de origem carrega suas referências de forma assíncrona. Uma seleção posterior invalida a atualização visual da escolha anterior. Preserva seleção de IDs moleculares e arquivos. |
| `tests/test_editor_responsiveness.py` | Onze testes novos para reprodução controlada, identidade de controles, cache, concorrência, formulários e estresse. |
| `tests/test_pipeline_execution.py` | Teste de configuração/retomada seleciona pela identidade persistida do arquivo, em vez do rótulo visual antigo. O rótulo já incluía origem e alvo na revisão original. |
| `scripts/benchmark_editor_interactions.py` | Comparação reproduzível com mediana, p95, contagem de nós e leituras de CSVs. |

A invalidação por identidade do objeto é intencional: quando restauração/merge substitui uma etapa por um dicionário equivalente, seus callbacks devem passar a editar a nova instância. Os controles de saída de arraste permanecem estáveis por UUID. Remoção limpa nós, snapshots visuais e controles de saída; undo/redo continua reconstruindo corretamente as instâncias restauradas.

O caminho legado DOCK6 `base_vina_path` continua considerado no pré-carregamento, antes da migração já existente no formulário. Não foram criadas representações persistidas novas.

A estrutura permanece incremental: não houve reescrita de framework, layout, backend ou execução. Arestas ainda são recalculadas no arraste; snapshots e rótulos ainda são comparados nas alterações estruturais. Não se alega atualização estritamente O(1) do editor.

## Medições

Os números completos, incluindo medianas e p95 por cenário, estão em [editor_responsiveness_2026-10-20.json](editor_responsiveness_2026-10-20.json).

O script carrega `app.py` e `flow_canvas.py` originais pelo Git e compara com as versões de trabalho no mesmo processo. Os demais módulos são compartilhados entre as variantes; as escolhas assíncronas de origem não entram no cenário cronometrado de formulário vazio. Usa SQLite real em diretório temporário, controles reais Flet e fachada de página sem socket/Flutter. Como nos testes existentes, as chamadas de armazenamento são despachadas diretamente: tempos de executor, rede, rasterização e resposta de clique no navegador não são medidos.

São 21 amostras por cenário, com uma rodada de aquecimento excluída para cada tamanho. p95 usa o método nearest rank. O cenário histórico usa cinco CSVs de mil registros, sem executar docking. Percentis em uma amostra curta podem oscilar; não foram substituídos por estimativas.

| Blocos após inclusão | Inclusão original, mediana / p95 (ms) | Atual, mediana / p95 (ms) | Nós construídos, original → atual |
|---:|---:|---:|---:|
| 1 | 1.512 / 1.672 | 1.515 / 1.727 | 1 → 1 |
| 10 | 6.432 / 6.995 | 3.149 / 29.715 | 10 → 1 |
| 25 | 14.748 / 15.202 | 6.140 / 6.567 | 25 → 1 |
| 50 | 29.444 / 49.430 | 11.911 / 12.497 | 50 → 1 |
| 100 | 64.705 / 91.008 | 27.295 / 28.014 | 100 → 1 |

| Cenário | Original, mediana / p95 (ms) | Atual, mediana / p95 (ms) |
|---|---:|---:|
| Abrir Vina sem histórico, 100 blocos | 6.021 / 31.662 | 6.110 / 6.470 |
| Mover bloco inferior, 100 blocos | 15.860 / 16.309 | 1.123 / 1.187 |
| Abrir Vina sem conexão, cinco CSVs históricos | 129.503 / 151.330 | 6.110 / 6.839 |

O ganho não foi uniforme em todos os percentis: nesta amostra, o p95 da inclusão em dez blocos aumentou de 6,995 ms para 29,715 ms, apesar da redução da mediana. Esse pico não foi atribuído a uma causa medida. Abrir um formulário sem histórico já era rápido e não apresentou redução relevante na mediana.

No cenário histórico, as leituras de referência caíram de cinco por abertura para zero. A adição constrói somente o novo bloco quando a seleção anterior está vazia; se a seleção muda de um bloco existente, esse bloco também pode precisar atualizar sua borda. Os valores para todos os tamanhos e interações estão no JSON.

Não foram medidos FPS, tempo até pintura/interatividade no navegador, CPU percentual, RSS ou crescimento de heap em uma sessão gráfica de horas. O estresse automatizado verifica limites e identidade de objetos Python, não prova ausência de todo vazamento de memória. A ausência de leituras históricas não necessárias e a contagem de construções são verificadas diretamente, independentemente de ruído nos tempos.

## Validação

- Inclusão de Vina seguida imediatamente de abertura do formulário.
- Abertura/fechamento de todos os 15 tipos de blocos, duas vezes seguidas.
- Histórico de compostos não lido para docking sem origem; escolha posterior carrega IDs reais; origem já configurada recupera esses IDs ao reabrir.
- Invalidação de CSV alterado, limite do cache e descarte de leitura obsoleta.
- Solicitações concorrentes publicam somente o formulário mais recente.
- Restauração de objetos equivalentes reatribui os callbacks de movimento ao estado vivo.
- Remoção/undo libera/restaura os controles esperados.
- Dez ciclos de conexão/desconexão em cada um dos tamanhos 10, 25, 50 e 100, preservando controles não afetados; renomear produtor atualiza o rótulo do consumidor.
- Estresse com 100 blocos: 50 ciclos de remoção/inclusão e 150 deslocamentos, com contagens limitadas e histórico de até 50 snapshots.
- Preservação de controles de origem e aceitação de arraste atrasado, conexão/desconexão, descoberta de resultados persistidos e reutilização de produtor executado.
- Edição/persistência de parâmetros, seletores por identidade, entradas diretas, merge de salvamento concorrente e proteção contra polling obsoleto nos testes existentes.
- Execução de dependências, configuração de etapa pendente e retomada pelo formulário, com workers controlados e armazenamento real.
- `compileall` e `git diff --check`.

A rodada direcionada passou com **100 testes** em 32,524 s; o módulo de execução passou separadamente com **22 testes** em 20,230 s. Esses números incluem testes sobrepostos à suíte geral e não devem ser somados como casos distintos.

A rodada final da suíte completa, com `PATH` do ambiente Conda e cache Matplotlib em `/tmp`, terminou com:

```text
Ran 895 tests in 333.269s

FAILED (failures=3, errors=1, skipped=2)
```

A rodada isolada final do módulo novo passou com **11 testes, todos aprovados**, em 7,009 s, incluindo rótulos de produtores e estresse de conexões. Esses casos se sobrepõem à suíte geral. A asserção adicional de espaço disponível no arraste também passou isoladamente.

As quatro ocorrências remanescentes foram reproduzidas **também na revisão original**, extraída para diretório temporário, com as mesmas dependências:

| Teste | Resultado também observado na revisão original |
|---|---|
| `test_network_refinement.RefinedViewerTests.test_mcs_persisted_result_translates_with_dynamic_counts` | `KeyError: 'ids'` em `molecular_inspector.show_result()`. |
| `test_preparation_settings.CandidatePreparationTests.test_exported_candidates_are_reused_without_conversion` | Esperava `ValueError` para o formato solicitado, mas a chamada não o lançou. |
| `test_ranking_scopes.CentralityFeasibilityTests.test_betweenness_sampling_is_repeatable_and_covers_nontrivial_components` | Divergência numérica na comparação da centralidade amostrada; o valor exato varia entre execuções. |
| `test_workspace.WorkspaceTests.test_catalog_matches_all_operations` | A expectativa inclui o parâmetro legado `base_vina_path`, ausente nos campos DOCK6 do catálogo. |

Não se alteraram algoritmos científicos nem módulos alheios ao editor para ocultar essas falhas. Elas permanecem documentadas; a suíte global não está totalmente aprovada.

Na primeira execução ampla, sem o `bin` Conda no `PATH`, houve 886 testes, seis falhas, três erros e dois ignorados. Além das quatro ocorrências acima, Open Babel não era encontrado e um teste de retomada dependia de um rótulo antigo. Os dois problemas também foram reproduzidos no checkout original. Com `PATH` corrigido, os testes de consenso/DOCK6 afetados passaram. O teste de retomada foi atualizado para selecionar pela identidade do arquivo e passou até concluir a execução, fechar o formulário e conferir a configuração persistida.


As execuções controladas não equivalem a executar um protocolo científico Vina/DOCK6 real nem validar gestos em um navegador. A validação final de fluidez na sessão gráfica e com seus dados históricos continua necessária; não se declara docking científico ou ausência de todo bloqueio persistente como validado por esses testes.

## Reprodução

Com o ambiente Conda existente ativado (incluindo seu diretório `bin` no `PATH`):

```bash
PYTHONPATH=src python scripts/benchmark_editor_interactions.py \
  --baseline 4ed01455a4d7b2e7c0109443183b764e973a38dd --repeats 21
PYTHONPATH=src:tests python -m unittest test_editor_responsiveness \
  test_pipeline_editor_regressions test_canvas_drag_regression \
  test_pipeline_execution
PYTHONPATH=src python -m unittest discover -s tests
```

O benchmark só cria diretórios temporários e não altera pipelines do usuário.
