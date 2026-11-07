# Auditoria de I/O de resultados — 2026-11-06

## Escopo e diagnóstico

Inspeção do paginador compartilhado (CompoundTables/ResultTables/DockingResults), catálogo de resultados publicados, ResultFiles, WorkspaceStore e publicação atômica em storage.py. O fluxo combina SQLite para metadados e CSV/JSON para artefatos científicos. Os índices runs_project e assets_project já existiam; não se atribui a lentidão a uma ausência generalizada de índices.

| Gargalo encontrado | Ação corretiva | Impacto esperado |
| --- | --- | --- |
| compound_tables.py, CompoundTables.page: mesmo com CSV em memória, cada página recriava registros/links para todos os resultados correspondentes e repetia a ordenação | Uma única visão em memória, limitada ao snapshot existente (até 8 MiB/50 mil linhas), identificada pelo conteúdo, busca e ordenação; construir respostas e links somente para a página | Para navegação na mesma visão: trabalho sobre linhas de O(n) ou O(n log n) para O(p), com p até 100. Primeiro filtro/ordenação mantém O(n)/O(n log n). Espaço adicional O(n) de referências; somente uma visão retida |
| published_results.py, verify: conexão e SELECT por asset referenciado | Uma conexão e consultas parametrizadas em lotes de até 500 IDs, projetando somente id e path e restringindo por projeto | De a conexões/consultas para uma conexão e ceil(a/500) consultas. Hashes continuam obrigatórios; espaço adicional O(a) |
| workspace.py, inicialização: consulta de execução ativa usava índice de projeto e precisava filtrar histórico | Índice parcial runs_active_project, criado também em bancos existentes | Busca restrita ao subconjunto ativo; evita percorrer todo o histórico do projeto. Pequeno custo de manutenção do índice em transições de status |

## Integridade e recursos

Autorização, validação do caminho, assinatura stat do CSV e consulta de permissão de edição continuam em cada chamada. Alterações no CSV invalidam o snapshot; respostas recebem cópias dos dicionários para evitar corromper o cache. Arquivos maiores mantêm leitura em fluxo e não são retidos pela nova visão. Ordenação global de arquivos grandes ainda pode exigir O(n) de memória. Hashes de integridade científica não foram substituídos por cache persistente. Conexões usam o context manager existente, que fecha o SQLite em finally. Formatos CSV/JSON, APIs, publicação atômica, fsync e rollback foram preservados.

## Medição

Comparação local entre CompoundTables da revisão HEAD anterior à alteração e a implementação atual, no mesmo processo, usando WorkspaceStore temporário, CSV de 40.000 linhas (ID, SMILES e score), páginas de 25 linhas ordenadas por score. Uma carga inicial e mediana de dez navegações subsequentes. Tempos incluem autorização e acesso SQLite; excluem renderização Flet e rede.

| Operação | Antes | Depois |
| --- | ---: | ---: |
| Primeira página ordenada | 185,6 ms | 68,0 ms |
| Mediana de navegação ordenada | 140,8 ms | 5,0 ms |

Redução observada de aproximadamente 96% na navegação (28 vezes). É um microbenchmark sintético, não uma medição de desempenho completo em projetos reais. Mudança de busca, ordenação ou conteúdo exige reconstruir a visão.

## Alterações

- src/biomolexplorer/compound_tables.py: inicialização de _view e processamento de página (linhas 43 e 135–168).
- src/biomolexplorer/published_results.py: consultas de assets em lotes (linhas 92–108).
- src/biomolexplorer/workspace.py: índice parcial (linhas 51–52).
- tests/test_result_io_performance.py: regressões de reutilização de ordenação, criação de links, isolamento de resposta, invalidação, autorização, leitura em fluxo, plano SQLite e batching.

## Limites e trabalho posterior

Ainda há varreduras e hashes no catálogo de resultados publicados, além de parsing completo em alguns visualizadores ADMET e de poses. São candidatos a perfilamento em um projeto real. Não foi introduzido paralelismo na escrita científica: a durabilidade das publicações atômicas e das transações de curadoria exige preservar a ordem existente.

## Validação executada

- 16 testes de cache/paginação, ordenação e resultados publicados: aprovados (`test_result_io_performance`, versão inicial com três testes, `test_published_results`, `test_result_sorting`).
- Suíte final de cinco regressões específicas de I/O: aprovada, incluindo lotes de 501 assets e rejeição de assets ausentes.
- Execução ampliada de apresentação/persistência e testes de tabelas herdados: 72 testes, 71 aprovados e uma falha de expectativa do nome exibido (`test_native_table_buttons_show_generated_views_and_enforce_read_only`). A mesma falha foi reproduzida isoladamente substituindo CompoundTables pela versão HEAD anterior à mudança.
- `git diff --check`: aprovado.

A validação não inclui navegação renderizada em projeto científico real nem benchmark de escrita em armazenamento remoto. Exclusões de arquivos temporários `.insyncdl` já estavam presentes no início e não foram alteradas por este trabalho.
