# Gerenciamento de artefatos: implementação e verificação

## Diagnóstico da arquitetura original

O executor já separava diretórios por execução, bloco e lote individual. Conexões
admitiam múltiplas referências `stage`/`asset`, e seletores qualificados evitavam
ambiguidades dentro de uma origem. Faltava um catálogo explícito que relacionasse
cada arquivo ao produtor, à execução e aos arquivos efetivamente consumidos.
A janela de handoff renderizava um checkbox por arquivo. As configurações podiam
crescer com uma linha por seleção, e tabelas independentes de docking podiam ser
reduzidas a uma única opção na passagem para outro bloco.

## Identidade e rastreabilidade

`src/biomolexplorer/artifacts.py` mantém a representação compartilhada:

- `artifact_id`: SHA-256 do identificador do produtor e do caminho lógico da saída,
  relativo ao produtor nos diretórios nativos. Não depende do nome apresentado,
  do diretório da execução ou da identificação transitória do worker.
- `occurrence_id`: identifica o artefato e a execução que o produziu. Reutilização
  conserva a execução de origem.
- Metadados: nome/identificador/operação do produtor, execução, nome original,
  seletor legível, localizador lógico, caminho físico, alvo/consulta e entradas.

O seletor pode ganhar uma qualificação adicional se aparecer outro arquivo com o
mesmo nome; isso não altera a identidade de uma saída nativa existente. A resolução
usa o ID quando presente, verificando a disponibilidade, e continua aceitando
seletores antigos. Ausência, ambiguidade e identidade inválida impedem o consumo.

## Persistência e consumo

O índice `artifact_records` é persistido junto das etapas no JSON das execuções
já armazenado em SQLite, sem migração de tabelas. Os bindings existentes recebem
um campo opcional `artifact_id`; `stage`, `asset`, `selector` e `compound_id`
continuam compatíveis. Uploads conservam seus identificadores de assets.

Cada saída registra as entradas efetivas identificadas por `input_files`.
No processamento individual, o índice usa os arquivos e a configuração do lote
correspondente. Em merge, registra todas as origens que participaram da saída.
Uploads diretos são associados pelo caminho autorizado e pelo ID do asset.
O cache preserva o registro de procedência dos resultados reutilizados.

Nomes, conteúdo e colunas de CSVs científicos permanecem intactos. As adaptações
atingem descoberta, seleção, resolução e metadados; não alteram os algoritmos
científicos. A resolução compartilhada também atende à identificação dos arquivos
no pareamento de entradas de DOCK6.

## Interface

`ui/artifact_picker.py` é reutilizado por `FileSelection` e `InputEditor`, tanto
na configuração guiada quanto nas pausas de execução. Oferece origem conectada,
pesquisa, arquivos disponíveis, Add File, Add All, Add All Sources, remoção e Clear
Selection. Seleções de outras origens são preservadas. As ações em lote abrangem
toda a origem, independentemente da busca ou da página atual.

Listas disponíveis e selecionadas exibem 25 itens por página; os menus têm altura
limitada. A inclusão em lote atualiza a página uma vez, e itens inválidos continuam
visíveis para serem removidos ou substituídos. O modo individual/merge e as escolhas
opcionais de composto para docking são preservados.

O canvas identifica o produtor nas entradas. Tabelas gerais de arquivos mostram
procedência e filtro por origem. Seletores de tabelas de compostos e docking usam
rótulos com produtor, contexto e caminho qualificado, mantendo o nome original
separado para downloads e integrações. A descoberta de docking reaproveita
`result_tables` para preservar cada resumo independente e omitir cópias redundantes
por receptor.

## Verificação

A suíte ampla de execução, grafos, configuração, seleção, cache, portabilidade,
docking, resultados, localização e imports passou com 170 testes. Uma rodada final
com 119 testes passou; os nove testes do registro de artefatos passaram novamente
após o ajuste de preservação do contexto efetivo dos lotes na redescoberta e no cache.

`tests/test_artifact_registry.py` cobre identidades entre produtores/execuções,
colisão de nomes adicionada posteriormente, IDs inválidos, arquivos ausentes,
processamento individual e merge, persistência, troca de conexão e múltiplas tabelas
de docking. O seletor foi exercitado com 202 arquivos físicos e com 5.000 itens,
renderizando apenas 25 itens selecionados e emitindo uma atualização no lote.

A renderização Flet foi inspecionada no Chrome em 1280 × 1000 e 390 × 844, sem
exceções JavaScript. A documentação bilíngue foi regenerada; o validador confirmou
32 páginas HTML com arquivos, âncoras e IDs válidos. `compileall` e `git diff --check`
também passaram.

As integrações usam workers controlados para verificar execução e resolução dos
arquivos. Não repetem cálculos nativos de docking nem dependem das APIs externas.
Os fixtures antigos de pausa com seleção já explícita e de preparação com caminhos
inexistentes foram atualizados para os contratos já vigentes, sem modificar esses
protocolos científicos.

## Compatibilidade e limites

- A procedência é registrada por arquivo/lote, não por molécula ou linha de CSV.
- Resultados históricos sem índice são descobertos a partir do manifest existente.
  Sua procedência detalhada depende de terem conservado `input_files` e lotes.
- Arquivos legados fora de uma árvore de produtor conhecida usam localização
  qualificada; se forem movidos manualmente, uma nova seleção pode ser necessária.
- O catálogo completo permanece em memória, enquanto a renderização é paginada.
  Pesquisa percorre esse catálogo e validação de CSVs continua sujeita ao custo
  de descoberta dos arquivos.
- A assinatura de implementação existente pode invalidar caches após atualizar o
  software; os resultados históricos permanecem disponíveis. A reutilização e a
  migração são verificadas com a mesma implementação instalada.
- Versões antigas da aplicação não entendem bindings novos com `artifact_id`;
  projetos antigos continuam aceitos pela versão atual.
