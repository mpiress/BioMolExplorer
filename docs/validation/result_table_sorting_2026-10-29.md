# Ordenação de resultados — 29/10/2026

## Diagnóstico e correção

Os logs de frontend registraram `each visible DataRow must contain exactly as many
visible DataCells as there are visible DataColumns (5)`, no caminho
`DockingResultsTable.load → CompoundTableViewer.load → page.update`.
A implementação de Docking chamava o carregamento da tabela genérica, que publicava
linhas de quatro células enquanto o componente já tinha cinco colunas de Docking.
Isso quebrava a atualização do Flet e comprometia a visualização. Não havia um link
de cabeçalho nem um submit de formulário responsável pelo erro.

O carregamento agora delega a renderização ao tipo concreto da tabela. Colunas e
linhas são construídas sem publicação intermediária, antes de uma atualização
consistente. O mesmo caminho atende Vina, DOCK6, consenso e a inspeção genérica de
CSVs. Os testes executam `DataTable.before_update()` em cada atualização, em vez de
simular `page.update()` como uma operação vazia; assim detectam a falha original.

## Padrão compartilhado

`table_sort.py` centraliza comparação estável, numérica ou textual sem distinção
de maiúsculas, com valores ausentes/não finitos no final em ambas as direções.
`ui/table_sort.py` centraliza o clique, a alternância e o indicador de direção.
As colunas de ações permanecem sem ordenação, pois contêm botões, não valores.

| Visualização de resultados | Ordenação |
|---|---|
| Vina / DOCK6 / consenso | Todas as colunas de dados, antes da paginação |
| Compostos ChEMBL / PubChem / expansão | Código e SMILES; demais campos pela inspeção completa |
| Redocking | Target, PDB, ligante, resíduo, cadeia, estado e RMSD |
| Arquivos publicados | Nome exibido e tamanho numérico, preservando filtro de origem |
| Laboratório de grafos | Colunas de nós e arestas, mantendo inspeção molecular, filtros e exportação |
| Contatos geométricos | Resíduo, cadeia, distâncias e pares de átomos, sem recarregar a cena |
| Relatório de preparação | Campos do relatório, com paginação |
| CSVs de ADMET, alertas, fingerprints, similaridades e demais etapas | Inspeção somente leitura pelo botão Ver tabela |

Resultados tabulares são ordenados antes de selecionar a página. O índice original
de cada linha e a versão do arquivo são preservados para ações de remoção, poses,
contatos e footprint. Filtros, seleção do arquivo, projeto e aba não são reiniciados.
A inspeção genérica recusa remoção também no serviço, não apenas na interface.
Permite consultar CSVs publicados de execuções concluídas ou interrompidas, inclusive
relatórios de falha, mantendo a autorização por projeto, execução e etapa.

O serviço mantém um snapshot por instância para CSVs de até 8 MiB e 50.000 linhas.
O cache é invalidado por caminho/dispositivo/inode/tamanho/mtime/ctime. Arquivos
maiores mantêm leitura em streaming e a ordenação global necessária, sem limitar
incorretamente a ordenação à página visível. A autorização é revalidada nas consultas.
Não são chamados provedores remotos nem refeitos cálculos científicos ao ordenar.

## Validação e integridade

- Suíte de regressão: **199 testes aprovados em 148,820 s**, incluindo persistência,
  reutilização, portabilidade, recuperação, grafos, resultados e ações de docking.
- Verificação final: **30 testes aprovados em 15,156 s**, incluindo schemas genéricos,
  CSV de etapa com falha, leitura somente, invalidação do snapshot e traduções.
  As suítes têm sobreposição; os totais não são testes distintos.
- Chrome real, largura 1000 px: quatro cliques em `score`, alternando ascendente e
  descendente. O mesmo projeto e a aba Execuções permaneceram ativos. O primeiro
  composto alternou entre M40 (-40) e M01 (-1), sem erros JavaScript. A contagem de
  leituras completas/hashes do CSV permaneceu em uma e seu SHA-256 não mudou.
- Ordenação de RMSD e demais campos sem consulta adicional de simulações; filtro
  textual e projeto preservados. Filtro de origem preservado na tabela de arquivos.
- Assinaturas científicas comparadas com os contratos de compatibilidade existentes:
  nenhuma alteração. A mudança não invalida os caches de processamento.
- Dois testes antigos de exclusão de falhas tiveram a montagem de fixtures adaptada
  à confirmação exigida pelo protocolo atual, sem relaxar a proteção de produção.
- Nenhuma biblioteca instalada; testes destrutivos somente em projetos temporários.
  Resultados originais do usuário não foram removidos ou recalculados.

[Logs, relatório do navegador e captura](result_table_sorting_2026-10-29/).
