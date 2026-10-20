# Exploração dos grafos e recursos do Gephi

O botão **Restaurar grafo completo** remove simultaneamente os filtros de
comunidade, componente, identificador, grau, centralidade, escore, vizinhança e
seleção molecular. Também desliga a agregação e ajusta a câmera. A troca entre
Full Graph e MCC usa a mesma limpeza. **Limpar filtros / voltar à rede** mantém
o universo escolhido (Full ou MCC), mas remove todas as suas restrições.

Redes acima de 3.000 nós ou 12.000 arestas usam um fundo raster e uma camada
vetorial para seleção/vizinhos. Todos os nós têm coordenadas e permanecem no
índice de cliques. A cena tem até 4.096 px por dimensão; zoom, busca e inspeção
continuam disponíveis. A cena e a imagem são atualizadas antes das chamadas de
câmera ao retornar de um subgrafo. Coordenadas ausentes de um layout de outro
universo são completadas pelo layout integral salvo.

## CSVs por componente e entradas de docking

Ao executar **Filtrar por grafos / Filter by graphs**, cada análise passa a
produzir tabelas moleculares independentes na pasta `Molecules`:

| Arquivo | Compostos incluídos |
| --- | --- |
| `<análise>__full_graph.csv` | Grafo completo, incluindo nós isolados |
| `<análise>__component_0001_mcc.csv` | MCC / LCC: maior componente conectado, identificado como componente 1 |
| `<análise>__component_0002.csv`, etc. | Somente os compostos do componente indicado, inclusive componentes unitários |

A numeração corresponde à visualização daquela análise. Cada CSV conserva
`molecule_chembl_id`, `canonical_smiles` e as demais propriedades moleculares,
acrescentando `degree`, `graph_analysis_id`, `graph_component`, `graph_is_mcc`
e `graph_scope`. A visualização salva relaciona cada subconjunto ao seu arquivo
em `compound_exports`. Análises de entradas diferentes continuam separadas.

Os arquivos também aparecem em **Arquivos desta análise**. Nas configurações de
Vina ou DOCK6, conecte a saída do bloco de grafos à entrada de compostos e escolha
os arquivos desejados no seletor. Os rótulos distinguem **Grafo completo**,
**MCC / LCC · Componente 1** e **Componente N**, mantendo a análise/origem e o
caminho do arquivo. É possível escolher um ou vários subconjuntos; a seleção
salva conserva a identidade dos arquivos ao reabrir o pipeline.

Tabelas de arestas (`source,target,value`) continuam disponíveis nos resultados,
mas não são apresentadas como saídas moleculares do bloco. Um grafo calculado
somente com arestas externas precisa da tabela de compostos com SMILES para
que seus subconjuntos possam ser usados em docking.

Os CSVs anteriores `<análise>.csv` e `molecules.csv` são preservados para
compatibilidade: o primeiro contém o MCC daquela análise; o segundo é a união
dos MCCs, **não** o grafo completo. Seus rótulos deixam isso explícito.
Para obter os novos subconjuntos em projetos antigos, execute novamente o bloco
de grafos; resultados já salvos não são reescritos retroativamente.

## Recursos disponíveis

Os recursos foram comparados com o [guia oficial do Gephi](https://gephi.org/quickstart/)
e o [catálogo de layouts](https://docs.gephi.org/Plugins/layout).

| Recurso | BioMolExplorer |
| --- | --- |
| Overview / navegação | Full, MCC, componentes, comunidades, ego de 1–3 passos, busca, zoom, ajuste, tela ampliada, investigação 2D/3D |
| ForceAtlas2 | Implementação CPU determinística, atração linear ou LinLog, massas grau + 1, gravidade e velocidade adaptativa; repulsão exata; 80 iterações |
| Fruchterman–Reingold | Rede até 1.500 nós / 20.000 arestas, ou cálculo por componente com empacotamento |
| Outras organizações | Circular/anéis, concêntrica por grau, Random (semente 42), Kamada–Kawai, BFS exploratória e grupos Louvain |
| Aparência / ranking | Cor por métrica ou grau; tamanho por métrica ou uniforme; Viridis, verde–amarelo e azul–vermelho |
| Aparência / partições | Cor por comunidade, componente ou origem, IDs e contagens na legenda paginada |
| Preview | Exibição de códigos, arestas e contornos dos componentes; seleção destaca vizinhos |
| Filtros | Identificador, componente, comunidade, intervalo de grau, centralidade (incluindo k-core), escore, ego e propriedades moleculares |
| Estatísticas e ranking | Grau, centralidade de grau, betweenness, closeness, eigenvector, PageRank, k-core, clustering e triângulos |
| Estatísticas globais | Densidade, componentes, isolados, grau médio, clustering médio, modularidade, assortatividade, diâmetro e caminho médio do LCC |
| Data Laboratory | Tabelas de nós e arestas, busca, ordenação, 50 registros por página, investigação por composto e CSV completo da seleção pesquisada |
| Exportação | CSV, relatório JSON, SVG e PNG com legenda; GEXF e GraphML com propriedades e métricas calculadas |
| Comunidades | Louvain sem pesos, resolução 1, semente 42; agregação preserva todos os grupos e abre os compostos ao clicar |

ForceAtlas2 segue as leis de força descritas por
[Jacomy et al. (2014)](https://doi.org/10.1371/journal.pone.0098679).
Esta implementação tem execução finita e não usa o motor interativo nem a
aceleração Barnes–Hut do Gephi. Os layouts por forças trabalham por componente
até 1.500 nós / 20.000 arestas, com orçamento total de 8 milhões de pares por
iteração. Componentes acima desse orçamento usam anéis e exibem um aviso. Nenhum
composto é amostrado ou retirado. As posições não serão idênticas às do Gephi.

Os cálculos de centralidade mantêm os limites informados no painel e na ajuda
do ranking. Louvain é permitido até 100.000 nós / 500.000 arestas, inclusive nas
redes do experimento Alzheimer que excediam o limite anterior. Triângulos e
clustering por nó exigem soma dos graus ao quadrado até 10 milhões. Diâmetro e
caminho médio descrevem explicitamente o LCC, até 1.500 nós / 20.000 arestas.
Métricas indisponíveis aparecem como `—` e recebem cor cinza.

## Interpretar cores e exportar

A barra contínua mostra zero, ponto médio e máximo da métrica **no universo
antes dos filtros**, para conservar a comparação ao abrir comunidades. Na
paleta verde–amarelo, verde representa valores baixos e amarelo o máximo. Na
Viridis, os valores baixos são roxos. Tamanho e cor têm seletores independentes.
O contorno azul indica pertencimento ao MCC. Cores de partições são categorias,
sem ordenação de importância; quando a paleta se repete, os IDs distinguem os
grupos. Para mais de 90 categorias, a imagem exportada mostra as primeiras 90;
a interface tem paginação e os dados completos ficam nos relatórios/arquivos.

SVG e PNG exportam toda a rede desenhada, mesmo quando a interação usa fundo
raster. GEXF e GraphML exportam o subgrafo molecular filtrado, preservam o escore
como atributo e não o convertem em distância ou peso topológico. GEXF inclui
posição, cor e tamanho no [formato de visualização do Gephi](https://docs.gephi.org/desktop/User_Manual/Import/GEXF_File_Format/).
GraphML inclui os mesmos dados como atributos; a aplicação importadora pode
exigir selecionar esses atributos para aplicar a aparência. Para exportar
compostos a partir da agregação, abra primeiro a representação individual.

Recursos que dependem de dados ou motores ausentes não são simulados: linha do
tempo exige observações temporais reais; layouts geográficos exigem coordenadas;
OpenOrd/Yifan Hu, edição destrutiva da rede, plugins Java e projetos `.gephi`
continuam no Gephi. GEXF permite continuar a investigação nessa aplicação.
