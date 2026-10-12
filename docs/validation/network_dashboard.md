# Filter by Graphs: diagnóstico, dashboard e validação

## Diagnóstico comprovado

A área cinza foi reproduzida no Chrome com um resultado real. Separando os
controles em quatro regiões, o problema ficou restrito à barra superior: o
seletor Full Graph/LCC estava dentro de uma `Row(wrap=True)`. O canvas, os
seletores em linhas delimitadas e os demais painéis renderizavam normalmente.
A correção colocou o seletor em sua própria linha com largura delimitada, separou
os botões e manteve entradas expansíveis dentro de linhas com restrições de
largura. O painel molecular usa a mesma regra. As capturas finais mostram Full
Graph, LCC, seleção de um composto e a apresentação em celular sem a área cinza.

Também foram identificados caminhos de custo desnecessário:

- Layout de força recalculado na abertura de resultados legados.
- Busca MCS sem limite de tempo iniciada implicitamente em arquivos sem fragmento.
- Cópias de subgrafos e repetição de contagens durante consultas à vista.
- Construção da representação nativa no caminho síncrono da interface.
- Arquivos maiores que 32 MB impedindo a descoberta das outras análises do bloco.
- Propriedades ChEMBL armazenadas em `molecule_properties`, frequentemente como
  texto de um dicionário, não acessíveis como descritores numéricos planos.

Não foi identificado vazamento do contexto GPU. A implementação mantém o canvas
Flet existente; não acrescenta um motor WebGL ou dependências de renderização.

## Arquitetura e desempenho

`network_analysis.py` concentra a topologia e os caches por universo. Layouts,
centralidades e comunidades são calculados sob demanda por `ui.call`, fora da
thread da interface. Coordenadas nativas válidas são reaproveitadas; layouts
legados recebem anéis de custo baixo. Não há cálculo automático de todos os
algoritmos ao abrir o dashboard.

Filtros usam vistas do grafo e uma chave de cache, sem copiar suas arestas.
Seleção e hover usam o índice espacial existente; selecionar um nó atualiza os
destaques e os detalhes, sem executar layouts ou reconstruir a topologia.
Layouts alternativos evitam pontos coincidentes e a área de desenho é limitada
a 4.096 px por dimensão, com aviso de ajuste. Isso impede que uma hierarquia
longa produza um canvas com dezenas de milhares de pixels, sem omitir vértices.
Contagens e métricas calculadas são reutilizadas. A construção do visualizador
inline também usa `ui.call`; resultados assíncronos obsoletos são descartados.
O inspetor carrega relações em grupos de 50, com ação explícita para carregar
mais. Buscas em componentes ou egos acima do orçamento continuam permitindo
inspecionar o composto e expandir sua vizinhança com filtros.
Ao trocar de análise ou fechar, shapes, painéis e caches derivados são liberados.

A descoberta identifica arquivos acima de 32 MB e mantém as demais análises
operacionais. Esses arquivos continuam disponíveis para download. Limites de
renderização não modificam o conteúdo científico nem produzem amostragem silenciosa.
Em uma busca que exceda o orçamento, a interface oferece ego de um salto;
se essa vizinhança também for excessiva, exibe somente o composto selecionado,
com aviso e controle explícito para ampliar a vista.

## Funcionalidades

- Full Graph/LCC, zoom, pan, busca, seleção, vizinhos e identificação por tooltip.
- Layouts circular, force-directed, Kamada–Kawai, comunidades e hierarquia BFS.
  BFS usa a maior conectividade como raiz por componente e não representa causalidade.
- Grau, Degree Centrality, Betweenness, Closeness, Eigenvector e PageRank; tamanho,
  intensidade, ranking e distribuições por métrica.
- Louvain, grupos, filtros de comunidades e rede agregada com contagens de
  compostos, relações internas e relações entre grupos.
- Componentes, isolados, densidade, grau médio, agrupamento, modularidade,
  distribuição de graus e tamanho de componentes.
- Filtros de identificador, grau, centralidade, escore e ego networks de 1–3 saltos.
- Inspeção 2D/3D existente, propriedades moleculares planas e aninhadas,
  comparação 2D com vizinho conectado.
- Histogramas moleculares, perfis por comunidade e LCC versus outros componentes,
  dispersão entre descritores ou descritor–centralidade, Pearson exploratório,
  seleção de ponto e filtragem de intervalo vinculadas à rede.
- Exportação de métricas CSV, relatório JSON com fonte, métodos e seleção,
  vista SVG; PNG e CSV do MCC existentes preservados.

Métricas topológicas ignoram escores como distâncias. Para fingerprints, a origem
registra o método; para entradas externas, `value` é apresentado como escore cuja
interpretação depende da fonte. Comunidades não são classes químicas, centralidade
não indica atividade biológica e associações não implicam causalidade.

## Integração molecular e procedência

O dashboard usa inicialmente atributos presentes no JSON do experimento. Foram
confirmados descritores em `molecule_properties`, incluindo `full_mwt`, `alogp`,
`psa`, `hba`, `hbd`, `rtb`, `aromatic_rings` e fórmula molecular. Objetos e textos
JSON/literais de dicionários são normalizados na camada de apresentação, sem
alterar os arquivos. Não são recalculados descritores de todos os compostos.

O usuário pode consultar CSVs do manifest da mesma execução e selecionar uma
única fonte adicional. `artifact_sources.py` usa metadados persistidos quando
existem; para manifests legados deriva a identidade do produtor, execução de
origem e localização lógica qualificada. Nomes iguais não identificam a fonte.

A correspondência usa o ID do composto e exige concordância dos SMILES canônicos
quando ambos existem. Duplicatas divergentes e estruturas conflitantes são
excluídas e contadas. Atributos do grafo têm prioridade; a fonte complementa
campos ausentes. O painel distingue campos do grafo e da fonte selecionada.
O relatório inclui SHA-256 da fonte. A associação verifica também o SHA-256 do
grafo aberto e rejeita resultados modificados enquanto o dashboard estava aberto.
Não existe junção automática de outras execuções ou alvos. O relatório exportado
preserva a identidade e contexto da fonte e a regra de correspondência.

Os serviços validam projeto, execução, manifest e identidade antes de ler a fonte.
Tabelas sem identificadores compatíveis ou removidas são rejeitadas. Dados ausentes
não são imputados; IDs numéricos, como CID, não são tratados como descritores.

## Validação com resultados reais

Foram encontrados nove arquivos de grafos na instalação local: oito abaixo de
32 MB e um com 82,94 MB. Os oito modelos foram carregados sem alterar os arquivos.
As medições usam `perf_counter` e `tracemalloc`, abrangendo leitura, validação,
construção e montagem dos controles; não representam o RSS total nem memória GPU.
O primeiro caso inclui importações frias.

| Compostos | Relações | LCC | Montagem (s) | Pico Python (MiB) | 100 consultas à vista (s) |
|---:|---:|---:|---:|---:|---:|
| 213 | 1.611 | 71 | 1,255 | 46,924 | 0,0014 |
| 271 | 2.167 | 95 | 0,196 | 4,414 | 0,0016 |
| 58 | 243 | 21 | 0,039 | 0,758 | 0,0007 |
| 6.446 | 24.495 | 1.134 | 1,239 | 31,748 | 0,0005 |
| 101 | 292 | 21 | 0,054 | 0,998 | 0,0009 |

O caso de 6.446 compostos mantém a topologia e informa que a vista individual
excede o orçamento; não desenha um subconjunto silencioso. Os outros três casos
de 101 compostos levaram aproximadamente 0,050–0,053 s.
As medições completas estão em [performance.json](network_dashboard/performance.json).

No navegador, o resultado real de 213 compostos/1.611 relações foi alternado para
LCC de 71 compostos/907 relações. Selecionar `CHEMBL3771360` destacou seus vizinhos
e exibiu sua estrutura 2D e atributos. Foram conferidas larguras de 1280 e 390 px.

A análise molecular real identificou 12 atributos numéricos e 213 pares de massa
molecular versus grau. Pearson foi 0,23658, apenas uma observação exploratória,
sem teste causal ou de significância. A associação a um CSV de outro produtor da
mesma execução encontrou 213 correspondências, sem conflitos de estrutura.
A referência e as contagens estão em [molecular_data.json](network_dashboard/molecular_data.json).

Capturas:

- [Desktop / Full Graph](network_dashboard/desktop.png)
- [Celular](network_dashboard/mobile.png)
- [LCC](network_dashboard/lcc.png)
- [Composto e estrutura 2D](network_dashboard/compound.png)
- [Histograma molecular e comparação de comunidades](network_dashboard/molecular.png)

## Testes automatizados e documentação

A rodada ampla passou com 112 testes de análise, integração, autorização, execução
científica do módulo, resultados, localização, zoom e layout. Após a extração dos
descritores aninhados e os ajustes finais, 77 testes direcionados passaram.
Rodadas posteriores de geometria (26 testes), navegação/preview (42),
tradução/legendas (31) e uma regressão específica de fechamento também passaram.
Os ajustes de identidade de fontes reutilizadas e hash passaram em 54 testes
direcionados, incluindo rejeição de grafos alterados e deslocamento do diretório.
As últimas verificações passaram com 14 testes de junção/normalização e 32 testes
de busca em redes extensas, layout e localização. Não se somam as rodadas, pois
existem testes em comum.
Também passaram regressões de carregamento progressivo de relações, foco
explícito em um composto de um ego excessivo e navegação até seu vizinho.
Os novos testes verificam grafos conectados, desconectados, isolados, comunidades,
limites de custo, filtros não destrutivos, métricas por universo, layouts sem
alteração de topologia, fontes persistidas/autorizadas, conflitos, dados ausentes,
seleção vinculada e limpeza ao fechar. Os testes antigos de tamanho uniforme foram
atualizados para tamanhos limitados por centralidade; a figura MCC exportada conserva
marcadores uniformes. As relações são descritas por escore, sem pressupor similaridade.

Documentação bilíngue regenerada; 32 páginas HTML passaram na verificação de arquivos,
âncoras e IDs. `compileall` e `git diff --check` passaram. Dados reais foram usados
na validação adicional, sem repetir cálculos de docking ou consultar APIs externas.

## Limites e compatibilidade

- Vista individual: até 3.000 nós/12.000 arestas; agregação: até 500 comunidades;
  dispersão: até 2.000 pares. Reduções são identificadas explicitamente.
- Kamada–Kawai: até 400 nós/8.000 arestas; force-directed: 1.500/20.000.
- Betweenness exata até 400 nós; acima disso, 64 pivôs com semente 42, até
  5.000 nós/50.000 arestas. Closeness: 1.500/20.000. Eigenvector exige conexão
  e até 50.000 arestas; PageRank: 20.000 nós/50.000 arestas. Falta de convergência
  é informada. Louvain: até 10.000 nós/50.000 arestas.
- Agrupamento médio: até 5.000 nós/20.000 arestas; métricas indisponíveis não são
  substituídas por números inventados. Pearson exige três pares e variância não nula.
- O grafo completo permanece em memória dentro do limite de arquivo. O sistema
  não oferece streaming de topologia ou WebGL. Limites podem ser revistos com
  benchmarks representativos, sem prometer consumo constante para qualquer rede.
- Fontes adicionais são restritas ao manifest da execução selecionada. Mistura
  entre execuções exigiria um fluxo explícito de comparação, não implementado.
- Procedência de atributos já incorporados a arquivos históricos depende do
  metadata original. O dashboard não inventa rastreabilidade por molécula perdida.
- Apenas propriedades registradas ou normalizadas são analisadas. Fórmulas e
  estruturas ausentes permanecem indisponíveis; não há cálculo químico em lote.
- Fragmentos antigos ausentes não são recalculados implicitamente. Fragmentos
  registrados e algoritmos científicos de construção do grafo permanecem intactos.
