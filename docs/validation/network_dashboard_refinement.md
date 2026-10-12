# Validação do refinamento do Molecular Network Analytics

Data: 12/10/2026. Implementação local em Flet, preservando `NetworkAnalysis`,
os resultados científicos do pipeline e a identidade visual do aplicativo.
A referência fornecida orientou hierarquia e proporções; os valores das capturas
provêm de um resultado ChEMBL real.

## Alterações

- Abas Rede, Indicadores, Análise molecular e Ranking, preservando a montagem e
  a câmera do visualizador entre as abas. Sem reconstrução da arquitetura científica.
- Área do grafo com altura responsiva (320–660 px), margem no enquadramento,
  ajustes no carregamento, Full Graph/LCC, layout, filtros que substituem vértices
  e redimensionamento. Modo expandido e recolhimento independente do inspector.
- Clique vazio e Esc removem seleção sem alterar filtros/câmera. Zoom, arraste,
  botões e o modal preservam a seleção; Esc respeita o modal aberto.
- Comunidades no seletor paginado e nos gráficos, sem lista textual redundante.
- Cartões de métricas e distribuições em grid, eixos e detalhes por hover.
  Reorganização visual não recalcula métricas.
- Pesquisa de fontes por nome/trecho/prefixo, sem diferenciar caixa; até 50 opções
  por página mais a seleção preservada quando fora da página. Scroll interno e
  tooltip com nome/origem. IDs de artefatos e fonte carregada permanecem intactos.
- Ranking com posição, ID, valor, componente/comunidade, barras proporcionais,
  ordenação, busca e limite de resultados. Seleção visível conserva a câmera;
  localização explícita ou composto fora da representação usa navegação apropriada.
- Selected Compound compacto com molécula 2D, métricas e scaffold Bemis–Murcko.
  Propriedades completas, identificadores, estruturas e relações ficam no modal
  responsivo, incluindo carregamento progressivo das relações e ação 3D existente.
- MCS sob demanda para conjunto explícito, com cache, timeout e callback de
  interrupção. Consulta SMARTS e desenho da referência permanecem distintos;
  átomos **e ligações da consulta** são destacados no composto pertencente ao conjunto.
  Trocar apenas a seleção individual não inicia outra busca MCS.

## Correção científica e limites

Bemis–Murcko é individual; moléculas acíclicas têm scaffold vazio informado.
MCS exige todas as estruturas do conjunto escolhido, cobertura 1.0, comparação
por elementos e ordem de ligações, anéis com anéis e anéis completos. Não exige
valências iguais. O desenho do fragmento vem de uma estrutura de referência,
identificada para as novas buscas; o SMARTS contém a consulta de correspondência.

Padrões interativos: 50 compostos / 10 s. Configuração permitida: 2–200 compostos,
1–30 s de busca, até 256 átomos por molécula. Nenhum composto é amostrado para
satisfazer esses limites. SMILES inválidos/ausentes impedem uma alegação de
fragmento comum. Timeout conserva um resultado parcial, sem confirmar maximalidade.
Cancelamento ocorre no próximo callback do RDKit; a interface bloqueia buscas
sobrepostas até a anterior terminar. O cache é limitado a oito buscas e 32
representações individuais. Resultados persistidos do LCC são reutilizados quando
IDs e contagem concordam; seus critérios originais são apresentados, sem substituí-los
pelos critérios da nova análise. JSON exportado inclui conjunto, estado, critérios,
SMARTS e procedência disponível. Resultados científicos originais não são regravados.

## Resultados reais

Dataset: 213 compostos, 1.611 relações, 22 componentes; LCC de 71 compostos / 907 relações.
Descritores registrados de ChEMBL aparecem na interface, incluindo massa molecular,
ALogP, PSA e contagens de doadores/aceptores. Sem imputação de dados ausentes.

| Análise | Conjunto | Resultado | Tempo observado |
|---|---:|---|---:|
| Full Graph, limite padrão 50 | 213 | Bloqueio explícito, nenhuma amostragem | 0.079 s |
| MCS do LCC, limite 100 / timeout 2 s | 71 | 14 átomos, completo, cobertura 100% | 1.036 s |
| MCS de comunidade, timeout 2 s | 10 | 18 átomos, completo, cobertura 100% | 0.012 s |

Tempos são medições locais desta execução, incluindo preparação e desenho;
não constituem garantia para outros conjuntos. O fragmento persistido exibido nas
capturas tem 17 átomos e critérios anteriores do pipeline. O novo MCS de 14 átomos
acima exige anéis completos; são análises distintas, com critérios identificados.
[IDs, SMARTS, parâmetros e medições](network_dashboard_refinement/real-validation.json).

## Verificações

Rodada final: **124 testes passaram**, cobrindo escopos Full/LCC, caches, filtros,
proveniência, limites, geometria, métricas, estruturas, exportações, navegação e
localização. Novos cenários incluem margem real do enquadramento, ausência de
reset em seleção visível, pesquisa/paginação sem perda de identidade, scaffold
acíclico, cobertura real do MCS, timeout parcial, interrupção durante callback,
invalidar gráficos ao trocar universo, impedir resultados de buscas canceladas
e preservar as relações progressivas dentro do novo modal.

Comando:

```bash
PYTHONPATH=src:tests MPLCONFIGDIR=/tmp/bme-matplotlib python -m unittest \
  test_network_refinement test_network_analysis test_molecular_network_dashboard \
  test_graph_explorer test_graph_layout_fragments test_zoom_visualizations \
  test_frontend_navigation test_localization
```

Chrome local com Flet: 1280×1000, 768×1024, 390×844 e 320×568. Sem exceções JS
ou erros de renderização nas capturas finais. Exercitados abas, ranking,
seleção, zoom, arraste, Fit, modo expandido, clique vazio, Full/LCC e modal
redimensionado no celular. O reconhecedor de toque não limpa seleção durante
arraste. O teste de fontes autorizado continua coberto com manifests reais de
fixtures; a sessão visual usa um workspace temporário e leitura do dataset real,
sem copiar credenciais nem alterar o projeto original. O modal nessa sessão
informa execução não registrada porque o dataset foi aberto fora do manifesto.

Compilação Python, `git diff --check` e geração dos guias estáticos concluídos.

## Capturas

- [Rede desktop](network_dashboard_refinement/desktop.png)
- [Seleção e scaffold](network_dashboard_refinement/selection.png)
- [Indicadores](network_dashboard_refinement/indicators.png)
- [Análise molecular](network_dashboard_refinement/molecular.png)
- [Ranking](network_dashboard_refinement/ranking.png)
- [Modal completo](network_dashboard_refinement/details.png)
- [Modal no celular](network_dashboard_refinement/mobile-details.png)
- [Modo expandido](network_dashboard_refinement/focus.png)
- [Seleção preservada após zoom](network_dashboard_refinement/zoom-selected.png)
- [Seleção preservada após arraste](network_dashboard_refinement/drag-selected.png)
- [Seleção limpa por clique vazio](network_dashboard_refinement/cleared.png)
- [LCC](network_dashboard_refinement/lcc.png)
- [Tablet](network_dashboard_refinement/tablet.png)
- [Celular](network_dashboard_refinement/phone.png)
- [Tela de 320 px](network_dashboard_refinement/small-phone.png)


## Correção de espaçamento, idioma e dados iniciais

- Parâmetros MCS em linhas separadas, largura disponível e espaçamento de 16 px.
  Desenho e resumo do resultado persistido lado a lado. Textos estáticos e
  dinâmicos traduzidos; timeout legado zero identificado como ilimitado.
- Comunidades calculadas na inicialização e troca de universo, sob os limites
  existentes, independentemente da opção de colorir comunidades. O dataset real
  mostrou **24 comunidades no Full Graph**.
- Análise molecular inicial no MCC: **71/71 valores de massa molecular** e
  Pearson grau–massa **−0,08654**, sem imputação. O universo pode ser alterado
  explicitamente; selecionar intervalo respeita o conjunto escolhido.
- Fonte inicial: CSV `Molecules` do mesmo produtor e identificador exato da
  análise, autorizado pelo manifesto. CSV de arestas não é escolhido por apenas
  ter o mesmo nome. Ausência/invalidade dessa correspondência usa os dados já
  registrados no grafo. A opção permanece disponível e a fonte pode ser trocada.
- Descritores, eixo e gráficos inicializados automaticamente. Opções sem dados
  têm estado explícito; comparar vizinho requer seleção de um composto.
- Rodada de regressão: **98 testes passaram**. Três verificações dirigidas finais
  confirmaram MCC/intervalo, tradução reversível e seleção do CSV molecular
  autorizado. Navegador local em inglês, desktop e celular, sem exceções JS.
  Validação com cópias do grafo e CSV reais em um workspace temporário.

Capturas:

- [MCS em inglês](network_dashboard_refinement/defaults/mcs-en.png)
- [Contagem de comunidades](network_dashboard_refinement/defaults/communities-en.png)
- [Fonte correlata e dados MCC iniciais](network_dashboard_refinement/defaults/molecular-mcc-en.png)
- [Celular em inglês](network_dashboard_refinement/defaults/phone-en.png)
