# ADMET orientado à hipótese de pesquisa

[Documentação](README.md) · Português · [English](en/admet_research.md)

Revisão e implementação: 10/10/2026. O objetivo é permitir exploração reproduzível de compostos, sem transformar um conjunto de regras físico-químicas em uma afirmação de segurança ou eficácia.

## O que foi encontrado na implementação anterior

O módulo calculava seis descritores RDKit, classificava BBB/HIA por elipses locais e penalidades e estimava P-gp por uma contagem de características. PAINS e pequenos conjuntos SMARTS próprios causavam exclusão automática, descrita como toxicidade. Não havia modelo de metabolismo, excreção ou toxicidade experimental, nem escolha de hipótese ou relatório dos critérios que rejeitavam cada molécula.

A nova implementação separa descritores calculados, estimativas QSPR, regras de seleção, alertas estruturais e modelos treinados opcionais. Resultados antigos não são recalculados automaticamente. Configurações antigas, sem perfil, continuam no modo `legacy`; novos blocos começam em exploração ampla, com saída descritiva.

## Fundamentos e escolhas metodológicas

| Método | Papel na exploração | Implementação e interpretação |
|---|---|---|
| Wildman–Crippen e TPSA | Lipofilicidade e polaridade | RDKit calcula WLOGP/TPSA a partir da estrutura submetida. WLOGP não é logD em um pH específico. |
| Lipinski | Triagem de pequenas moléculas com hipótese oral | Contagem de MW > 500, WLOGP > 5, HBD > 5 e HBA > 10. Não é uma regra universal para outras vias, produtos naturais ou macrociclos. |
| Veber | Flexibilidade e polaridade em seleção oral | Violações de RB > 10 e TPSA > 140 Å². A evidência original envolve biodisponibilidade em ratos; não equivale a absorção humana medida. [Veber et al., 2002](https://pubmed.ncbi.nlm.nih.gov/12036371/). |
| Regra de três | Bibliotecas para descoberta por fragmentos | MW ≤ 300, WLOGP ≤ 3, HBD/HBA ≤ 3, com extensões RB ≤ 3 e TPSA ≤ 60. A extensão é identificada no nome da regra e pode ser substituída por intervalos próprios. [Congreve et al., 2003](https://pubmed.ncbi.nlm.nih.gov/14554012/). |
| QED | Caracterização contínua de semelhança com espaço de fármacos | Valor de 0 a 1, sem corte obrigatório. O QED RDKit pode diferir da implementação da publicação. Não mede potência, exposição ou toxicidade. [Bickerton et al., 2012](https://www.nature.com/articles/nchem.1243), [documentação RDKit](https://www.rdkit.org/docs/source/rdkit.Chem.QED.html). |
| ESOL adaptado ao RDKit | Estimativa de solubilidade aquosa | `ESOL_logS`, em log10 mol/L, usa o refit RDKit documentado por Pat Walters, com MW, MolLogP, RB e proporção aromática. Não é a equação original com CLOGP nem o ESOL SwissADME. A forma cristalina, pH e formulação não são modelados. [Implementação e referência de Delaney](https://github.com/PatWalters/solubility). |
| Potts–Guy adaptado | Hipóteses cutâneas e comparação de permeação | `logKp = −6.3 + 0.71 × WLOGP − 0.0061 × MW`, em log10 cm/s. A substituição do logP pelo RDKit é explícita. Veículo, estado da pele e ionização não são modelados. [Potts & Guy, 1992](https://pubmed.ncbi.nlm.nih.gov/1608900/). |
| PAINS/BRENK/NIH | Revisão de interferência, reatividade e motivos estruturais | Os nomes dos padrões encontrados ficam no CSV. O padrão é anotar; rejeitar exige escolha do pesquisador. Um match não comprova toxicidade. [RDKit: catálogos e referências](https://www.rdkit.org/docs/GettingStartedInPython.html#filtering-unwanted-substructures). |

Também são registrados átomos pesados, anéis, anéis aromáticos, fração de carbonos sp³, carga formal, fragmentos desconectados e avisos sobre a estrutura submetida. A interpretação de sais, misturas e espécies ionizadas deve considerar que esses descritores não representam automaticamente a espécie predominante no ensaio.

### BBB, HIA e o EGG

O BOILED-Egg publicado organiza a permeação passiva intestinal e cerebral no espaço WLOGP–TPSA. Ele não estabelece exposição livre no cérebro, transporte ativo ou biodisponibilidade oral completa. [Daina & Zoete, 2016](https://pubmed.ncbi.nlm.nih.gov/27218427/).

O EGG do BioMolExplorer conserva a geometria histórica local: HIA usa centro (75, 2), raios (75, 3); BBB usa centro (42, 2.3), raios (47, 2.2), e o classificador BBB inclui penalidades adicionais. Portanto, um ponto dentro da gema ainda pode ser BBB− pelo classificador. A legenda descreve regiões e classes separadamente. Essas fórmulas não foram validadas aqui como uma reprodução do modelo publicado.

No SwissADME original, a cor dos pontos informa a predição de substrato de P-gp; em nossa ferramenta, laranja/círculo e azul/quadrado informam BBB+/BBB−. A distinção aparece na interface. Não se atribui P-gp por descritores nas análises novas; `PGP=not_evaluated`. O perfil legado conserva a antiga heurística e identifica seu método. [SwissADME, Daina et al., 2017](https://www.nature.com/articles/srep42717).

## Bibliotecas avaliadas

| Recurso | Adequação ao BioMolExplorer | Decisão |
|---|---|---|
| RDKit | Cálculos determinísticos locais e catálogos interpretáveis, já disponível no ambiente científico | Base obrigatória desta ampliação; não é apresentado como preditor universal ADMET. |
| ADMET-AI | Inferência local com modelos Chemprop treinados em conjuntos TDC; disponível como API Python | Adaptador opcional para a linha 2.x, com extra fixado em 2.0.1. Sem chamadas a serviço remoto ou percentis DrugBank. [Repositório dos autores](https://github.com/swansonk14/admet_ai), [distribuição 2.0.1](https://pypi.org/project/admet-ai/2.0.1/). |
| Therapeutics Data Commons | Definições de endpoints, conjuntos de dados e avaliação para treinamento/comparação | Referência para significado dos endpoints; importar a biblioteca não cria, por si só, um modelo validado. [TDC](https://tdcommons.ai/single_pred_tasks/overview/). |
| DeepChem | Framework para treinamento e avaliação de modelos específicos | Adequado a uma futura hipótese com dados próprios; exige escolha de dados, splits, arquitetura e validação externa. Não substituído por um modelo genérico inventado. [Documentação oficial](https://deepchem.readthedocs.io/en/stable/index.html). |
| OpenADMET | Ferramentas para construir e usar modelos de propriedades com conjuntos experimentais | Alternativa para desenvolvimento e validação de modelos por projeto; exige integração e validação por endpoint. [Projeto](https://www.openadmet.org/), [exemplo de inferência](https://demos.openadmet.org/en/latest/demos/05_Ensemble_Model_Inference/05_Model_Ensemble_Inference.html). |
| SwissADME | Referência científica e ferramenta web | Usado para fundamentar interpretação, não como dependência de scraping ou promessa de equivalência numérica. |

ADMET-AI v2 e v1 foram treinados separadamente e podem produzir valores diferentes. O relatório registra a versão usada. Os endpoints, nomes, tipos de tarefa, espécies e unidades são lidos dos metadados que acompanham a instalação. Não se aplica corte 0.5 universal nem se substitui a classe BBB local por uma saída de rede neural.

### Cobertura que depende de modelos treinados

A API opcional acrescenta as saídas disponíveis naquela versão, prefixadas por `ML_`. Elas abrangem absorção/permeabilidade e solubilidade; distribuição, ligação a proteínas e volume de distribuição; metabolismo por CYP; meia-vida e clearance; e endpoints como hERG, AMES e DILI. A presença de uma saída não garante sua confiabilidade para a série química estudada. Os significados e escalas devem ser conferidos por endpoint. [TDC ADME](https://tdcommons.ai/single_pred_tasks/adme/), [TDC toxicidade](https://tdcommons.ai/single_pred_tasks/tox/).

Um exemplo crítico é `ML_Pgp_Broccatelli`: trata-se de **inibição de P-gp**, não de classificação de substrato. Não é usado para colorir o EGG ou inferir efluxo. `ML_BBB_Martins` também permanece distinto de `BBB`, calculado pela heurística local. Os limites personalizados podem selecionar qualquer saída `ML_` efetivamente retornada. Valor ausente não atende a um limite exigido.

O preditor pode ser instalado em um ambiente isolado da aplicação, sem substituir as dependências científicas existentes:

```bash
python scripts/setup_admet_ai.py
```

O instalador usa o Python que o invoca para preparar `.admet-env` com PyTorch CPU e ADMET-AI 2.0.1, verifica uma inferência e só então marca o ambiente como pronto. O worker descobre esse ambiente automaticamente. Para outra localização, configure `BIOMOL_ADMET_PYTHON` com o caminho do Python correspondente. Também é possível instalar o extra `admet-ml` diretamente em um worker dedicado. A configuração mostra a disponibilidade antes de executar. Use um Python compatível; o ambiente verificado utiliza Python 3.12. Consulte [instalação e diagnóstico do preditor](installation.md#preditor-admet-ai-opcional).

O adaptador processa SMILES únicos em lotes limitados, verifica identidade/ordem e valores finitos e reutiliza a predição para estruturas repetidas. Nesta implementação foi executada inferência real para etanol e ácido acético, retornando 41 endpoints via ambiente isolado. Essa verificação funcional não constitui validação científica dos modelos para uma biblioteca de pesquisa; domínio de aplicabilidade e incerteza continuam não avaliados.

## Como configurar uma hipótese

| Perfil | Ponto de partida | Cuidados na interpretação |
|---|---|---|
| Exploração ampla | Todos os compostos válidos, anotação de alertas | Não perder diversidade antes de conhecer a série. |
| Candidatos orais | HIA+, Lipinski com até uma violação, Veber sem violações | Absorção passiva não inclui metabolismo de primeira passagem. |
| Acesso ao SNC | BBB+ local | Investigar efluxo, pKa/logD, ligação e exposição cerebral livre. |
| Restrição periférica | BBB− local | Não demonstra ausência de exposição cerebral. |
| Fragmentos | Regra de três ampliada, sem violações | Ajustar para estratégia de ligação e crescimento de fragmentos. |
| Produtos naturais / além da regra de cinco | Sem corte oral obrigatório | Macrociclos, transporte ativo e outras vias podem exigir modelos próprios. |
| Personalizada | Combinação livre de critérios | Documentar racional e unidades de cada limite. |
| Legado | Exclusão histórica por padrões | Apenas compatibilidade; os padrões não constituem toxicidade validada. |

O perfil oferece valores iniciais. BBB e HIA aceitam não restringir, positivo ou negativo. Regras, tolerância de violações, catálogos e tratamento dos alertas são editáveis. Intervalos são inclusivos, e todos os critérios são combinados com **E**. Um mínimo ou máximo vazio não restringe aquele lado. Não há aprovação por média que oculte uma falha relevante.

`Anotar todos os compostos` envia todos os avaliados ao próximo bloco e apenas marca a hipótese. `Filtrar pela hipótese` filtra o CSV principal. Ambos preservam uma tabela `_evaluated.csv`, `ADMET_selected` e `ADMET_reasons`; análises novas também exportam `_selected.csv`, `_rejected.csv`, BBB± e HIA±. O relatório `_admet_report.json` registra configurações, versões, métodos, quantidades e limitações. Registros moleculares inválidos continuam tratados pelo relatório de qualidade molecular do pipeline.

Exemplo: uma hipótese periférica oral pode combinar BBB−, HIA+, regras orais e `ESOL_logS ≥ −4`. Esse último corte é uma escolha do pesquisador, não uma recomendação clínica. A exploração de produtos naturais pode manter todos os compostos, registrar QED e solubilidade e revisar alertas sem exclusão automática.

## Visualização e interpretação

O EGG interativo mede seu próprio painel e redesenha eixos, regiões e pontos na largura disponível. O painel usa altura de 560 px; detalhes e downloads passam para baixo em larguras menores. O ajuste não altera os descritores, não desloca pontos por jitter e inclui outliers com margem nos limites. Hover identifica código, TPSA, WLOGP e classificações; clique permite consultar estrutura e propriedades. Pontos coincidentes continuam disponíveis pelo seletor já existente.

Laranja/azul com círculo/quadrado reduz dependência da percepção de cor. Grade discreta, eixos com unidades, legenda explícita e contraste moderado favorecem leitura. O PNG usa a mesma geometria e codificação. Esta aplicação dos princípios visuais é uma decisão de interface, apoiada na discussão de acessibilidade e transparência das figuras científicas, sem reivindicar avaliação de usabilidade experimental. [Crameri et al., 2020](https://www.nature.com/articles/s41467-020-19160-7), [Nature Methods, Points of View, 2026](https://www.nature.com/articles/s41592-026-03143-5).

[Veja um EGG exportado](../imgs/examples/admet-egg.png) e suas [estruturas de entrada e propriedades calculadas](../imgs/examples/admet-example.json). O exemplo ilustra a triagem local; não é uma predição SwissADME validada.

## Limites atuais e evolução científica necessária

Nenhuma regra, alerta ou saída neural constitui medição experimental. Não há ainda avaliação validada de domínio de aplicabilidade por composto, intervalos de predição, pKa/logD, permeação de espécies dependentes de pH, agregação coloidal experimental, exposição livre cerebral ou integração de dados específicos do ensaio. O relatório marca domínio/incerteza como não avaliados, em vez de gerar uma confiança fictícia.

Para uma hipótese específica, o próximo avanço científico é reunir dados relevantes à espécie, tecido, via e ensaio; definir splits que evitem vazamento de séries; comparar modelos e baselines; avaliar calibração, extrapolação e desempenho externo; e confirmar os compostos selecionados em ensaios apropriados. Os perfis implementados facilitam essa exploração, mas não substituem essa validação.

## Organização da configuração e gestão dos resultados

Os campos são agrupados por objetivo/saída e critérios de seleção, em duas colunas quando há largura suficiente e uma coluna em telas estreitas. Alertas e intervalos ficam em seções recolhíveis; as explicações extensas ficam na ajuda contextual. A organização aplica proximidade entre campos relacionados, espaçamento consistente e divulgação progressiva, conforme [agrupamento e espaço em formulários](https://www.nngroup.com/articles/form-design-white-space/) e [divulgação progressiva](https://www.nngroup.com/articles/progressive-disclosure/).

A execução respeita seletores de arquivos explícitos e nomes de arquivos salvos nos blocos; solicita curadoria quando alguma entrada conectada ainda usa seleção automática. Cada etapa dispõe de “Remover todos os resultados”, com confirmação. A remoção inclui os arquivos auxiliares e invalida referências reaproveitadas; preserva registro e log e permite restauração pelo histórico. Execuções ativas bloqueiam a remoção.
