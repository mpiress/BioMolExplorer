# Documentação científica e metodológica

Esta referência descreve o comportamento científico do BioMolExplorer 2.2.0 no estado atual do código. Fórmulas e parâmetros abaixo são contratos implementados, não recomendações clínicas nem validação prospectiva. O arquivo de origem de cada método é indicado para permitir auditoria.

## Escopo, proveniência e interpretação

O pipeline separa descritores, heurísticas, predições treinadas, escores de docking e análises topológicas. Eles não são probabilidades intercambiáveis. A interface conserva arquivos de origem, parâmetros e artefatos por etapa; resultados históricos podem ter sido produzidos por versões anteriores e devem ser avaliados junto ao manifesto do projeto.

Não há, no código atual, cálculo de PCA, t-SNE, UMAP, AUC-ROC ou Enrichment Factor. Também não há classificador de ponte salina no diagrama de interações. Esses métodos não são atribuídos à plataforma.

## ADMET, propriedades e ADMET-AI

### Descritores RDKit

`src/caad/admet.py::MoleculeEvaluator.calculate_properties` calcula a partir do SMILES sanitizado: massa molecular `MW`, área de superfície polar topológica `TPSA`, Wildman-Crippen `WLOGP`, doadores e aceptores de hidrogênio (`HBD`, `HBA`), ligações rotacionáveis `RB`, QED, fração de carbonos sp3, átomos pesados, anéis totais/aromáticos, carga formal, fragmentos e proporção aromática.

A proporção aromática implementada é

```text
AromaticProportion = N_atomos_aromaticos / max(1, N_atomos_pesados)
```

QED é o método do RDKit baseado na desirability ponderada de Bickerton et al.; o BioMolExplorer não reimplementa nem altera seus pesos.

### Heurísticas e regras

`MoleculeEvaluator.screening_properties` implementa:

```text
ESOL_logS = 0.26121066137801696
            - 0.0066138847738667125 MW
            - 0.7416739523408995 WLOGP
            + 0.003451545565957996 RB
            - 0.42624840441316975 AromaticProportion

logKp = -6.3 + 0.71 WLOGP - 0.0061 MW
```

`ESOL_logS` é um ajuste local compatível com a implementação RDKit documentada e tem unidade `log10(mol/L)`. `logKp` é uma heurística Potts-Guy com logP RDKit, em `log10(cm/s)`; nenhuma delas é uma medição experimental.

As violações são contagens booleanas exatas:

```text
Lipinski = [MW > 500] + [WLOGP > 5] + [HBD > 5] + [HBA > 10]
Veber    = [RB > 10] + [TPSA > 140]
Fragment = [MW > 300] + [WLOGP > 3] + [HBD > 3]
           + [HBA > 3] + [RB > 3] + [TPSA > 60]
```

O perfil oral aceita por padrão no máximo uma violação de Lipinski e zero de Veber. O perfil de fragmentos aceita zero violações. Os perfis `exploratory`, `oral`, `cns`, `peripheral`, `fragment`, `natural_products`, `custom` e `legacy`, bem como `annotate`/`filter`, estão em `src/biomolexplorer/admet_config.py`. PAINS, BRENK, NIH, ZINC e ChEMBL são catálogos estruturais do RDKit; correspondência é alerta de revisão, não prova de toxicidade.

### BOILED-Egg, BBB, HIA e P-gp

As classificações locais em `src/caad/admet.py` são heurísticas explícitas:

```text
egg_BBB = ((TPSA - 42)^2 / 47^2) + ((WLOGP - 2.3)^2 / 2.2^2)
BBB+ se egg_BBB <= 1 e penalidades <= 1

penalidades = [TPSA > 90] + [MW > 500] + [RB > 10]
              + [HBD > 2] + [WLOGP < 0] + [WLOGP > 6]

egg_HIA = ((TPSA - 75)^2 / 75^2) + ((WLOGP - 2)^2 / 3^2)
HIA+ se egg_HIA <= 1

Pgp_score = [MW > 400] + [TPSA > 75] + [HBD >= 2] + [WLOGP > 3.5]
PGP+ se Pgp_score >= 2
```

Esses rótulos não modelam ionização, exposição livre ou domínio de aplicabilidade e não devem ser lidos como probabilidades.

### ADMET-AI 2.x

`src/biomolexplorer/admet_predictors.py` integra `admet-ai==2.0.1`, localmente ou em ambiente isolado definido por `BIOMOL_ADMET_PYTHON`. O modelo é instanciado como `ADMETModel(include_physchem=False, drugbank_path=None, num_workers=0)` e processado em lotes de até 256 SMILES únicos. As colunas são preservadas com prefixo `ML_`, nomes/unidades/metadados do provedor e sem percentis DrugBank. O BioMolExplorer não aplica limiar probabilístico genérico, não inventa incerteza e registra `domain_assessment=not_assessed` e `uncertainty=not_provided`.

## Assinaturas moleculares e similaridade

### Morgan/ECFP

`src/kernel/descriptors.py::Descriptors.morgan_worker` chama `GetMorganFingerprintAsBitVect` com padrão `radius=2`, `nBits=2048`, `useChirality=True`, `useFeatures=True` e `useBondTypes=False`. Assim, o raio equivale ao ECFP4 por diâmetro, mas os invariantes são farmacofóricos/funcionais do RDKit, incluem quiralidade e deliberadamente não diferenciam tipos de ligação durante a atualização. `morgan_n_bits` e `radius` podem ser alterados na operação `fingerprints`.

### MACCS e farmacóforo 2D

MACCS usa `rdkit.Chem.MACCSkeys.GenMACCSKeys`. O RDKit devolve vetor de 167 posições, com a posição zero não utilizada e 166 chaves públicas/SMARTS; não se trata de uma reprodução das chaves privadas originais.

O farmacóforo usa `Gobbi_Pharm2D.factory` e `Generate.Gen2DFingerprint`. A fábrica RDKit codifica pares e trios de famílias de características químicas em faixas de distância topológica 2D. Portanto, a implementação não calcula uma matriz euclidiana 3D própria nem expõe um conjunto manual de bins no BioMolExplorer; famílias e bins pertencem à fábrica Gobbi do RDKit usada pela versão instalada.

### Métricas binárias e busca aproximada

Para bits ativos `A` e `B`, seja `a=|A|`, `b=|B|`, `c=|A interseção B|` e `n` o tamanho do vetor. As métricas mais usuais expostas são:

```text
Tanimoto = c / (a + b - c)
Dice     = 2c / (a + b)
Cosine   = c / sqrt(ab)
Sokal    = c / (2a + 2b - 3c)
Russel   = c / n
Kulczynski = c(a + b) / (2ab)
```

Também são encaminhadas diretamente às implementações RDKit `RogotGoldberg`, `AllBit`, `McConnaughey`, `Asymmetric` e `BraunBlanquet`. `MolSimilarity.get_similarity` arredonda a três casas antes de aplicar o threshold. Com `approximate=True`, candidatos são pré-selecionados por MinHash LSH (`num_perm=256`); a métrica RDKit exata é calculada somente nos candidatos. Com `approximate=False`, todos os vetores são comparados. O threshold da operação é convertido de porcentagem para fração pelo wrapper.

### MCS e redes

O MCS usa `rdFMCS.FindMCS`; opções de anel e timeout vêm da configuração da operação de grafos. Redes de similaridade são não direcionadas e, para centralidade/comunidades, não ponderadas: os valores de similaridade não são tratados como distâncias. `src/biomolexplorer/network_analysis.py` oferece grau, centralidade de grau, betweenness, closeness com correção Wasserman-Faust, eigenvector em universo conectado, PageRank com `alpha=0.85`, k-core, clustering, triângulos e Louvain sem pesos com semente 42. Há orçamentos computacionais explícitos em `network_centrality.py`; aproximação ocorre apenas no betweenness grande e é informada.

## Docking e consenso Vina + DOCK 6

### Escores brutos

Vina e DOCK 6 mantêm seus escores brutos e poses. O consenso não compara diretamente energias entre motores/receptores. Para cada motor, alvo e composto, `src/biomolexplorer/consensus_model.py` conserva o menor score; empates usam deterministicamente o caminho do conformero.

### Quantil com empates e ECR

Em um grupo ordenado de tamanho `N`, uma faixa empatada que ocupa índices zero-based de `s` a `e-1` recebe:

```text
q = ((s + e - 1) / 2) / max(N - 1, 1)
r = exp(-q / tau), com tau = 0.1
```

Logo, melhor rank tem `q=0` e relevância 1; pior rank tende a `q=1`. Para `T` alvos, ranks ausentes recebem `q=1` na penalização, e a cobertura observada multiplica a média geométrica:

```text
coverage_engine = observacoes_engine / T
factor_engine = coverage_engine * exp(-(sum_t q_engine,t) / (T tau))
consensus = (factor_vina + factor_dock6) / 2
```

Somente compostos com pelo menos um par alvo-composto observado nos dois motores entram no resultado. A cobertura parcial permanece explícita. `repulsion_weight` é aceito por compatibilidade, mas é atualmente um argumento sem efeito. O modelo é registrado como `target-balanced-ecr-v1`. O resultado é ordenado por consenso decrescente e código molecular crescente.

Esse rank-consensus reduz dependência de escalas incompatíveis, mas não é probabilidade calibrada de ligação e não garante redução de falsos positivos; essa hipótese precisa ser validada em conjuntos independentes.

## Footprints energéticos do DOCK 6

`src/biomolexplorer/resources/dock6/footprint.template` executa `footprint_similarity_score_primary` usando uma molécula de referência MOL2. O DOCK 6 produz por resíduo componentes de van der Waals e eletrostática e o escore de similaridade do footprint. O BioMolExplorer não recalcula o potencial: `src/biomolexplorer/footprints.py` valida e lê a saída finita do motor, enquanto `src/biomolexplorer/ui/footprint.py` apresenta os perfis.

Conceitualmente, um footprint é o vetor de contribuições por resíduo, por exemplo `F=(E_vdW,1, E_es,1, ..., E_vdW,n, E_es,n)`, comparado ao vetor de referência conforme a função FPS do DOCK 6. As ligações de hidrogênio podem contribuir ao scoring do motor/protocolo, mas o parser atual não publica uma terceira série independente de energia de H-bond por resíduo. Portanto, a documentação e a interface não devem inventar essa coluna.

Antes da avaliação, `footprint_receptor` ordena átomos pelo identificador de subestrutura MOL2 e renumera átomos/ligações para manter resíduos contíguos. Isso é preparação estrutural, não normalização energética.

## AMBER Score e saneamento estrutural

O fluxo `dynamics_dock6` usa o Amber Score do DOCK 6 e AmberTools via executor central (`src/biomolexplorer/preparation_tools.py`). Os parâmetros públicos incluem região móvel, limiar RMS gradient, ciclos de minimização antes/depois da dinâmica, passos MD, modelo GB, cutoff não ligado, temperatura e seleções móveis. Perfis de resíduos e PTMs são descritos em `amber_score_structural_support.md`; componentes desconhecidos não são silenciosamente convertidos em aminoácidos padrão. A política de resíduos, arquivos preparados e diagnósticos devem acompanhar o resultado.

Amber Score e dinâmica curta são refinamentos do protocolo configurado, não estimativas rigorosas de energia livre. Metais, crosslinks, estados de protonação/tautomeria, loops ausentes e parâmetros não cobertos exigem revisão especializada.

## Validação física das poses

`src/biomolexplorer/pose_geometry.py` implementa a bateria de critérios identificada como `posebusters_inspired_v1`, inspirada em Buttenschoen, Morris e Deane (2024). A implementação é própria e não executa o pacote PoseBusters.

| Critério | Medida implementada |
| --- | --- |
| Química | sanitização RDKit da molécula lida, sem sanitização automática na leitura |
| Conectividade | `len(GetMolFrags(molécula)) == 1` |
| Coordenadas | posições finitas e distintas |
| Geometria interna | maior desvio relativo de comprimento de ligação e maior desvio absoluto de ângulo contra o confôrmero relaxado; padrões `0,25` e `25°` |
| Planaridade | maior distância ao plano ajustado de anéis aromáticos e ligações duplas fora de anel; padrão `0,25 Å` |
| Energia interna | `ΔE = E(pose) − E(relaxado)` no mesmo campo de força; padrão `100 kcal/mol` |
| Choque estérico | menor razão `d / (r_vdW(ligante) + r_vdW(receptor))` entre átomos pesados; padrão `0,75` |

A referência interna é uma cópia da mesma conectividade minimizada com MMFF94 e, sem parametrização completa, com UFF; o método efetivo é registrado por pose. Raios de van der Waals seguem Bondi (1964), com `1,70 Å` para elementos não listados. Solvente e o ligante cristalográfico do próprio sítio são excluídos do teste estérico.

Poses lidas por proximidade, como PDBQT, não recebem os critérios que dependem de ordens de ligação: sem elas, o campo de força e o teste de planaridade mediriam outra molécula. Sem receptor associado, o critério estérico fica sem medida. Em ambos os casos o relatório registra a condição, em vez de reprovar a pose.

ΔE é tensão interna conformacional no campo de força usado. Não é energia de ligação, energia livre nem probabilidade de acerto do modo de ligação. Aprovação indica plausibilidade física da geometria entregue. Consulte [validação de poses](pose_validation.md).

## Interações ligante-receptor e visualização 3D

`src/biomolexplorer/docking_interactions.py` é um classificador geométrico conservador inspirado nos thresholds do PLIP; não é o motor PLIP. Ele exige topologia química coerente, coordenadas finitas e hidrogênios explícitos para H-bonds. Se essas condições falham, a capacidade é marcada indisponível em vez de inferir química.

| Interação | Critério implementado |
| --- | --- |
| H-bond | distância doador-aceptor `0.5 < d <= 4.1 Å`; H ligado a `0.6 < d(D,H) <= 1.25 Å`; ângulo D-H-A `> 100 graus`; `d(H,A) < d(D,A)` |
| Hidrofóbica | par C-C elegível mais próximo por resíduo, `0.5 < d <= 4.0 Å`; carbono proteico não pode estar a menos de `1.9 Å` de N/O/S |
| pi paralelo | centros aromáticos `0.5 < d <= 5.5 Å`, offset `<= 2 Å`, ângulo entre normais `<= 30 graus` |
| pi-T | mesmos limites de centro/offset e ângulo entre normais `>= 60 graus` |
| van der Waals | gap `d - r_vdW(proteina) - r_vdW(ligante)` entre `-0.4` e `0.5 Å` |

A superfície do sítio de ligação do visualizador usa a escala de hidrofobicidade de Kyte e Doolittle (1982) por resíduo, interpolada entre −4,5 e +4,5, sobre a superfície de exclusão de solvente da estrutura recebida; estados de protonação equivalentes (HID/HIE/HIP, CYX, ASH, GLH, LYN, MSE) são mapeados para o resíduo padrão, e componentes fora da escala ficam neutros. É uma representação, não um cálculo de energia, área enterrada ou druggability.

Anéis proteicos reconhecidos: Phe, Tyr, Trp e estados His/HID/HIE/HIP. A vizinhança de resíduos é pré-limitada a 7.5 Å. O detector não classifica ponte salina, cátion-pi, halogênio ou coordenação metálica. O contato van der Waals é geométrico, não energia. A interface renderiza linhas, cores, resíduos, distâncias e ângulos junto à pose/receptor efetivamente selecionados.

## Referências fundamentais

- Rogers D, Hahn M. Extended-Connectivity Fingerprints. Journal of Chemical Information and Modeling. 2010;50:742-754. DOI 10.1021/ci100050t.
- Durant JL et al. Reoptimization of MDL Keys for Use in Drug Discovery. Journal of Chemical Information and Computer Sciences. 2002;42:1273-1280. DOI 10.1021/ci010132r.
- Gobbi A, Poppinger D. Genetic optimization of combinatorial libraries. Biotechnology and Bioengineering. 1998;61:47-54.
- Bajusz D, Racz A, Heberger K. Why is Tanimoto index an appropriate choice for fingerprint-based similarity calculations? Journal of Cheminformatics. 2015;7:20.
- Lipinski CA et al. Experimental and computational approaches to estimate solubility and permeability. Advanced Drug Delivery Reviews. 2001;46:3-26.
- Veber DF et al. Molecular properties that influence oral bioavailability. Journal of Medicinal Chemistry. 2002;45:2615-2623.
- Delaney JS. ESOL: A Simple Model for Aqueous Solubility. Journal of Chemical Information and Computer Sciences. 2004;44:1000-1005.
- Daina A, Zoete V. A BOILED-Egg to predict gastrointestinal absorption and brain penetration. ChemMedChem. 2016;11:1117-1121.
- Bickerton GR et al. Quantifying the chemical beauty of drugs. Nature Chemistry. 2012;4:90-98.
- Trott O, Olson AJ. AutoDock Vina. Journal of Computational Chemistry. 2010;31:455-461.
- Allen WJ et al. DOCK 6: Impact of new features and current docking performance. Journal of Computational Chemistry. 2015;36:1132-1156.
- Balius TE et al. Testing inhomogeneous solvation theory in structure-based ligand discovery. Proceedings of the National Academy of Sciences. 2017;114:E6839-E6846.
- Salentin S et al. PLIP: fully automated protein-ligand interaction profiler. Nucleic Acids Research. 2015;43:W443-W447.
- Buttenschoen M, Morris GM, Deane CM. PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences. Chemical Science. 2024;15:3130-3139.
- Bondi A. van der Waals Volumes and Radii. The Journal of Physical Chemistry. 1964;68:441-451.
- Halgren TA. Merck molecular force field. I. Basis, form, scope, parameterization, and performance of MMFF94. Journal of Computational Chemistry. 1996;17:490-519.
- Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. Journal of Molecular Biology. 1982;157:105-132.
- Freeman LC. A set of measures of centrality based on betweenness. Sociometry. 1977;40:35-41.
- Blondel VD et al. Fast unfolding of communities in large networks. Journal of Statistical Mechanics. 2008;P10008.

## Rastreabilidade de implementação

Os principais pontos auditáveis são `src/caad/admet.py`, `src/biomolexplorer/admet_config.py`, `src/biomolexplorer/admet_predictors.py`, `src/kernel/descriptors.py`, `src/biomolexplorer/consensus_model.py`, `src/biomolexplorer/footprints.py`, `src/biomolexplorer/docking_interactions.py`, `src/biomolexplorer/network_analysis.py`, `src/biomolexplorer/network_centrality.py`, `src/biomolexplorer/pose_geometry.py`, `src/biomolexplorer/molecular_surface.py`, `src/biomolexplorer/operations.py` e os templates em `src/biomolexplorer/resources/dock6/`.
