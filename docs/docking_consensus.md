# Consenso de docking: modelo matemático e interpretação

[Portal de documentos](../index_doc.html) · [English](en/docking_consensus.md) · [Parâmetros da API](backend_usage.md#api-consensus)

O bloco **Docking consensus** transforma avaliações do AutoDock Vina e DOCK6 em uma lista única de prioridades **por composto**. O modelo implementado é `target-balanced-ecr-v1`, definido por `consensus_model.py`, com `TAU = 0.1`. Este guia explica exatamente o cálculo atual, com exemplos reproduzíveis, sem alterar o algoritmo.

**Leitura rápida:** o score bruto indica energia segundo um motor; o rank indica a posição relativa do composto; a relevância transforma esse rank para uma escala adimensional; o consenso combina as relevâncias dos dois motores. **Menor score bruto é melhor; maior relevância/consenso é melhor.** O valor final não é energia, probabilidade de ligação, p-valor nem confiança estatística calibrada.

## 1. O problema que o modelo resolve

Vina e DOCK6 usam funções de pontuação diferentes. Seus valores não devem ser somados ou promediados diretamente como se tivessem a mesma calibração. Mesmo dentro de um motor, alvos distintos podem produzir distribuições diferentes. O BioMolExplorer compara posições relativas **separadamente em cada combinação motor/alvo**, depois agrega essas informações.

Se 1.000 compostos forem avaliados nos mesmos três alvos por ambos os motores, há 3.000 avaliações Vina e 3.000 DOCK6, mas o consenso retorna **1.000 linhas**. Isso exige que todos os compostos tenham pelo menos um par composto/alvo comum aos dois motores. Sem qualquer par comum, a operação informa que não há interseção; não inventa resultados.

O bloco não realiza docking, preparação de estruturas, minimização ou geração de footprint. Ele usa os scores e as referências estruturais já produzidos.

## 2. Notação e dados necessários

| Símbolo | Significado |
| --- | --- |
| i | Identificador do composto, em `molecule_chembl_id` |
| e | Motor: Vina ou DOCK6 |
| t | Identificador do alvo, em `receptor_id` |
| x_iet | Menor score bruto observado para o composto i no motor e e alvo t |
| N_et | Número de compostos distintos observados naquele motor/alvo |
| r_iet | Rank de x_iet, começando em 1 e usando média em empates |
| q_iet | Rank normalizado entre 0 e 1 |
| p_iet | Relevância exponencial do composto naquele motor/alvo |
| T | Quantidade de alvos na união das entradas dos dois motores |
| c_ie | Fração dos T alvos com observação para o composto i no motor e |
| R_ie | Relevância global do composto no motor e, após agregar os alvos |
| C_i | Consenso final entre R_i,Vina e R_i,DOCK6 |

As linhas precisam identificar corretamente composto, alvo e score finito. Os identificadores são correspondências textuais: o modelo não descobre que `AChE`, `ACHE` e `4M0E_A` representam o mesmo alvo. Uma identificação ausente/ambígua compromete a comparação; um rótulo vazio não permite interpretar biologicamente aquela categoria.

Os SMILES isoméricos não vazios são canonizados para conferir a identidade: o mesmo código com estruturas ou estereoisômeros diferentes causa erro. Se o SMILES não foi fornecido, o modelo não reconstrói a estereoquímica nem comprova a identidade estrutural apenas pelo código. Códigos diferentes não são fundidos automaticamente.

### Exemplo mínimo do formato de entrada

Para reproduzir o exemplo 1, as tabelas podem estar em `Vina/docking_results.csv` e `Dock6/docking_results.csv` sob a pasta de entrada. Os nomes de motor são `vina` e `dock6`.

Tabela Vina:

```csv
molecule_chembl_id,canonical_smiles,receptor_id,engine,score
A,CCO,T,vina,-9
B,CCC,T,vina,-8
C,CCN,T,vina,-7
```

Tabela DOCK6:

```csv
molecule_chembl_id,canonical_smiles,receptor_id,engine,score
A,CCO,T,dock6,-30
B,CCC,T,dock6,-40
C,CCN,T,dock6,-20
```

Esses arquivos demonstram o cálculo dos fatores. Não incluem anexos estruturais: para abrir 3D/footprint, as entradas reais precisam conservar os respectivos arquivos e referências de origem. **Ausência de observação significa ausência da linha do par**, não preencher `score` com texto vazio, NaN ou infinito; linhas com score não numérico/não finito são inválidas.


## 3. Seleção de observações e população de referência

Antes de normalizar, o algoritmo conserva o menor score de cada chave `(motor, alvo, composto)`. Repetir uma pose ou importar o mesmo resultado em mais de um arquivo não aumenta seu peso. Scores iguais para a mesma chave usam o caminho estrutural como desempate determinístico, sem mudar o valor numérico.

A normalização usa todos os compostos observados em cada motor/alvo, **inclusive aqueles que depois não forem elegíveis para a saída final**. A elegibilidade é aplicada à lista final: é necessário pelo menos um alvo em que o composto tenha sido avaliado pelos dois motores. Assim, populações diferentes entre motores podem gerar ranks diferentes mesmo para um composto compartilhado.

T é a união dos alvos encontrados nas entradas, não apenas a interseção dos alvos completamente avaliados. Um alvo presente somente em um motor participa da política de cobertura. Por isso, confira a correspondência e a completude dos conjuntos antes de interpretar o ranking.

## 4. Normalização ordinal por motor e alvo

Os scores x são ordenados em ordem crescente. O melhor recebe rank 1. Quando os valores empatam, usa-se a média das posições ocupadas pelo empate:

```text
r_iet = média das posições de x_iet na ordenação crescente
q_iet = (r_iet - 1) / max(N_et - 1, 1)
```

Com N>1, q=0 indica a primeira posição e q=1 a última. Se quatro compostos têm scores `[-10, -8, -8, -6]`, os ranks são `[1, 2.5, 2.5, 4]` e q é `[0, 0.5, 0.5, 1]`. O identificador não resolve artificialmente um empate energético nesta etapa.

Casos especiais:

- **Um composto:** N=1, r=1 e q=0. Isso indica somente que não há concorrentes observados.
- **Todos empatados com N>1:** r=(N+1)/2 e q=0.5 para todos. A distribuição não permite discriminação entre compostos.
- **Score positivo:** permanece positivo na ordenação. Não é substituído por zero nem convertido para valor absoluto.
- **Score ausente:** não recebe energia inventada; seu tratamento ocorre na agregação descrita adiante.
- **NaN/infinito:** são recusados como scores não finitos.

A transformação usa a ordem, não a distância energética. Uma transformação estritamente crescente dos scores preserva os ranks. Isso reduz dependência de unidade/escala e da magnitude dos extremos, mas descarta a informação de quanto duas energias diferem. Acrescentar ou remover concorrentes pode mudar N e os ranks normalizados; relevâncias de execuções com bibliotecas diferentes não são diretamente comparáveis.

## 5. Relevância exponencial e relação com ECR

O ECR publicado por Palacio-Rodríguez et al. usa contribuições exponenciais de ranks e sua soma entre funções de pontuação. Na expressão original, cada contribuição tem a forma `(1/σ) × exp(-r/σ)`. O BioMolExplorer **adapta essa ideia**, usando rank reescalado e uma relevância relativa adimensional; não afirma reproduzir literalmente a densidade original. A agregação multialvo abaixo é uma política adicional da plataforma. [Fonte primária do ECR](https://doi.org/10.1038/s41598-019-41594-3).

```text
τ = 0.1
p_iet = exp(-q_iet / τ)
```

Para N>1, a mesma expressão pode ser escrita como `p = exp(-(r-1)/(τ × (N-1)))`: o melhor rank recebe 1 e não é usado o prefator de densidade do ECR original.

τ é fixo na implementação atual; não há controle desse parâmetro no formulário. Não é um corte que remove automaticamente 90% das moléculas. Todos os compostos elegíveis continuam na saída; o decaimento concentra a relevância relativa nos melhores ranks.

| q | Relevância p = exp(-q/0.1) |
| --- | --- |
| 0.0 | 1.000000 |
| 0.1 | 0.367879 |
| 0.2 | 0.135335 |
| 0.5 | 0.006738 |
| 1.0 | 0.000045 |

As relevâncias não somam necessariamente 1. Logo, p não é uma probabilidade. Um singleton recebe p=1 sem demonstrar afinidade; uma distribuição constante com mais de um composto recebe `exp(-5) ≈ 0.006738` para todos, não zero. O [estudo de combinação de poses e ranks](https://pubs.rsc.org/en/content/articlehtml/2021/ra/d1ra05785e) exemplifica o uso de uma escala exponencial de 10% da biblioteca; isso não demonstra que a adaptação atual seja universalmente ótima.

## 6. Agregação de alvos e tratamento de dados ausentes

Para cada composto/motor, contamos os alvos em que há observação. Um alvo ausente recebe **q*=1**, correspondente à pior posição relativa, sem preencher um score bruto:

```text
q*_iet = q_iet, se houver observação
q*_iet = 1,     se não houver observação

c_ie = número de alvos observados / T
R_ie = c_ie × exp( -sum_t(q*_iet) / (T × τ) )
```

A forma equivalente torna explícita a média geométrica:

```text
p*_iet = exp(-q*_iet / τ)
R_ie = c_ie × (produto_t p*_iet)^(1/T)
```

O cálculo usa a soma de ranks no expoente, evitando multiplicar diretamente muitos fatores pequenos. Todos os alvos têm o mesmo peso. Com cobertura completa, c=1 e R é a média geométrica das relevâncias por alvo. Com ausência, há duas consequências explícitas: o rank pior é incluído na média e o fator de cobertura reduz o resultado.

**Importante:** esta fórmula depende da média de q*, não de uma variância ou de um termo específico para o pior alvo. Ela não garante um desempenho mínimo em cada alvo. Dois perfis com a mesma média de ranks e cobertura produzem o mesmo R, mesmo que um seja mais desigual. Esse limite é demonstrado no exemplo 4.

## 7. Consenso entre Vina e DOCK6

As relevâncias globais dos dois motores são combinadas por média aritmética de pesos iguais:

```text
C_i = (R_i,Vina + R_i,DOCK6) / 2
```

Maior C implica maior prioridade segundo a política implementada. Não existe uma segunda normalização, teste de hipótese ou penalidade explícita de discordância entre motores. Um motor pode compensar parcialmente o outro: `(0.9 + 0.1)/2` e `(0.5 + 0.5)/2` produzem 0.5. Consenso não significa, por si só, concordância entre as poses ou entre as energias.

A ordenação final usa C decrescente e, em empate exato, o identificador do composto em ordem textual crescente. O código não calcula intervalo de confiança, erro-padrão ou significância estatística para C.

## 8. Exemplos numéricos completos

### Exemplo 1 — um alvo e três compostos

Entradas brutas, com menor valor melhor dentro de cada motor:

| Composto | Score Vina | Score DOCK6 |
| --- | --- | --- |
| A | -9 | -30 |
| B | -8 | -40 |
| C | -7 | -20 |

A população de cada motor tem N=3. Aplicando ranks, exponencial e média:

| Composto | q Vina | q DOCK6 | R Vina | R DOCK6 | Consenso C |
| --- | --- | --- | --- | --- | --- |
| A | 0.0 | 0.5 | 1.000000 | 0.006738 | 0.503369 |
| B | 0.5 | 0.0 | 0.006738 | 1.000000 | 0.503369 |
| C | 1.0 | 1.0 | 0.000045 | 0.000045 | 0.000045 |

A e B empatam no fator; a ordenação textual coloca A antes de B. O resultado tem três linhas. Na UI, C aparece como `0.000` por arredondamento, embora seu valor interno seja positivo. Não compare `-9` diretamente com `-30`: seus ranks foram calculados em populações e motores distintos.

### Exemplo 2 — AChE, BChE e BACE1

Considere ranks já normalizados em bibliotecas de 11 compostos por motor/alvo, de modo que incrementos de 0.1 são possíveis. Para um composto X:

| Alvo | q Vina | p Vina | q DOCK6 | p DOCK6 |
| --- | --- | --- | --- | --- |
| AChE | 0.0 | 1.000000 | 0.1 | 0.367879 |
| BChE | 0.1 | 0.367879 | 0.2 | 0.135335 |
| BACE1 | 0.2 | 0.135335 | 0.3 | 0.049787 |

A cobertura é completa nos dois motores:

```text
R_X,Vina  = exp(-(0.0 + 0.1 + 0.2)/(3 × 0.1)) = exp(-1) = 0.367879
R_X,DOCK6 = exp(-(0.1 + 0.2 + 0.3)/(3 × 0.1)) = exp(-2) = 0.135335
C_X       = (0.3678794412 + 0.1353352832)/2 = 0.2516073622
```

A tabela exibe aproximadamente Vina `0.368`, DOCK6 `0.135` e Consenso `0.252`. As três avaliações de X se tornam **uma linha**, com detalhes por alvo preservados em `target_scores`.

### Exemplo 3 — um alvo ausente

Um composto tem q=0 em AChE/BChE, mas não foi avaliado em BACE1, em determinado motor:

```text
q* = [0, 0, 1]
c  = 2/3
R  = (2/3) × exp(-(0 + 0 + 1)/(3 × 0.1))
R  = 0.0237826622
```

Se isso ocorre nos dois motores, C também vale 0.0237826622. Se o outro motor tivesse cobertura completa e q=0 nos três alvos, seu R seria 1 e C seria 0.5118913311. Portanto, a ausência é penalizada, mas não equivale a exclusão automática quando existe pelo menos um par comum. Se nenhum alvo do composto for compartilhado pelos dois motores, ele não é elegível, ainda que cada motor tenha avaliado o composto em alvos diferentes.

### Exemplo 4 — limite da política de amplitude

| Perfil, com cobertura completa | q por alvo | Média de q | R |
| --- | --- | --- | --- |
| Equilibrado | 0.3, 0.3, 0.3 | 0.3 | 0.049787 |
| Desigual | 0.0, 0.0, 0.9 | 0.3 | 0.049787 |

Não há penalização adicional de dispersão: os dois valores são iguais. Para exigir relevância mínima em todos os alvos, é necessário avaliar os detalhes por alvo e estabelecer um critério científico separado; esse filtro não é aplicado silenciosamente pelo consenso.

## 9. Estrutura da tabela final e auditoria

| Coluna/campo | Significado |
| --- | --- |
| `molecule_chembl_id` | Uma referência por composto elegível |
| `vina` | R_i,Vina global, não a energia bruta |
| `dock6` | R_i,DOCK6 global, não a energia bruta |
| `normalized_score`, `consensus` | Mesmo C_i final |
| `consensus_model` | Identidade do modelo `target-balanced-ecr-v1` |
| `target_count` | Número T de alvos na união das entradas |
| `paired_target_count` | Quantos alvos desse composto têm avaliações nos dois motores |
| `vina_target_coverage`, `dock6_target_coverage` | Cobertura de cada motor, entre 0 e 1 |
| `target_scores` | JSON com alvo, raw scores, ranks normalizados e relevâncias de cada motor; ausências são null |
| `receptor_id`, `vina_raw`, `dock6_raw` | Contexto representativo salvo, não resultados globais |
| Campos de pose/receptor/footprint | Referências estruturais para as ações existentes |

O alvo representativo é escolhido entre os alvos compartilhados do composto, pela maior soma das relevâncias dos dois motores naquele alvo; o identificador do alvo resolve empates. Essa escolha fornece o contexto estrutural salvo e **não modifica** os fatores globais.

Vina relevance, DOCK6 relevance e Consensus mostram três casas decimais. O cálculo, a ordenação e o CSV conservam a precisão original. Os botões existentes são mantidos: com vários alvos, 3D Vina, 3D DOCK6 e footprint abrem uma combobox pesquisável de alvo; com um alvo, abrem diretamente. O arquivo escolhido deve corresponder ao composto e alvo. Ausência de anexo não autoriza substituir outro alvo nem executar novamente docking para criar o anexo.

## 10. Por que não z-score, min-max ou uma soma de energias?

O modelo atual foi definido para usar ordem relativa em vez da escala energética. Min-max depende das amplitudes/extremos; z-score depende da média e do desvio da população; uma soma direta exige calibração entre motores. Essas alternativas não são calculadas nesta implementação. Isso não prova que o ranking exponencial supere todos os métodos em todos os alvos: a adequação deve ser validada com dados do estudo.

`repulsion_weight` é um argumento legado sem efeito em `consensus_rows`; não muda os pesos dos motores nem aplica correção de repulsão ao fator atual. τ e os pesos 1/2 são fixos. O resultado não inclui bônus de estereoquímica, lipofilicidade, QED, toxicidade ou similaridade estrutural. Compatibilidade estereoquímica verifica identidade quando há SMILES, não uma pontuação extra.

## 11. Como interpretar e validar cientificamente

1. Confira se identificadores de composto/alvo, sinais dos scores e populações de referência estão corretos.
2. Inspecione `target_count`, coberturas e alvos compartilhados; não interprete ausência como medida de baixa afinidade experimental.
3. Use os valores não arredondados para seleção automática; a UI pode mostrar `0.000` para fatores positivos pequenos.
4. Inspecione `target_scores` e poses dos alvos de interesse. O fator global não substitui análise de seletividade, contatos ou qualidade geométrica.
5. Para avaliar utilidade preditiva, use ativos/decoys ou referências experimentais adequadas ao estudo e métricas de enriquecimento previamente definidas. O bloco não calcula essa validação automaticamente.
6. Evite comparar C entre execuções com bibliotecas/targets diferentes sem avaliar a mudança de população. Não traduza C em probabilidade de sucesso clínico ou afinidade absoluta.

## 12. Reproduzir a agregação com Python padrão

O exemplo abaixo recebe **q já normalizado**, não energias brutas. Reproduz o exemplo 2 usando somente a biblioteca padrão; não é uma nova API nem executa docking. Os vetores dos dois motores devem representar a mesma união de alvos, na mesma ordem; use `None` para observação ausente.

```python
from math import exp, fsum

def relevance(q_values):
    if not q_values:
        raise ValueError("A união de alvos deve ser não vazia")
    if any(q is not None and not 0 <= q <= 1 for q in q_values):
        raise ValueError("Informe ranks normalizados entre 0 e 1")
    coverage = sum(q is not None for q in q_values) / len(q_values)
    average_rank = fsum(1.0 if q is None else q for q in q_values) / len(q_values)
    return coverage * exp(-average_rank / 0.1)

vina = relevance([0.0, 0.1, 0.2])
dock6 = relevance([0.1, 0.2, 0.3])
consensus = (vina + dock6) / 2
print(vina, dock6, consensus)
```

Resultado aproximado: `0.3678794412 0.1353352832 0.2516073622`. Os quatro exemplos e suas variantes de cobertura foram conferidos contra a função atual `consensus_rows`, sem executar Vina, DOCK6 ou modificar o código da plataforma. Trata-se de verificação aritmética da documentação, não de validação preditiva ou experimental.

## 13. Referências e origem das decisões

- [ECR, Palacio-Rodríguez et al. (2019)](https://doi.org/10.1038/s41598-019-41594-3): fundamento de combinar ranks exponencialmente.
- [Combinação de consenso de poses e ranks (2021)](https://pubs.rsc.org/en/content/articlehtml/2021/ra/d1ra05785e): exemplo do uso da escala exponencial relativa à biblioteca.
- `consensus_model.py`, `consensus_io.py` e `wrappers/docking.py`: fonte da verdade para a adaptação, deduplicação, cobertura, exportação e modos de finalização implementados.

A normalização em [0,1], a agregação geométrica entre alvos e a penalização de cobertura são escolhas explícitas do BioMolExplorer. A documentação distingue essas escolhas dos resultados publicados na literatura; não apresenta o modelo inteiro como fórmula original do ECR.

## Persistência, desempenho e evidências anteriores

O modelo possui identidade científica própria: consensos legados não podem ser
reutilizados como resultados deste cálculo. Mudanças de entradas invalidam o
cache pelos contratos existentes. Resultados Vina/DOCK6 e suas preparações
permanecem reutilizáveis, pois seus algoritmos não foram alterados.

Os testes verificam fórmula fechada, sinais positivos/negativos, invariância de
escala, empates, singleton, cobertura ausente, duplicação, ordem de entrada,
estereoisômeros, contexto estrutural e cardinalidade de 1.000 × 3 alvos.
Essas verificações não substituem validação prospectiva com ativos e decoys.

As colunas visíveis de relevância Vina, relevância DOCK6 e consenso apresentam
três casas decimais. Os dados usados para ordenar, calcular e exportar o CSV
conservam sua precisão original; a formatação é aplicada apenas às células.

## Execução e finalização

O transporte das entradas conserva referências absolutas já autorizadas às
estruturas produtoras. Não copia todas as poses para cada entrada; somente as
poses representativas são materializadas na saída do consenso. A entrada Vina
e seu alias legado base_input_path compartilham uma única preparação. Um
manifesto de transporte registra os arquivos estruturais referenciados e suas
assinaturas de dispositivo, inode, tamanho, mtime e ctime, para detectar mudanças.

Quando existe docking_results.csv na raiz da entrada, essa tabela é a fonte
autoritativa, inclusive vazia: a leitura não percorre as árvores de poses.
Diretórios legados continuam usando a resolução anterior como alternativa.

A tabela é publicada atomicamente antes do gráfico. A figura utiliza o backend
Agg diretamente, sem fábrica de janelas; falhas de visualização registram um
consensus_report.json e não invalidam a tabela calculada. Gráficos incompletos
ou antigos são removidos quando essa geração falha.

O worker informa leitura, normalização, agregação, exportação e conclusão.
Atualizações intermediárias persistidas do consenso são limitadas a uma a cada
dez segundos. A confirmação terminal do worker segue diretamente para finalizar
o bloco, evitando regravar todo o histórico como running apenas para publicar
a última mensagem. As demais operações mantêm seu fluxo de atualização.

Validação local: um worker real concluiu com 1.729 compostos em aproximadamente
6,4 segundos, usando scores e anexos existentes, com equivalência dos fatores
calculados e hashes inalterados dos CSVs originais. Essa medição cobre o worker,
não a renderização completa da interface ou a publicação de todo o pipeline.

## Consulta do histórico e uso de memória

A busca de reuso em pipelines com consenso lê provas indexadas por execução,
mantendo no máximo duas execuções decodificadas em memória. O primeiro resultado
válido encerra a busca; somente seus descritores de apresentação são carregados.
A consulta inicial de resultados e a prévia de invalidação também usam essas
provas. As verificações de hashes, parâmetros, entradas, permissões e proveniência
continuam obrigatórias.

A resposta de aceitação do pipeline é capturada antes de iniciar o agendador,
evitando que a resposta da fila incorpore descritores volumosos de resultados
reutilizados. Durante a validação do consenso o bloco aparece como em execução,
com uma mensagem específica de verificação dos resultados disponíveis.

Na auditoria de pipelines com consenso, o JSON integral de cada execução é
preservado como conteúdo imutável em `.history/run-stages`, transferido do SQLite
em blocos de 1 MiB. Execuções inalteradas reutilizam seus conteúdos já arquivados.
O snapshot anterior é referenciado sem ser decodificado novamente. Snapshots
anteriores e arquivos científicos permanecem intactos. Restauração verifica o
SHA-256 e reconstrói os mesmos dados; a exportação explícita converte os snapshots
compactos para o formato portátil já existente. A primeira conversão do histórico
tem custo de leitura e armazenamento; execuções seguintes reaproveitam os blobs.

Validação adicional: um histórico real com 18 execuções e aproximadamente 3,4 GB
de JSON foi consultado sem escrita no banco original, recuperando o cache a partir
de uma execução. Um pipeline isolado com scores de 1.000 compostos em três alvos
produziu exatamente 1.000 linhas e terminou em aproximadamente 4,2 segundos; a
reexecução com cache levou 0,28 segundo, sem iniciar outro worker. Esses tempos
cobrem o pipeline de teste, incluindo auditoria, e não a renderização no navegador.
