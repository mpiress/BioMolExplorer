# Consenso de docking por composto

Modelo `target-balanced-ecr-v1`. Cada composto produz **uma única linha**, mesmo
quando avaliado em vários alvos. Os botões existentes da tabela são preservados.
O cálculo não executa novamente docking, preparação, minimização ou footprint.

## Fundamentação

Usamos uma adaptação normalizada do Exponential Consensus Ranking (ECR), robusta
a unidades, deslocamentos, escalas e valores extremos das energias. Ordenamos
sempre os scores brutos em ordem crescente, sem aplicar valor absoluto: energias
positivas grandes não podem ser promovidas artificialmente.

ECR tem suporte em benchmarks, mas **não existe normalização universalmente
superior para todos os alvos**. O parâmetro de cauda de 10% segue a prática ECR;
a média geométrica entre alvos e a penalização de cobertura abaixo são decisões
explícitas deste modelo para priorizar amplitude terapêutica. Não constituem uma
alegação de superioridade experimental ou uma probabilidade de ligação.

- Palacio-Rodríguez et al. (2019), [ECR original](https://doi.org/10.1038/s41598-019-41594-3).
- [Combinação de consenso de poses e ranks](https://pmc.ncbi.nlm.nih.gov/articles/PMC8965822/) (2022).
- [Comparação de normalizações e consensos](https://pmc.ncbi.nlm.nih.gov/articles/PMC9941253/) (2023).
- [DockM8](https://doi.org/10.1186/s13321-026-01287-2) (2026): resultados dependem do alvo e do protocolo.

## Normalização por motor e alvo

Sejam i o composto, e o motor (Vina ou DOCK6) e t o alvo. Primeiro selecionamos
o menor score bruto de cada par (i,t,e); repetir arquivos, poses ou execuções não
aumenta seu peso. Ordenamos os N_et compostos observados nesse motor/alvo.
Empates recebem rank médio r_iet, começando em 1:

```text
q_iet = (r_iet - 1) / max(N_et - 1, 1)
p_iet = exp(-q_iet / tau),   tau = 0,10
```

q varia de 0 (melhor) a 1 (pior); p varia de exp(-10) a 1. Trata-se de uma
versão deslocada e reescalada do ECR, não da sua densidade original. Todos os
scores entram no consenso exclusivamente após essa normalização. A transformação
é invariante a transformações estritamente crescentes dos scores.
Um singleton recebe q=0 e p=1, sem implicar evidência estatística de afinidade;
uma distribuição constante com N>1 recebe q=0,5 para todos. Os tamanhos das
bibliotecas e a cobertura devem ser considerados na interpretação do ranking.

## Agregação entre alvos e motores

T é a união dos identificadores de receptores presentes nas entradas. Todos os
alvos e ambos os motores têm pesos iguais, sem pesos inferidos a partir das
energias. Para cada motor, c_ie é a fração de alvos com score observado.
Atribuímos q*=1 aos scores ausentes, sem inventar energia bruta:

```text
R_ie = c_ie × exp(-sum_t(q*_iet) / (T × tau))
C_i  = (R_i,Vina + R_i,DOCK6) / 2
```

Com cobertura completa, R é a média geométrica das relevâncias por alvo. Isso
penaliza compostos fortes em apenas um alvo e fracos nos demais. A cobertura
evita favorecer compostos com avaliações ausentes. Em alvo único completo,
C=(p_Vina+p_DOCK6)/2: consenso ECR direto com pesos iguais. C está entre 0 e 1;
**maior é melhor**. Não se aplica outra normalização após agregar.

A lista inclui os compostos com ao menos um par alvo/composto avaliado pelos
dois motores, preservando o contrato de correspondência entre entradas.
Compostos sem qualquer par comum são excluídos; ausência total de pares mantém
o comportamento de bloco não executado com explicação. Coberturas parciais
continuam visíveis no CSV e não são confundidas com dados completos.
Para 1.000 compostos avaliados nos mesmos três alvos pelos dois motores,
as 6.000 avaliações geram exatamente **1.000 linhas**.

## Identidade e estereoquímica

O agrupamento usa o identificador molecular existente, sem fundir códigos
distintos pela semelhança dos nomes ou SMILES. Validamos o SMILES canônico
isomérico em todos os alvos e motores: o mesmo código associado a enantiômeros
ou estruturas distintas causa erro. Não inferimos configuração estereoquímica
ausente, nem convertemos concordância energética em qualidade geométrica.
RMSD de poses, atividade experimental, seletividade e toxicidade exigem
validação própria e não recebem pesos arbitrários neste consenso.

## Tabela e rastreabilidade

| Coluna visível | Campo CSV | Conteúdo |
| --- | --- | --- |
| Código do composto | molecule_chembl_id | Uma referência por composto |
| Relevância Vina | vina | R_i,Vina, normalizado e agregado nos alvos |
| Relevância DOCK6 | dock6 | R_i,DOCK6, normalizado e agregado nos alvos |
| Consenso | normalized_score | C_i, também exportado como consensus |
| Ações | Campos estruturais existentes | Mesmos botões e mecanismos de acesso |

O CSV inclui consensus_model, canonical_smiles, target_count,
paired_target_count, vina_target_coverage, dock6_target_coverage e target_scores.
Este último é JSON com energias brutas, ranks normalizados e relevâncias de
**todos** os alvos; valores ausentes são null. Não arredondamos o cálculo.
z-score e min-max deixam de ser produzidos pelo novo modelo; leitores de
resultados antigos continuam aceitando esses campos.

As ações 3D, resíduos e footprint referem-se ao **mesmo alvo representativo**:
o alvo com melhor consenso local entre os pares completos, com desempate pelo
identificador do receptor. receptor_id, vina_raw e dock6_raw descrevem esse
contexto, sem serem tratados como scores globais. Fontes estruturais originais
permanecem nas execuções produtoras; o consenso copia apenas o contexto
representativo, reutilizando anexos compartilhados. Ausência de anexos não
impede calcular scores. A remoção continua auditada pelos serviços existentes.

## Cache e validação

O modelo possui identidade científica própria: consensos legados não podem ser
reutilizados como resultados deste cálculo. Mudanças de entradas invalidam o
cache pelos contratos existentes. Resultados Vina/DOCK6 e suas preparações
permanecem reutilizáveis, pois seus algoritmos não foram alterados.

Os testes verificam fórmula fechada, sinais positivos/negativos, invariância de
escala, empates, singleton, cobertura ausente, duplicação, ordem de entrada,
estereoisômeros, contexto estrutural e cardinalidade de 1.000 × 3 alvos.
Essas verificações não substituem validação prospectiva com ativos e decoys.
