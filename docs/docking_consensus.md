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

Com múltiplos alvos, as ações 3D Vina, 3D DOCK6, resíduos e footprint abrem
um popup com combobox pesquisável antes de apresentar o resultado. A opção
Continuar é habilitada após selecionar um alvo. A lista não valida poses
antecipadamente; somente o resultado escolhido é validado ao abrir. Cada ação recupera a pose,
o receptor ou o PDF correspondente ao composto e alvo escolhidos nas tabelas
de origem, inclusive nas execuções anteriores do mesmo projeto que
produziram entradas reutilizadas, validando as referências e os hashes
disponíveis. Anexos ausentes ou alterados impedem abrir aquele resultado; nunca são substituídos pelo alvo representativo.
Com um único alvo, o acesso permanece direto.

O consenso mantém apenas uma linha por composto. receptor_id, vina_raw e
dock6_raw descrevem o contexto representativo salvo, sem serem tratados como
scores globais. A seleção utiliza resultados já existentes, sem recalcular o
docking nem alterar os arquivos persistidos. A remoção continua auditada pelos
serviços existentes.

## Cache e validação

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
