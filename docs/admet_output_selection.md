# Saídas ADMET por arquivo de entrada

Cada entrada é avaliada separadamente, preservando alvo, base e origem MOL/SIM.
O bloco publica as seguintes tabelas:

| Identificação | Conteúdo |
| --- | --- |
| BBB+ | Compostos classificados como BBB+ |
| BBB− | Compostos classificados como BBB− |
| HIA+ | Compostos classificados como HIA+ |
| HIA− | Compostos classificados como HIA− |
| Regras · Lipinski, Veber, Lipinski + Veber ou Fragmentos | Um arquivo da regra configurada, com os compostos que passam pelo limite de violações |
| Alertas químicos | Uma tabela de todos os compostos válidos, com flags `Alert_PAINS`, `Alert_BRENK` etc. em `true`/`false`, descrições dos alertas e os descritores avaliados |

BBB/HIA são partições independentes da mesma entrada: a combinação de hipóteses
não remove compostos de uma classificação diferente. `ADMET_selected` e
`ADMET_reasons` preservam o resultado da hipótese combinada em cada registro.
No modo combinado oral, Veber exige zero violações e o limite informado se
aplica a Lipinski. O arquivo de regras não depende de BBB/HIA nem dos alertas.

As quatro tabelas BBB/HIA são publicadas mesmo vazias, mantendo seus cabeçalhos.
Sem regra ou sem catálogos selecionados, o arquivo correspondente não é gerado.
Não são criadas cópias genéricas `compounds`, `evaluated`, `selected` ou
`rejected` nas novas execuções. Os compostos válidos e os descritores continuam
nas partições e na tabela de alertas; exclusões de estruturas inválidas seguem
no relatório de qualidade. O relatório JSON registra critérios, métodos,
versões, entrada e os caminhos exatos das saídas, inclusive após publicação.

A interface mostra BBB+/BBB− e HIA+/HIA−. Os nomes físicos curtos seguem o padrão
do projeto, por exemplo `B04__admet__bbb-pos.csv`,
`B04__admet__rules-lipinski.csv` e `B04__admet__alerts.csv`. Diretórios e
metadados distinguem as entradas. Downloads incluem a identificação da entrada;
os seletores preservam o alvo e MOL/SIM nos rótulos.

Em Resultados, **Ver alertas químicos** abre uma tabela paginada e pesquisável,
com código da molécula, uma coluna booleana por catálogo, os alertas encontrados
e um botão para visualizar a estrutura 2D. `true` significa correspondência ao
catálogo, não toxicidade comprovada. A revisão de alertas de uma partição BBB ou
HIA usa a tabela completa da mesma entrada, sem misturar outros alvos.

Os blocos seguintes listam somente saídas processadas. Até 100 arquivos ficam
no mesmo menu pesquisável; coleções maiores continuam paginadas. Por exemplo,
seis entradas de ACHE/BCHE/BACE1 com Lipinski e alertas produzem 36 CSVs, sendo
12 de BACE1. Arquivos e referências históricos continuam legíveis; a mudança
do contrato vale para novas execuções. Descritores compartilhados continuam
em cache, enquanto os critérios da configuração atual são reaplicados.
