# Conjuntos de saída e processamento

As seleções usam nomes públicos na convenção
`TG[ALVO]-SRC[ORIGEM]-OP[OPERAÇÃO]-TIPO-ID[hash].ext`. O identificador curto
distingue ocorrências; a identidade interna e a linhagem completa permanecem
nos metadados. `TG[NA]` indica alvo desconhecido ou não aplicável.
As publicações físicas anteriores `Bnn__operação__conteúdo.ext` e os nomes
operacionais dos motores continuam compatíveis, sem renomeação destrutiva.
Consulte a [convenção compartilhada](artifact_naming.md).

| Bloco | Identificação dos conjuntos |
| --- | --- |
| Importação | Nome original e origem do arquivo enviado |
| ChEMBL | MOL: moléculas recuperadas; SIM: novos similares, excluindo IDs de MOL |
| PubChem | SIM: novos similares, sem repetir os compostos de referência |
| Expansão de similares | MOL/SIM para conjuntos separados; nomes anteriores continuam aceitos |
| ZINC | MOL e arquivos de estruturas/conformações |
| PDB | PDB e identificador da estrutura; complexos e ligantes |
| ADMET | BBB+, BBB−, HIA+, HIA−; todos os avaliados, selecionados e fora dos critérios |
| Fingerprints | Morgan, MACCS ou farmacóforo, com o conjunto de entrada identificado |
| Similaridade | Métrica utilizada e conjunto de entrada |
| Grafos | MCC/LCC, componente, grafo completo e relações |
| Preparação | REC: receptor; LIG: ligante; centros de docking |
| Redocking | RMSD e arquivos do complexo correspondente |
| Vina/DOCK6 | SCORE: resultados; POSE: pose com identificador molecular |
| Consenso | RANK: consenso dos resultados Vina e DOCK6 |
| Amber Score · DOCK6 | Score, contribuições energéticas, topologias/coordenadas, protocolo e comparação inicial/final disponível; não há dinâmica convencional nem RMSF |

BBB e HIA são as classificações previstas pelo método ADMET configurado. Os
relatórios preservam os métodos e parâmetros usados.

## Escolher como processar

- **Individual:** analisa separadamente cada arquivo selecionado. Havendo vários
  campos de entrada, processa suas combinações.
- **Merge:** reúne as entradas selecionadas e analisa o conjunto completo.
- **Individual + conjunto completo (ALL):** mantém as análises individuais e
  acrescenta a análise conjunta quando existe mais de uma combinação. A seleção
  de saída identifica essa análise adicional com ALL.

Consenso e dinâmica têm regras próprias para reunir resultados ou pares e
reaproveitar seus cálculos; o formulário informa esse comportamento.

A recuperação ChEMBL por alvo publica MOL e SIM sem criar `_FULL` ou outra cópia
normalizada contendo a mesma união. Os registros originais do provedor e os
relatórios continuam disponíveis para rastrear a recuperação. Buscas diretas e
opções legadas de expansão mantêm compatibilidade com seus arquivos existentes.
Para analisar MOL e SIM juntos, selecione ambos no próximo bloco e escolha Merge
ou Individual + ALL.

Se uma configuração antiga aponta para um arquivo que deixou de existir (por
exemplo, FULL após uma nova recuperação), o pipeline solicita uma nova seleção.
As etapas produtoras já concluídas são mantidas; confirmar a seleção continua a
execução a partir do bloco pendente.

Lotes já calculados são reutilizados por conteúdo e configuração. ALL de
fingerprints é construído a partir dos vetores existentes, sem outro worker.
Descritores ADMET são armazenados por estrutura, versão do RDKit, implementação,
catálogos de alertas e uso do método legado; os critérios de seleção são
reaplicados a cada execução. Modelos ADMET-AI continuam com sua execução própria.
Similaridade e grafos do conjunto completo são análises próprias: podem incluir
relações entre moléculas de arquivos diferentes e não equivalem à união das
análises individuais.

## Excluir resultados

Após excluir os resultados, os cartões, ícones e detalhes indicam **Resultados
excluídos**, sem selos de conclusão ou reaproveitamento dos arquivos removidos.
O resumo diferencia os resultados disponíveis do histórico de etapas concluídas.
Entradas enviadas, configurações, logs e histórico permanecem preservados. Uma
nova execução não reutiliza os artefatos excluídos.
