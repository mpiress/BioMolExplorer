# Diagnóstico da execução Vina de 21–22/10/2026

## Resultado

A falha não foi causada por lentidão ou indisponibilidade do AutoDock Vina. O
job `4d3b0940689a4a35bda4bab8ebd7d15d` terminou porque o PDBQT preparado para
`PUBCHEM10367945` continha `nan` em coordenadas e cargas. A biblioteca nativa do
Vina recebeu esse arquivo e encerrou todo o worker com `internal_error`, sem
lançar uma exceção Python recuperável.

O problema foi reproduzido no arquivo real e o novo validador o rejeita antes da
chamada nativa com uma mensagem que identifica o composto. O cache também passa
a recusar entradas antigas com valores não finitos.

## Evidência observada

Os números abaixo foram calculados de `execution.log`, `events.jsonl` e
`diagnostic.json` do job; não são projeções:

| Medida | Valor observado |
| --- | ---: |
| Início / fim UTC | 21/10 05:33:37 / 22/10 18:43:56 |
| Tempo total | 133.819 s (37,17 h) |
| Preparações iniciadas | 733 |
| Preparações Chimera concluídas | 668 |
| Preparações publicadas | 666 |
| Falhas/avisos por composto | 67 |
| Workers de preparação | 2 |
| Soma dos tempos Chimera | 266.439,603 s |
| Chimera: média / mediana / p95 / máximo | 398,862 / 306,822 / 993,063 / 2.459,517 s |
| Construção dos mapas Vina | 1,431 s, uma vez |
| Dockings Vina concluídos | 6 |
| Vina: soma / média / mediana / máximo | 155,009 / 25,835 / 25,794 / 28,265 s |

A soma do tempo Chimera é aproximadamente duas vezes o tempo de parede da fase,
confirmando que os dois workers estavam ativos. A demora dominante era a
preparação molecular solicitada (cargas AM1 e minimização), não o docking. O job
usava pH 7,4, caixa 24 × 24 × 24 Å, `exhaustiveness=20` e 10 poses. Esses
parâmetros científicos não foram alterados.

## Correção

- PDBQT e MOL2 são validados ao serem produzidos, salvos no cache, restaurados do
  cache e imediatamente antes do Vina. Coordenadas/cargas não finitas, ausência
  de átomos e tipo atômico inválido são rejeitados.
- A preparação continua isolando falhas por composto. Moléculas válidas seguem
  para o docking e o relatório de preparação registra as excluídas.
- O lote Vina roda em processo filho com checkpoint atômico por configuração.
  Uma exceção afeta apenas o ligante corrente; até mesmo um `abort()` nativo é
  detectado pelo coordenador, que marca o composto e reinicia o lote restante.
- Receptor e mapas permanecem em uma única sessão por processo. Não há
  paralelização simultânea de dockings que competiria com as CPUs internas do
  próprio Vina.
- `vina_execution.json` reúne versão, composto, configuração, duração, falhas e
  tentativas. Eventos estruturados carregam job, etapa, configuração e composto.

O benchmark real, pequeno e controlado, da preparação permanece em
`vina_preparation_performance_2026-10-20.md`. Ele mediu 3,1490 s em série,
2,1755 s em paralelo e 0,5063 s com cache para dois compostos, com PDBQTs
idênticos. Ele valida o mecanismo, mas não deve ser extrapolado para as 733
moléculas desta execução. A campanha completa não foi reiniciada automaticamente.

## Critérios de retomada

Ao repetir a mesma etapa, preparações válidas com a mesma estrutura, pH, opções e
identidade das ferramentas são reaproveitadas pelo cache do projeto. A entrada
corrompida é recalculada; se continuar inválida, será registrada e ignorada sem
cancelar as demais. Alterações científicas invalidam a chave de cache.

