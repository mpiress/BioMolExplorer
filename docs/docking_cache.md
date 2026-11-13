# Reutilização de docking concluído

Vina e DOCK6 possuem dois níveis de reaproveitamento: o experimento publicado
completo e os checkpoints de preparação/cálculo. Se o primeiro nível é válido,
o agendador não materializa novas entradas, não inicia um worker e não chama
Chimera, Vina ou os programas de grid/minimização/docking do DOCK6.

## Causa da repetição

Os três caminhos de reutilização — resultados publicados, cache de estágio e
cache de lote — rejeitavam qualquer contagem positiva de `docking_failures` ou
`excluded_records`. Esses contadores também incluem moléculas excluídas durante
a preparação. Assim, resultados concluídos e íntegros eram tratados como
execuções incompletas, mesmo com entradas e parâmetros idênticos.

Na execução analisada, Vina tinha 640 compostos excluídos na preparação, em
três receptores: 1.920 pares registrados como `preparation_failed`, nenhum par
sem explicação e 5.451 pares com scores. DOCK6 tinha 376 exclusões distribuídas
em três lotes. Os relatórios já explicavam essas exclusões; a decisão de cache
não considerava essa informação.

## Critérios de aceitação

A validação comum em `docking_reuse.accounted_preparation_exclusions` permite
reaproveitar esse resultado somente quando:

- O estágio e seus lotes terminaram com sucesso, sem remoção, obsolescência ou
  marcação de resultados parciais.
- Os hashes dos arquivos correspondem ao manifesto original.
- Existem `docking_accounting.json`, `docking_unprocessed.csv` e
  `docking_results.csv` para cada lote validado.
- Os pares esperados equivalem à soma dos pares com scores e dos pares
  excluídos, com `unaccounted_pairs = 0`.
- Todo par sem score está explicitamente classificado como
  `preparation_failed`, com identificador, receptor e motivo.
- Scores são finitos, os pares não se repetem e as contagens por receptor e
  composto correspondem ao relatório.

As verificações já existentes de configuração, entradas selecionadas,
proveniência, escopo do projeto e integridade dos demais artefatos continuam
obrigatórias. Alterar entradas ou parâmetros invalida a correspondência.

Exclusões **não são convertidas em sucessos**. Contadores, relatórios, scores e
arquivos originais são mantidos, inclusive quando o resultado é reaproveitado
novamente após reiniciar o serviço.

Falhas do motor (`docking_failed`), pares sem explicação, relatórios ausentes,
arquivos alterados e execuções interrompidas continuam impedindo a reutilização
integral. Nesses casos, o fluxo anterior permanece e pode recuperar checkpoints
válidos de preparação e cálculos já terminados. Resultados sem exclusões mantêm
os critérios de cache já existentes.

O contador isolado de exclusões de qualidade não é interpretado como falha do
motor. A busca de cache e sua prévia de invalidação usam o histórico indexado
também em pipelines com Vina/DOCK6 sem consenso, evitando carregar simultaneamente
os JSONs completos de todas as execuções. A auditoria reutiliza o armazenamento
incremental de execuções já utilizado pelo consenso.

Para tentar novamente compostos excluídos, execute sem reaproveitar os resultados
ou altere as entradas/opções pertinentes. Não é necessário apagar scores válidos.

## Validação

`tests/test_docking_completed_reuse.py` simula um pipeline Vina → DOCK6 já
concluído e reinicia o serviço duas vezes. O teste impede tanto a criação de
`JobManager` quanto nova materialização de entradas, verifica a reutilização dos
dois blocos, conserva as exclusões e compara os bytes dos artefatos originais.
Também verifica os caches de estágio/lote e rejeita alterações de configuração,
entradas, relatórios, scores e falhas do motor.

Uma verificação somente de leitura, usando os resultados e as configurações
atuais do projeto Alzheimer, confirmou a reutilização dos dois blocos, sem
iniciar workers ou regenerar entradas. Ela não executou novos motores científicos.

```sh
PYTHONPATH=src:tests python -m unittest test_docking_completed_reuse \
  test_docking_incremental test_docking_accounting test_published_results \
  test_preparation_failures test_vina_reuse test_dock6_operation_modes \
  test_redocking_incremental
```
