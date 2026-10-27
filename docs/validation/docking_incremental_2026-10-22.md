# Reuso e execução incremental de redocking/docking

## Causas corrigidas

O cache de etapas completas não cobria a mudança parcial de um lote. No Vina,
a simples existência de uma pasta ou pose podia encerrar a execução ou pular um
experimento sem conferir suas entradas e parâmetros. Preparações de referência
e checkpoints de redocking dependiam do sucesso do run completo. Alguns
resultados parciais também podiam ser aceitos como lote concluído pelo scheduler.

O reuso passa a combinar o atalho de etapa/lote já existente com checkpoints
por preparação, sítio e cálculo. Arquivos históricos não são usados apenas por
nome. A tabela atual contém somente a seleção solicitada e os cálculos
concluídos dessa requisição; o histórico anterior é preservado.

## Implementação

| Componente | Unidade e dependências de reuso |
| --- | --- |
| `calculation_cache.py` | SQLite por projeto; experimentos, entradas identificáveis, artefatos e objetos SHA-256; transações e locks |
| `structure_cache.py` | Receptor e referência independentes; PDB, cadeia, opções, cofatores, pH, templates e ferramentas |
| `redocking_cache.py` | Par completo, pose e RMSD; publicação incremental independente do estado final do run |
| `vina_batch.py` | Conteúdo do receptor/ligante, configuração efetiva e identidade do Vina; validação numérica e de score |
| `dock6_cache.py` | Sítio/grids compartilhados; minimização, docking e footprint por composto; templates e binários/parâmetros DOCK6 |
| Scheduler | Cache integral e por lote, entradas efetivamente materializadas e identidade das ferramentas; repetição dos itens que falharam |

O índice fica em `.biomolexplorer/cache/calculations/index.sqlite3`. As tabelas
`experiments`, `subjects`, `artifacts` e `objects` permitem localizar a origem,
o protocolo, o composto/receptor e os arquivos de um experimento. Os objetos
são armazenados uma vez por conteúdo no cache. Reflinks, quando disponíveis,
reduzem o custo de materialização e conservam a independência das saídas
editáveis; em outros sistemas de arquivos usa-se cópia independente. Há no
máximo 64 arquivos de lock, em vez de um arquivo por experimento.

Exemplo de consulta somente de leitura, para manutenção:

```sql
SELECT e.namespace, s.role, s.identifier, a.role AS artifact_role, a.sha256
FROM subjects AS s
JOIN experiments AS e ON e.namespace = s.namespace AND e.key = s.experiment_key
JOIN artifacts AS a ON a.namespace = e.namespace AND a.experiment_key = e.key
WHERE s.identifier = 'PUBCHEM10367945';
```

Os hashes de grids e executáveis são memorizados por assinatura do arquivo
durante o lote. O Vina verifica hits antes de iniciar o worker nativo. Uma
requisição totalmente recuperável não importa o motor, não gera mapas e não
inicia subprocessos Vina. A opção existente de executar novamente desativa as
leituras dos checkpoints; entradas preparadas fornecidas pelo usuário continuam
sendo entradas do protocolo, sem uma preparação científica implícita adicional.

## Validação

Os testes usam o store/scheduler reais e arquivos temporários. Os adaptadores
Chimera/DOCK6 são substituídos na fronteira dos comandos externos; validação,
índice, materialização, seleção e montagem das tabelas são executados de fato.

A suíte focada final concluiu **190 testes, todos aprovados**, em 122,296 s.
Ela cobre os módulos incrementais, preparação, interoperabilidade, cache de
etapas, execução do pipeline, processamento de entradas, progresso, diagnósticos,
cancelamento e localização. A validação estática confirmou sintaxe e ausência
de erros de whitespace no diff. O verificador da documentação validou arquivos,
âncoras e IDs de **32 páginas HTML**.

| Cenário | Comportamento verificado |
| --- | --- |
| Vina A+B → A+B → A+B+C | Somente A e B na primeira chamada, nenhum novo cálculo na repetição, somente C na ampliação |
| Mudança de um ligante | Recalcula somente o composto alterado |
| Mudança de parâmetros/receptor | Invalida os cálculos dependentes, sem reduzir parâmetros científicos |
| DOCK6 A+B → A+B+C | Mantém sítio, grids, minimizações e dockings já concluídos; processa C |
| DOCK6 flex → rigid | Reutiliza sítio e minimizações; refaz o docking |
| Novo PDB/par de referência | Prepara e calcula o novo item; preserva os resultados selecionados anteriores |
| Falha de B/interrupção após A | Mantém o checkpoint de A e tenta B novamente |
| Requisições concorrentes iguais | Uma execução por chave, com saídas independentes |
| Objeto corrompido/arquivo ausente | Recusa o checkpoint afetado |
| Reuso desativado | Executa novamente os cálculos solicitados |
| Remoção de seleção/resultado falho | Mantém os pares restantes sem recolocar itens removidos na tabela atual |

Também foi executado Vina real, via processo filho, para dois ligantes pequenos
com tipos atômicos distintos, CPU=1, seed=73, exhaustiveness=1 e uma pose. A
primeira chamada concluiu ambos; a repetição recuperou poses byte a byte
idênticas sem iniciar outro worker. Esse teste verifica integração e reuso,
sem extrapolar tempos para uma campanha real ou afirmar equivalência científica
entre protocolos diferentes. Não foi iniciada uma campanha extensa.

Os manuais e a arquitetura foram atualizados em português/inglês, com regeneração
das páginas HTML e atualização do `index.html`.
