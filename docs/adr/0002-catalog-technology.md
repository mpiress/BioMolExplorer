# ADR-002 — Tecnologia do catálogo persistente de metadados

- **Status:** Aceita
- **Data:** 2026-11-28
- **Fase:** PHASE-01 (#8)
- **Requisitos:** REQ-ARCH-002, REQ-ARCH-003, REQ-CAT-002, REQ-CAT-004, AC-01-03
- **Autores e revisores:** Michel

## Contexto

O catálogo já existe e já é SQLite: `workspace.py` abre
`workspace.sqlite3`, cria as tabelas do workspace e do índice de artefatos, e
`result_catalog.create_schema` acrescenta o modelo de leitura dos resultados. A
decisão desta fase não é escolher uma tecnologia nova, e sim registrar por que a
atual é mantida, o que ela não entrega, e como o esquema passa a evoluir.

O que motiva o registro:

- **O esquema não tinha versão.** Todas as tabelas eram criadas com
  `CREATE TABLE IF NOT EXISTS`, que responde se uma tabela existe, nunca qual
  versão dela o banco contém. As correções de esquema viviam em
  `artifact_registry.migrate_schema`, disparadas pela inspeção de
  `PRAGMA table_info` — um mecanismo que funciona para renomear uma coluna uma
  vez, mas não ordena mudanças nem é idempotente por construção (REQ-CAT-002).
- **A plataforma é local por requisito.** 01.4 de #6 exige funcionar sem
  servidor de banco, serviço em nuvem ou qualquer infraestrutura obrigatória.
- **A concorrência real é modesta.** Um host, um processo de interface e
  trabalhadores em subprocesso. A escrita já é serializada por
  `BEGIN IMMEDIATE`, e o WAL, habilitado uma vez por arquivo em
  `WorkspaceStore._enable_concurrent_reads`, permite que a interface leia o
  andamento de uma execução sem esperar por ela.

## Alternativas consideradas

| Alternativa | Vantagens | Desvantagens | Evidência |
| --- | --- | --- | --- |
| A — SQLite com WAL e migrações ordenadas *(selecionada)* | Sem infraestrutura; transacional; índices e planos verificáveis; já em produção no projeto | Um único escritor por vez; sem fila distribuída; estatísticas de planejamento ausentes sem `ANALYZE` | Planos em `docs/validation/phase_01_identity_catalog_2026-11-28/README.md` |
| B — PostgreSQL | Concorrência de escrita real; tipos ricos; particionamento | Exige serviço instalado e configurado, contra REQ-ARCH-003; um projeto deixaria de ser uma pasta copiável | 01.4 de #6 |
| C — DuckDB | Excelente para varredura analítica colunar de tabelas grandes | Modelo de concorrência inadequado para escritas pequenas e frequentes de metadados durante uma execução; duplicaria o catálogo em vez de generalizá-lo (REQ-ARCH-002) | — |
| D — Armazenamento chave-valor embutido (LMDB, SQLite sem SQL) | Escrita rápida de pares | Sem consultas por índice composto nem planos verificáveis; AC-01-03 e REQ-CAT-004 pedem exatamente isso | — |
| E — Arquivos JSON com índice em memória | Nenhuma dependência | É o estado que F-04 e F-07 mediram como custo dominante; o índice se perde a cada abertura | `docs/validation/artifact_architecture_assessment_2026-11-26/README.md` |

Para os **conjuntos de dados analíticos** — tabelas de scores com milhões de
linhas, não metadados — a escolha não é esta. O custo de varredura colunar é o
problema de ADR-006 (PHASE-05), e DuckDB ou Parquet continuam candidatos lá.
Esta decisão cobre o catálogo de metadados.

## Solução selecionada

Manter SQLite e dar-lhe um mecanismo explícito de evolução.

### Migrações ordenadas, idempotentes e não destrutivas

`catalog_schema.MIGRATIONS` declara cada mudança como uma sequência de
instruções, com um número permanente e um nome. `catalog_schema.apply` cria
`schema_migrations(version, name, applied)`, aplica em ordem o que falta e
grava a linha de versão **por último**, no mesmo commit das instruções que
descreve.

Daí decorrem as três propriedades exigidas:

- **Ordenada.** As migrações são aplicadas na ordem dos seus números, que são
  permanentes e não reutilizados.
- **Idempotente.** Toda instrução é segura para rodar duas vezes
  (`IF NOT EXISTS`), e uma migração já registrada não é reaplicada.
- **Retomável.** Uma interrupção deixa as instruções que já rodaram no lugar e
  nenhuma linha de versão. A abertura seguinte reaplica a migração inteira
  desde o início, sem destruir nada — e é por isso que a idempotência das
  instruções é a condição, não um detalhe de implementação.

Nenhuma migração remove ou reescreve uma linha. A única remoção que a primeira
migração executa é a de um índice redundante (`artifact_entities_lookup`),
substituído por um que o planejador consegue usar por inteiro.

As migrações rodam na abertura do workspace, logo depois de
`artifact_registry.migrate_schema` e antes de qualquer indexação, para que o
preenchimento do histórico escreva no esquema corrente.

### Limites aceitos e registrados

| Limite | Consequência prática |
| --- | --- |
| Um escritor por vez | Duas execuções simultâneas serializam suas escritas de índice. O WAL mantém a leitura livre; a fila de computação é assunto de PHASE-04 |
| Um host | Não há execução distribuída. 01.4 de #6 já estabelece isso |
| Limite de variáveis por consulta (999) | Consultas por lista de ocorrências são fatiadas em blocos de 400 (`artifact_catalog._instances`) |
| Sem `ANALYZE` | O planejador escolhe sem estatísticas, então um índice só é usado com confiança quando cobre as colunas devolvidas. Foi o que exigiu o índice de entidades de ADR-001 |
| `PRAGMA journal_mode=WAL` pode ser recusado | Sistemas de arquivos sem memória compartilhada mantêm o modo anterior, e o workspace abre de todo modo — comportamento já existente, preservado |

## Justificativa

A tecnologia atual atende aos requisitos verificáveis desta fase: os planos de
consulta usam índices em todos os caminhos de acesso representativos (AC-01-03),
a transação é real e o projeto continua sendo uma pasta que o usuário copia.
Trocar o banco traria concorrência de escrita que a plataforma não usa, ao custo
de uma dependência que 01.4 proíbe. O que faltava não era o banco, era a
disciplina de esquema — e é isso que esta decisão acrescenta.

## Compromissos

- O catálogo não escala horizontalmente, e não se pretende que escale.
- Sem `ANALYZE`, a escolha de índice depende de cobertura. Índices cobertos
  custam escrita e espaço; o relatório da fase mede esse custo.
- `schema_migrations` passa a ser um ponto de verdade: uma migração publicada
  nunca é editada, apenas sucedida por outra.

## Implicações científicas

Nenhuma. Nenhum dado científico é lido, escrito ou convertido pela decisão.
A disponibilidade do catálogo afeta o que a interface consegue listar, não o que
os algoritmos calculam.

## Implicações de migração

Um banco de uma versão anterior é atualizado na primeira abertura, de forma
aditiva. Um banco já atualizado que seja aberto por uma versão anterior da
aplicação continua funcionando: as tabelas novas são ignoradas, e as colunas que
a versão anterior lê permanecem gravadas.

Projetos exportados (`.bme.zip`) não mudam de formato nesta fase. A
compatibilidade de cache é preservada: nenhuma identidade se deslocou.

## Requisitos de validação

| Critério | Verificação |
| --- | --- |
| REQ-CAT-002 | `test_artifact_catalog.CatalogMigrationTests` — versão registrada, reaplicação sem efeito, interrupção concluída na abertura seguinte, linhas e tabelas antigas preservadas |
| REQ-CAT-004, AC-01-03 | `test_representative_queries_use_their_indexes`, com os planos no relatório da fase |
| REQ-MIG-001 | `CatalogBackfillTests` — preenchimento idempotente, execução ilegível registrada e ignorada |

## Desvios da especificação de referência

Nenhum. A sub-issue #8 pede explicitamente a manutenção do SQLite com WAL,
conexões curtas e o fallback existente, que é o que esta decisão registra.
