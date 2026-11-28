# Modelo de identidade e catálogo versionado: PHASE-01 — 28/11/2026

PHASE-01 de #6 ([#8](https://github.com/mpiress/BioMolExplorer/issues/8)), a
partir da [avaliação da PHASE-00](../artifact_architecture_assessment_2026-11-26/README.md)
e com a suíte de testes já verde após a
[triagem](../test_suite_triage_2026-11-27/README.md) e seus três passos.
Branch `feat/phase-01-identity-catalog`, a partir de `origin/test`.

A fase formaliza as quatro identidades exigidas por #6 e dá ao catálogo
existente um esquema versionado e migrável, com índices verificados. As duas
decisões estão registradas em
[ADR-001](../../adr/0001-artifact-identity-model.md) e
[ADR-002](../../adr/0002-catalog-technology.md), nos dois idiomas.

| Medida | Antes | Agora |
| --- | --- | --- |
| `unittest discover -s tests` | OK, 1.670 testes | **OK**, 1.692 testes |
| Testes da fase | — | **22** (`tests/test_artifact_catalog.py`) |
| Falhas esperadas rastreadas | 3 (I1 ×2, I2) | **3** (I1 ×2, I2) |
| Versão do esquema do catálogo | — | **1**, registrada em `schema_migrations` |
| Identidades de cache deslocadas | — | **0 de 26** |
| Planos de consulta com varredura de tabela | — | **0 de 10** |

## O que faltava

A plataforma já distinguia três identidades, todas preservadas: o artefato
lógico (`artifact_id`, `occurrence_id`), a entidade científica
(`artifact_entities`) e o conteúdo, este apenas como uma coluna. Faltavam a
quarta identidade e a consultabilidade das relações.

| Lacuna | Evidência da PHASE-00 |
| --- | --- |
| Dependências gravadas como JSON em `artifact_instances.dependencies`; nenhum índice responde *quem consumiu esta ocorrência* | F-04 |
| O conteúdo não é um objeto: sem tamanho, formato nem local, e contável várias vezes | F-02, F-03, F-14 |
| A operação não existe no catálogo: a assinatura do cache vive na chave, não em linha consultável | 01.2, ADR-004 |
| O esquema não tinha versão: `CREATE TABLE IF NOT EXISTS` responde se a tabela existe, nunca qual versão ela é | REQ-CAT-002 |

## O que passou a existir

### Tabelas

| Tabela | Papel |
| --- | --- |
| `content_objects(sha256, size, media_type, storage_ref, created)` | Conteúdo físico, identificado pelo SHA-256 dos bytes. `storage_ref` fica vazio até a PHASE-02 |
| `operations(...)` | Uma execução de um bloco, com assinatura provisória, digests de parâmetros e de software, estado e tempos |
| `operation_outputs(operation_id, occurrence_id, port)` | Arestas de saída |
| `operation_inputs(operation_id, port, position, occurrence_id, asset_id, role)` | Arestas de entrada, que podem referenciar uma ocorrência, um arquivo do projeto ou um arquivo externo |
| `catalog_indexed_runs(project_id, run_id, version, updated)` | Marca de execução já derivada, por versão do catálogo |
| `catalog_backfill_issues(project_id, run_id, stage_id, reason, details, recorded)` | O que não pôde ser derivado, registrado em vez de adivinhado |
| `schema_migrations(version, name, applied)` | Versão do esquema |

Os desvios em relação aos nomes de referência da #8 — `signature_version`,
`position`/`asset_id`/`role`, as duas tabelas de controle e a ausência de chave
estrangeira sobre `content_sha256` — estão justificados um a um na seção de
desvios de [ADR-001](../../adr/0001-artifact-identity-model.md).

### Módulos

| Arquivo | Conteúdo |
| --- | --- |
| `src/biomolexplorer/catalog_schema.py` | Migrações ordenadas, idempotentes e não destrutivas; `apply`, `schema_version`, `pending` |
| `src/biomolexplorer/artifact_catalog.py` | Derivação das operações e dos objetos de conteúdo; consultas por entidade, conteúdo, produtor, montante e jusante |

A derivação não tem caminho próprio de escrita: `artifact_registry.register_run`,
que já é o único ponto por onde o índice de uma execução é gravado, passa a
derivar as operações da execução que registra. O histórico usa a mesma
derivação, por `artifact_catalog.backfill_operations`, lendo `runs.stages`.

### Uma invariante que faltava

Reindexar uma execução apaga as operações que ela produziu, antes de
reconstruí-las. A primeira versão desta fase apagava todas elas, e isso tinha um
efeito que só apareceu ao revisar o próprio código: ao limpar os resultados de
uma execução produtora, a linhagem de outra que a reaproveitou era perdida — ou
não, dependendo da ordem em que `clear_results` reindexa as execuções do
projeto, que a consulta `SELECT id,stages FROM runs WHERE project_id=?` não
especifica.

A regra passou a ser explícita: **uma operação sobrevive enquanto alguma
ocorrência indexada por outra execução ainda a declarar como produtora.** Limpar
os resultados de uma execução não desfaz o fato de ela ter sido executada, e o
resultado não depende mais de ordem alguma.

A guarda roda a cada mudança de resultado, então seu custo precisa seguir o que
a execução publicou, e não o tamanho do projeto. O planejador reordenava o join
e varria as instâncias do projeto; `CROSS JOIN` fixa a ordem, e `indexed_run_id`
entrou no índice `artifacts_occurrence` para que a pergunta *outra execução
ainda guarda esta ocorrência* seja respondida sem ler a linha. O plano está na
tabela abaixo e é verificado pelo mesmo teste dos demais.

| Teste | Verifica |
| --- | --- |
| `test_clearing_a_producer_keeps_the_lineage_of_a_run_that_reused_it` | A linhagem sobrevive à limpeza do produtor |
| `test_forgetting_a_run_drops_the_operations_nothing_claims` | A guarda não retém o que ninguém reivindica |

### O que não é adivinhado

| Motivo registrado | Situação |
| --- | --- |
| `missing_artifact_records` | A etapa publicou arquivos sem descritor |
| `missing_occurrence` | Descritor sem `occurrence_id` |
| `missing_software_identity` | A execução não registrou versão de implementação; a operação existe, sem assinatura |
| `unreadable_run` | A coluna `stages` não é JSON válido |

O terceiro caso é o das execuções anteriores ao digest de implementação. Um
digest da ausência seria uma assinatura incapaz de distinguir duas
implementações, então a assinatura fica `NULL` e a lacuna é registrada.

## Planos de consulta

Medidos com `EXPLAIN QUERY PLAN` sobre um banco real, criado pelo próprio
`WorkspaceStore`. Nenhum plano contém `SCAN`, e o teste
`test_representative_queries_use_their_indexes` verifica, para cada caminho de
acesso, a ausência de varredura **e** o índice esperado.

| Caminho de acesso | Plano |
| --- | --- |
| Artefatos de uma entidade científica | `SEARCH artifact_entities USING COVERING INDEX artifact_entities_entity (project_id=? AND entity_type=? AND entity_id=?)` |
| Ocorrência para artefato | `SEARCH artifact_instances USING INDEX artifacts_occurrence (project_id=? AND occurrence_id=?)` |
| Guarda de `forget_run` | `SEARCH o USING INDEX operations_project_run (project_id=? AND run_id=?)` + `SEARCH e ... (operation_id=?)` + `SEARCH i USING COVERING INDEX artifacts_occurrence (project_id=? AND occurrence_id=?)` |
| Artefatos de um conteúdo | `SEARCH artifact_instances USING INDEX artifacts_content (project_id=? AND content_sha256=?)` |
| Artefatos de um produtor | `SEARCH artifact_instances USING INDEX artifacts_project_stage (project_id=? AND stage_id=? AND indexed_run_id=?)` |
| Operações de uma execução | `SEARCH operations USING INDEX operations_project_run (project_id=? AND run_id=?)` |
| Operações por assinatura | `SEARCH operations USING INDEX operations_signature (project_id=? AND signature=?)` |
| Operação produtora de uma ocorrência | `SEARCH e USING INDEX operation_outputs_occurrence (occurrence_id=?)` + `SEARCH o USING INDEX sqlite_autoindex_operations_1 (operation_id=?)` |
| Montante | `SEARCH e USING INDEX operation_outputs_occurrence (occurrence_id=?)` + `SEARCH o ... (operation_id=?)` + `SEARCH i USING INDEX sqlite_autoindex_operation_inputs_1 (operation_id=?)` |
| Jusante | `SEARCH i USING INDEX operation_inputs_occurrence (occurrence_id=?)` + `SEARCH o ... (operation_id=?)` + `SEARCH e USING COVERING INDEX sqlite_autoindex_operation_outputs_1 (operation_id=?)` |

**Um índice precisou ser substituído.** A consulta por entidade usava
`artifact_entities_lookup(project_id, entity_type, entity_id)`, e o planejador a
ignorava: preferia a chave primária com apenas `project_id`, isto é, lia todas
as entidades do projeto. Sem `ANALYZE`, o SQLite escolhe sem estatísticas e
privilegia um índice que cobre as colunas devolvidas. Incluir `occurrence_id` e
`qualifiers` no índice (`artifact_entities_entity`) fez o planejador usar os
três termos de igualdade; o índice anterior tornou-se redundante e foi removido
pela migração. Como `workspace.py` recriava o índice antigo a cada abertura, a
sua declaração saiu de lá: os índices do catálogo passam a ser declarados apenas
em `catalog_schema.py`.

A independência do sistema de arquivos é verificada diretamente:
`test_indexed_queries_never_read_a_directory` substitui `Path.glob`,
`Path.rglob`, `Path.iterdir`, `os.walk`, `os.listdir` e `os.scandir` por funções
que falham, e executa as oito consultas.

## Custo das linhas novas

Medido com 2.000 poses publicadas por uma execução de docking, comparando o
tamanho do banco com e sem as tabelas desta fase (`VACUUM` nos dois casos):

| Medida | Valor |
| --- | --- |
| Banco completo | 9.277.440 B |
| Sem as tabelas da PHASE-01 | 7.938.048 B |
| Custo da fase | 1.339.392 B, ou **670 B por artefato** (~14% do catálogo) |
| Tempo de `register_run` para 2.000 artefatos | 0,19 s |

O custo é dominado pelos identificadores de 64 caracteres hexadecimais, que são
texto. Nenhum `stat` adicional foi introduzido: o tamanho do objeto de conteúdo
vem do mesmo `stat` que já decidia se o arquivo estava disponível.

**O custo não cresce com o reuso.** Uma execução que reaproveita 200 poses de
outra acrescenta às tabelas desta fase exatamente **uma** linha, a marca da
execução:

| Tabela | Linhas acrescentadas pela execução que reaproveita |
| --- | --- |
| `operations`, `operation_outputs`, `operation_inputs`, `content_objects` | 0 |
| `catalog_indexed_runs` | 1 |
| `artifact_instances` | 200 |

As 200 linhas de `artifact_instances` são a duplicação já medida em F-02 e F-14,
anterior a esta fase e endereçada pelas PHASE-02 e PHASE-05. As operações
continuam apontando para a execução que produziu os resultados, de modo que a
procedência sobrevive ao reuso (REQ-PROV-004).

## Identidades de cache: nenhuma deslocada

`artifact_registry.py` e `workspace.py` foram alterados, e os dois podem entrar
no fecho de importação de `stage_cache.operation_digest`. As 26 identidades que
uma edição de código-fonte poderia mover foram medidas antes e depois, com o
mesmo procedimento do
[passo 3 da triagem](../test_suite_port_contract_2026-11-28/README.md):
**nenhuma mudou**. Nenhum alias de compatibilidade foi necessário, e
preparações, redockings e docagens já calculados seguem reaproveitáveis.

## Migrações

| Propriedade | Como é garantida | Teste |
| --- | --- | --- |
| Ordenada | Números permanentes em `catalog_schema.MIGRATIONS`, aplicados em ordem | `test_the_schema_records_its_version_and_has_nothing_pending` |
| Idempotente | Toda instrução é `IF NOT EXISTS`; uma migração registrada não é reaplicada | `test_applying_the_migrations_again_changes_nothing`, `test_reopening_the_workspace_keeps_one_row_per_migration` |
| Retomável | A linha de versão é gravada por último, no mesmo commit das instruções | `test_an_interrupted_migration_completes_on_the_next_open` |
| Não destrutiva | Nenhuma migração remove ou reescreve uma linha; a única remoção é de um índice redundante | `test_legacy_rows_and_tables_survive_the_migration` |
| Preenchimento idempotente | Uma execução é processada uma vez por versão do catálogo | `test_the_backfill_derives_operations_and_runs_only_once` |

As migrações rodam na abertura do workspace, depois de
`artifact_registry.migrate_schema` e antes de qualquer indexação.

## Compatibilidade

- **Projetos anteriores.** O esquema é atualizado de forma aditiva na primeira
  abertura, e as operações do histórico são derivadas de `runs.stages`. Nenhuma
  linha é reescrita; `artifact_id`, `occurrence_id` e `artifact_type` já
  persistidos permanecem idênticos.
- **Versões anteriores da aplicação.** Um banco já atualizado continua legível:
  as tabelas novas são ignoradas e a coluna `dependencies` continua gravada.
- **Pacotes `.bme.zip`.** Formato inalterado nesta fase.

## Estado dos requisitos

| Identificador | Estado | Evidência |
| --- | --- | --- |
| AC-01-01 | `IMPLEMENTED` | `ContentIdentityTests` |
| AC-01-02 | `IMPLEMENTED` | `test_entities_resolve_to_artifacts_through_the_index` |
| AC-01-03 | `IMPLEMENTED` | `test_indexed_queries_never_read_a_directory`, `test_representative_queries_use_their_indexes` e os planos acima |
| AC-01-04 | `IMPLEMENTED` | ADR-001 e ADR-002 nos dois idiomas |
| REQ-ID-001, REQ-ID-002, REQ-ID-003 | `IMPLEMENTED` | `content_objects`; identidades lógicas inalteradas; `artifact_entities` consultável por índice |
| REQ-ID-004 | `IN_PROGRESS` | Identidade de operação registrada; a assinatura é provisória até ADR-004 (PHASE-03) |
| REQ-CAT-001 | `IMPLEMENTED` | Conteúdo, operações e arestas como linhas |
| REQ-CAT-002 | `IMPLEMENTED` | `CatalogMigrationTests` |
| REQ-CAT-003, REQ-CAT-004 | `IMPLEMENTED` | Planos acima, sem varredura |
| REQ-CAT-005 | `IN_PROGRESS` | Overhead medido e independente do reuso nas tabelas desta fase; a duplicação de `artifact_instances` segue para a PHASE-02 |
| REQ-DEP-001 | `IN_PROGRESS` | Arestas registradas para operações que publicam artefatos; execuções sem saída são da PHASE-04 |
| REQ-ARCH-004 | `IN_PROGRESS` | Dois dos doze ADRs aceitos |
| REQ-TEST-001 | `IMPLEMENTED` nesta fase | Suíte verde, sem falha nova e sem sucesso inesperado nos três itens rastreados |
| REQ-TEST-002 | `IMPLEMENTED` nesta fase | 22 testes novos |

## Limites desta fase

- A granularidade da operação é a etapa, não o lote. Decidido em ADR-004.
- Uma etapa que não publica arquivo não gera operação: a derivação parte dos
  artefatos publicados. Registrar execuções sem saída é tarefa da PHASE-04.
- A assinatura é provisória e não deve ser usada como chave de cache antes de
  ADR-004.
- `content_objects.storage_ref` fica vazio, e `size` fica `NULL` para artefatos
  preenchidos pelo histórico, que não são medidos no disco na abertura do
  workspace. A PHASE-02 preenche os dois ao trazer o conteúdo para o
  armazenamento gerenciado.
- As dependências ficam gravadas duas vezes, como JSON e como linhas, até que os
  leitores antigos sejam migrados (PHASE-05).
- Fora de #6, segue aberta a investigação do passo 1 da triagem: topologia ou
  tipos atômicos MOL2 respondem por 694 compostos excluídos no projeto
  Alzheimer.

## Reprodução

```sh
conda activate BioMolExplorer
export PYTHONPATH=src:tests PYTHONDONTWRITEBYTECODE=1
python -m unittest discover -s tests
python -m unittest test_artifact_catalog -v
```

O ambiente precisa estar ativado, e não apenas o seu interpretador invocado por
caminho: `test_molecular_core` e `test_ligand_preparation_performance` chamam
`obabel` por nome.

A medição das identidades de cache usou apenas leitura, por
`stage_cache.implementation_digest`, `stage_cache.operation_digest` e os digests
dos arquivos do adaptador de `ligand_cache.preparation_identity`.
