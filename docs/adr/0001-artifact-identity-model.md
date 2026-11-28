# ADR-001 — Modelo de identidade dos artefatos científicos

- **Status:** Aceita
- **Data:** 2026-11-28
- **Fase:** PHASE-01 (#8)
- **Requisitos:** REQ-ID-001, REQ-ID-002, REQ-ID-003, REQ-ID-004, REQ-CAT-001, REQ-DEP-001, REQ-PROV-004, AC-01-01, AC-01-02
- **Autores e revisores:** Michel

## Contexto

A plataforma já distingue identidades de artefatos, mas apenas três delas são
de primeira classe e nenhuma operação tem identidade própria.

| Identidade | Onde vive hoje | O que responde |
| --- | --- | --- |
| Artefato lógico | `artifacts.slot_id` e `artifacts.records` gravam `artifact_id` e `occurrence_id` em `artifact_instances` | *Qual saída de qual bloco é esta, e de qual execução* |
| Entidade científica | `artifact_registry.index_entities` e a tabela `artifact_entities` | *Qual composto, estrutura, receptor ou pose está neste arquivo* |
| Conteúdo | coluna `artifact_instances.content_sha256`, preenchida pelo `artifact_manifest` da etapa | *Quais bytes foram publicados* |
| Operação | — | *Qual execução de qual bloco produziu isto, a partir de quê* |

O `artifact_id` é o par produtor + saída qualificada (`slot_id(stage_id,
locator)`); a `occurrence_id` acrescenta a execução que produziu o arquivo. Nem
um nem outro depende do caminho físico, e a
[avaliação da PHASE-00](../validation/artifact_architecture_assessment_2026-11-26/README.md)
confirmou que a relocação de um projeto preserva os dois. Esse não é o problema.

O que falta está medido no mesmo relatório:

- **F-04 — dependências como JSON.** `artifact_instances.dependencies` guarda
  uma lista de `artifact_id` serializada e `provenance` guarda um objeto. Não há
  índice que responda *quem consumiu esta ocorrência*: a pergunta exige ler e
  desserializar todas as linhas do projeto.
- **O conteúdo não é um objeto.** Sendo uma coluna, o SHA-256 não tem tamanho,
  formato nem local de armazenamento associados, e nada impede que o mesmo
  conteúdo seja contado como vários artefatos físicos. F-02, F-03 e F-14
  mediram justamente a cópia repetida do mesmo conteúdo entre execuções.
- **A operação não existe no catálogo.** O cache de etapas
  (`stage_cache.stage_key`, `operation_digest`) já calcula uma assinatura por
  execução, mas ela vive na chave do cache, não em uma linha consultável. Não é
  possível perguntar *quais operações compartilham esta assinatura* nem *o que
  esta operação consumiu*, que são as perguntas das PHASE-03 e PHASE-04.

## Alternativas consideradas

| Alternativa | Vantagens | Desvantagens | Evidência |
| --- | --- | --- | --- |
| A — Identidade única por conteúdo: o SHA-256 passa a ser *a* identidade do artefato | Deduplicação trivial; um só identificador | Dois experimentos distintos que produzem bytes idênticos deixam de ser distinguíveis, violando REQ-PROV-004; a procedência de uma pose deixaria de dizer de qual redocking ela veio | `artifact_instances` já registra ocorrências distintas com o mesmo `content_sha256` |
| B — Manter as dependências em JSON e consultar com `json_extract` | Nenhuma migração de esquema | Nenhum índice é utilizável sobre o conteúdo do JSON; o custo da travessia cresce com o número de artefatos do projeto, que é exatamente o que F-04 aponta | Plano de consulta sobre `json_each` percorre a tabela inteira |
| C — Quatro identidades separadas, com arestas como linhas *(selecionada)* | Cada pergunta tem um índice; reuso é referência, não cópia; a procedência sobrevive à deduplicação | Mais tabelas e uma derivação a manter | Planos de consulta em `docs/validation/phase_01_identity_catalog_2026-11-28/README.md` |
| D — Modelo de procedência completo (W3C PROV / RDF) | Vocabulário padronizado e interoperável | Exige dependência externa e um armazenamento adicional, contra REQ-ARCH-003; o vocabulário não cobre as entidades científicas do domínio sem extensão própria | 01.4 de #6 |

## Solução selecionada

Quatro identidades, com papéis disjuntos e nenhuma delas derivada de um caminho
de arquivo.

### 1. Conteúdo

`content_objects(sha256, size, media_type, storage_ref, created)`. A identidade
é o SHA-256 dos bytes (REQ-ID-001). `media_type` é o formato dos bytes
(`chemical/x-pdb`, `text/csv`), que não se confunde com o `artifact_type`
científico da ocorrência. `storage_ref` permanece vazio nesta fase: o
armazenamento gerenciado é decidido em ADR-003 (PHASE-02), e um caminho é
localização, nunca identidade.

Vários artefatos lógicos referenciam um objeto de conteúdo. É o que a consulta
`artifact_catalog.artifacts_by_content` retorna, e o que AC-01-01 verifica.

### 2. Artefato lógico

`artifact_id` e `occurrence_id` permanecem **inalterados**, com o mesmo cálculo
e as mesmas colunas. Nenhuma identidade já persistida é recalculada: esta fase
não reescreve nenhuma linha de `artifact_instances`.

### 3. Entidade científica

`artifact_entities` permanece como está, com `entity_type` ∈ {`compound`,
`structure`, `receptor`, `pose`, `community`, `molecular_pair`} e os
qualificadores (`chain_id`, `assembly_id`, `ligand_id`) em JSON. A consulta
`artifacts_by_entity` resolve entidade → ocorrências → artefatos (AC-01-02).

O índice da tabela passa a incluir as colunas devolvidas pela consulta
(`artifact_entities_entity`), porque sem elas o planejador preferia a chave
primária usando apenas `project_id` — isto é, lia todas as entidades do
projeto. O índice anterior, `artifact_entities_lookup`, torna-se redundante e é
removido pela migração; remover um índice não remove nenhuma linha.

### 4. Operação

`operations` registra **uma execução de um bloco**:

```
operations(operation_id, project_id, operation_type, stage_id, block_id, run_id,
           signature, signature_version, parameters_digest, software_digest,
           state, started, finished, updated)
```

- `operation_id = SHA-256(['operation', project_id, stage_id, run_id])`, onde
  `run_id` é a execução que **produziu** os artefatos. É determinístico, de modo
  que a derivação pode ser repetida sem criar linhas novas.
- `run_id` não é necessariamente a execução que indexou os artefatos. Quando uma
  execução reaproveita um resultado em cache, suas ocorrências continuam
  apontando para a operação que as produziu. É assim que o reuso aparece como
  referência, e não como cópia (REQ-CAT-005).
- `parameters_digest` cobre os parâmetros efetivos da etapa;
  `software_digest` cobre `processing_version` (o digest de implementação),
  `producer_software` (versões das ferramentas) e `runtime_id`.
- `signature` combina tipo de operação, parâmetros, software e os digests de
  conteúdo das entradas, prefixada por `signature_version`. **É provisória:**
  ADR-004 (PHASE-03) fixa o modelo definitivo, e a coluna de versão existe
  para que uma assinatura antiga nunca seja confundida com uma nova.
- `software_digest` é `NULL`, e com ele a assinatura, quando a execução não
  registrou identidade de software — o caso de toda execução anterior ao digest
  de implementação. Um digest da ausência seria uma assinatura incapaz de
  distinguir duas implementações, então a lacuna é registrada em
  `catalog_backfill_issues` em vez de preenchida.

### Arestas de dependência

```
operation_outputs(operation_id, occurrence_id, port)
operation_inputs(operation_id, port, position, occurrence_id, asset_id, role)
```

Uma entrada é a ocorrência de outro artefato, um arquivo do projeto ou um
arquivo externo, então `role` ∈ {`artifact`, `asset`, `external`} diz qual
referência foi registrada e as duas colunas de referência são anuláveis.
`position` ordena as entradas de uma porta e mantém a chave primária estável
quando a derivação é repetida.

A travessia da linhagem é feita sobre essas linhas:
`upstream` vai de uma ocorrência à operação que a publicou e daí às entradas
dela; `downstream` faz o inverso. A coluna JSON `dependencies` continua gravada,
porque é o que a derivação lê e o que versões anteriores da aplicação leem.

### Onde a derivação acontece

`artifact_registry.register_run` é o único caminho de escrita do índice, e
passa a derivar as operações da execução que registra. A mesma derivação atende
o histórico: `artifact_catalog.backfill_operations` lê `runs.stages`, que já
guarda os descritores dos artefatos. Uma execução é processada uma vez por
versão do catálogo, registrada em `catalog_indexed_runs` — a presença de linhas
em `operations` não serve como marca, justamente porque o reuso dá a uma
execução operações que pertencem a outras.

Reindexar uma execução apaga as operações que ela produziu antes de
reconstruí-las, sob uma invariante explícita: **uma operação sobrevive enquanto
alguma ocorrência indexada por outra execução ainda a declarar como produtora.**
Limpar os resultados de uma execução não desfaz o fato de ela ter sido
executada, e a linhagem de uma execução que reaproveitou aqueles resultados não
pode depender da ordem em que as execuções do projeto são reindexadas.

O que a derivação não consegue estabelecer é registrado, nunca adivinhado,
seguindo a política de `artifact_migration_issues`:

| Motivo | Situação |
| --- | --- |
| `missing_artifact_records` | A etapa publicou arquivos sem descritor; sem ocorrência não há aresta |
| `missing_occurrence` | Um descritor sem `occurrence_id` |
| `missing_software_identity` | A execução não registrou versão de implementação; a operação existe, sem assinatura |
| `unreadable_run` | A coluna `stages` da execução não é JSON válido |

## Justificativa

As quatro perguntas são diferentes e cada uma precisa de um índice próprio. A
alternativa A responde à pergunta do conteúdo e perde a da procedência: duas
docagens com parâmetros distintos que convergem para a mesma pose passariam a
ser o mesmo artefato, e a plataforma deixaria de poder dizer de qual protocolo
cada resultado veio — um requisito científico, não uma conveniência. A
alternativa B mantém a estrutura atual e com ela o custo que F-04 mediu. A
alternativa D traria um vocabulário que ainda precisaria ser estendido para
compostos, cadeias e pares cristalográficos, ao custo de infraestrutura que
01.4 proíbe.

Separar a operação, e não apenas as arestas, é o que permite às fases seguintes
perguntar *esta operação já foi executada com esta assinatura* (PHASE-03) e
*quais operações dependem desta ocorrência* (PHASE-04) sem reabrir o esquema.

## Compromissos

- **Granularidade da operação é a etapa, não o lote.** Uma etapa em modo
  individual processa vários lotes com versões e tempos próprios; o catálogo
  registra uma operação por etapa e execução, e consolida os digests de software
  dos lotes em um só. A granularidade por lote é decidida em ADR-004, com o
  modelo de assinatura.
- **Uma etapa que não publica arquivo não gera operação.** A derivação parte dos
  artefatos publicados. Registrar execuções sem saída é tarefa da PHASE-04, que
  observa a execução em vez de reconstruí-la.
- **A assinatura desta fase é provisória** e não deve ser usada como chave de
  cache antes de ADR-004.
- **Duplicação temporária.** As dependências ficam gravadas duas vezes, como
  JSON e como linhas, até que os leitores antigos sejam migrados (PHASE-05).
- **Um `stat` por artefato publicado** passa a alimentar o tamanho do objeto de
  conteúdo. Não é um custo novo: é o mesmo `stat` que já decidia se o arquivo
  estava disponível.

## Implicações científicas

Nenhuma. Nenhum algoritmo, protocolo, tolerância ou arredondamento é tocado.
Em particular:

- `artifact_id`, `occurrence_id` e `artifact_type` das linhas já persistidas
  permanecem idênticos — inclusive o `molecular_poses` das poses de redocking,
  preservado deliberadamente pela decisão C5 da
  [triagem da suíte](../validation/test_suite_port_contract_2026-11-28/README.md).
- Nenhuma identidade de cache se desloca. As 26 identidades que um edição de
  código-fonte poderia mover foram medidas antes e depois: nenhuma mudou, de
  modo que preparações, redockings e docagens já calculados seguem
  reaproveitáveis.
- A rastreabilidade aumenta: a operação produtora de uma ocorrência passa a ser
  uma linha consultável, e não uma inferência sobre nomes de arquivo.

## Implicações de migração

A migração é aditiva: cria tabelas e índices, não altera nem remove nenhuma
linha. Projetos de versões anteriores continuam legíveis, e o preenchimento do
histórico é idempotente e retomável — uma passagem interrompida deixa marcadas
as execuções que concluiu.

Uma versão anterior da aplicação que abra o mesmo banco ignora as tabelas novas
e continua lendo a coluna `dependencies`, que permanece gravada. Pacotes
`.bme.zip` não mudam de formato nesta fase.

## Requisitos de validação

| Critério | Verificação |
| --- | --- |
| AC-01-01 | `test_artifact_catalog.ContentIdentityTests` — dois artefatos lógicos com os mesmos bytes compartilham um objeto de conteúdo, com tamanho e tipo de mídia |
| AC-01-02 | `test_artifact_catalog.CatalogQueryTests.test_entities_resolve_to_artifacts_through_the_index` |
| REQ-DEP-001 | `OperationIdentityTests.test_each_block_execution_becomes_one_operation_with_its_edges` |
| REQ-PROV-004 | `test_a_reused_result_points_to_the_operation_that_produced_it` |
| Política de não adivinhar | `test_a_published_file_without_a_descriptor_is_recorded_not_guessed`, `test_a_run_without_software_identity_has_no_signature` |

Relatório completo, com os planos de consulta e a medição das identidades de
cache: `docs/validation/phase_01_identity_catalog_2026-11-28/README.md`.

## Desvios da especificação de referência

O Apêndice A.2 e a sub-issue #8 dão nomes de referência para as tabelas. Os
desvios, todos documentados aqui, são:

| Desvio | Razão |
| --- | --- |
| `operations` ganha `signature_version` e `updated` | Sem a versão, uma assinatura provisória desta fase seria indistinguível da definitiva de ADR-004 |
| `operation_inputs` ganha `position`, `asset_id` e `role` em lugar de `(operation_id, occurrence_id, port, role)` | Uma entrada pode não ser um artefato; `position` dá uma chave primária estável a entradas repetidas na mesma porta |
| Tabelas adicionais `catalog_indexed_runs` e `catalog_backfill_issues` | Idempotência do preenchimento e registro do que não pode ser derivado |
| `artifact_instances.content_sha256` referencia `content_objects` sem chave estrangeira declarada | Linhas antigas podem ter conteúdo não registrado; uma restrição as tornaria ilegíveis, contra REQ-MIG-001 |
