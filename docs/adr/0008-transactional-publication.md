# ADR-008 — Publicação transacional e recuperação

- **Status:** Aceita
- **Data:** 2026-11-28
- **Fase:** PHASE-02 (#9)
- **Requisitos:** REQ-STO-004, REQ-STO-005, REQ-INT-002, REQ-INT-003, REQ-TEST-003, AC-02-03
- **Autores e revisores:** Michel

## Contexto

Publicar um artefato toca três coisas que precisam concordar: os **bytes** no
armazenamento, as **linhas** do catálogo e o **estado da execução**. Cada uma
já tinha a sua garantia isolada:

| Parte | Garantia existente |
| --- | --- |
| Bytes | `storage.py` escreve em temporário com `fsync` e `rename`; `content_store.put` acrescenta a verificação do digest antes do rename |
| Linhas e estado da execução | Uma só transação SQLite: `artifact_registry.register_run` roda dentro da transação do chamador, que é a mesma que grava `runs.stages` |
| Movimentos de arquivo | `result_transactions.ResultTransaction` — journal em disco mais uma linha marcadora em `result_file_transactions`, e `recover` na abertura do workspace |

O que faltava era a relação entre elas. Uma interrupção pode cair de qualquer
lado do commit, e sem uma ordem declarada os dois estados ruins são possíveis:
uma linha que aponta para bytes ausentes, ou bytes que ninguém referencia.

Desses dois, apenas um é aceitável. Uma referência pendente quebra a leitura de
um resultado científico; bytes não referenciados são inertes.

## Alternativas consideradas

| Alternativa | Vantagens | Desvantagens | Evidência |
| --- | --- | --- | --- |
| A — Ordem declarada (bytes antes das linhas) mais journal do conteúdo introduzido *(selecionada)* | Referência pendente torna-se impossível por construção; a recuperação só precisa tratar o caso inerte | Uma interrupção pode deixar bytes a limpar | Testes de injeção de falha do relatório |
| B — Commit em duas fases próprio, com log de intenções sobre as três partes | Simetria completa | Reimplementa o que a transação do SQLite já dá para duas das três partes, e acrescenta um log a recuperar | — |
| C — Linhas antes dos bytes | A falha deixa o armazenamento limpo | Deixa exatamente o estado inaceitável: o catálogo prometendo conteúdo ausente | 01.3 |
| D — Nada; confiar na reexecução | Simples | Uma reexecução científica custa horas, e REQ-INT-003 exige retomar de saídas validadas | A.3 |

## Solução selecionada

### A invariante

> **Os bytes chegam ao armazenamento antes de qualquer linha os referenciar.**

Daí decorre que nenhuma linha do catálogo pode apontar para conteúdo ausente.
O único resíduo possível de uma interrupção é conteúdo que ninguém referencia.

### O journal

`ResultTransaction.record_objects` nomeia o conteúdo que a transação vai
introduzir, **antes** de ele ser escrito, no mesmo journal que já registrava os
movimentos de arquivo. O journal passa a ter uma chave `objects`.

A linha marcadora em `result_file_transactions` continua sendo o que distingue
uma transação concluída de uma interrompida: ela é inserida na transação do
chamador, então ela existe se e somente se aquela transação commitou.

`ArtifactManager.publish` é a sequência completa: journal, bytes, linhas, commit
— e as linhas commitam com o estado da execução, porque são a mesma transação.

### A recuperação

`result_transactions.recover`, já chamado na abertura do workspace:

| Interrupção | Estado encontrado | O que a recuperação faz |
| --- | --- | --- |
| Antes do commit | Journal presente, marcador ausente | Restaura os movimentos e descarta o conteúdo nomeado no journal **que nada referencia** |
| Depois do commit, antes da limpeza | Journal e marcador presentes | Conclui a limpeza e **mantém** o conteúdo: as linhas o referenciam |
| Durante a transação do catálogo | O SQLite desfaz as linhas | Igual ao primeiro caso; os bytes verificados permanecem no armazenamento |
| Depois de tudo | Nem journal nem marcador | Nada a fazer |

### A recuperação é conservadora

Descarta-se um objeto apenas quando **nenhum** artefato lógico do projeto
declara aquele conteúdo e **nenhum** objeto de conteúdo registrado aponta para
ele. REQ-STO-005 proíbe apagar conteúdo científico automaticamente, e a
consistência do catálogo nunca é motivo suficiente para apagar bytes.

A consequência honesta: uma publicação interrompida cujo conteúdo pertence a um
artefato ainda publicado deixa bytes gerenciados no armazenamento, sem a linha
que os registra. Eles são conteúdo correto, a publicação seguinte deduplica
sobre eles, e a varredura de órfãos em escala é da PHASE-07 (REQ-INT-001 e
REQ-INT-002). `ContentStore.sweep` já existe para quando ela chegar.

## Justificativa

A assimetria entre os dois estados ruins é o que decide. Como só um deles é
tolerável, basta ordenar as escritas para tornar o outro impossível — e isso é
mais forte do que qualquer protocolo que tente consertar os dois depois. O
journal não existe para garantir consistência; existe para que a limpeza do
caso tolerável seja possível em vez de ficar para sempre.

A alternativa B reimplementaria, com um log novo a recuperar, o que a transação
do SQLite já entrega para as linhas e o estado da execução.

## Compromissos

- Uma interrupção pode deixar bytes a limpar, e a limpeza é conservadora por
  decisão.
- O journal é por projeto e por transação; ele não substitui uma verificação de
  integridade do armazenamento inteiro, que é da PHASE-07.
- `publish` abre a sua própria transação. Uma publicação feita dentro de uma
  transação já aberta pelo chamador passa `db` e participa dela, mas então o
  journal é responsabilidade do chamador.

## Implicações científicas

Nenhuma. Nenhum resultado é calculado, convertido ou descartado pela decisão.
O que ela garante é que um resultado já calculado não se torne ilegível por uma
interrupção, o que é um requisito de rastreabilidade (01.3), não de cálculo.

## Implicações de migração

A chave `objects` é aditiva no journal:

- um journal escrito por esta versão e lido por uma anterior é lido
  normalmente, porque a versão anterior lê `moves` e ignora chaves que não
  conhece;
- um journal escrito por uma versão anterior e lido por esta não tem `objects`,
  e a recuperação trata a lista como vazia.

Nenhuma tabela nova, nenhuma linha reescrita.

## Requisitos de validação

| Critério | Verificação |
| --- | --- |
| Interrupção antes do commit | `test_content_store.TransactionalPublicationTests.test_an_interruption_before_the_commit_leaves_no_dangling_reference` |
| Limpeza do conteúdo não reivindicado | `test_an_interruption_before_the_commit_discards_unclaimed_content` |
| Interrupção depois do commit | `test_an_interruption_after_the_commit_keeps_the_published_content` |
| Transação do catálogo desfeita | `test_a_rolled_back_catalog_transaction_leaves_the_store_consistent` |
| Publicação completa não deixa resíduo | `test_a_complete_publication_leaves_no_journal_and_no_marker` |
| REQ-STO-004 | `test_publication_leaves_no_partial_object_when_content_changes` |

## Desvios da especificação de referência

Nenhum. A #9 pede que a escrita do conteúdo, as linhas do catálogo e o commit
do estado da execução se recuperem de forma consistente nos dois lados do
commit, e é isso que a invariante mais o journal entregam. A varredura de
órfãos em escala permanece onde a própria #6 a coloca, na PHASE-07.
