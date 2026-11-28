# ADR-003 — Estratégia de armazenamento endereçado por conteúdo

- **Status:** Aceita
- **Data:** 2026-11-28
- **Fase:** PHASE-02 (#9)
- **Requisitos:** REQ-STO-001, REQ-STO-002, REQ-STO-003, REQ-STO-005, REQ-ARCH-002, REQ-ARCH-003, REQ-CAT-005, AC-02-01, AC-02-02, AC-02-04
- **Autores e revisores:** Michel

## Contexto

A plataforma já armazenava conteúdo pelo seu digest, em dois lugares:

| Mecanismo | O que faz | O que falta |
| --- | --- | --- |
| `calculation_cache` | Tabela e diretório `objects`, `copy_artifact` com copy-on-write, publicação transacional depois da validação, travas por chave | Objetos graváveis; verificação dos bytes escritos só pelo chamador; nenhuma política de entrega ao consumidor |
| `preparation_cache` | Entradas imutáveis publicadas por `rename` atômico após validar cada arquivo, com todo digest verificado na leitura | Layout próprio, por sufixo e não por conteúdo; migra na PHASE-03 |

Peças relacionadas: `storage.py` (escrita atômica e `fsync` de diretório),
`result_transactions.py` (movimentos de arquivo recuperados em torno de um
commit) e os nomes públicos de `artifact_naming.py`.

Três fatos medidos nesta máquina decidem a política, e estão no
[relatório da fase](../validation/phase_02_content_store_2026-11-28/README.md):

1. **O sistema de arquivos é ext4 e não suporta reflink.** `FICLONE` falha com
   `EOPNOTSUPP`, então toda materialização por cópia é uma cópia integral. A
   política "reflink, depois cópia" do Apêndice A degrada para cópia aqui.
2. **Hard links funcionam.** Um link para um objeto custa zero bytes.
3. **Um objeto em modo 0444 recusa escrita.** O núcleo nega o `open` para
   escrita em vez de deixar corromper conteúdo que outras execuções referenciam.

Juntos, (2) e (3) dão o que (1) nega: uma entrega sem custo para um consumidor
somente-leitura, que se protege sozinha.

## Alternativas consideradas

| Alternativa | Vantagens | Desvantagens | Evidência |
| --- | --- | --- | --- |
| A — Um armazenamento, no layout que `calculation_cache` já usa, com objetos somente-leitura *(selecionada)* | Caches existentes continuam válidos sem serem movidos; uma só implementação de publicação, verificação e entrega | Mantém dois índices (o do projeto e o do catálogo) sobre um só armazenamento | Medições do relatório |
| B — Um terceiro armazenamento, novo, e migração dos dois atuais | Modelo limpo | Contraria REQ-ARCH-002 (generalizar, não duplicar) e invalidaria caches já calculados do usuário | A.4, F-08 |
| C — Objetos graváveis, como hoje | Nenhuma mudança de permissão | Um motor que reescreve a própria entrada corrompe conteúdo compartilhado; é exatamente o que REQ-STO-003 proíbe | — |
| D — Compressão dos objetos | Menos bytes | PDB, MOL2 e PDBQT são lidos por ferramentas nativas que esperam o arquivo; descomprimir a cada materialização troca espaço por tempo em cada reuso, e o ganho não foi medido como necessário | 01.2 |
| E — Hard link para todo consumidor | Custo zero sempre | Um motor que escreve em lugar corrompe o objeto; só é seguro para consumidores comprovadamente somente-leitura | A.2 |

## Solução selecionada

### Escopo e raiz

O armazenamento é **do projeto**: `<projeto>/.biomolexplorer/store/objects`.
Um projeto continua sendo uma pasta que o usuário copia, que é o que sustenta a
portabilidade verificada na PHASE-00 e a exportação `.bme.zip`. O
compartilhamento entre projetos permanece como já era, por
`BIOMOL_PREPARATION_CACHE` e `BIOMOL_CACHE_DIR`, que definem a raiz do
`calculation_cache`; 01.4 exclui cache *distribuído*, não um diretório local
compartilhado que o usuário escolhe.

### Layout

`objects/<dois primeiros caracteres do digest>/<digest>` — exatamente o layout
que `calculation_cache` já usava. É por isso que um cache existente continua
funcionando sem ser movido: a implementação mudou, o endereço não.

### Permissões

Um objeto é gravado em modo **0444**. Substituí-lo é governado pelo diretório,
não pelo arquivo, então a publicação atômica continua possível; escrever nele,
não.

### Publicação

`ContentStore.put` escreve em um temporário no diretório de destino, faz
`fsync`, **verifica o digest dos bytes escritos**, aplica o modo somente-leitura
e só então renomeia, sincronizando o diretório. Uma interrupção em qualquer
ponto deixa ou nada, ou um objeto completo — nunca um parcial. Uma publicação
concorrente do mesmo conteúdo renomeia um objeto idêntico sobre outro idêntico.

`validate` recebe os bytes que estão sendo publicados e recusa-os levantando
exceção: **conteúdo que falha no contrato do seu tipo nunca é publicado**.

### Verificação

A verificação completa acontece **na leitura e antes de todo reuso**
(REQ-STO-002): `open`, `materialize` e `ArtifactManager.resolve` recomputam o
digest e recusam uma divergência em vez de tratá-la como acerto. Na publicação,
um objeto já presente é confiado — o `calculation_cache` passa a sua própria
função de verificação, com o cache de estampa do sistema de arquivos que ele já
mantinha, de modo que o comportamento dele não muda.

### Materialização

| Modo | Como | Para quem |
| --- | --- | --- |
| `shared` | Hard link para o objeto somente-leitura | Consumidor que apenas lê. Custa zero bytes e o núcleo recusa quem escrever em lugar |
| `isolated` | Cópia independente e gravável (copy-on-write quando houver), verificada depois | Motor que reescreve a própria entrada (REQ-STO-003) |

Quando o link não é possível (outro dispositivo, sistema de arquivos sem
links), o modo `shared` faz uma cópia e **mantém-na somente-leitura**: o
contrato do modo vale, o custo deixa de valer.

### Compressão

Nenhuma, por padrão. Os formatos científicos são lidos por ferramentas nativas
que esperam o arquivo; comprimir trocaria espaço por descompressão em cada
materialização. A decisão é revisível se uma medição mostrar necessidade.

### Por que `calculation_cache` é generalizado, e não substituído

Porque os seus checkpoints são resultados científicos já calculados do usuário.
Mover o armazenamento invalidaria todos. A mudança é de implementação:
`copy_artifact` delega ao armazenamento, `object_path` devolve o mesmo caminho,
as tabelas e os leitores seguem idênticos. O custo dessa decisão aparece na
seção de identidades de cache do relatório: `calculation_cache.py` é lido
diretamente por `stage_cache.implementation_digest`, de modo que 16 das 26
identidades se deslocaram e foram registradas como equivalência auditada.

## Justificativa

A alternativa B seria mais limpa no papel e destrutiva na prática. A C mantém o
risco que REQ-STO-003 existe para eliminar. A E é a tentação óbvia num sistema
sem reflink, e é segura apenas para o consumidor que não escreve — por isso o
modo existe, explícito, em vez de ser o padrão.

O que a decisão acrescenta ao que já havia não é um armazenamento, é o
fechamento de três lacunas: os bytes escritos passam a ser verificados *pelo
armazenamento*, o objeto passa a recusar escrita, e a entrega ao consumidor
passa a ser uma escolha declarada entre compartilhar e isolar.

## Compromissos

- **Dois índices, um armazenamento.** O índice do projeto
  (`calculation_cache.objects`) e o catálogo do workspace (`content_objects`)
  descrevem conteúdo do mesmo armazenamento. Unificar os índices é da PHASE-05.
- **Resultados publicados continuam nas pastas das execuções.** Nesta fase o
  armazenamento gerenciado contém os checkpoints de cálculo e o que o
  gerenciador registra explicitamente. Trazer todo resultado publicado para o
  armazenamento duplicaria bytes num sistema sem reflink; fazê-lo por link
  exige provar que todo consumidor é somente-leitura, o que a própria #9
  condiciona. A decisão fica para ADR-010 (PHASE-07), com a retenção.
- **O contrato de conteúdo é tão forte quanto a declaração do tipo.** `scores`
  não declara esquema porque as suas colunas variam por motor, então o seu
  contrato verifica presença e cabeçalho, não colunas nomeadas.
- **A verificação na leitura custa uma leitura.** Um reuso de 15 mil poses
  relê os bytes para verificá-los. É o preço de nunca reusar conteúdo
  divergente, e é o que `calculation_cache` já pagava.
- **O catálogo não é reescrito.** Uma ocorrência indexada por uma versão
  anterior, sem digest gravado, tem o digest **computado** na resolução, não
  escrito de volta.

## Implicações científicas

Nenhum algoritmo, protocolo, tolerância ou byte científico muda. Em particular:

- Os digests, os caminhos dos objetos, as chaves e as tabelas do
  `calculation_cache` são os mesmos, e as 26 identidades de cache voltaram ao
  valor anterior após o registro da equivalência auditada: preparações,
  redockings e docagens já calculados continuam reaproveitáveis.
- A entrega passa a ser mais segura, não diferente: um motor que antes poderia
  reescrever um objeto compartilhado agora recebe uma cópia isolada ou é
  recusado pelo núcleo.
- Conteúdo que falha no contrato do seu tipo deixa de poder ser publicado ou
  reusado. É uma recusa nova, e nunca um resultado diferente.

## Implicações de migração

Aditiva. Caches existentes continuam no mesmo lugar, com os mesmos endereços, e
são lidos pela nova implementação. Objetos já gravados permanecem graváveis —
só os novos nascem somente-leitura; nada é reescrito para mudar o modo dos
antigos.

Uma versão anterior da aplicação que abra o mesmo projeto lê os objetos
normalmente, porque o layout não mudou. Pacotes `.bme.zip` não mudam de formato
nesta fase; a inclusão do armazenamento no pacote é da PHASE-07.

## Requisitos de validação

| Critério | Verificação |
| --- | --- |
| AC-02-01 | `test_content_store.StoreTests.test_identical_content_from_two_producers_is_stored_once`, `ArtifactManagerTests.test_two_runs_publishing_the_same_bytes_share_one_managed_object` |
| AC-02-02 | `test_a_writable_materialization_cannot_reach_the_stored_object`, `test_materialization_is_shared_by_default_and_isolated_on_request` |
| AC-02-04 | `test_usage_reports_the_store_and_what_the_catalog_claims` e a medição `R_storage` do relatório |
| REQ-STO-002 | `test_a_tampered_object_is_refused_on_read_and_on_materialization`, `test_resolve_returns_the_managed_object_and_detects_tampering` |
| REQ-STO-004 | `test_publication_leaves_no_partial_object_when_content_changes`, `test_concurrent_publication_of_the_same_content_is_safe` |
| Contratos de conteúdo | `ValidationContractTests` |

## Desvios da especificação de referência

| Desvio | Razão |
| --- | --- |
| A política de materialização não é "reflink, depois cópia", e sim "link para quem só lê, cópia isolada para quem escreve" | Medido: não há reflink neste sistema de arquivos. A tentativa de reflink permanece na cópia, para onde houver |
| O armazenamento não recebe os resultados publicados nesta fase | Sem reflink, duplicaria bytes; por link, exigiria provar todo consumidor somente-leitura, o que a própria #9 condiciona |
| `preparation_cache` não é roteado aqui | A própria #9 o adia para a PHASE-03 |
