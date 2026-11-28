# Armazenamento endereçado por conteúdo e gerenciador de artefatos: PHASE-02 — 28/11/2026

PHASE-02 de #6 ([#9](https://github.com/mpiress/BioMolExplorer/issues/9)), sobre
o catálogo da [PHASE-01](../phase_01_identity_catalog_2026-11-28/README.md).
Branch `feat/phase-02-content-store`, a partir de `origin/test`.

A fase dá à plataforma um só armazenamento imutável endereçado por conteúdo e um
gerenciador que registra, deduplica, valida, resolve e materializa artefatos,
generalizando os mecanismos existentes em vez de acrescentar um terceiro. As
decisões estão em [ADR-003](../../adr/0003-content-addressable-storage.md) e
[ADR-008](../../adr/0008-transactional-publication.md), nos dois idiomas.

| Medida | Antes | Agora |
| --- | --- | --- |
| `unittest discover -s tests` | OK, 1.692 testes | **OK**, 1.720 testes |
| Testes da fase | — | **28** (`tests/test_content_store.py`) |
| Falhas esperadas rastreadas | 3 (I1 ×2, I2) | **3** (I1 ×2, I2) |
| Identidades de cache deslocadas | — | **0 de 26**, após equivalência auditada |
| Aliases de compatibilidade | 68 | **79** |

## Três fatos medidos antes de decidir

A política de materialização do Apêndice A — "reflink, depois cópia" — pressupõe
copy-on-write. Nesta máquina ele não existe:

| Fato | Medição |
| --- | --- |
| Sistema de arquivos | ext4 (`/dev/nvme0n1p2`) |
| `FICLONE` (reflink) | **Não suportado** — `Operation not supported` |
| Hard link | Suportado |
| Escrita em objeto modo 0444 | **Recusada** pelo núcleo |

Os dois últimos dão o que o segundo nega: uma entrega sem custo para um
consumidor somente-leitura, que se protege sozinha. É por isso que a política
desta fase é **link para quem só lê, cópia isolada para quem escreve**, e não
"reflink, depois cópia" — o desvio está registrado em ADR-003.

## O que passou a existir

| Arquivo | Conteúdo |
| --- | --- |
| `src/biomolexplorer/content_store.py` | `ContentStore`: `put`, `has`, `open`, `verify`, `materialize`, `size`, `discard`, `stored`, `sweep`, `usage`. Publica em temporário, verifica o digest dos bytes escritos, aplica modo 0444 e renomeia atomicamente |
| `src/biomolexplorer/artifact_manager.py` | `ArtifactManager`: `register`, `register_many`, `publish`, `resolve`, `resolve_entity`, `materialize`, `provenance`, `usage` |
| `src/biomolexplorer/artifact_validation.py` | Contratos de conteúdo por tipo; `InvalidContent` |
| `artifact_contracts.CONTENT_CONTRACTS` | Declaração do contrato de cada tipo, ao lado dos tipos, por nome e sem importar nada |
| `result_transactions.py` | `record_objects` e `discard_unreferenced`: a recuperação do journal passa a cobrir o conteúdo |

`calculation_cache` passou a delegar: `copy_artifact` chama
`content_store.clone_or_copy`, `object_path` devolve o caminho do armazenamento,
`remember` publica por `store.put` e `restore` materializa por
`store.materialize`. As tabelas, os leitores, as chaves e os endereços são os
mesmos — e é justamente por isso que caches já calculados continuam válidos.

## Deduplicação e custo: a medição

Cenário real do próximo passo: **um receptor preparado publicado por três
blocos** (redocking, Vina e DOCK6), 316 KB cada.

| Medida | Valor |
| --- | --- |
| Bytes nas pastas das execuções | 948.000 B (três cópias) |
| Objetos no armazenamento gerenciado | **1** |
| Bytes no armazenamento gerenciado | 316.000 B |
| **`R_storage`** (gerenciado ÷ publicado) | **0,333** |
| Referências lógicas ao mesmo conteúdo | 3 |
| Tempo de `publish` das três | 0,081 s |

E a entrega aos jobs, que é onde o reuso aparece:

| Medida | Valor |
| --- | --- |
| Tamanho aparente das três materializações | 948.000 B |
| Hard links entre elas e o objeto | **3 de 3** |
| **Blocos extras ocupados no disco** | **0 B** |

Três blocos recebem o mesmo receptor preparado e o disco não cresce. Num sistema
sem reflink, é o hard link para um objeto somente-leitura que entrega isso — e o
modo 0444 é o que torna a partilha segura, porque o núcleo recusa um motor que
escreva em lugar.

## Identidades de cache: 16 deslocadas, 0 ao final

`stage_cache.implementation_digest` **lê `calculation_cache.py` diretamente**, e
por ele o arquivo alcança os fechos de importação das operações. Rotear o
armazenamento por ele deslocou 16 das 26 identidades.

A contribuição de cada mudança foi isolada, revertendo apenas
`calculation_cache.py` e medindo de novo:

| Mudança | Identidades deslocadas |
| --- | --- |
| Módulos novos, `CONTENT_CONTRACTS`, `result_transactions`, testes | **0** |
| `calculation_cache.py` roteado pelo armazenamento | **16** — `implementation_digest` e 15 operações |

As 16 reduzem-se a **11 pares distintos** de digests, registrados como
equivalência auditada em `resources/processing_compatibility.json`
(`audit_2026_11_28_content_store`, 68 → 79 aliases). O fundamento é verificável
item a item:

- `copy_artifact` delega com a mesma tentativa de copy-on-write, o mesmo
  fallback de cópia independente e o mesmo `replace` atômico;
- `object_path` devolve o mesmo `objects/<prefixo>/<digest>`;
- `remember` publica o mesmo digest dos mesmos bytes, confiando num objeto já
  presente pelo mesmo cache de estampa daquele índice;
- `restore` materializa a mesma cópia isolada e mantém a mesma verificação
  posterior;
- o que é novo — verificar os bytes escritos antes do rename e gravar o objeto
  somente-leitura — protege conteúdo e não calcula nada.

Medido depois do registro: **as 26 identidades voltaram ao valor anterior**.
Preparações, redockings e docagens já calculados continuam reaproveitáveis, e os
87 testes dos módulos de fronteira que dependem dessa compatibilidade
(`test_ligand_preparation_performance`, `test_chimera_only_preparation`,
`test_docking_performance`, `test_unified_results_review`, `test_stage_cache`,
`test_restart_reuse`, `test_project_portability`, `test_artifact_publication`,
`test_preparation_resources`) passam.

## Publicação transacional: os quatro pontos de interrupção

A invariante é a ordem: **os bytes chegam ao armazenamento antes de qualquer
linha os referenciar**. Uma referência pendente torna-se impossível; o único
resíduo possível é conteúdo que ninguém referencia.

| Interrupção | Estado | Recuperação | Teste |
| --- | --- | --- | --- |
| Antes do commit | Journal sem marcador | Descarta o conteúdo do journal que nada referencia | `test_an_interruption_before_the_commit_leaves_no_dangling_reference` |
| Antes do commit, conteúdo não reivindicado | Idem, e o artefato já não existe | Descarta o objeto | `test_an_interruption_before_the_commit_discards_unclaimed_content` |
| Depois do commit | Journal e marcador | Conclui a limpeza e mantém o conteúdo | `test_an_interruption_after_the_commit_keeps_the_published_content` |
| Transação do catálogo desfeita | SQLite desfaz as linhas | Bytes verificados permanecem; nenhuma linha os registra | `test_a_rolled_back_catalog_transaction_leaves_the_store_consistent` |

A recuperação é conservadora por decisão: só descarta bytes quando nenhum
artefato lógico do projeto declara aquele conteúdo e nenhum objeto registrado
aponta para ele. REQ-STO-005 proíbe apagar conteúdo científico
automaticamente, e consistência de catálogo não é motivo para apagar bytes.

## Contratos de conteúdo

| Contrato | Tipos | O que verifica |
| --- | --- | --- |
| `table` | `compounds`, `chembl`, `fingerprints`, `similarity`, `scores` | Cabeçalho presente e as colunas que o tipo declara |
| `structure` | `structures`, `receptor_pdb` | Pelo menos um átomo com coordenadas finitas; `_atom_site` em mmCIF; seção de átomos em MOL2 |
| `prepared` | `prepared_structures` | `structure` mais o validador que o motor exige (`ligand_cache.validate_prepared_ligand`) |
| `pose` | `molecular_poses`, `vina`, `dock6` | `structure` sobre o primeiro modelo |
| `download_list` | `zinc_urls` | A lista é interpretável por `zinc_retrieval.read_download_list` |
| `present` | `dynamics`, `other` | Arquivo presente e não vazio |

Dois limites merecem registro. **`scores` não declara esquema** em
`artifact_contracts.TYPES`, porque as suas colunas variam por motor — então o
contrato dele verifica cabeçalho, não colunas nomeadas; um contrato é tão forte
quanto a declaração do tipo. E **o tipo `other` verifica apenas presença**:
um arquivo cujo tipo o usuário não declarou não tem esquema contra o qual ser
verificado, o que é a mesma razão pela qual a decisão C5 recusou `other` na
recuperação ZINC.

As amostras moleculares são examinadas em até 4.096 linhas: a pergunta é se o
arquivo é um artefato molecular, e publicar quinze mil poses não deve custar a
leitura integral delas duas vezes.

## Uma mudança de teste que a proteção exigiu

Três testes corrompiam um objeto do cache escrevendo nele, para provar que
corrupção não é acerto. Com o objeto em modo 0444 o núcleo recusa a escrita, e
os três passaram a falhar — a proteção funcionando. Eles agora adulteram como um
ator externo precisaria: desbloqueando o arquivo, escrevendo e bloqueando de
novo, o que preserva o inode e, com ele, a intenção exata de cada teste
(inclusive a do que restaura o `mtime`). A intenção científica — corrupção
detectada no restore — está intacta.

## Localização: 19 mensagens novas, e duas armadilhas

A auditoria de localização exige que todo literal de interface em `src/` tenha
tradução nos dois sentidos. Os três módulos novos introduziram 19 mensagens, e
a suíte completa acusou 20 erros que a execução dos módulos de fronteira não
mostrava.

**Reserializar os catálogos apagaria linhas.** `en.json` contém chaves
duplicadas — `"Visualização"` aparece com `"View"` e com `"Visualization"`. O
parser colapsa duplicatas no último valor, então gravar o dicionário de volta
deletaria silenciosamente as linhas sombreadas: o runtime não mudaria, porque
ele já só vê a última, mas o arquivo perderia informação. As entradas passaram
a ser inseridas como texto antes da chave de fechamento, e o diff é aditivo —
as quatro linhas "removidas" são as últimas entradas de cada arquivo ganhando
uma vírgula.

**Mensagem com valor interpolado não vai no catálogo de frases.** A auditoria
substitui `{0}` por um marcador antes de traduzir, e quem resolve isso é o
dicionário de padrões. Das 19 mensagens, 5 são frases fixas (`en.json`,
`pt.json`) e 14 têm marcador (`patterns.json`, `patterns_pt.json`). Com a
separação correta, a auditoria passou de 20 erros a zero.

## Estado dos requisitos

| Identificador | Estado | Evidência |
| --- | --- | --- |
| AC-02-01 | `IMPLEMENTED` | `test_identical_content_from_two_producers_is_stored_once`, `test_two_runs_publishing_the_same_bytes_share_one_managed_object`, `R_storage` = 0,333 |
| AC-02-02 | `IMPLEMENTED` | `test_a_writable_materialization_cannot_reach_the_stored_object`, `test_materialization_is_shared_by_default_and_isolated_on_request` |
| AC-02-03 | `IMPLEMENTED` | `test_provenance_reads_the_producer_and_its_inputs_from_the_catalog` |
| AC-02-04 | `IMPLEMENTED` | `test_usage_reports_the_store_and_what_the_catalog_claims` e as medições acima |
| REQ-STO-001 | `IMPLEMENTED` | Conteúdo idêntico é um objeto no armazenamento gerenciado |
| REQ-STO-002 | `IMPLEMENTED` | Verificação na leitura e antes de todo reuso; adulteração recusada |
| REQ-STO-003 | `IMPLEMENTED` | `isolated` para quem escreve; o objeto recusa escrita |
| REQ-STO-004 | `IMPLEMENTED` | Temporário, `fsync`, verificação, `rename`; publicação concorrente segura |
| REQ-STO-005 | `IMPLEMENTED` nesta fase | A recuperação nunca apaga conteúdo referenciado |
| REQ-STO-006 | `IN_PROGRESS` | Resultados publicados continuam nas pastas das execuções; ADR-010 (PHASE-07) decide a mudança |
| REQ-INT-002 | `IN_PROGRESS` | Publicações interrompidas recuperam-se; reconstruir o catálogo é da PHASE-07 |
| REQ-CAT-005 | `IN_PROGRESS` | `R_storage` medido; a duplicação de `artifact_instances` segue para a PHASE-05 |
| REQ-ARCH-002, REQ-ARCH-004 | `IN_PROGRESS` | Um armazenamento, dois índices; quatro dos doze ADRs aceitos |
| REQ-MIG-004 | `IN_PROGRESS` | Equivalência auditada registrada e verificada; compatibilidade de projetos é da PHASE-07 |
| REQ-TEST-003 | `IMPLEMENTED` nesta fase | Injeção de falha nos quatro pontos e publicação concorrente |

## Limites desta fase

- **Resultados publicados não entram no armazenamento.** Sem reflink,
  duplicariam bytes; por link, exigiriam provar que todo consumidor é
  somente-leitura, o que a própria #9 condiciona. ADR-010 (PHASE-07) decide,
  junto com a retenção.
- **`preparation_cache` não foi roteado.** A própria #9 o adia para a PHASE-03.
- **Dois índices sobre um armazenamento**, o do projeto e o do workspace.
  Unificá-los é da PHASE-05.
- **A varredura de órfãos em escala** é da PHASE-07; `ContentStore.sweep` já
  existe para ela.
- **Sem compressão.** Decisão revisável por medição.
- Fora de #6, segue aberta a investigação do passo 1 da triagem: topologia ou
  tipos atômicos MOL2 respondem por 694 compostos excluídos no projeto
  Alzheimer.

## Reprodução

```sh
conda activate BioMolExplorer
export PYTHONPATH=src:tests PYTHONDONTWRITEBYTECODE=1
python -m unittest discover -s tests
python -m unittest test_content_store -v
```

O ambiente precisa estar ativado, e não apenas o seu interpretador invocado por
caminho: `test_molecular_core` e `test_ligand_preparation_performance` chamam
`obabel` por nome.

Qualquer nova edição em `calculation_cache.py` desloca as mesmas 16 identidades
e exige repetir a medição e o registro.
