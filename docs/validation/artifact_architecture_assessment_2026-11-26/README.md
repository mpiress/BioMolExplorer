# Avaliação da arquitetura de artefatos científicos — PHASE-00 (#7)

Data: 26/11/2026 · Commit avaliado: `2b3d910` (branch `test`) · Issue principal: #6

Esta é a avaliação inicial exigida pela seção 02.4 de #6: avaliação da arquitetura,
inventário de artefatos, modelo de dependências, medições de referência e plano de
implementação. **Nenhum código de produção foi alterado.** As únicas adições são
dois scripts somente de leitura em `scripts/`, esta pasta de evidências e o
processo de ADRs em `docs/`.

| Documento | Conteúdo |
| --- | --- |
| `README.md` (este) | Arquitetura atual, inventários, medições, achados, invariantes, riscos e plano |
| [`requirements.md`](requirements.md) | Registro proposto de requisitos `REQ-*` e critérios `AC-*`, para revisão em #6 |
| [`dependency_model.md`](dependency_model.md) | Grafo de contratos entre operações e grafo real do projeto Alzheimer |
| `audit.json` | Saída completa de `scripts/audit_artifact_architecture.py` |
| `baseline.json` | Saída completa de `scripts/benchmark_artifact_baseline.py` |
| `tests.json` | Resultado da suíte por módulo, com cada teste que falhou e a última linha do erro |

## 1. Método

- **Somente leitura.** O banco do workspace foi aberto com `mode=ro`. Os projetos
  foram percorridos e lidos, nunca escritos. Nenhuma aplicação, worker, migração
  ou motor científico foi iniciado. Os índices descartáveis dos benchmarks foram
  criados numa pasta temporária fora dos projetos.
- **Dados reais.** Foram usados todos os projetos presentes nesta máquina:

| Rótulo | Projeto | Situação | Execuções | Tamanho em disco |
| --- | --- | --- | --- | --- |
| small | `5187a2bd…` (pasta `5b53ff40…`) | excluído (`deleted=1`) | 2 | 0,199 GB |
| medium-history | `ba933e6a…` | excluído | 14 | 0,386 GB |
| medium | `44701cb5…` (pasta `9cd0748b…`) | excluído | 6 | 2,027 GB |
| large | `947c3e18…` Alzheimer (`~/Downloads/Alzheimer`) | ativo | 23 | 234,394 GB |

  Os três projetos menores estão marcados como excluídos, mas seus dados continuam
  no disco e no banco (F-09). Eles foram usados porque são os únicos dados reais
  de tamanho pequeno e médio disponíveis. Não existe projeto pequeno ativo.
- **Reprodução:**

```sh
PY=~/progs/miniconda3/envs/BioMolExplorer/bin/python
DB=~/.local/share/biomolexplorer/workspace.sqlite3
P=~/.local/share/biomolexplorer/projects
$PY scripts/audit_artifact_architecture.py --database $DB \
  --project small=$P/5b53ff40bf9c466eac16655c6aea5c57 \
  --project medium=$P/9cd0748b39e9430fbf74e9c9e77c65d5 \
  --project medium-history=$P/ba933e6af6dc48ecabc566977d4a9e90 \
  --project large=~/Downloads/Alzheimer \
  --hash-duplicates --snapshot-sample 2 --table-sizes --output audit.json
$PY scripts/benchmark_artifact_baseline.py --database $DB <mesmos --project> \
  --repeats 5 --scratch /tmp/bme-baseline --output baseline.json
# A suíte completa trava (F-17); execute por módulo, com limite de tempo:
ls tests/test_*.py | xargs -n1 basename | sed 's/\.py$//' | xargs -P 6 -I{} \
  sh -c 'PYTHONPATH=src:tests timeout 600 '"$PY"' -m unittest {} > {}.log 2>&1; echo "{} $?"'
```

## 2. Ambiente

| Item | Valor |
| --- | --- |
| Sistema | Zorin OS 18.1 (base Ubuntu 24.04), Linux 7.0.0-44, x86_64 |
| Hardware | 16 CPUs lógicas, 69 GiB de RAM, NVMe |
| Sistema de arquivos dos projetos | ext4 (`/dev/nvme0n1p2`). **Reflink (`FICLONE`) não é suportado** (sondado no mesmo dispositivo) |
| Python científico | 3.12.3 (`envs/BioMolExplorer`), RDKit 2023.09.6, pandas 2.1.4, NumPy 1.26.4, NetworkX 3.3 |
| Motores | Vina e Open Babel no ambiente conda; Chimera em `/usr/local/bin`; DOCK6 6.11 e Antechamber em `~/progs/dock6/bin` |
| Implantação | Backend local Linux com interface Flet (desktop ou servida); SQLite em WAL; um supervisor por workspace (`docs/architecture.md`) |

## 3. Arquitetura atual por camada (seção 03.2)

| Camada | Módulos atuais |
| --- | --- |
| Motores computacionais | `caad/` (`docking.py` ativo, `redocking.py` alternativo), `kernel/` (descritores, fingerprints, LSH, filtros), `wrappers/`, `crawlers/`, `vina_batch.py`, `vina_session.py`, `dock6_*`, `md_engine.py`, `chimera_*`, `molecular_core.py`, `molecular_adapters.py`, `admet_*` |
| Integração do pipeline | `pipeline.py` (`PipelineService`), `flow.py`, `execution_order.py`, `jobs.py`, `worker.py`, `operations.py`, `artifact_contracts.py`, `bindings.py`, `semantic_references.py`, `input_validation.py` |
| Serviços de dados científicos | `artifacts.py`, `artifact_registry.py`, `artifact_naming.py`, `stage_cache.py`, `published_results.py`, `calculation_cache.py`, `preparation_cache.py`, `ligand_cache.py`, `structure_cache.py`, `dock6_cache.py`, `redocking_cache.py`, `md_cache.py`, `molecular_cache.py`, `legacy_cache.py`, `docking_reuse.py`, `redocking_reuse.py` |
| Persistência | `workspace.py` (SQLite do workspace), `jobs.sqlite3`, `storage.py`, `result_transactions.py`, `project_state.py` (histórico), `project_archive.py` (`.bme.zip`) |
| Acesso unificado a dados | `result_catalog.py` (projeções `result_read_*`), `result_index.py`, `result_files.py`, `result_tables.py`, `compound_tables.py`, `docking_results.py`, `selection_views.py` |
| Apresentação | `ui/` (`app.py`, `results_viewer.py`, `graph_explorer.py`, `molecular_inspector.py`, `pose_settings.py`…), `visualizations.py`, `docking_scene.py`, `pdb_view.py`, `scientific_report.py` |

Fronteiras violadas hoje: leitores e telas ainda descobrem resultados percorrendo
pastas (F-07), e `docking_data.py` mistura preparação, descoberta e leitura.

## 4. Inventário de persistência

| Armazenamento | Local | Chave e validação | Publicação | Escopo |
| --- | --- | --- | --- | --- |
| Estado das execuções (autoritativo) | `workspace.sqlite3`, `runs.stages` (JSON) | id da execução; triggers incrementam `result_read_versions` | linha inteira reescrita a cada mudança de resultado | workspace |
| Catálogo de artefatos | `execution_manifests`, `artifact_instances`, `artifact_entities`, `artifact_entity_indexes`, `artifact_migration_issues` | `(project, indexed_run, occurrence)`; SHA-256 por conteúdo; dependências e procedência como JSON | transacional, na publicação | workspace |
| Projeções de leitura | `result_read_versions`, `result_read_catalog`, `result_read_stages` | revisão + `schema_version`; zlib + JSON | preenchidas sob demanda na leitura | workspace |
| Assinatura de etapa | `stage_cache.stage_key` (dentro de `runs.stages`) | SHA-256 de JSON canônico: operação, parâmetros efetivos, vínculos, hashes de entrada, `processing_version`, runtime, templates | — | projeto |
| Checkpoints de cálculo | `CalculationCache` (`experiments`, `objects`, `artifacts`, `subjects`, `content_checks`) em `.input-cache/.calculation-cache` ou `BIOMOL_PREPARATION_CACHE/../calculations` | namespace + chave; objetos por SHA-256 conferidos antes do uso | transacional, após validação, com locks por chave (≤ 64) | projeto ou compartilhado |
| Preparo nativo | `preparation_cache` (`.biomolexplorer/cache/preparation`) | `receptor_key` + `tool_manifest` | renomeação atômica após validar todos os arquivos | projeto |
| Caches especializados | `ligand_cache`, `structure_cache`, `dock6_cache`, `redocking_cache`, `md_cache`, `vina_batch` | chaves próprias de cada módulo | checkpoints por par, composto ou lote | projeto |
| Fingerprints | `molecular_cache` (SQLite próprio, `cache/calculations/molecular`) | chave por molécula e protocolo | memoização | projeto |
| Índices de tabelas | `.biomolexplorer/cache/result-index/*.sqlite3` | dispositivo/inode/tamanho/mtime/ctime + SHA-256 | descartável | projeto |
| Histórico | `.history/snapshots/*.json`, `.history/blobs/<sha256>`, `.history/run-stages/<sha256>.json` | blobs por SHA-256; snapshot = estado completo antes/depois do evento | cópia completa a cada evento | projeto |
| Visões de entrada | `.input-cache/<hash>/…` e `runs/<run>/<stage>/inputs/…` | hash das seleções | cópia materializada | projeto |
| Tarefas | `.biomolexplorer/jobs.sqlite3`; pastas `jobs/<uuid>/artifacts/` | UUID da tarefa | as saídas publicadas continuam dentro da pasta do job | execução |

## 5. Inventário de artefatos

São 17 tipos de bloco: 16 operações (`operations.OPERATIONS`) mais a importação.
Os tipos e portas tipados estão em `dependency_model.md`. Distribuição física no
projeto grande, por extensão, em `runs/`:

| Extensão | Arquivos | Volume | Origem principal |
| --- | --- | --- | --- |
| `.pdbqt` | 113.030 | 2,58 GB | ligantes e receptores preparados e poses do Vina |
| `.csv` | 4.789 | 1,31 GB | tabelas de compostos, scores, consenso e similaridade |
| `.txt` | 29.090 | 0,97 GB | logs nativos e saídas dos motores |
| `.nrg` | 14 | 0,69 GB | grids de energia do DOCK6 |
| `.pdf` | 25.410 | 0,60 GB | relatórios e figuras |
| `.json` | 6.558 | 0,51 GB | manifests, relatórios e visões |
| `.mol2` | 66.242 | 0,34 GB | ligantes e poses do DOCK6 |

No catálogo do projeto grande, **63% das linhas (256.160 de 405.919) têm tipo
`other`**: a maior parte dos arquivos publicados não tem contrato científico nem
entidade associada (F-10). Os tipos com entidade são `molecular_poses` (142.751),
`compounds` (4.919), `receptor_pdb` (879), `structures` (546), `similarity` e
`fingerprints` (308 cada), `dynamics` (36) e `scores` (12).

## 6. Medições de referência

### 6.1 Armazenamento

| Projeto | Total | Histórico (`.history`) | Resultados (`runs/`, incluindo pastas de jobs) | Cópias de entrada | Bytes duplicados idênticos |
| --- | --- | --- | --- | --- | --- |
| small | 0,199 GB | 0,091 GB (46%) | 0,101 GB | 0,007 GB | **0,108 GB (54%)** |
| medium-history | 0,386 GB | 0,238 GB (62%) | 0,148 GB | 0,000 GB | **0,160 GB (41%)** |
| medium | 2,027 GB | 0,877 GB (43%) | 1,025 GB | 0,125 GB | **1,154 GB (57%)** |
| large | 234,394 GB | **224,185 GB (95,6%)** | 8,189 GB | 1,996 GB | **18,689 GB (8,0%)** |

Projeto grande, detalhado:

| Área | Volume | Arquivos |
| --- | --- | --- |
| `.history/snapshots` | 211,688 GB | 463 |
| `.history/run-stages` | 9,769 GB | 32 |
| `runs` (pastas `.biomolexplorer/jobs/<uuid>`) | 5,741 GB | 240.954 |
| `.history/blobs` | 2,728 GB | 116.387 |
| `runs` (pastas das etapas) | 2,448 GB | 69.855 |
| `.input-cache` | 1,975 GB | 73.338 |
| `.input-cache/.calculation-cache` | 0,021 GB | 37 |

- **Conteúdo científico versus metadados:** o conteúdo científico (resultados e
  cópias de entrada) soma 10,2 GB, e o histórico mais o estado das execuções somam
  224 GB. Fora do projeto, o banco do workspace tem mais 9,16 GB.
- **Snapshots:** 463 arquivos. 39 têm mais de 1 GiB, o maior tem 6,92 GB e a
  mediana é 27,4 MB. Só 44 usam `_compact_runs` (4,64 GB) e 6 são apenas
  referências. Num snapshot amostrado de 1,57 GB, **1.523,8 MB (97%) são o JSON
  das execuções** e 46,3 MB são a lista de 141.799 arquivos.
- **Bytes duplicados idênticos no projeto grande (18,689 GB):**
  - 8,540 GB são snapshots idênticos byte a byte (o "depois" de um evento é igual
    ao "antes" do evento seguinte; o maior par tem 2,85 GB);
  - 10,15 GB são conteúdo copiado entre `.history/blobs`, `.input-cache`, as pastas
    das etapas e as pastas dos jobs (`.json` 9,05 GB, `.pdbqt` 4,05 GB, `.csv`
    1,64 GB, `.txt` 0,97 GB, `.nrg` 0,69 GB).
- **Receptores preparados:** 55 nomes em 332 arquivos e 57 conteúdos distintos
  (217,6 MB). Exemplos: `4M0E_A.dockprep.mol2`, `4M0E_A.dockprep.pdbqt`,
  `2WSL_A.dockprep.*` e `2QMG_A.dockprep.*` aparecem em **36 a 37 cópias cada,
  todas com o mesmo conteúdo**.

### 6.2 Banco do workspace e catálogo

- **Banco:** 9,158 GB, dos quais 1,351 GB (14,8%) são páginas livres. Por tabela:
  `runs` 6.382,5 MB, `artifact_instances` 702,0 MB (+ 69,1 MB de chave primária e
  101 MB de índices), `artifact_entities` 147,1 MB (+ 146 MB de índices),
  `result_read_stages` 108,6 MB, `result_read_catalog` 74,6 MB,
  `artifact_entity_indexes` 30,8 MB e `execution_manifests` 30,2 MB.
- **Projeções de leitura:** `result_read_catalog` tem 23 linhas para 52
  execuções. Só parte das execuções tem a projeção pronta (F-14).
- **JSON das execuções do projeto grande:** 23 execuções somam 6,368 GB, dos quais
  **5,834 GB (91,6%) são `artifact_records`**. A execução mais recente tem
  675,7 MB, mas todas as suas etapas foram reaproveitadas e ela durou cerca de
  5 minutos. Ali, `docking_vina` tem 5.460 registros (253,7 MB, cerca de 47 KB
  cada) e `dynamics_dock6` tem 66 registros (133,4 MB, cerca de 2 MB cada).
  Decodificar esse JSON leva 1,8 a 1,9 s.
- **Catálogo do projeto grande:**

| Medida | Valor |
| --- | --- |
| Linhas em `artifact_instances` | 405.919 |
| Ocorrências distintas | 122.572 |
| Artefatos lógicos distintos | 36.873 |
| Conteúdos distintos (SHA-256) | 34.271 |
| Linhas sem hash de conteúdo | 5.138 |
| Linhas reindexadas a partir de outra execução | **283.347 (69,8%)** |
| Conteúdos com mais de um caminho físico | 25.154 (83.208 caminhos extras) |
| JSON de dependências / procedência / metadados | 262,4 MB / 17,9 MB / 10,0 MB |
| Entidades indexadas | 118.125 ocorrências: pose 256.301, compound 222.286, molecular_pair 17.880, community 3.372, structure 524, receptor 352 |
| Problemas de migração registrados | 0 neste projeto (1 no workspace) |

### 6.3 Computação e reuso (histórico real)

Dados extraídos de `runs.stages` (estados e horários gravados pela aplicação).
Cada valor é uma observação única, não uma repetição controlada.

| Projeto | Etapas concluídas | Calculadas | Reaproveitadas | Acerto no nível da etapa | Recálculos redundantes | ACR histórico |
| --- | --- | --- | --- | --- | --- | --- |
| small | 4 | 4 | 0 | 0% | 0 | — |
| medium | 13 | 13 | 0 | 0% | 0 | — |
| medium-history | 42 | 22 | 20 | 47,6% | 0 | 100% |
| large | 295 | 51 | 244 | **82,7%** | **13 (3.310 s)** | **94,9%** |

Um recálculo é **redundante** quando a mesma `cache_key` já tinha sido concluída e a
etapa foi calculada de novo. A métrica não separa recusas do usuário de falhas
de correspondência.

Tempo calculado no projeto grande, por operação: `docking_vina` 10,21 h,
`docking_dock6` 3,12 h, `redocking` 2,39 h, `retrieve_pubchem` 0,83 h,
`retrieve_compounds` 0,60 h. As demais operações ficaram abaixo de 0,2 h cada.

Frio, morno e quente, a partir de execuções reais do projeto grande:

| Operação | Primeiro cálculo | Recalculada com checkpoints internos | Reuso do resultado publicado |
| --- | --- | --- | --- |
| `docking_vina` (mesma chave `ea80b750e5`) | 18.377 s | 337 – 502 s (4 vezes) | 0 – 0,1 s |
| `redocking` (mesma chave `502f099600`) | 465 s (5.773 s na primeira configuração) | 181 – 233 s (5 vezes) | 0 s |
| `docking_dock6` | 8.226 s | 977 – 1.022 s (chaves diferentes) | 1,6 – 2 s |

Os checkpoints internos (`vina_batch`, `redocking_cache`, `dock6_cache`) evitam
até 97% do custo quando o reuso no nível da etapa não acontece. Mas o estado
gravado não diz por que a etapa não foi reaproveitada (F-11).

### 6.4 Recuperação de dados (latência)

Medida com `scripts/benchmark_artifact_baseline.py`: 5 repetições por caso, mediana em
segundos (desvio-padrão em `baseline.json`). **Frio** significa que os arquivos
envolvidos foram despejados do cache de páginas antes de cada repetição. **Morno**
significa uma chamada de preparo seguida das repetições. Todas as repetições, frias e
mornas, deram saída idêntica. Cada projeto foi medido na execução concluída mais
recente.

| Caso | small (frio / morno) | medium (frio / morno) | large (frio / morno) |
| --- | --- | --- | --- |
| `os.walk` de `runs/` (só morno) | — / 0,002 (1.915 arq.) | — / 0,004 (4.393) | — / **0,143 (310.809)** |
| `rglob('*.csv')` em `runs/` (só morno) | — / 0,005 | — / 0,008 | — / 0,300 |
| Catálogo: artefatos da maior etapa | 0,003 / 0,0003 (120) | 0,003 / 0,0002 (113) | 0,390 / 0,016 (15.258) |
| Catálogo: 50 consultas por entidade | sem entidades | sem entidades | 0,120 / 0,003 |
| Catálogo: 50 consultas por conteúdo | 0,004 / 0,0004 | 0,009 / 0,0004 | 0,070 / 0,002 |
| JSON completo da execução (`json.loads`, caminho de `get_run`) | 0,011 / 0,002 | 0,008 / 0,001 | **6,83 / 2,65** (675,7 MB) |
| Projeção sem cache (`_project_blob`, lê descritor a descritor) | 0,012 / 0,002 | 0,008 / 0,002 | **13,34 / 9,82** |
| Projeção em cache (`read_run`) | sem cache | sem cache | 5,09 / 0,62 |
| Provas em cache (`read_proofs`) | sem cache | sem cache | 0,24 / 0,13 |
| Snapshot de histórico de tamanho mediano | 0,0001 (2 KB) | 0,0001 (2 KB) | 0,17 / 0,20 (27,4 MB) |

O projeto medium-history teve resultados equivalentes aos do medium e está em
`baseline.json`.

Tabelas CSV registradas (as três maiores da execução mais recente de cada projeto):

| Projeto | Tabela (linhas, tamanho) | Leitura completa | Índice: construção / tamanho | Página de 25 linhas (frio / morno) |
| --- | --- | --- | --- | --- |
| small | `matches.csv` (306.090 linhas, 23,3 MB) | 0,31 s | 1,71 s / 82,4 MB | 0,022 / 0,008 |
| medium-history | `matches.csv` (519.562 linhas, 39,1 MB) | 0,54 s | 2,81 s / 139,0 MB | 0,036 / 0,013 |
| medium | `morgan_selected_compounds.csv` (6.446 linhas, 40,7 MB) | 0,16 s | 0,60 s / 86,1 MB | 0,045 / 0,001 |
| large | `…_SIMS_MIX_evaluated.csv` (1.705 linhas, 9,3 MB) | 0,05 s | 0,36 s / 24,0 MB | 0,033 / 0,002 |

Leituras:
- O **custo dominante é o descritor das execuções, não a descoberta.** No projeto
  grande, abrir uma execução custa 2,7 a 6,8 s, ou 9,8 a 13,3 s pelo caminho
  incremental sem cache. Uma consulta indexada ao catálogo da mesma etapa custa 0,016
  a 0,39 s. As consultas usam índice:
  - `artifacts_project_stage`, `artifact_entities_lookup` e `artifacts_content` são
    todos `SEARCH … USING INDEX`.
- A projeção sem cache é cerca de 2× mais lenta que `json.loads` do mesmo JSON. Ela só
  compensa quando é gravada uma vez e reutilizada: em cache, `read_run` leva 0,62 s
  morno.
- Com metadados de diretório em cache, a varredura de `runs/` custa só 0,14 a 0,30 s
  no projeto grande. O custo de F-07 aparece quando a varredura se repete (por
  exemplo, a cada evento de histórico) e quando ela também faz `stat` e hash, como na
  auditoria completa (11,6 s).
- O índice de tabelas ocupa de 2,1× a 4,2× o tamanho do CSV de origem (por exemplo,
  23,3 MB → 82,4 MB). Esse sobrecusto entra em `IO`.

### 6.5 Suíte de testes

A suíte completa (`unittest discover -s tests`) **travou**: parou de produzir saída
depois de cerca de 10 minutos e ficou mais de 11 minutos esperando num laço
assíncrono. Por isso, cada um dos 169 módulos foi executado separadamente, com
limite de 600 s e 6 processos em paralelo. Os 45 módulos com problema foram
repetidos **em série**, com o mesmo resultado, então nenhuma falha vem do
paralelismo.

| Medida | Valor |
| --- | --- |
| Módulos | 169 (124 sem falhas, 45 com falhas, erros ou travamento) |
| Testes executados | 1.663 |
| Falhas / erros | 50 / 38 |
| Ignorados | 1 |
| Travamento | `test_run_reuse_ui.ReusePromptTests.test_late_result_lookup_cannot_open_prompt_after_leaving` (não termina; os outros 7 testes do módulo dão erro por fixture desatualizada: `'WorkspaceUI' object has no attribute 'store'`, `KeyError: 'pipeline'`) |

Causas mais frequentes:
- 15 erros `ValueError: Confirme a invalidação dos resultados antes de salvar esta
  alteração`: testes anteriores à regra de confirmação de invalidação.
- Contagens de arquivos ADMET esperadas que não batem com as atuais.
- Rótulos e traduções divergentes (`test_localization*`, `test_pipeline_guidance`).
- `ImportError: MINIMIZATION` em `test_native_structures`.
- `test_chembl_rest`, que depende de rede.

Módulos com falhas em áreas que a #6 vai alterar:
- `test_stage_cache`, `test_restart_reuse` e `test_result_persistence`, incluindo
  `CrashRecoveryTests` antes e depois do commit;
- `test_pipeline_failure_recovery`, `test_ligand_preparation_performance`,
  `test_input_processing` e `test_global_artifact_names`.

A lista completa, com a última linha de cada erro, está em `tests.json`.
O código de produção não foi alterado nesta avaliação. A única mudança que toca
o que os testes leem é a documentação, e `test_installation_guide` passa.

## 7. Achados, ordenados por impacto medido

| # | Achado | Impacto medido (projeto grande, salvo indicação) | Requisitos |
| --- | --- | --- | --- |
| F-01 | Cada evento de histórico grava o estado completo do projeto, antes e depois, incluindo o JSON completo das execuções (`project_state.record`/`snapshot`; `compact_runs` só em 44 de 463 arquivos). | 211,7 GB em 463 snapshots (90% do projeto); 8,54 GB de snapshots idênticos | REQ-STO-006 |
| F-02 | Uma execução que reaproveita resultados copia todos os `artifact_records` da execução de origem para o seu próprio `runs.stages`. | 5,83 GB de registros em 6,37 GB de JSON; 675,7 MB numa execução de cerca de 5 min | REQ-CAT-005 |
| F-03 | O catálogo reindexa por execução os artefatos reaproveitados, e `.history/run-stages` guarda outra cópia do JSON de cada execução. | 283.347 linhas reindexadas (69,8%); 9,77 GB em `run-stages`; banco com 9,16 GB | REQ-CAT-005 |
| F-05 | `.history/blobs` guarda cópias completas do conteúdo dos resultados (sem reflink no ext4). | Maior fonte de duplicação nos três projetos menores (41 – 57% do tamanho) | REQ-STO-006 |
| F-06 | As saídas publicadas ficam dentro das pastas de trabalho dos jobs. Entradas são materializadas por cópia em `.input-cache` e em `runs/…/inputs/`. | 5,74 GB e 240.954 arquivos em pastas de jobs; receptores com 36 – 37 cópias idênticas; 10,15 GB de conteúdo duplicado | REQ-STO-006, REQ-STO-003 |
| F-04 | As dependências existem só como JSON (`artifact_instances.dependencies`, `artifact_records`). Não há tabela de operações nem de arestas. | 262,4 MB de JSON de dependências; a pergunta "quem depende de X" exige decodificar as execuções | REQ-CAT-001, REQ-DEP-001 |
| F-07 | Varreduras recursivas em caminhos de rotina, por exemplo `project_state.snapshot` (todo o projeto a cada evento), `docking_data.py:307,325`, `pose_sources.py:214,218`, `pipeline.py:1168,1584,1734`, `md_inputs.py:133`, `ui/app.py:907`, `clear_results.py`, `result_sync.py:135`. | Percorrer o projeto grande leva 11,6 s com metadados já em cache (501.082 arquivos); ver 6.4 | REQ-VIS-001, REQ-CAT-003 |
| F-08 | Há sete caches especializados com raízes e formatos de chave próprios, mais `stage_key` e `CalculationCache`. | Políticas de chave e validação diferentes por módulo; o mesmo receptor é materializado várias vezes | REQ-ARCH-002, REQ-CACHE-002 |
| F-11 | Recálculos redundantes com a mesma `cache_key` não registram o motivo da recusa do reuso. | 13 recálculos (0,92 h); `redocking` 6× e `docking_vina` 5× com a mesma chave | REQ-REUSE-003, REQ-CACHE-001 |
| F-09 | Projetos excluídos mantêm dados no disco e linhas no banco (exclusão lógica sem política de retenção). | 4 projetos excluídos, 29 execuções, cerca de 2,6 GB em disco | REQ-MIG-003 |
| F-10 | A maior parte dos artefatos publicados não tem tipo científico (`other`). | 63% das linhas do catálogo | REQ-ID-003, REQ-API-002 |
| F-14 | Projeções de leitura só existem para parte das execuções, e o banco tem páginas livres não recuperadas. | 23 de 52 execuções com projeção; 1,35 GB livres | REQ-CAT-005, REQ-PERF-002 |
| F-12 | O sistema de arquivos não oferece reflink, então toda cópia custa os bytes inteiros. | Sondagem `FICLONE`: "Operation not supported" | ADR-003 |
| F-15 | Há uma localização de projeto que aponta para uma pasta inexistente (`f55932e9…` → `~/Downloads/BioMol`). | 1 projeto | REQ-INT-001 |
| F-16 | A suíte de testes de referência não passa: 45 de 169 módulos falham, incluindo cache, reuso e recuperação após falha. | 50 falhas e 38 erros em 1.663 testes | REQ-TEST-001, REQ-ENG-008 |
| F-17 | Um teste de interface trava indefinidamente e impede a execução completa da suíte. | `test_late_result_lookup_cannot_open_prompt_after_leaving` | REQ-TEST-001 |
| F-13 | Arquivos `.pyc` em `src/**/__pycache__` são versionados, e a suíte de testes os reescreve. | 10 arquivos versionados alterados por uma execução da suíte (restaurados depois) | REQ-TEST-001 (higiene) |

## 8. Invariantes científicos (seção 01.3)

| Invariante | Onde é implementado | Proteção atual contra reuso incompatível |
| --- | --- | --- |
| Algoritmos científicos existentes | `caad/`, `kernel/`, `wrappers/`, módulos científicos de `biomolexplorer/` | `stage_cache.operation_digest` / `implementation_digest`; equivalências exatas e auditadas em `resources/processing_compatibility.json` |
| Protocolos de preparação molecular | `chimera_preparation.py`, `docking_preparation.py`, `molecular_core.py`, `ligand_mol2.py`, `redocking_preparation.py`, `native_preparation.py` | `preparation_cache.receptor_key` + `tool_manifest`; `ligand_cache.preparation_identity`; `PREPARATION_SCHEMA_VERSION` na chave do redocking |
| Padronização molecular | `molecule_quality.py` (remoção auditável), `structure_sanitation.py`, `molecular_parsing.py`, `kernel/filters.py` | Cópias de execução; arquivos de coordenadas nunca são editados no lugar |
| Docking e redocking | `caad/docking.py` (ativo), `vina_batch.py`, `vina_session.py`, `dock6_protocol.py`, `dock6_redocking.py`, `redocking_dock6.py`, `redocking_config.py` | Identidade por conteúdo de receptor e ligante, opções efetivas e binário (`vina_batch`, `dock6_cache`); `docking_reuse.accounted_preparation_exclusions` |
| Semântica das predições ADMET | `caad/admet.py`, `wrappers/admet.py`, `admet_config.py`, `admet_outputs.py`, `admet_predictors.py`, `admet_ai_runner.py` | `admet_output_schema` na chave da etapa |
| Fingerprints | `kernel/descriptors.py` (gerador Morgan, MinHash/LSH), `molecular_cache.py` | Chave por molécula e protocolo |
| Similaridade e redes | `kernel/descriptors.py` (Tanimoto), `caad/complex_network.py`, `caad/graph_results.py`, `network_analysis.py`, `network_centrality.py`, `network_scaffolds.py`, `molecular_network.py`, `graph_layouts.py` | Assinatura da etapa; entradas de grafo por SHA-256 |
| Formulação matemática do ECR | `consensus_model.py` ([modelo](../../docking_consensus.md)) | Assinatura da etapa; **fora do escopo de mudança (01.4)** |
| Agregação multialvo e cobertura | `consensus_model.py`, `consensus_targets.py` | Assinatura da etapa |
| Precisão numérica | Testes de regressão (por exemplo `test_consensus_model`, `test_molecular_core`, `test_docking_reference_regressions`) | O reuso só aceita bytes idênticos (SHA-256) |
| Reprodutibilidade | `stage_cache.runtime_signature` / `runtime_manifest` (Python, pacotes, ferramentas), binários e templates nas chaves | Mudanças de ambiente invalidam o reuso; equivalências só por lista auditada |
| Rastreabilidade experimental | `artifacts.py` (`artifact_id`, `occurrence_id`), `artifact_registry.py`, `artifact_naming.py`, `execution_manifests` | Procedência por arquivo ou lote (não por molécula ou linha) |

Regras que todas as fases precisam manter:
- uma exclusão nunca é convertida em sucesso ([reuso de docking](../../docking_cache.md));
- uma correspondência inválida é sempre um erro de cache, nunca um acerto;
- nada é reaproveitado entre ambientes diferentes sem uma equivalência auditada.

## 9. Riscos

| Risco | Efeito | Mitigação no plano |
| --- | --- | --- |
| Mudar o formato do histórico (F-01) | Perda da capacidade de restauração ou de snapshots antigos | Leitura dupla (antigo e novo) e migração não destrutiva com teste de ida e volta na PHASE-07; nenhum snapshot antigo é apagado automaticamente |
| Deduplicar sem reflink (F-12) | Hard links compartilham inode, então uma escrita no lugar corromperia o original | Objetos somente leitura; cópia isolada para quem escreve; hard link só para consumidores comprovadamente somente leitura (ADR-003) |
| Unificar assinaturas (F-08) | Invalidar caches válidos ou, pior, aceitar incompatíveis | Mapear cada chave antiga para o modelo novo com teste por dimensão de incompatibilidade; leitura das chaves antigas preservada (PHASE-03) |
| Referenciar em vez de copiar `artifact_records` (F-02) | Leitores atuais esperam registros completos em cada execução | Resolver referências na camada de acesso (PHASE-05) antes de mudar a gravação; testes de equivalência dos leitores |
| Retomar execuções interrompidas | Muda um comportamento documentado | Aprovação explícita na PHASE-04 |
| Projetos de tamanho pequeno e médio são projetos excluídos | O perfil medido pode não representar projetos ativos | Recriar um projeto de referência pequeno no início da PHASE-03 e medir de novo |
| Volume do projeto grande | Migração, auditoria e testes de ida e volta levam minutos a horas | Indexação retomável e incremental; benchmarks com amostragem documentada |

## 10. Plano de implementação

O plano segue as sub-issues #8 a #15, com estes ajustes motivados pelas medições:

0. **Pré-requisito da PHASE-01 (F-16, F-17):** a REQ-ENG-008 exige a suíte verde
   depois de cada fase, mas ela já falha hoje. Antes de alterar código, é preciso
   triar as 88 falhas e erros, corrigir ou atualizar os testes desatualizados,
   eliminar o travamento e tirar os `.pyc` do controle de versão (F-13). Isso
   deve ser feito numa tarefa à parte, fora do escopo de #6, para não misturar
   correções com a nova arquitetura. Até lá, o critério por fase é **nenhuma
   falha nova** em relação a `tests.json`.
1. **PHASE-01 (#8):** além do esquema de identidade e catálogo, definir como uma
   ocorrência reaproveitada é **referenciada** por outra execução, sem copiar a
   linha (F-03), e como `runs.stages` passa a apontar para registros do catálogo
   (F-02). As mudanças de gravação acontecem nas fases seguintes, com leitura
   dupla.
2. **PHASE-02 (#9):** o armazenamento por conteúdo também atende
   `.history/blobs` e as visões de entrada (`.input-cache`, `inputs/`) (F-05,
   F-06). As saídas passam a ser publicadas fora das pastas dos jobs. Os
   snapshots passam a registrar o conteúdo uma única vez e a referenciar o estado
   das execuções; esse é o maior ganho medido (F-01). Por envolver restauração,
   **a mudança de gravação do histórico é antecipada da PHASE-07 para a
   PHASE-02**, com leitura dos dois formatos. A migração dos snapshots existentes
   continua na PHASE-07.
3. **PHASE-03 (#10):** piloto do receptor, como planejado. Ponto de partida
   medido: 332 arquivos de receptor preparado e 57 conteúdos. Começar recriando
   um projeto de referência pequeno e ativo e medindo o preparo de receptor com
   motor real, a frio e a quente.
4. **PHASE-04 (#11):** registrar o motivo de cada decisão de reuso (F-11) antes do
   planejador incremental.
5. **PHASE-05 (#12) e PHASE-06 (#13):** prioridade para os caminhos de F-07 e para
   a projeção de leitura de todas as execuções (F-14).
6. **PHASE-07 (#14):** migrar os snapshots existentes, aplicar política de
   retenção aos projetos excluídos (F-09) e verificar localizações inválidas
   (F-15).
7. **PHASE-08 (#15):** integração restante e medições finais com os mesmos scripts.

### Limiares de aceitação propostos (para revisão — REQ-PERF-003)

| Métrica | Linha de base | Limiar proposto |
| --- | --- | --- |
| `R_storage`, projeto grande, mantendo todo o histórico restaurável | 234,4 GB | ≥ 80% (≤ 46,9 GB) |
| `R_storage`, projetos small e medium | 0,199 / 2,027 GB | ≥ 40% |
| `IO` (catálogo + índices ÷ projeto) | ver 6.2 | não crescer com o número de reusos; reportado por projeto |
| JSON gravado por uma execução totalmente reaproveitada | 675,7 MB | ≤ 5% da execução de origem |
| Abrir essa execução (projeto grande) | 6,83 s frio / 2,65 s morno (`json.loads`); 0,62 s morno com projeção | ≤ 10% da linha de base, frio e morno |
| `CHR` no nível da etapa, cenário de reexecução idêntica | 82,7% (histórico) | 100% das etapas elegíveis |
| `ACR`, cenário de reexecução idêntica | 94,9% (histórico) | 100%, com motivo registrado para toda recusa |
| Regressão em cargas sem reuso | ver 6.4 | ≤ 5% ou dentro de 2 desvios-padrão da linha de base, o que for maior |
| Descoberta de resultados nas telas prioritárias | ver 6.4 | sem varredura recursiva; latência ≤ a da linha de base com cache frio |

## 11. Limitações desta avaliação

- Os tempos de computação da seção 6.3 são observações históricas únicas, não
  repetições controladas. Repetir o pipeline grande levaria várias horas por
  execução. O benchmark controlado de preparo de receptor com motor real fica
  para o início da PHASE-03.
- O modo "frio" despeja do cache de páginas os arquivos envolvidos
  (`posix_fadvise`). Ele não despeja os metadados de diretórios, o que exigiria
  root. Por isso, as varreduras foram medidas só no modo morno.
- A latência foi medida no backend de leitura. Tempo de renderização Flet e FPS
  não foram medidos.
- O critério de recálculo redundante não separa recusas do usuário de falhas de
  correspondência, porque o estado gravado não diferencia os dois casos (F-11).
- `scripts/benchmark_molecular_preparation.py` usa um backend `python` que não
  corresponde à preparação atual (só Chimera) e não foi executado.
