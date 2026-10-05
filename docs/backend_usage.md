# Uso da camada de aplicação

[Documentação](README.md) · Português · [English](en/backend_usage.md)

Para utilizar a interface pronta com autenticação e projetos, consulte o
[guia do workspace Flet](frontend.md). Este documento descreve a CLI e os serviços
Python utilizados pela interface.

Antes de executar a CLI ou iniciar a interface, siga o [guia de instalação e configuração](installation.md): baixe o [repositório do GitHub](https://github.com/mpiress/BioMolExplorer) e instale **UCSF Chimera 1.17, DOCK6 6.11 e DMS** no computador dos cálculos. Configure seus executáveis no `PATH`. Sem essas ferramentas, as etapas que dependem delas falharão.

Crie/ative o ambiente Conda e instale o projeto a partir da raiz do código baixado:

```bash
conda env create -f environment.yml
conda activate BioMolExplorer
python -m pip install -e . --no-deps --no-build-isolation
```

As dependências científicas são gerenciadas pelos arquivos Conda. O pacote
Python fornece a camada de aplicação, os módulos científicos e os recursos;
a instalação do pacote por si só não instala Chimera 1.17, DOCK6 6.11 ou DMS.
Open Babel e Vina são declarados no ambiente Conda. Para um ambiente já existente,
atualize-o conforme necessário e verifique também esses executáveis.

## CLI

```bash
biomolexplorer retrieve_compounds \
  --parameters examples/retrieve_compounds.json \
  --output ./datasets
```

A forma equivalente é `python -m biomolexplorer ...`. A CLI executa a operação
até terminar e devolve JSON com os artefatos. Sem instalar o pacote, use
`PYTHONPATH=src python -m biomolexplorer ...` a partir da raiz do projeto.
Os scripts em `workflow/` continuam disponíveis e encontram `src` pelo local
do próprio arquivo. `/datasets` mantém o significado de pasta relativa ao
workspace; `/tmp/estudo`, por exemplo, é um caminho absoluto normal.

Os filtros ChEMBL ficam em `src/biomolexplorer/resources/crawlers`.
Também é possível fornecer um dicionário `chembl_filters` em parâmetros, com
chaves `target`, `bioactivity`, `molecules` e `similars`. Cada chave fornecida
substitui o conjunto de filtros daquela etapa; etapas omitidas usam os padrões.
Essa opção permite configurar execuções na interface sem alterar arquivos globais.

## Tarefas para Flet

```python
from pathlib import Path
from biomolexplorer.config import AppConfig
from biomolexplorer.jobs import JobManager

config = AppConfig(
    workspace=Path('/tmp/meu-estudo'),
    max_jobs=1,
    cpu_workers=2,
    max_queued_jobs=100,
    job_timeout=86400,
)
manager = JobManager(config)
job = manager.submit('retrieve_compounds', {
    'search_term': 'CHEMBL220',
    'include_pubchem': True,
    'pubchem_threshold': 75,
    'pubchem_max_records': 1000,
})
# Esta chamada retorna sem aguardar a recuperação científica.
status = manager.get(job['id'])
# manager.cancel(job['id'])
# Ao encerrar a aplicação: manager.close()
```

Use uma única instância duradoura de `JobManager` por workspace. Um bloco `with`
fecha o gerenciador ao sair e cancela tarefas pendentes; não crie um gerenciador
novo a cada clique da interface.

`examples/flet_controller.py` fornece um controlador sem dependência de widgets.
Seus métodos usam `asyncio.to_thread` para as chamadas de gerenciamento e
`updates` fornece snapshots para atualizar os controles no loop da interface:

```python
job = await controller.submit('admet', {
    'base_input_path': '/tmp/meu-estudo/datasets/compounds/CHEMBL220',
    'input_file': 'compounds.csv',
})
async for snapshot in controller.updates(job['id']):
    # Atualize controles Flet neste loop; snapshots são dicionários serializáveis.
    print(snapshot['status'], snapshot['error'])
```

Estados terminais: `succeeded`, `failed`, `cancelled`, `interrupted`.
`result.artifacts` fornece os caminhos dos arquivos; `log_path` aponta para o
log do worker. Há estado de execução, mas ainda não há percentual de progresso
por etapa. A interface deve tratar erro de validação/submissão e o estado
`failed` de uma tarefa já aceita.

Se a interface usar outro ambiente Python, configure
`AppConfig(worker_python=Path('/caminho/conda/envs/BioMolExplorer/bin/python'), ...)`.
O processo científico usa esse interpretador. Nenhuma dependência Flet é
necessária no worker. O controlador não constrói uma interface visual.

## Operações disponíveis

Consulte `biomolexplorer.operations.OPERATIONS` para os parâmetros completos.

| Operação | Entradas principais | Saídas |
| --- | --- | --- |
| `retrieve_compounds` | `search_term`, opções PubChem e filtros ChEMBL | Dados por fonte e `compounds/<alvo>/compounds.csv` |
| `expand_similar_compounds` | `search_term`, `base_input_path` contendo `ChEMBL/molecules` e `ChEMBL/similars` | Compostos novos, relações e conjunto consolidado |
| `retrieve_structures` | `target`, filtros PDB | PDBs e `pdb_codes.csv` |
| `retrieve_zinc` | `filename`, `base_input_path` contendo o arquivo de URIs | CSV ZINC |
| `admet` | `base_input_path`, opcional `input_file` | CSVs de avaliação, subconjuntos e figura |
| `fingerprints` | `base_input_path`, algoritmos, `chunk_size` | CSVs de fingerprints |
| `similarity` | `base_input_path` dos fingerprints, métrica, limiar e `approximate` | CSVs de arestas na subpasta `Similarity` |
| `graphs` | `similarity_path`; `base_input_path` opcional dos compostos; opções MCS | Grafos/MCC independentes, fragmentos e CSVs |
| `prepare_structures` | PDBs próprios, `target`, registros PDB, pH e cargas | Estruturas preparadas, centros e `pdb_codes.csv` |
| `redocking` | Diretório PDB, `target`, parâmetros de preparação | Cópias de trabalho, resultados Vina e RMSD |
| `docking_vina` | Complexos preparados, compostos selecionados, `mol_filename` | PDBQT e resultados Vina |
| `docking_dock6` | Complexos preparados, compostos, `base_vina_path`, `pdb_code`, instalação DOCK6 | Conformações, scores e footprints |
| `consensus` | Diretório contendo `Vina` e `Dock6`, ou entradas separadas `base_vina_path` e `base_dock6_path`, `target` | CSV e figura de consenso |

Enums aceitam nomes ou valores, como `TanimotoSimilarity`/`Tanimoto` e
`Morgan`/`morgan`. Listas de filtros PDB também usam strings serializáveis.
O parâmetro Vina `pdb_code` segue a função existente: uma lista de registros
`[PDB_CODE, LIGAND, RESNUM, CHAIN]`; DOCK6 recebe um único registro desses quatro
campos. Redocking pode ler os registros de `pdb_codes.csv`.

Cada tarefa recebe um novo diretório de saída. Para encadear tarefas, use os
artefatos retornados pela etapa anterior como entrada da próxima. A operação
`graphs` recebe somente similaridades prontas em `similarity_path` (CSV ou pasta).
Cada arquivo gera uma análise independente. `base_input_path` aceita uma tabela
ou pasta de compostos complementar. Os CSVs de arestas usam `source,target,value`,
com pesos entre 0 e 1; não exigem prefixos de métrica no nome. Tabelas com
`molecule_chembl_id,canonical_smiles` preservam nós isolados e permitem visualizar
estruturas e encontrar o fragmento comum. Sem elas, apenas a topologia é apresentada.

Métrica, fingerprint e limiar são definidos na operação `similarity`. A operação
`graphs` não recalcula nem filtra os pesos. `mcs_timeout` tem padrão 30 segundos,
`mcs_ring_matches_ring_only` é true e `mcs_complete_rings_only` é false.
Veja o [guia do workspace](frontend.md) para seleção, interação e downloads.

Chamadas avançadas podem fornecer `graph_inputs`, uma lista de dicionários com
`kind: "similarity"`, `file`, `compound_files` (lista, vazia se não houver estruturas)
e metadados opcionais `label`, `metric` e `fingerprint`. Estes últimos descrevem
a origem e não configuram cálculos. Entradas de tipo `fingerprints` são rejeitadas.
No workspace, os caminhos precisam pertencer ao projeto. O pipeline resolve os
SMILES pela cadeia anterior e mantém cada análise separada, em `plots/`,
`Molecules/`, `data/maxcomp/` e `centroids/`. A consulta MCS continua como SMARTS interno; o modelo e a interface também oferecem
o SMILES do fragmento extraído da molécula de referência e o estado da busca.

## Testes

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

Execute no ambiente científico. Veja [a revisão da arquitetura](architecture.md)
para os comportamentos corrigidos, as verificações e os limites atuais.

## Projetos portáteis e entradas múltiplas

`WorkspaceStore.create_project(..., directory="/pasta/nova")` define o diretório
de arquivos do projeto. `project_dir(id)` resolve tanto projetos antigos quanto
novas pastas selecionadas. Use `export_project(token,id)` e
`import_project(token,arquivo,directory)` para transferir configurações, resultados
e versões com verificação de integridade. Contas e sessões não entram no pacote.
`history(token,id)` lista alterações; `rollback(token,id,event_id)` restaura o estado
anterior ao evento e exige proprietário e ausência de execução ativa.

Um binding aceita uma referência simples ou `{"sources":[...referências...]}`.
Cada referência usa `stage` ou `asset` e `selector` opcional. CSVs selecionados são
normalizados em cópias próprias. `input_processing="individual"` separa arquivos;
`input_processing="merge"` combina as entradas escolhidas. `provided_results` contém
`kind` e `asset_ids` para bypass do cálculo com resultados validados; no ADMET,
exija código, SMILES, TPSA e WLOGP. Consulte [projetos e versões](projects.md).


## Seleção, retomada e idioma da interface

`PipelineService.submit(token, project_id, reuse_results=True)` reaproveita resultados compatíveis. A interface consulta `existing_results` e pede a decisão ao usuário; `reuse_results=False` força o recálculo. Etapas que precisam de dados ficam em `awaiting_input`; `resume(token, run_id, configuration)` confirma referências explícitas e modo. Resultados completos persistem entre inicializações. A CLI executa uma operação e não possui os popups do workspace.

`biomolexplorer-ui --language pt` inicia em português; `--language en` (padrão) inicia em inglês. A tela de login permite trocar por sessão. `ui/localization.py` aplica catálogos em `resources/i18n/` apenas à apresentação, preservando parâmetros, identificadores e dados editáveis. Consulte [manual](user_manual.md) e [validação](pipeline_validation.md).


## Seleção, retomada e idioma da interface

`PipelineService.submit(token, project_id, reuse_results=True)` reaproveita resultados compatíveis. A interface consulta `existing_results` e pede a decisão ao usuário; `reuse_results=False` força o recálculo. Etapas que precisam de dados ficam em `awaiting_input`; `resume(token, run_id, configuration)` confirma referências explícitas e modo. Resultados completos persistem entre inicializações. A CLI executa uma operação e não possui os popups do workspace.

`biomolexplorer-ui --language pt` inicia em português; `--language en` (padrão) inicia em inglês. A tela de login permite trocar por sessão. `ui/localization.py` aplica catálogos em `resources/i18n/` apenas à apresentação, preservando parâmetros, identificadores e dados editáveis. Consulte [manual](user_manual.md) e [validação](pipeline_validation.md).
