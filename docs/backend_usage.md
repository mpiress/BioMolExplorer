# Uso da camada de aplicação

[Documentação](README.md) · Português · [English](en/backend_usage.md)

Para utilizar a interface pronta com autenticação e projetos, consulte o
[guia do workspace Flet](frontend.md). Este documento descreve a CLI e os serviços
Python utilizados pela interface.

Antes de executar a CLI ou iniciar a interface, siga o [guia de instalação e configuração](installation.md): baixe o [repositório do GitHub](https://github.com/mpiress/BioMolExplorer) e instale **UCSF Chimera 1.17 e DOCK6 6.11** no computador dos cálculos. Configure seus executáveis no `PATH`. Sem essas ferramentas, as etapas que dependem delas falharão.

Crie/ative o ambiente Conda e instale o projeto a partir da raiz do código baixado:

```bash
conda env create -f requirements.yml
conda activate BioMolExplorer
python -m pip install -e . --no-deps --no-build-isolation
```

As dependências científicas são gerenciadas pelos arquivos Conda. O pacote
Python fornece a camada de aplicação, os módulos científicos e os recursos;
a instalação do pacote por si só não instala Chimera 1.17 ou DOCK6 6.11.
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
| `retrieve_structures` | Texto/IDs/UniProt/ligantes ou filtros; `target` opcional | PDBs, `pdb_codes.csv`, relatório |
| `retrieve_zinc` | `base_input_path`, `filename` do arquivo de tranches | `compounds.csv`, MOL2 individuais e relatório |
| `admet` | `base_input_path`, opcional `input_file` | CSVs de avaliação, subconjuntos e figura |
| `fingerprints` | `base_input_path`, algoritmos, `chunk_size` | CSVs de fingerprints |
| `similarity` | `base_input_path` dos fingerprints, métrica, limiar e `approximate` | CSVs de arestas na subpasta `Similarity` |
| `graphs` | `similarity_path`; `base_input_path` opcional dos compostos; opções MCS | Grafos/MCC independentes, fragmentos e CSVs |
| `prepare_structures` | Receptor PDB bruto ou preparado, compostos, `target`, registros PDB e opções por papel | Receptores, centros, `pdb_codes.csv` e candidatos nos formatos selecionados |
| `redocking` | Diretório PDB, `target`, parâmetros de preparação | Cópias de trabalho, resultados Vina e RMSD |
| `docking_vina` | Complexos preparados, compostos selecionados, `mol_filename` | PDBQT e resultados Vina |
| `docking_dock6` | Complexos preparados, compostos, `base_vina_path`, `pdb_code`, instalação DOCK6 | Conformações, scores e footprints |
| `consensus` | Diretório contendo `Vina` e `Dock6`, ou entradas separadas `base_vina_path` e `base_dock6_path`, `target` | CSV e figura de consenso |

Enums aceitam nomes ou valores, como `TanimotoSimilarity`/`Tanimoto` e
`Morgan`/`morgan`. Listas de filtros PDB também usam strings serializáveis.
O parâmetro Vina `pdb_code` segue a função existente: uma lista de registros
`[PDB_CODE, LIGAND, RESNUM, CHAIN]`; DOCK6 recebe um único registro desses quatro
campos. Redocking exige `pdb_codes` explicitamente selecionados: uma lista de
`[PDB_CODE, LIGAND, RESNUM, CHAIN]`, com resolução opcional no quinto elemento,
complementada pelos metadados de `pdb_codes.csv`. `preparation_pairs` guarda as
opções por chave `PDB|LIGAND|RESNUM|CHAIN`; receptor e ligante usam a mesma cadeia.
O parâmetro legado `verbose` é ignorado; a verbosidade permanece em zero.
Consulte [configuração do redocking](redocking_configuration.md).

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

## Exemplo de redocking

O exemplo exige uma entrada existente em `/tmp/study/PDB/Estruturas/4M0E.pdb`, contendo o resíduo 1YL / 604 na cadeia A. Adapte os caminhos e o par à sua coleção; esta chamada não recupera o PDB.

```python
from biomolexplorer.operations import execute_operation

parameters = {
    "base_input_path": "/tmp/study/PDB",
    "target": "Estruturas",
    "pdb_codes": [["4M0E", "1YL", 604, "A", 2.0]],
    "prepare_complex": True,
    "pH": 7.4,
    "sizeof_box": [24, 24, 24],
    "exhaustiveness": 20,
    "num_modes": 10,
    "charge_type": "gas",
    "preparation_pairs": {
        "4M0E|1YL|604|A": {
            "cofactors": [],
            "receptor": {
                "remove_solvent": True,
                "remove_hydrogens": True,
                "add_hydrogens": True,
                "minimize": True,
                "charge_type": "gas",
            },
            "ligand": {
                "remove_solvent": True,
                "remove_hydrogens": True,
                "add_hydrogens": True,
                "minimize": True,
                "charge_type": "gas",
            },
        },
    },
}
result = execute_operation("redocking", parameters, "/tmp/study/redocking-output")
```

Sem opções por par, remoção de solvente e de hidrogênios existentes, adição de hidrogênios e minimização usam `True`; cargas do receptor usam `gas` e as do ligante herdam `charge_type` do bloco quando omitidas. Não há cofator padrão. Ao usar `prepare_complex=False`, mantenha os receptores `<PDB>_<CHAIN>.dockprep.pdbqt`, os ligantes de referência `<PDB>_<LIGAND>_<RESNUM><CHAIN>.lig.pdbqt` e três coordenadas finitas por complexo em `Prepared/centers.csv`.

`execute_operation` prepara uma cópia em `structures/<target>/`, preservando a entrada original. Os metadados com `RMSD` ficam em `structures/<target>/pdb_codes.csv`, os preparados em `structures/<target>/Prepared/` e as poses em `<target>/` dentro da saída. `result.artifacts` lista os arquivos exportáveis. As chamadas diretas a `wrappers.redocking.perform_redocking` não têm o isolamento fornecido pela camada de aplicação. A CLI recebe o mesmo dicionário por `--parameters`, sem as janelas de consulta da interface.

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


## Parâmetros de recuperação flexível

Consulte [recuperação da informação](retrieval.md) para os modos ChEMBL `search_mode`, exemplos sem alvo/EC obrigatório, limites, filtros e relatórios. Os nomes das operações e contratos CSV das etapas seguintes permanecem disponíveis.

Consulte [Logs e diagnóstico](logging.md) para o formato comum, contexto por execução, códigos de falha, resumo do job e o comando `python -m biomolexplorer.log_report`.


## Contrato interoperável de docking

**Docking com Vina** e **Docking com DOCK6** recebem compostos de qualquer bloco que publique moléculas: ChEMBL, PubChem, ZINC, ADMET, fingerprints, grafos, consenso ou importação do usuário. Selecione a tabela na entrada de compostos e conecte separadamente o receptor preparado. Não é obrigatório passar por grafos, ADMET ou pelo outro motor. PubMed fornece referências bibliográficas; moléculas obtidas dessas referências devem ser importadas com código e SMILES ou estrutura válida.

Para reutilizar conformações, conecte **Vina → DOCK6** ou **DOCK6 → Vina** na entrada de compostos, escolhendo `docking_results.csv` ou a pose desejada. A identidade e o SMILES são preservados; a conversão utiliza a conformação selecionada, sem gerar outra a partir do SMILES. Sem uma pose de entrada, a preparação gera um conformero 3D. DOCK6 usa o centro do sítio do receptor para selecionar esferas e posicionar ligantes recém-gerados; a entrada opcional de poses Vina continua disponível para refinamento de candidatos correspondentes.

Confira a caixa, pH e esforço do Vina e as cargas, superfície, distância, raio, busca flexível/rígida e footprint do DOCK6. Ambos exigem receptores preparados, metadados e centros; DOCK6 também exige MOL2 do receptor e PDB sem hidrogênios. Seus resultados apresentam código molecular, receptor, score, SMILES, **3D** da pose calculada e **Remover**. O botão abre a estrutura obtida no docking.

**Consenso de docking** reúne os resultados selecionados de cada motor, incluindo vários lotes e ramificações. Calcula somente a interseção por código molecular e receptor, usando o menor score quando há poses repetidas. Resultados do mesmo código com estruturas diferentes são recusados. Não havendo interseção, o bloco é marcado como não executado e explica que os motores não avaliaram compostos em comum para o mesmo receptor.

A tabela de consenso apresenta os scores Vina e DOCK6, SMILES, receptor e normalizações z-score/min-max, com botões **3D Vina**, **3D DOCK6** e **Remover** em cada linha. O score DOCK6 usado no consenso é `min(0, Grid_Score + repulsion_weight × Internal_energy_repulsive)`; a repulsão ausente vale zero. Para um composto ou scores constantes, as normalizações são zero. A remoção exige permissão de edição, fica registrada e afeta próximas entradas; não recalcula análises já concluídas. As tabelas principais de cada lote são usadas sem reintroduzir linhas removidas a partir das tabelas auxiliares.

As saídas nativas usam `docking_results.csv`: `molecule_chembl_id,canonical_smiles,receptor_id,engine,score,conformer_file`. Caminhos de pose são relativos à tabela. Importe também os arquivos de conformação. No backend direto, `base_selected_mols` é a pasta da tabela e `mol_filename` é seu nome sem `.csv`; `base_vina_path` do DOCK6 é opcional. O consenso recebe pastas em `base_vina_path` e `base_dock6_path`, exporta as duas poses em `poses/` e retorna `skipped_reason` quando não há interseção.


A operação `prepare_structures` aceita `preparation_options` com `receptor`, `ligand` e `cofactors`. Cada papel aceita `remove_solvent`, `remove_hydrogens`, `add_hydrogens`, `minimize` (booleanos) e `charge_type` (`gas` ou `am1`). O wrapper aplica essas opções aos registros explícitos ou inferidos de `pdb_codes.csv`. O pH permanece comum aos dois papéis. Chamadas legadas sem esse parâmetro conservam a preparação anterior; a interface usa as opções por papel e migra configurações antigas sem alterar templates personalizados.


`prepare_structures` também recebe `base_selected_mols`, `mol_filename` (padrão `compounds`), `receptor_prepared` e `docking_engines` (`vina`, `dock6` ou `both`). O pipeline identifica automaticamente receptores preparados nas entradas do redocking. Receptores prontos são copiados com seus arquivos complementares e centros, sem repetir o preparo; PDBs brutos seguem o preparo usado pelo redocking. Os candidatos externos são preparados em `Target/Compounds/compounds.csv`, com códigos e SMILES preservados e colunas `prepared_pdbqt` e/ou `prepared_mol2`. Vina e DOCK6 reutilizam esses arquivos sem outra minimização. O manifesto informa os formatos disponíveis; uma ferramenta cujo formato não foi exportado não pode usar essa saída. O CSV e todos os arquivos referenciados devem acompanhar a entrada. Chamadas legadas sem `base_selected_mols` continuam preparando somente estruturas.

Use um bloco separado para cada modo de receptor: não combine PDBs brutos com receptores preparados no mesmo bloco. O pH continua disponível para o preparo dos candidatos quando o receptor é reutilizado. Antes de preparar os compostos, o bloco verifica os centros do sítio (três coordenadas finitas) e os arquivos do receptor exigidos pela saída escolhida. Vina requer `.dockprep.pdbqt`; DOCK6 requer também `.dockprep.mol2` e `.noH.pdb`. O PDBQT permanece como arquivo de seleção do receptor nos dois casos. Arquivos ausentes interrompem o processo com o nome do arquivo necessário, sem refazer o preparo do receptor.

## Importação por arquivo e tranches ZINC

A etapa `import_results` aceita `asset_types`, um mapa `{asset_id: tipo}` cobrindo exatamente os `asset_ids` selecionados. Um bloco pode reunir CSVs de compostos, estruturas e listas ZINC, publicando seus respectivos tipos. Sem `asset_types`, o parâmetro legado `kind` continua valendo para todos os arquivos. O pipeline revalida arquivos e conjuntos antes da execução; receptores preparados continuam exigindo os metadados. Remover um arquivo da tabela do formulário altera apenas a seleção do bloco.

`download_workers` configura os downloads paralelos de `retrieve_zinc` e `wrappers.crawlers.load_zinc`: inteiro de 1 a 16, padrão 4; 1 executa sequencialmente. Cada tarefa usa sua própria sessão HTTP. A fila de antecipação é limitada ao número de threads; a normalização molecular ocorre na ordem da lista. A configuração é validada antes dos downloads e o relatório registra `download_workers` efetivo, limitado ao número de links.

`retrieve_zinc` chama `wrappers.crawlers.load_zinc`. `base_input_path` indica a pasta da lista e `filename` seu nome (padrão `zinc_urls.txt`); na interface, o pipeline resolve ambos a partir dos arquivos selecionados. A saída tem nome estável `compounds.csv`, independentemente do nome da lista. MOL2s individuais ficam em `Conformers/`, com caminhos autorizados no manifesto. `conformer_origin=library` identifica conformações que ainda precisam de posicionamento no sítio; o bloco de preparo preserva essa origem em `prepared_origin`.

```python
from biomolexplorer.operations import execute_operation

execute_operation('retrieve_zinc', {
    'base_input_path': '/path/to/lists',
    'filename': 'zinc-download.uri',
}, '/path/to/output')
```

Consulte [tranches 2D e 3D](retrieval.md#tranches-zinc-2d-e-3d) para formatos e organização dos artefatos.

## Composto específico e entradas pré-configuradas do docking

Referências de `base_selected_mols` em `docking_vina` e `docking_dock6` aceitam `compound_id` opcional, junto de `stage`/`asset` e `selector`. Sem essa chave, todos os compostos do arquivo participam. O filtro ocorre por referência antes da deduplicação e preserva conformações e arquivos preparados. Um código ausente interrompe a etapa com mensagem de validação. A seleção integra a configuração, o cache e a proveniência dos arquivos materializados.

```json
{"stage": "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", "selector": "compounds.csv", "compound_id": "CHEMBL1"}
```

`requires_curation` dispensa uma nova confirmação de Vina/DOCK6 quando receptor e compostos estão configurados e todas as referências identificam arquivos (`asset` ou `selector` explícito). Referências `stage` com `selector=auto` continuam exigindo seleção. A configuração é revalidada antes de executar. Conexões incompatíveis são recusadas com a mesma mensagem padrão no canvas e em `validate_pipeline`.
