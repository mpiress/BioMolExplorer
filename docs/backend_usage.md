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
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
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

A forma equivalente é `PYTHONPATH=src python -m biomolexplorer ...`. A CLI executa a operação
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
    job_timeout=None,
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
| `retrieve_pubchem` | Referência manual SMILES/CID/nome ou compostos, limiar e máximo | Tabela de similares e relatório por referência |
| `expand_similar_compounds` | `search_term`, `base_input_path` contendo `ChEMBL/molecules` e `ChEMBL/similars` | Compostos novos, relações e conjunto consolidado |
| `retrieve_structures` | Texto/IDs/UniProt/ligantes ou filtros; `target` opcional | PDBs, `pdb_codes.csv`, relatório |
| `retrieve_zinc` | `base_input_path`, `filename` do arquivo de tranches | `compounds.csv`, MOL2 individuais e relatório |
| `admet` | `base_input_path`, opcional `input_file` | CSVs de avaliação, subconjuntos e figura |
| `fingerprints` | `base_input_path`, algoritmos, `chunk_size` | CSVs de fingerprints |
| `similarity` | `base_input_path` dos fingerprints, métrica, limiar e `approximate` | CSVs de arestas na subpasta `Similarity` |
| `graphs` | `similarity_path`; `base_input_path` opcional dos compostos; opções MCS | Grafos/MCC independentes, fragmentos e CSVs |
| `prepare_structures` | Receptor PDB bruto ou preparado, compostos, `target`, registros PDB e opções por papel | Receptores, centros, `pdb_codes.csv` e candidatos nos formatos selecionados |
| `redocking` | Diretório PDB, `target`, parâmetros de preparação | Cópias de trabalho, resultados Vina ou DOCK6 6.11 e RMSD |
| `docking_vina` | Complexos preparados, compostos selecionados, `mol_filename` | PDBQT e resultados Vina |
| `docking_dock6` | Complexos preparados, compostos, `base_vina_path`, `pdb_code`, instalação DOCK6 | Conformações, scores e footprints |
| `dynamics_dock6` | `base_input_path`, `base_selected_mols`, seleção de pares e DOCK6 6.11 | Amber Score, termos energéticos, protocolo, estruturas e análise inicial/final disponível |
| `consensus` | Diretório contendo `Vina` e `Dock6`, ou entradas separadas `base_vina_path` e `base_dock6_path`, `target` | CSV e figura de consenso |

Enums aceitam nomes ou valores, como `TanimotoSimilarity`/`Tanimoto` e
`Morgan`/`morgan`. Listas de filtros PDB também usam strings serializáveis.
O parâmetro Vina `pdb_code` segue a função existente: uma lista de registros
`[PDB_CODE, LIGAND, RESNUM, CHAIN]`; DOCK6 recebe um único registro desses quatro
campos. Redocking exige `pdb_codes` explicitamente selecionados: uma lista de
`[PDB_CODE, LIGAND, RESNUM, CHAIN]`, com resolução opcional no quinto elemento,
complementada pelos metadados de `pdb_codes.csv`. `preparation_options` define o preparo comum de todos os receptores e ligantes e tem precedência. Para compatibilidade, `preparation_pairs` guarda as
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
`graphs` não recalcula nem filtra os pesos. `mcs_timeout` tem padrão 0 (sem limite de tempo),
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

Escolha `redocking_engine="vina"` (padrão) ou `redocking_engine="dock6"`. DOCK6 usa `dock6_app_path` da instalação 6.11, `conformer_search_type`, `density`, `radius`, `distance`, `max_orientations`, `num_scored_conformers` e `dock6_seed`. Com entradas preparadas, exige receptor `.dockprep.mol2` e `.noH.pdb`, ligante `.lig.mol2` e centros. Vina usa PDBQT. O exemplo legado a seguir continua válido; clientes novos podem substituir `preparation_pairs` por um único `preparation_options={"receptor": {...}, "ligand": {...}, "cofactors": []}`.

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

Consulte [Logs e diagnóstico](logging.md) para o formato comum, contexto por execução, códigos de falha, resumo do job e o comando `PYTHONPATH=src python -m biomolexplorer.log_report`.


## Contrato interoperável de docking

**Docking com Vina** e **Docking com DOCK6** recebem compostos de qualquer bloco que publique moléculas: ChEMBL, PubChem, ZINC, ADMET, fingerprints, grafos, consenso ou importação do usuário. Selecione a tabela na entrada de compostos e conecte separadamente o receptor preparado. Não é obrigatório passar por grafos, ADMET ou pelo outro motor.

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

`requires_curation` dispensa uma nova confirmação de entradas quando os arquivos necessários estão configurados e todas as referências identificam arquivos (`asset` ou `selector` explícito). Referências `stage` com `selector=auto` continuam exigindo seleção. A configuração é revalidada antes de executar. Conexões incompatíveis são recusadas com a mesma mensagem padrão no canvas e em `validate_pipeline`.

## Duração das execuções

Os cálculos científicos e as etapas do pipeline não possuem limite automático de duração. Chimera, Open Babel, Vina, DOCK6 e os demais comandos aguardam conclusão ou erro. O cancelamento do usuário encerra o grupo de processos, incluindo comandos filhos. `job_timeout=None` é o padrão; um limite positivo só é aplicado se o chamador da API o definir explicitamente. `run_command` também aceita um timeout explícito para usos específicos, mas não usa mais `BIOMOL_COMMAND_TIMEOUT`. Na busca MCS, 0 permite concluir a busca; valores positivos são um orçamento opcional definido no bloco. Blocos salvos anteriormente podem conservar esse orçamento. Timeouts de comunicação HTTP continuam específicos por requisição, sem limitar a duração total da etapa.

## Configuração ADMET e proveniência

Os parâmetros de pesquisa são `admet_profile`, `admet_mode`, `prediction_backend`, `bbb_target`, `hia_target`, `rule_set`, `max_rule_violations`, `alert_catalogs`, `alert_policy` e `property_limits`. A ausência de perfil mantém `legacy`; declare o perfil para uma análise nova. Exemplo de parâmetros de seleção, a combinar com as entradas da operação:

```json
{
  "admet_profile": "custom",
  "admet_mode": "filter",
  "prediction_backend": "rdkit",
  "bbb_target": "BBB-",
  "hia_target": "HIA+",
  "rule_set": "oral",
  "max_rule_violations": 1,
  "alert_policy": "annotate",
  "alert_catalogs": ["PAINS", "BRENK"],
  "property_limits": {"ESOL_logS": [-4, null], "QED": [0.3, null]}
}
```

Os limites são inclusivos e combinados com E; `null` deixa um lado aberto. Os cortes do exemplo são escolhas de pesquisa. Chaves `ML_` só são aceitas com backend `admet_ai` e devem corresponder às saídas disponíveis. Valores ausentes não satisfazem critérios exigidos. Configuração, versões, métodos e contagens ficam em `_admet_report.json`; tabelas avaliadas/selecionadas/rejeitadas preservam motivos. Veja [guia ADMET](admet_research.md) e [instalação do preditor](installation.md#preditor-admet-ai-opcional).

## Docking em lotes e relatório de compostos não utilizados

Receptores preparados são validados e reutilizados. A API Python Vina reutiliza mapas por receptor/caixa/configuração dentro do lote; templates não suportados ou API ausente usam CLI, registrado no log. DOCK6 compartilha superfície, esferas e grid por receptor/sítio, valida os MOL2 e executa minimização/docking/footprint por composto. Falhas individuais excluem somente o candidato; falhas globais interrompem a etapa. Veja [reutilização Vina](user_manual.md#reutilizacao-do-receptor-no-vina) e [validação DOCK6](user_manual.md#validacao-e-exclusoes-no-dock6).

`preparation_failures.csv` é o nome preservado do relatório também para falhas posteriores do DOCK6. Colunas: `compound_identifier`, `submitted_smiles`, `input_file`, `input_row`, `preparation_step`, `preparation_context`, `reason`. Campos `preparation_step/context` identificam a operação concreta, mesmo para docking/footprint. Os detalhes da operação incluem `preparation_failures`, contagem e caminho do relatório quando disponíveis. Não inclua esses compostos no denominador de resultados utilizáveis sem registrar as exclusões.

## Persistência e reutilização incremental

Os resultados publicados pertencem a execuções duráveis no SQLite. Cada lote
concluído registra configuração, entradas e hashes, manifesto de saídas, versão
científica e ambiente do worker. Lotes completos de execuções posteriormente
interrompidas, canceladas ou com falhas podem ser reutilizados; lotes incompletos,
ativos, removidos, modificados ou com protocolos incompatíveis não são aceitos.
Os arquivos são referenciados pela nova execução, sem copiar os resultados antigos.

As chaves consideram somente entradas realmente selecionadas. Mudanças em um
conjunto não invalidam outros conjuntos independentes. Rótulos de apresentação de
grafos não fazem parte da identidade científica. Operações de merge e métricas
que dependem do conjunto completo são recalculadas quando esse conjunto muda;
não são aproximadas por resultados individuais incompatíveis.

Caches do projeto ficam em `.biomolexplorer/cache`: receptores preparados são
compartilhados entre redocking e docking quando PDB, cadeia, opções do receptor,
cofatores, pH, templates e ferramentas correspondem. A extração do complexo para
preparar o ligante continua própria de cada chamada, preservando seu inventário
original. Fingerprints também são reutilizados por molécula, algoritmo, opções e
versão do RDKit, em SQLite com payload comprimido e checksum. Desabilitar a
reutilização ignora esses caches computacionais.

`ResultFiles.remove_project_results(token, pid)` remove resultados e caches,
preservando todas as configurações. `ResultFiles.clear_project_results(token, pid)`
restaura os parâmetros padrão e limpa seleções e templates, mantendo a topologia.
Ambas preservam originais da biblioteca e histórico e recusam execuções ativas.
A interface sempre confirma a exclusão antes de chamar esses métodos.

`WorkspaceStore.result_removal_plan` compara requisições científicas normalizadas e
hashes de conteúdo persistidos, prevendo arquivos e blocos afetados. Alterações
destrutivas são recusadas pelo backend até a confirmação. Após confirmar, salve com
`save_pipeline(..., prune_results=True, confirmation=proposal['confirmation'])`.
A confirmação é vinculada ao estado previsto e revalidada sob o lock SQLite.
Mudanças concorrentes exigem nova confirmação. Cancelar não salva o rascunho.
Uma fila de propagação e índices de consumidores por arquivo percorrem o DAG sem
varrer repetidamente todos os descendentes. Lotes independentes e arquivos
compartilhados são preservados; operações agregadas afetadas são invalidadas.
Movimentos, metadados e auditoria são revertidos se a transação falhar.

Não existe uma ação manual de sincronização: salvamento, `submit` e `resume`
verificam a equivalência. Curadorias confirmadas também são persistidas no projeto,
evitando que uma escolha feita em execução pareça uma alteração posterior.
Requisições normalizadas e hashes por assinatura de arquivo possuem caches
limitados; o hash de conteúdo é refeito quando inode, tamanho, mtime ou ctime mudam.
Fingerprints possuem unidades independentes por entrada e algoritmo. A preparação
conjunta mantém os dois formatos na mesma conformação nativa. Índices parciais de redocking preservam pares válidos
e receptores compartilhados, sem ressuscitar pares removidos. Resultados antigos
combinando algoritmos podem fornecer componentes verificáveis sem repetir jobs.

Resultados antigos permanecem legíveis. Quando a versão científica ou a
proveniência antiga não puderem ser verificadas, pode ser necessário recalcular o
artefato; ele não é apagado por essa invalidação. A versão desconhecida de um
produtor antigo não é apresentada como se fosse a versão atual.

O cache versão 1 possui migração conservadora para builds e ambientes científicos
explicitamente auditados em `resources/legacy_cache_compatibility.json`. A migração
confere a requisição e os hashes originais, preserva a proveniência e pode aproveitar
lotes e pares anteriores mesmo quando novas entradas são adicionadas. A compatibilidade
não é presumida para versões desconhecidas nem para pesos antigos do ADMET-AI.

As exclusões de resultados usam um journal durável de movimentações, coordenado com
a transação SQLite. Na abertura do workspace, uma interrupção antes do commit restaura
os arquivos; após o commit, conclui a limpeza. Recuperações conflitantes preservam os
arquivos e registram o problema, sem sobrescrever dados existentes.

## Referência exata dos parâmetros ativos

Inventário extraído do registro atual, sem executar cálculos. Obrigatório significa exigido pela API; a interface pode preencher esse campo. Opcional não significa válido para toda combinação: os validadores verificam tipos, domínios e dependências. O CLI expõe 15 operações científicas; `import_results` é o 16º tipo de bloco, gerenciado pelo workspace e não pelo CLI. Nomes abaixo são chaves da API, não rótulos traduzidos da interface.

### API dynamics_dock6

`wrappers.dynamics.perform_dynamics_dock6`

| Obrigatórios | Opcionais |
| --- | --- |
| `base_input_path`, `base_selected_mols` | `dock6_app_path`, `md_inputs`, `md_seed_path`, `selected_pairs`, `excluded_pairs`, `amber_score_type`, `amber_score_movable_region`, `amber_score_minimization_rmsgrad`, `amber_score_before_md_minimization_cycles`, `amber_score_md_steps`, `amber_score_after_md_minimization_cycles`, `amber_score_gb_model`, `amber_score_nonbonded_cutoff`, `amber_score_temperature`, `amber_score_movable_distance_cutoff`, `amber_score_receptor_movable_atom_expr`, `amber_score_complex_movable_atom_expr`, `amber_score_abort_on_unprepped_ligand`, `use_existing_ligand_charges` |

### API retrieve_pubchem

`biomolexplorer.pubchem_retrieval.retrieve_pubchem`

| Obrigatórios | Opcionais |
| --- | --- |
| — | `target`, `reference_source`, `reference_type`, `reference`, `selection_mode`, `compound_id`, `threshold`, `max_records`, `base_input_path`, `input_file` |

### API retrieve_compounds

`wrappers.crawlers.retrieve_compounds`

| Obrigatórios | Opcionais |
| --- | --- |
| `search_term` | `include_pubchem`, `pubchem_threshold`, `pubchem_max_records`, `chembl_filters`, `search_mode`, `max_targets`, `max_records`, `expand_chembl`, `similarity_threshold` |

### API expand_similar_compounds

`wrappers.crawlers.expand_similar_compounds`

| Obrigatórios | Opcionais |
| --- | --- |
| `search_term` | `threshold`, `max_records`, `base_input_path`, `input_file` |

### API retrieve_structures

`wrappers.crawlers.load_pdb`

| Obrigatórios | Opcionais |
| --- | --- |
| — | `target`, `pdb_query`, `pdb_ids`, `uniprot_ids`, `ligand_ids`, `max_records`, `pdb_ec`, `organism`, `PolymerEntityTypeID`, `ExperimentalMethodID`, `max_resolution`, `must_have_ligand`, `ligand_selection`, `excluded_pdb_ids` |

### API retrieve_zinc

`wrappers.crawlers.load_zinc`

| Obrigatórios | Opcionais |
| --- | --- |
| `base_input_path` | `filename`, `verbose`, `download_workers` |

### API prepare_structures

`wrappers.redocking.prepare_structures`

| Obrigatórios | Opcionais |
| --- | --- |
| `base_input_path`, `target` | `pdb_codes`, `pH`, `charge_type`, `preparation_options`, `base_selected_mols`, `mol_filename`, `receptor_prepared`, `docking_engines` |

### API admet

`wrappers.admet.ADMETWrapper`

| Obrigatórios | Opcionais |
| --- | --- |
| `base_input_path` | `input_file`, `verbose`, `admet_profile`, `admet_mode`, `prediction_backend`, `bbb_target`, `hia_target`, `rule_set`, `max_rule_violations`, `alert_policy`, `alert_catalogs`, `property_limits`, `dataset_origin` |

### API fingerprints

`wrappers.molecular_analyzer.generate_fingerprints`

| Obrigatórios | Opcionais |
| --- | --- |
| `base_input_path` | `morgan_n_bits`, `radius`, `files`, `morgan`, `maccs`, `pharmacophore`, `chunk_size` |

### API similarity

`wrappers.molecular_analyzer.compute_similarity`

| Obrigatórios | Opcionais |
| --- | --- |
| `base_input_path` | `metric`, `fingerprint`, `filename`, `threshold`, `approximate` |

### API graphs

`wrappers.molecular_analyzer.analyze_graphs`

| Obrigatórios | Opcionais |
| --- | --- |
| — | `similarity_path`, `base_input_path`, `mcs_timeout`, `mcs_ring_matches_ring_only`, `mcs_complete_rings_only`, `graph_inputs` |

### API redocking

`wrappers.redocking.perform_redocking`

| Obrigatórios | Opcionais |
| --- | --- |
| `base_input_path`, `target` | `pdb_codes`, `pH`, `sizeof_box`, `exhaustiveness`, `num_modes`, `prepare_complex`, `charge_type`, `preparation_pairs`, `preparation_options`, `redocking_engine`, `dock6_app_path`, `density`, `radius`, `distance`, `conformer_search_type`, `max_orientations`, `num_scored_conformers`, `dock6_seed` |

### API docking_vina

`wrappers.docking.perform_docking_vina`

| Obrigatórios | Opcionais |
| --- | --- |
| `base_input_path`, `target`, `base_selected_mols`, `mol_filename` | `pdb_code`, `pH`, `sizeof_box`, `exhaustiveness`, `num_modes`, `preparation_options`, `receptor_prepared` |

### API docking_dock6

`wrappers.docking.perform_docking_dock6`

| Obrigatórios | Opcionais |
| --- | --- |
| `base_input_path`, `target`, `base_selected_mols`, `dock6_app_path`, `charge_type`, `mol_filename` | `pdb_code`, `base_vina_path`, `density`, `radius`, `distance`, `conformer_search_type`, `plot_max_residues`, `pH`, `preparation_options`, `receptor_prepared` |

### API consensus

`wrappers.docking.generate_consensus`

| Obrigatórios | Opcionais |
| --- | --- |
| `base_input_path`, `target` | `repulsion_weight`, `base_vina_path`, `base_dock6_path` |

## Amber Score nativo: requisitos, configuração e resultados

`dynamics_dock6` exige receptor estrutural e ligante com coordenadas 3D compatíveis; uma tabela apenas com SMILES não substitui uma pose. A seleção inicia vazia ou com escolhas salvas. Separe receptores e ligantes nas tabelas próprias, confira os pares e execute somente os escolhidos. Redocking fornece exclusivamente receptores aos estágios posteriores, não seus ligantes cristalográficos.

O verificador exige DOCK6 que se identifique como **6.11**, `bin/` e `parameters/leap`, além de `prepare_amber.pl`, `amberize_receptor`, `amberize_ligand`, `amberize_complex`, `tleap`, `teLeap`, `antechamber` e `parmchk` na distribuição DOCK6. Esses utilitários AMBER internos não constituem um ambiente AmberTools externo configurável. A ausência de um componente necessário bloqueia o protocolo; não há substituição por outro score.

| Parâmetro | Padrão implementado | Interpretação |
| --- | --- | --- |
| `amber_score_type` | `dock6` | Único método ativo |
| `amber_score_movable_region` | `ligand` | ligand, everything, nothing, distance ou nab_atom_expression |
| `amber_score_minimization_rmsgrad` | 0.01 | Critério RMS do gradiente, positivo |
| `amber_score_before_md_minimization_cycles` | 100 | Ciclos antes do refinamento |
| `amber_score_md_steps` | 3000 | Passos internos; passo fixo 0.001 ps, não uma trajetória de equilíbrio |
| `amber_score_after_md_minimization_cycles` | 100 | Ciclos após o refinamento |
| `amber_score_gb_model` | 5 | 1, 2 ou 5; solvente implícito |
| `amber_score_nonbonded_cutoff` | 18.0 | Corte em Å |
| `amber_score_temperature` | 300.0 | Temperatura em K |
| `amber_score_movable_distance_cutoff` | 10.0 | Seleção móvel por distância em Å |
| `amber_score_receptor_movable_atom_expr` | `::` | Expressão NAB do receptor, não máscara cpptraj |
| `amber_score_complex_movable_atom_expr` | `::` | Expressão NAB do complexo |
| `amber_score_abort_on_unprepped_ligand` | true | Não pontuar ligante sem preparo AMBER válido |
| `use_existing_ligand_charges` | false | true preserva cargas MOL2; false usa a preparação AM1-BCC do bundle DOCK6 |

O resultado contém score e contribuições de complexo/receptor/ligante, protocolo, estruturas e arquivos de diagnóstico. Havendo estrutura final compatível, a análise compara RMSD, raio de giro e contagem de hidrogênios inicial/final. Contagem de hidrogênios não mede ligações H. Não se oferece RMSF ou ocupância temporal sem trajetória; o modo `nothing` pode não fornecer estrutura final. Não há controle de caixa explícita, pressão ou protocolo convencional de equilíbrio. Defeitos de tipagem, metais/cofatores não parametrizados e inconsistências de carga devem ser tratados explicitamente; a documentação não garante que todo par seja viável.
