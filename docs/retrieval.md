# Recuperação flexível de estruturas e compostos

[Documentação](README.md) · Português · [English](en/retrieval.md)

As etapas de recuperação podem iniciar ramificações independentes. Escolha primeiro o que deseja obter: estruturas PDB, compostos ligados a evidências de bioatividade, compostos encontrados diretamente ou expansão por similaridade. Revise os resultados antes de encaminhá-los ao próximo bloco.

## Buscar estruturas PDB

EC e nome da enzima não são obrigatórios. O nome da coleção é opcional e organiza a saída em `PDB/<coleção>/`; seu padrão é `Estruturas`. Informe pelo menos um critério real de pesquisa:

| Critério | Exemplo | Uso |
| --- | --- | --- |
| Texto | `acetylcholinesterase`, `EGFR` | Anotações textuais da RCSB, incluindo nomes e descrições |
| IDs PDB | `1CRN, 4HHB` | Seleção explícita de estruturas |
| Acessos UniProt | `P00533` | Estruturas associadas à proteína de referência |
| Códigos de ligantes | `ATP, ADP` | Estruturas que contêm os componentes químicos informados |
| EC | `3.1.1.7` | Classificação enzimática, quando conhecida |
| Organismo, polímero, método ou resolução | `Homo sapiens`, `Protein`, `2.5` Å | Pesquisa ou refinamento por atributos |

Critérios distintos são combinados com **E**. Valores dentro de uma lista são alternativas (**OU**). Assim, IDs PDB e organismo preenchidos juntos restringem a seleção à interseção; limpar campos amplia a busca. O texto livre segue a semântica da RCSB: termos separados podem ser alternativas; use aspas para pesquisar uma expressão. [Referência oficial de busca RCSB](https://search.rcsb.org/).

Para chamadas legadas sem outros critérios, um nome de coleção diferente do padrão também é usado como texto de busca.

A opção **Exigir ligante** permanece ativada por padrão. Desative-a para recuperar estruturas sem ligantes; elas continuam úteis para inspeção, mas precisam de uma referência apropriada para o redocking implementado. Componentes não aquosos detectados no PDB podem incluir íons/aditivos, portanto confira o ligante antes de usá-lo como referência farmacológica.

O máximo padrão é 100 estruturas. A busca é paginada e limitada, com timeout e repetição de falhas transitórias. Os downloads são concorrentes, com até quatro tarefas. Cada resposta deve conter coordenadas e ser lida pelo parser antes da publicação atômica do `.pdb`. Como antes, registros LINK/SSBOND são removidos para o fluxo de preparação existente.

Arquivos produzidos:

- `<PDB>.pdb`: estruturas baixadas com sucesso.
- `pdb_codes.csv`: contratos existentes `PDB_CODE,LIGAND,RESNUM,CHAIN,RESOLUTION`, com os ligantes detectados; estruturas sem ligantes não geram linhas de complexo.
- `retrieval_report.json`: filtros, consulta, total encontrado quando informado pela API, limite aplicado, IDs selecionados, downloads bem-sucedidos e erros individuais.

Uma falha individual não perde os downloads válidos. Se nenhum arquivo puder ser baixado, a etapa falha com orientação para o relatório. Uma pesquisa sem resultados também falha explicitamente. Entradas disponíveis somente em mmCIF não são convertidas automaticamente: o pipeline atual consome PDB e o relatório identifica a falha de obtenção desse formato.

## Revisar estruturas e ligantes na interface

A lista de resultados mostra os PDBs e oculta `pdb_codes.csv` e `retrieval_report.json`; os dois arquivos continuam armazenados e acompanham os contratos do pipeline. Cada estrutura oferece três ações:

- **Ligantes**: lista os registros detectados. Remova os indesejados ou adicione um código de ligante, número do resíduo e cadeia existentes no PDB baixado. Clique em **Salvar ligantes** para confirmar.
- **Visualizar estrutura 3D**: abre o explorador WebGL em uma aba do navegador, com rotação pelo mouse, zoom pela roda e deslocamento pelo botão direito. Há fitas, bastões, esferas e linhas, cores por cadeia/elemento/sequência, seleção de cadeia e modelos, ligantes e água opcionais, tela inteira e exportação PNG. O arquivo completo continua disponível para download.
- **Abrir no RCSB PDB**: abre a página oficial da estrutura pelo código PDB.

As ações também estão disponíveis na seleção dos PDBs antes do próximo bloco. Leitores consultam os registros; editores/proprietários podem alterá-los quando não há execução calculando ou na pausa de seleção de arquivos. Salvamentos entram no histórico, preservam os outros PDBs e atualizam a integridade dos resultados. A entrada modificada invalida o reaproveitamento dos consumidores afetados. Uma alteração concorrente exige reabrir a lista. Excluir todos os ligantes de uma estrutura deixa-a sem referência de complexo; para preparar/redockar, adicione uma referência válida ou deixe esse PDB fora da seleção. Uma nova recuperação pode detectar os candidatos novamente.

## Buscar na ChEMBL

Selecione **Buscar na ChEMBL por** antes de preencher a consulta:

| Modo (`search_mode`) | Consulta | Resultado |
| --- | --- | --- |
| `target` | Nome, CHEMBL220 ou P00533 | Identificação automática do alvo; preserva o fluxo anterior |
| `target_name` | Nome parcial | Alvos cujo nome contém o texto |
| `target_id` | CHEMBL220, CHEMBL240 | Um ou vários alvos explícitos |
| `uniprot` | P00533, P22303 | Alvos associados a acessos UniProt |
| `target_text` | Texto livre | Busca textual de alvos na API |
| `molecule_id` | CHEMBL25, CHEMBL50 | Compostos explícitos, sem exigir bioatividades |
| `molecule_name` | aspirin | Compostos pelo nome preferencial parcial |
| `similarity` | ID de composto ou SMILES | Compostos semelhantes à referência |
| `substructure` | c1ccccc1 | Compostos que contêm o fragmento SMILES |

IDs ChEMBL de alvo e de composto pertencem a recursos diferentes: escolha o modo correspondente. A identificação automática de UniProt valida o formato do acesso, evitando classificar qualquer texto com seis caracteres e um dígito como proteína.

Os modos de alvo executam **alvos → bioatividades → moléculas**, com expansão ChEMBL opcional. Os modos diretos consultam moléculas sem executar bioatividades. Similaridade/subestrutura não demonstram atividade contra um alvo; o relatório marca explicitamente a ausência dessa evidência. As rotas correspondem aos recursos documentados da [API ChEMBL](https://www.ebi.ac.uk/chembl/api/data/docs).

A interface mostra campos e filtros relevantes ao modo selecionado. **Ampliar com similares ChEMBL** controla a expansão após a busca por alvos. **Ampliar a seleção com PubChem** é independente e funciona também após buscas diretas. No modo `similarity`, o limiar ChEMBL padrão é 70%; a expansão de bioatividades mantém o limiar definido em seus filtros, inicialmente 65%.

## Filtros, limites e evidência

As novas configurações não restringem silenciosamente a alvos humanos, produtos naturais ou exclusão de produtos naturais. Organismo e tipo de alvo podem ficar vazios. **Produto natural → Qualquer** remove a restrição; outros campos opcionais também podem ser limpos. Configurações com overrides salvos continuam usando seus filtros. Configurações antigas sem overrides passam a usar os recursos padrão atualizados, que a interface exibe para revisão.

Os filtros padrão de bioatividade permanecem Ki/IC50, nM, ensaio B, pChEMBL preenchido e valor máximo 5000. Eles são editáveis. Para não aplicar um limite de valor, limpe **Atividade máxima**; para buscar sem unidade fixa, limpe também **Unidade de atividade**. Um limite numérico exige unidade padrão explícita, evitando comparar valores em unidades incompatíveis. Registros sem valor numérico podem permanecer quando não há limite. Relações (`=`, `<`, `>`) e outras evidências continuam nos dados exportados; não trate todas as medidas como valores exatos intercambiáveis.

O tipo de molécula pertence aos filtros de moléculas, não aos de alvo/bioatividade. Os grupos backend `chembl_filters.target`, `bioactivity`, `molecules` e `similars` substituem o respectivo grupo padrão quando fornecidos. `{}` remove as restrições daquele grupo. Na interface, não é necessário editar JSON para remover os filtros disponíveis.

Limites iniciais: 25 alvos; 1000 atividades por alvo; 1000 registros por consulta direta ou por referência de expansão ChEMBL. Os limites restringem os registros retornados antes de alguns filtros locais; um conjunto final menor não significa que todo o banco foi examinado. Aumente os limites ou refine a busca conforme o objetivo. Expansões por múltiplas referências podem gerar muitas consultas e muitos compostos; ambas podem ser desativadas.

Alvos/bioatividades são consultados novamente para evitar reutilizar filtros antigos. O fluxo mantém o cache de registros moleculares e o cache PubChem. Cada tarefa do sistema usa uma saída isolada; ao usar wrappers diretamente, utilize uma pasta nova para cada consulta/filtro, pois exportações antigas de moléculas podem continuar naquela pasta.

A saída consolidada continua em `compounds/<coleção>/compounds.csv`, com `molecule_chembl_id,canonical_smiles,molecule_properties,source`, compatível com ADMET, fingerprints e grafos. Estruturas são validadas e deduplicadas pelo contrato existente. Consultas SMILES ou textos com sintaxe especial recebem uma pasta `consulta_<hash>` para não transformar barras em caminhos. O texto original fica em `retrieval_report.json`.

## PubChem, ZINC e arquivos próprios

**Expandir similares** continua aceitando downloads ChEMBL ou uma tabela curada escolhida pelo usuário. A expansão PubChem conserva limiar/máximo por referência, cache, controle de frequência, repetição de falhas e exportação de relações. Compostos adicionais passam pela validação e deduplicação da consolidação. Eles também não constituem evidência de atividade no alvo.

**Recuperar ZINC** recebe listas de downloads TXT/URI e scripts exportados pelas tranches. Consulte [Tranches ZINC 2D e 3D](#tranches-zinc-2d-e-3d) para entradas, formatos e saídas padronizadas.

## Exemplos de backend

Parâmetros JSON para a CLI `biomolexplorer <operação> --parameters arquivo.json --output pasta`:

```json
{
  "pdb_ids": ["1CRN", "4HHB"],
  "must_have_ligand": false,
  "max_records": 10
}
```

Use esse arquivo com `retrieve_structures`; `target` e EC podem ser omitidos. Para compostos sem alvo:

```json
{
  "search_term": "CHEMBL25, CHEMBL50",
  "search_mode": "molecule_id",
  "max_records": 100,
  "include_pubchem": false
}
```

## Verificação realizada

Testes automatizados sem rede cobrem validação, modos, IDs, paginação limitada, timeouts, rotas estruturais codificadas, saída compatível, falha parcial e controles da interface. Consultas reais pequenas confirmaram os atributos RCSB de ID PDB/UniProt/ligante e os endpoints ChEMBL de nome, busca textual, UniProt, similaridade e subestrutura. Os wrappers também baixaram e leram 1CRN e produziram tabelas consolidadas por ID, nome, similaridade, subestrutura e alvo ChEMBL.

Esses testes confirmam contratos e casos de referência; disponibilidade e cobertura dos provedores variam. Uma busca limitada não é uma revisão exaustiva das bases.

```bash
PYTHONPATH=src:tests python -m unittest test_flexible_retrieval test_chembl_rest test_compound_retrieval test_guided_ui
```

## Seletores ChEMBL

**Organismo** oferece sugestões por nome científico, incluindo `Homo sapiens`, `Mus musculus`, `Rattus norvegicus`, `Escherichia coli`, `Saccharomyces cerevisiae` e `Danio rerio`. Há busca e digitação de outros nomes; não existe uma lista fechada de organismos. Deixe vazio ou selecione **Qualquer** para remover esse filtro. O exemplo e a orientação explicam que se trata do organismo do alvo.

**Tipo de ensaio** apresenta os códigos oficiais com descrição: B — Ligação, F — Funcional, A — ADMET, T — Toxicidade, P — Físico-químico e U — Não atribuído, além de **Qualquer**. Os códigos enviados à API permanecem em inglês, independentemente do idioma da interface. [Vocabulário oficial ChEMBL](https://chembl.gitbook.io/chembl-data-deposition-guide/file-structure/field-names-and-data-types-minimal-data-submission/assay.tsv).

**Medidas de atividade** utiliza a lista de 6.433 valores distintos de `standard_type` encontrados em consulta pública à ChEMBL em 05/10/2026. A lista é incluída no aplicativo e funciona sem consulta adicional ao abrir o formulário. As opções aparecem em até três colunas, com 24 itens por página e busca por nome. A seleção permanece ao pesquisar ou mudar de página; **Limpar seleção** aceita qualquer medida. **Outra medida** permite informar o nome exato de tipos novos ou específicos não presentes nessa versão da lista. Os valores científicos, inclusive maiúsculas/minúsculas e espaços, são preservados. [Fonte do catálogo público](https://www.ebi.ac.uk/chembl/elk/es/chembl_activity/_search).

O vocabulário de medidas cresce com a base. A lista representa a consulta registrada em `resources/crawlers/activity_types.json`, não uma enumeração imutável. Medidas como Inhibition, atividade percentual ou parâmetros cinéticos podem exigir remover o filtro nM, ajustar o limite numérico e permitir ausência de pChEMBL. Selecionar uma medida não modifica esses outros filtros automaticamente. A API e o pipeline continuam sujeitos aos limites de registros configurados.

Nomes que contêm vírgulas literais, como `K(p,uu,brain)`, são consultados com filtros exatos separados, preservando o limite total de atividades por alvo.

Na janela de ligantes, **Ligantes presentes na estrutura** sugere resíduos não aquosos do próprio arquivo e preenche código, resíduo e cadeia ao adicionar. Também é possível preencher uma linha manualmente.

## Visualizador PDB com interação 3D

Clique no ícone **Visualizar estrutura 3D** da estrutura desejada, nos resultados ou na seleção de arquivos. O visualizador abre em uma aba do navegador tanto no desktop quanto no modo web. Essa integração também atende Linux, onde o WebView do Flet não oferece suporte nativo.

- Arraste com o botão esquerdo para rotacionar a estrutura; use a roda do mouse para aproximar ou afastar.
- Arraste com o botão direito para deslocar a estrutura. **Recentrar** ajusta a câmera à seleção.
- Escolha a representação: fitas com ligantes, bastões, esferas ou linhas. Selecione a coloração por cadeia, elemento ou sequência.
- Filtre uma cadeia ou selecione outro modelo quando o PDB contiver vários modelos. O primeiro modelo é mostrado inicialmente.
- Ligue/desligue ligantes e água. Clique em um átomo para identificar resíduo, cadeia, nome do átomo e elemento.
- Use **Tela inteira** para ampliar a área e **Salvar imagem** para baixar o enquadramento atual em PNG.

O motor [3Dmol.js 2.5.5](https://github.com/3dmol/3Dmol.js/releases/tag/2.5.5) está incluído nos recursos do aplicativo, com licença e hashes de integridade. A visualização não depende de CDN, não consulta RCSB e não envia a estrutura a provedores externos. A biblioteca utiliza WebGL e perspectiva no próprio navegador. Um navegador com WebGL disponível é necessário; a página orienta o usuário se a inicialização gráfica falhar.

No modo web, as rotas do visualizador compartilham a origem/porta do Flet. No desktop, um servidor temporário da própria aplicação atende somente em `127.0.0.1` e encerra com a aplicação. Os links opacos são temporários e verificam sessão e permissões ao servir a página e o PDB. Logout, revogação de acesso e expiração impedem novas requisições; o conteúdo já carregado em uma aba continua sendo uma cópia local, como um download. Páginas e estruturas usam `Cache-Control: no-store`; pastas dos projetos não são publicadas como assets. O limite de abertura continua em 32 MB. O visualizador não altera coordenadas, curadoria de ligantes ou resultados do pipeline.

A biblioteca JavaScript é distribuída junto com a aplicação e não requer instalação de outro programa ou pacote Python. A hospedagem web utiliza FastAPI/Uvicorn já fornecidos pelo extra `ui` (Flet web).

### Verificar o visualizador

Os testes de autorização e integração estão em `tests/test_pdb_viewer.py`. Para validar WebGL, gestos e exportação em um Chrome/Chromium instalado, execute no ambiente da interface:

```bash
PYTHONPATH=src python scripts/validate_pdb_viewer.py --chrome /usr/bin/google-chrome
```

O script usa apenas um servidor temporário em loopback e um perfil de navegador temporário, confirma rotação, zoom, representações e PNG e recusa chamadas a serviços externos. A captura visual padrão fica em `/tmp/biomol-pdb-modern.png`; `--screenshot` permite outro destino.


## PubChem como bloco independente

**Recuperar ChEMBL** e **Recuperar PubChem** são blocos separados na interface. PubChem aceita referência manual por SMILES, CID ou nome, ou compostos conectados/enviados por CSV. Com uma tabela, escolha todas as referências ou um código específico. O limiar padrão é 75% e o limite padrão é 1.000 similares por referência. A saída contém somente novos similares, sem repetir as referências; preserva códigos `PUBCHEM<CID>` e SMILES. Conecte essa saída diretamente a Vina ou DOCK6, a ADMET ou a fingerprints. A operação backend é `retrieve_pubchem`; a API legada ChEMBL ainda aceita ampliação integrada, mas essa opção não aparece no bloco ChEMBL da interface.

## Tranches ZINC 2D e 3D

No navegador de [tranches ZINC20](https://zinc20.docking.org/tranches/home/), selecione 2D/SMI ou 3D/MOL2 e exporte a lista de downloads. O bloco aceita listas TXT/URI/URLS e arquivos de comandos cURL, wget ou PowerShell. São extraídos somente os links; os comandos não são executados. Links HTTP dos servidores oficiais são convertidos para HTTPS, inclusive nos redirecionamentos. Os formatos aceitos são `.smi`, `.mol2`, `.smi.gz`, `.mol2.gz` e as versões `.bz2`.

Abra **Recuperar ZINC**, use **Enviar meus arquivos** em **Lista de downloads das tranches ZINC** e escolha o arquivo exportado. Alternativamente, importe-o com o tipo **Lista de downloads ZINC** e conecte a importação. Selecione uma ou mais listas; merge reúne seus links sem repetições. Um [exemplo pequeno](../examples/zinc_tranches.uri) está disponível para conferir o fluxo. A distribuição [2D oficial](https://cache.docking.org/2D/) descreve os arquivos SMI, e o [guia oficial de screening](https://wiki.docking.org/index.php?title=ZINC15:examples:screening) descreve as exportações 3D comprimidas.

| Saída | Conteúdo e uso |
| --- | --- |
| `compounds.csv` | `molecule_chembl_id` com código ZINC, `canonical_smiles`, `source`, `source_url`, `structure_format`, `conformer_file` e `conformer_origin` |
| `Conformers/<ZINC>.mol2` | Estrutura 3D original, preservando coordenadas, tipos atômicos e cargas; acompanha a tabela |
| `retrieval_report.json` | URLs, formato, número de registros/compostos, duplicatas e SHA-256 de cada download |

Em **Downloads simultâneos**, escolha de 1 a 16 threads (padrão: 4). Use 1 para execução sequencial; aumente conforme a conexão e a capacidade dos servidores. O ganho depende da rede e dos limites do ZINC. Cada download mantém as validações de endereço e redirecionamento e as tentativas automáticas para falhas HTTP temporárias. O número de arquivos baixados antecipadamente é limitado às threads configuradas; os arquivos temporários são removidos após o processamento. A consolidação segue a ordem da lista, preservando a deduplicação e a escolha do primeiro conformero. `retrieval_report.json` registra o número de downloads simultâneos efetivamente usado.

A leitura SMI aceita cabeçalho opcional, espaços e tabs. Um identificador numérico recebe o prefixo ZINC; códigos ZINC existentes são preservados. Os blocos de um MOL2 com várias moléculas são separados por código. Duplicatas idênticas são reunidas; se há 2D e 3D para o mesmo código, a conformação 3D acompanha o registro. O primeiro conformero válido é mantido. Códigos com estruturas divergentes, SMILES inválidos, arquivos vazios e formatos incompatíveis interrompem a etapa.

Conecte `compounds.csv` a ADMET, fingerprints, similaridade/grafos por seu fluxo habitual, **Preparar para docking**, Vina ou DOCK6. Em 2D, o preparo gera uma conformação; em 3D, parte do MOL2 fornecido. O manifesto marca as estruturas 3D como conformações de biblioteca (`conformer_origin=library`), para posicioná-las no centro do sítio ao preparar a entrada DOCK6. Poses já calculadas conservam seu posicionamento. Não se aplica o tipo **Resultados DOCK6** ao MOL2 da biblioteca, pois ele ainda não contém scores de docking. Mantenha a tabela e os arquivos referenciados no mesmo conjunto ao exportar/importar resultados.
