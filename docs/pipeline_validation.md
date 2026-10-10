# Validação das etapas e conexões do pipeline

[Documentação](README.md) · Português · [English](en/pipeline_validation.md)

A revisão cobre o catálogo de operações, contratos de arquivos, relações entre blocos, seleção e materialização de entradas, processamento individual/merge, execução de ramificações, persistência e reaproveitamento. O [manual do usuário](user_manual.md) ensina o percurso operacional; esta página registra como verificá-lo e os limites da evidência.

## Matriz de contratos

| Operação | Entradas e condições | Saída para outras etapas |
| --- | --- | --- |
| Importar arquivos | Assets autorizados, validados conforme o tipo de cada arquivo | Um ou mais tipos publicados pelo conjunto selecionado |
| Recuperar compostos | Alvo/filtros; consulta ChEMBL e ampliação PubChem opcional | Compostos e contexto ChEMBL original |
| Recuperar PubChem | Referência manual, arquivo ou compostos conectados | Novos similares PubChem com código e SMILES |
| Expandir similares | Recuperação com os downloads ChEMBL originais | Compostos |
| Recuperar estruturas | Critérios PDB e consulta externa | Estruturas brutas e metadados |
| Recuperar ZINC | Lista de tranches do tipo `zinc_urls` ou `other` | `compounds.csv`, conformações MOL2 e relatório |
| Preparar para docking | Receptor PDB bruto ou do redocking; uma ou mais fontes de compostos | Receptores, centros, metadados e candidatos preparados para Vina, DOCK6 ou ambos |
| ADMET | Compostos | Compostos avaliados/filtrados e visualizações |
| Fingerprints | Compostos | Fingerprints com tipo e largura determinados |
| Similaridade | Fingerprints compatíveis | Relações `source,target,value` e linhagem molecular |
| Grafos | Relações de blocos de similaridade; dados externos pelo formulário | Compostos do MCC e visualizações |
| Redocking | Brutos se `prepare_complex=true`; preparados se `false` | Estruturas, metadados de avaliação e preparados |
| Vina | Preparados e compostos | Poses identificadas pela referência completa e composto |
| DOCK6 | Preparados e candidatos de qualquer origem molecular; poses Vina opcionais | Scores MOL2 e demais resultados do protocolo |
| Consenso | Poses Vina e DOCK6 com identificadores correspondentes | Tabela de scores e comparação |

`tests/test_connection_matrix.py` percorre todos os produtores nativos, saídas prontas e tipos de importação contra todas as portas, inclusive as duas configurações de redocking. Ligações permitidas são aplicadas e submetidas à ordenação do pipeline. Ligações recusadas não alteram a configuração. Também são verificadas fontes desativadas, autorrelações e o alias interno usado pelo consenso. Os testes de fluxo existentes verificam ciclos, exclusão, duplicação e normalização de projetos antigos.

A compatibilidade de tipos não substitui o contrato dos arquivos: a execução verifica cabeçalhos, conteúdos, caminhos e dependências. Grafos só recebem uma conexão de bloco **Calcular similaridade**, para manter a linhagem; CSVs externos podem ser selecionados diretamente no formulário. A expansão ChEMBL exige os downloads originais; resultados prontos genéricos não publicam esse contexto.

## Correções e regressões cobertas

| Caso | Comportamento verificado |
| --- | --- |
| Seleção repetida | Configurar a etapa recebe escolhas do popup; confirmar conserva o modo e as referências |
| Individual versus merge | Cálculos reais mantêm datasets separados ou geram um conjunto combinado conforme a decisão |
| Várias origens Vina no consenso | O alias acompanha todo o grupo escolhido, sem perder origens anteriores |
| Nome de resultado repetido | Um seletor ambíguo falha explicitamente; caminhos distinguem arquivos |
| Preparação → docking | Metadados, centros nativos/legados e auxiliares MOL2/noH acompanham o receptor correto; candidatos entram separadamente, sem ligantes de referência |
| Redocking preparado | Registros de quatro campos são normalizados; a porta aceita preparados somente com preparação desativada |
| Vina com várias referências | PDB/ligante/resíduo/cadeia distinguem poses; resultados de uma referência não fazem outra ser ignorada |
| DOCK6 individual | Apenas combinações com receptores, candidatos e poses correspondentes são executadas |
| Seleção de compostos para DOCK6 | A tabela restringe os conformeros passados ao motor |
| Consenso individual | Todos os resultados selecionados são reunidos antes da interseção por receptor e composto |
| Scores de consenso | Inteiros e notação científica são lidos; ausentes/não finitos são recusados |
| Um composto ou scores constantes | Normalizações são finitas e zero; a saída nativa pode ser importada novamente |
| Reinicialização/reaproveitamento | Entradas/configurações compatíveis retomam resultados; dados alterados ou hashes adulterados impedem adoção |
| Grafos/fragmentos | Componentes e vértices separados; SMILES válido extraído da referência, com estereoquímica de cortes verificada |
| Idioma | Sessões independentes, troca no login, traduções e preservação de parâmetros/dados do usuário |

Casos específicos de docking estão em `test_docking_handoffs.py`; os motores são substituídos na fronteira de comando, mas cópias, metadados, correspondências e leitura/normalização de scores são reais. `test_input_processing.py` executa fingerprints → similaridade → grafos em workers reais, com dois datasets e ambos os modos. Os demais testes cobrem ADMET real, cancelamento, timeout, isolamento, permissões, importação/exportação, histórico e cache.

## Reproduzir a verificação

Execute na raiz, usando o Python do ambiente científico e com a interface instalada:

```bash
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/biomol-mpl PYTHONPATH=src:tests \
  python -m unittest discover -s tests -q
python scripts/build_docs.py
```

Para revisar apenas os contratos acrescentados:

```bash
PYTHONPATH=src:tests python -m unittest \
  test_connection_matrix test_docking_handoffs test_localization -v
```

O teste de conexões é exaustivo para o catálogo declarado. Isso não significa executar todas as combinações possíveis de parâmetros científicos nem comprovar disponibilidade dos provedores ou funcionamento de cada versão de motor externo.

## Validar o protocolo no ambiente de uso

1. Confirme o Python do worker e a instalação dos executáveis utilizados.
2. Escolha um complexo de referência e um pequeno conjunto de candidatos conhecidos.
3. Execute recuperação/importação, preparação e redocking; confira ligante, cadeia, centros, RMSD e logs.
4. Execute Vina e DOCK6 com os mesmos receptores/códigos; confira que candidatos e poses correspondem aos arquivos selecionados.
5. Gere consenso e compare scores brutos com os arquivos dos motores; examine o efeito do peso de repulsão e das normalizações.
6. Repita com dois datasets em individual e merge; confirme a composição dos conjuntos em cada etapa.
7. Reinicie, escolha reaproveitar e compare com uma nova execução quando precisar validar mudanças de ambiente.

Consultas externas dependem de rede, limites e disponibilidade de ChEMBL, PubChem, PDB e ZINC. Os testes determinísticos não fazem uma campanha de consultas ao vivo. Preparação/docking precisam de validação com motores e dados reais na instalação de destino. O reaproveitamento após mudanças de ambiente é uma escolha explícita do pesquisador; para recalcular com a instalação atual, escolha **Não, executar novamente**. Regras ADMET e scores de docking são resultados computacionais, não comprovação experimental.

## Regressões de redocking e idioma

A revisão inclui registros estruturados do pandas (que não aceitam fatiamento como listas), seleção pré-configurada sem nova solicitação, ausência de pares, resíduos inválidos, ferramentas ausentes, caminhos com espaços, cofatores e centros preparados. Erros de preparação são propagados com sua causa em vez de retornar um conjunto vazio.

Os testes percorrem formulários de todas as operações em inglês, verificam mensagens aninhadas de workers e preservam valores científicos e nomes personalizados. Rótulos dos templates, nomes padrão das etapas e seletores também acompanham o idioma da sessão. Mensagens originalmente em inglês são traduzidas para português quando apresentadas nessa sessão.

Os testes automatizados de integração substituem ferramentas externas na fronteira de execução. Após a análise dos logs de 6 de outubro de 2026, também foi executado o caso real 4M0E / 1YL / 604 / A com Chimera, Open Babel, Vina e PyMOL: com as opções originais, incluindo minimização do receptor e do ligante, o redocking concluiu e produziu RMSD de aproximadamente 0,149 Å. O teste complementar sem minimização produziu aproximadamente 0,215 Å. As execuções utilizaram cópias temporárias da entrada; esses resultados verificam o fluxo desse caso, sem estabelecer tolerâncias científicas para outros complexos.

## Resultados e visualização 3D

`tests/test_redocking_results.py` cobre agrupamento por simulação, coleções importadas com identificadores iguais, artefatos duplicados, valores RMSD inválidos, nomes repetidos entre referência e poses no ZIP, permissões, sessões revogadas e ações do popup em português e inglês. `tests/test_pdb_viewer.py` valida formatos e acesso ao visualizador.

Também foram verificados no navegador um ligante MOL2 (35 átomos), um receptor PDBQT (4.142 átomos) e poses PDBQT (21 átomos), incluindo geometria visível, rotação, zoom, quatro representações e exportação PNG, sem erros JavaScript nem solicitações externas. Essas verificações complementam os testes automatizados; não substituem a avaliação científica das estruturas. Para repetir com um arquivo próprio, execute:

```bash
python scripts/validate_pdb_viewer.py --structure /path/to/structure.mol2
```

Consulte [Logs e diagnóstico](logging.md) para o formato comum, contexto por execução, códigos de falha, resumo do job e o comando `python -m biomolexplorer.log_report`.


## Docking interoperável: validação de 7 de outubro de 2026

O teste real `scripts/validate_docking.py` utilizou o receptor 1ABE e o ligante carregado do tutorial `ligand_sampling_demo` distribuído com DOCK6, além de etanol fornecido por SMILES. Vina e DOCK6 executaram separadamente, seguidos de Vina → DOCK6 e DOCK6 → Vina. As quatro execuções conservaram `REF` e `ETHANOL`, produziram scores finitos e exportaram as poses. Os três consensos (independente e ambos os encadeamentos) encontraram os dois compostos para `1ABE_A`; o caso sem interseção retornou `skipped_reason` sem tabela de scores. Os valores estão no [relatório JSON](validation/docking_2026-10-07.json).

| Composto | Vina independente | DOCK6 independente | Vina → DOCK6 | DOCK6 → Vina |
| --- | --- | --- | --- | --- |
| REF | -6.511 | -12.755546 | -23.418610 | -6.509 |
| ETHANOL | -2.673 | -13.771717 | -11.354554 | -2.673 |

O teste usa esforço Vina 1, duas poses e busca DOCK6 rígida. Verifica execução e contratos de software desse caso, sem definir critérios de qualidade científica ou comparar numericamente motores diferentes. No consenso, o score DOCK6 inclui o peso de repulsão e pode diferir do score bruto dessa tabela.

Para repetir, ative o ambiente científico e informe uma pasta de saída nova:

```bash
PYTHONPATH=src python scripts/validate_docking.py --dock6-root /path/to/dock6 --output /tmp/docking-check
PYTHONPATH=src:tests python -m unittest test_docking_interoperability test_docking_handoffs test_connection_matrix test_results_display test_localization
```

As regressões cobrem entradas isoladas e aliases CSV, poses nos dois sentidos, identidade de uma pose individual do consenso, múltiplos lotes, receptores distintos, códigos com estruturas divergentes, remoção auditada e autorização das poses. Tabelas principais vazias permanecem válidas para o consenso; seus compostos não reaparecem a partir de auxiliares. Os cálculos DOCK6 usam cópias em caminhos temporários curtos e exportam artefatos para o projeto, evitando truncamento nos auxiliares legados; falhas dos subprocessos preservam seu diagnóstico ao atravessar o executor paralelo.

As poses PDBQT e MOL2 do consenso também foram verificadas no navegador: rotação, zoom, quatro representações, exportação PNG e ausência de solicitações externas ou erros JavaScript. A leitura de PDBQT remove registros de torção que seriam confundidos com limites de modelos, preserva os átomos de cada pose e identifica arquivos com scores Vina como ligantes mesmo após renomeação pelo consenso.

A suíte geral final executou 442 testes, sem falhas ou erros; duas verificações de servidor localhost foram ignoradas pelas restrições do sandbox. Os testes de navegador foram executados com acesso local autorizado e verificaram a pose Vina com todos os 14 átomos e a pose DOCK6 com 20 átomos. Um cálculo complementar DOCK6 **flexível**, reutilizando as grades e a preparação do caso independente, produziu scores -38,876862 para REF e -15,546004 para ETHANOL.

Para validar os links locais e âncoras após alterações na documentação, execute `python scripts/build_docs.py` e `python scripts/validate_docs.py`. Foram verificadas 27 páginas HTML, incluindo `index.html`, 24 guias e dois portais de idioma.

## Preparar para docking: validação de 8 de outubro de 2026

As regressões de `test_preparation_settings` e `test_docking_preparation` verificam receptores brutos e reutilizados, seleção sem ligantes de referência, configurações independentes, fontes ChEMBL/PubChem/ZINC combinadas, exportação para cada motor a partir de um único preparo e reaproveitamento dos arquivos sem conversão adicional. Centros inválidos, formatos ausentes e mistura de receptores brutos/preparados são rejeitados. As conexões respeitam os motores selecionados, e uma conexão rejeitada não altera a configuração do receptor.

Para repetir os testes focados e a validação com ferramentas reais, ative o ambiente científico e execute:

```bash
PYTHONPATH=src:tests python -m unittest test_preparation_settings test_docking_preparation test_docking_handoffs test_connection_matrix test_docking_interoperability
PYTHONPATH=src python scripts/validate_docking.py --dock6-root /path/to/dock6 --prepare-inputs --output /tmp/preparation-check
```

O segundo comando reutiliza o receptor do tutorial, prepara os candidatos com Chimera e Open Babel e executa Vina, DOCK6 e os dois encadeamentos. Os testes de seleção e fontes múltiplas usam dados locais; não consultam os provedores ao vivo.

Os 48 testes focados passaram. A execução real com `--prepare-inputs` concluiu o preparo de REF e ETHANOL, preservou o receptor byte a byte e validou Vina e DOCK6 independentes, ambos os encadeamentos e três consensos. Consulte o [relatório JSON](validation/preparation_2026-10-08.json). O caso sem compostos em comum foi sinalizado sem tabela de scores. Essa execução valida os contratos desse caso; não estabelece qualidade científica para outros complexos.

A suíte geral final executou 465 testes sem falhas ou erros, com duas verificações de servidor local ignoradas pelas restrições do sandbox. A revisão confirmou os links locais dos 26 documentos Markdown e os arquivos, IDs e âncoras das 27 páginas HTML.

## Importação e tranches ZINC: revisão de 8 de outubro de 2026

Os 75 testes focados passaram. `test_import_files` cobre seleção no disco/projeto, tabela, remoção sem apagar o original, tipos por arquivo, validação ao alterar o tipo, compatibilidade com receptores e compostos e modo de leitura. `test_zinc_tranches` cobre listas URI e scripts, HTTPS, formatos comprimidos, SMI sem cabeçalho, ordem das colunas, MOL2 individual, duplicatas, conflitos, redirecionamentos e materialização de múltiplas listas. Também verifica posicionamento das conformações de biblioteca no sítio DOCK6 e preservação dos originais.

Um download real da tranche pública `AAIA.smi` produziu 11 compostos com identificadores ZINC e SMILES válidos. O [relatório JSON](validation/zinc_2026-10-08.json) registra URL e hash. A leitura MOL2 comprimida e os contratos de conformação foram verificados com respostas determinísticas; este teste ao vivo não baixou uma biblioteca 3D completa.

Para repetir os testes focados:

```bash
PYTHONPATH=src:tests python -m unittest test_import_files test_zinc_tranches test_preparation_settings test_localization test_docking_preparation test_docking_handoffs test_guided_ui test_connection_matrix
```

A suíte geral desta revisão executou 482 testes sem falhas ou erros; duas verificações de servidor local foram ignoradas pelas restrições do sandbox. Também foram validados os links dos 26 documentos Markdown e os arquivos, IDs e âncoras das 27 páginas HTML.

## Ajuda e configuração do docking: revisão de 8 de outubro de 2026

`test_pipeline_guidance` cobre a recusa de conexões ChEMBL/PubChem → similaridade pela ação do canvas, preservação das conexões e do histórico e criação de uma conexão válida após a recusa. Também verifica a ajuda de todos os blocos para leitores e a tradução das explicações.

Os testes da seleção de compostos verificam os dois motores: filtro por identificador, padrão de todos os compostos, seleção por arquivo, troca de arquivo, preservação do estado do popup, cópia de arquivos moleculares preparados, materializações distintas para seleções distintas e mensagem para código ausente. Na execução real do serviço de pipeline, o supervisor científico é substituído por uma resposta determinística: entradas explícitas executam sem novo popup; conexões automáticas continuam aguardando seleção.

```bash
PYTHONPATH=src:tests python -m unittest test_pipeline_guidance test_connection_matrix test_visual_flow test_file_selection test_pipeline_execution test_localization
```

A suíte geral desta revisão executou 495 testes sem falhas ou erros; duas verificações de servidor local foram ignoradas pelas restrições do sandbox. Os oito testes novos passaram. As 27 páginas HTML foram verificadas quanto a arquivos locais, IDs e âncoras.

## Revisão do projeto PDB em 9 de outubro de 2026

Os logs do projeto PDB registraram **O arquivo escolhido não está nos resultados
da etapa** antes de qualquer cálculo de redocking. O seletor de `4M0E.pdb`
continha o identificador de um worker anterior; nova recuperação ou cache mudavam
essa pasta. Os seletores agora preservam lote e caminho lógico, removendo somente
a pasta temporária do job. Seletores antigos são resolvidos pelo mesmo contrato,
inclusive nas entradas de compostos, receptores e resultados de Vina/DOCK6.
Arquivos ausentes e seleções ambíguas continuam sendo recusados.

Em cópias temporárias, o par selecionado **4M0E / NO3 / 608 / B** concluiu o
redocking com as escolhas de preparação do projeto, incluindo minimização. Com
busca reduzida (`exhaustiveness=1`, `num_modes=2`), produziu RMSD de **2,888 Å** e
score Vina de **−2,887**. Sem minimização, produziu RMSD de **3,245 Å**. Esses testes
verificam execução, preparação, pose e cálculo de RMSD; não reproduzem a busca
completa configurada pelo usuário (`20` e `10`) nem estabelecem aceitação
científica do resultado. O [relatório de PDB/redocking](validation/pdb_redocking_2026-10-09.json)
registra causa, correção e limites.

Vina e DOCK6 também concluíram a validação real com os compostos REF e ETHANOL e
receptor 1ABE_A do exemplo de DOCK6: preparação para ambos os motores, execuções
independentes, reutilização de poses Vina → DOCK6 e DOCK6 → Vina, consenso e
tratamento da interseção vazia. Consulte o
[relatório de docking](validation/docking_2026-10-09.json).

## Visualizador e interações: revisão de 10 de outubro de 2026

Os testes focados de geometria, cenas/autorização, visualizador, resultados de
redocking e visualizações executaram 40 verificações com sucesso; duas verificações
de servidor local foram ignoradas no sandbox. A validação Chrome separada abriu
cenas sintéticas equivalentes de Vina e DOCK6 e um conformero isolado, verificando
estilos independentes, H, rotação, zoom, PNG, traçados, filtros e visibilidade.
Movimentos reais do mouse sobre os segmentos confirmaram o tooltip de cada tipo
presente e seu fechamento ao sair, sem erros JavaScript ou requisições externas.

Essa cobertura valida o fluxo e a geometria controlada; não substitui validação
científica em complexos experimentais nem cobre todos os tipos de interação de
outros classificadores. Consulte o
[relatório do visualizador](validation/viewer_controls_2026-10-10.md) e o
[guia de uso](molecular_viewer.md).

```bash
PYTHONPATH=src:tests python -m unittest test_docking_interactions test_docking_scene test_pdb_viewer test_redocking_results test_result_visualizations
python scripts/build_docs.py
python scripts/validate_docs.py
```
