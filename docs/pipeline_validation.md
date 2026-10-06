# Validação das etapas e conexões do pipeline

[Documentação](README.md) · Português · [English](en/pipeline_validation.md)

A revisão cobre o catálogo das 14 operações, contratos de arquivos, relações entre blocos, seleção e materialização de entradas, processamento individual/merge, execução de ramificações, persistência e reaproveitamento. O [manual do usuário](user_manual.md) ensina o percurso operacional; esta página registra como verificá-lo e os limites da evidência.

## Matriz de contratos

| Operação | Entradas e condições | Saída para outras etapas |
| --- | --- | --- |
| Importar arquivos | Assets autorizados, validados conforme o tipo selecionado | Tipo explicitamente publicado |
| Recuperar compostos | Alvo/filtros; consulta ChEMBL e ampliação PubChem opcional | Compostos e contexto ChEMBL original |
| Expandir similares | Recuperação com os downloads ChEMBL originais | Compostos |
| Recuperar estruturas | Critérios PDB e consulta externa | Estruturas brutas e metadados |
| Recuperar ZINC | Lista de endereços do tipo `other` | Compostos |
| Preparar estruturas | PDBs brutos e registros de referência | Receptores/ligantes preparados, centros e metadados |
| ADMET | Compostos | Compostos avaliados/filtrados e visualizações |
| Fingerprints | Compostos | Fingerprints com tipo e largura determinados |
| Similaridade | Fingerprints compatíveis | Relações `source,target,value` e linhagem molecular |
| Grafos | Relações de blocos de similaridade; dados externos pelo formulário | Compostos do MCC e visualizações |
| Redocking | Brutos se `prepare_complex=true`; preparados se `false` | Estruturas, metadados de avaliação e preparados |
| Vina | Preparados e compostos | Poses identificadas pela referência completa e composto |
| DOCK6 | Preparados, candidatos e poses Vina correspondentes | Scores MOL2 e demais resultados do protocolo |
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
| Preparação → docking | Metadados, centros nativos/legados, ligantes e auxiliares MOL2/noH acompanham o receptor correto |
| Redocking preparado | Registros de quatro campos são normalizados; a porta aceita preparados somente com preparação desativada |
| Vina com várias referências | PDB/ligante/resíduo/cadeia distinguem poses; resultados de uma referência não fazem outra ser ignorada |
| DOCK6 individual | Apenas combinações com receptores, candidatos e poses correspondentes são executadas |
| Seleção de compostos para DOCK6 | A tabela restringe os conformeros passados ao motor |
| Consenso individual | Resultados correspondentes são pareados, em vez de cruzar todos os arquivos |
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
