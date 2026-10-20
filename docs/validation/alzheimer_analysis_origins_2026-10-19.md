# Alzheimer: seis conjuntos e origem dos resultados — 2026-10-19

## Inspeção do projeto utilizado

O workspace existente foi consultado em modo SQLite somente leitura. O projeto **Alzheimer** contém as recuperações **ACHE / CHEMBL220**, **BCHE / CHEMBL1914** e **BACE1 / CHEMBL4822**. ADMET seleciona BBB+ e fornece três entradas para a expansão PubChem. Fingerprints recebe essas três saídas PubChem e, após a alteração do usuário, também as três saídas BBB+ diretas do ADMET.

Na execução consultada, fingerprints possuía **seis lotes**, mas a configuração salva da similaridade ainda selecionava somente os três arquivos `morgan_compounds.csv` originados da PubChem. Grafos mantinha os três seletores correspondentes de similaridade. A perda dos conjuntos ADMET acontecia na seleção de entrada da similaridade, antes da elaboração dos grafos. Não era uma substituição das moléculas ChEMBL pelas moléculas PubChem dentro do algoritmo de fingerprints.

A seleção salva anterior era uma seleção completa quando a origem produzia apenas três arquivos. Ela não acompanhou o crescimento para seis. O problema de apresentação escondia essa diferença: arquivos repetidos exibiam o nome genérico e o produtor imediato, sem o alvo e a procedência completa.

O arquivo JSON deste relatório registra as contagens e a matriz esperada. A fixture `tests/fixtures/alzheimer_analysis_pipeline.json` reproduz as oito etapas relevantes, suas conexões e seus parâmetros, removendo identificadores de artefatos que pertencem a arquivos privados e caminhos absolutos. Nenhum dado ou seleção foi escrito no projeto de produção durante a validação.

## Correção da execução

- Antes de verificar o cache, similaridade e grafos reconciliam a seleção com as saídas atuais de seus produtores.
- Uma seleção antiga que cobria todas as saídas anteriores pode acompanhar a ampliação de três para seis. A verificação usa os seletores e o histórico real, sem nomes de alvos fixos no algoritmo.
- A escolha efetiva é persistida e auditada, preservando mudanças concorrentes no projeto. O cache é consultado com os seis arquivos efetivamente selecionados.
- Seleções parciais antigas são preservadas. A opção **Acompanhar novos conjuntos da origem** permite controlar esse comportamento; desativá-la conserva apenas os arquivos explicitamente escolhidos.
- O processamento individual gera seis análises independentes. Cada similaridade mantém o seu arquivo de compostos; moléculas e relações de AChE, BChE, BACE1, ADMET e PubChem não são misturadas entre lotes.

Resultado esperado:

| Origem | ADMET BBB+ direto | PubChem a partir de ADMET BBB+ |
|---|---|---|
| ACHE · CHEMBL220 | 1 grafo | 1 grafo |
| BCHE · CHEMBL1914 | 1 grafo | 1 grafo |
| BACE1 · CHEMBL4822 | 1 grafo | 1 grafo |

## Correção da apresentação

A origem é reconstruída pela procedência transitiva dos arquivos efetivamente usados, inclusive quando produzidos numa execução anterior e reutilizados. As tabelas exibem o alvo/identificador e permitem filtrar a origem. Os seletores pesquisam o identificador e distinguem ADMET direto de PubChem. ADMET e grafos exibem um cabeçalho com a origem e o caminho de arquivos no pipeline.

A fonte molecular, a inspeção dos compostos e os relatórios JSON/GEXF/GraphML preservam essa identificação. Os títulos das apresentações do grafo incluem a origem. Arquivos auxiliares com o mesmo nome são limitados ao lote da análise selecionada, inclusive na escolha automática da tabela molecular correlata.

Resultados já existentes recebem os novos rótulos ao serem abertos, sem alterar seus CSVs ou recalcular o grafo para descobrir sua origem. Não há atribuição de alvo baseada somente na aparência de um nome de arquivo.

## Validação

1. A configuração e os metadados de execuções do projeto real foram inspecionados sem autenticação por tokens de sessão ou escrita no banco. O resolvedor reconstruiu as seis origens de fingerprints: três ADMET diretas e três via PubChem. A expansão em memória da configuração real passou de três para seis seletores de similaridade.
2. Uma réplica isolada executou inicialmente três conjuntos PubChem. Em seguida, recebeu as três conexões diretas de ADMET presentes no projeto do usuário. O agendador e o armazenamento reais produziram **seis lotes de fingerprints, seis de similaridade e seis de grafos**.
3. Fingerprints, similaridade e elaboração dos grafos foram executados pelos adaptadores científicos reais com RDKit e os componentes locais. Downloads e avaliação ADMET foram substituídos por pequenos conjuntos controlados: a validação não faz novas consultas às bases nem executa novamente o modelo ADMET-AI.
4. Cada um dos seis grafos possui dois compostos e uma relação no conjunto de teste. Os testes verificam o alvo e o ramo de cada ID, a ausência de mistura entre as entradas, a manutenção das fontes anteriores, a persistência das seis seleções e os seis rótulos distintos. Também confirmam que as recuperações e ADMET anteriores não são repetidos.
5. A pesquisa por CHEMBL220 na seleção de fingerprints encontra os dois ramos. Testes adicionais cobrem seleção parcial, procedência de resultados reutilizados quando o alvo atual mudou e proteção contra ciclos no histórico.

Bateria de apresentação/seleção/execução/cache:

```text
python -m unittest test_artifact_registry test_file_selection
  test_results_display test_graph_input_contract test_analysis_handoff
  test_stage_cache test_restart_reuse test_pipeline_guidance test_localization -q
```

**98 testes aprovados em 54,552 segundos.**

Bateria de seis grafos, análise molecular e regressões de redocking:

```text
python -m unittest test_analysis_origins test_molecular_network_dashboard
  test_graph_navigation_gephi test_graph_explorer test_redocking_incremental
  test_redocking_partial_failures -q
```

**103 testes aprovados em 61,035 segundos.** As duas baterias somam **201 testes aprovados**. Os avisos de depreciação de `fork` emitidos pela biblioteca durante a réplica não impediram a conclusão dos cálculos.

A documentação e as traduções foram atualizadas. Verificação de espaçamento: `git diff --check`.

## Retomada no projeto do usuário

Reinicie a aplicação para carregar a correção e execute o pipeline com reaproveitamento dos resultados. A seleção completa antiga é reconciliada com os novos conjuntos antes de decidir quais análises reutilizar. O projeto original não foi reexecutado integralmente nesta revisão; a produção de seis grafos foi comprovada na réplica isolada da sua topologia.
