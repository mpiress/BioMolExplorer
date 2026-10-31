# Docking with Vina: seleção múltipla de BBB+ — 31/10/2026

## Resultado e escopo

O gargalo do preparo de entradas foi reproduzido e corrigido. Foram usados todos os 12 BBB+ das saídas atuais de **ADMET** e **Evaluate ADMET**, com as três entradas de receptores da configuração existente: **6ZWI_A, 4M0E_A e 2IQG_A**. As 36 combinações de receptor/conjunto completaram a resolução e normalização das entradas em um projeto temporário.

Esses arquivos somam **47.536 linhas**: 183 em ADMET e 47.353 em Evaluate ADMET. Essa soma não significa compostos únicos entre todos os conjuntos. Os parâmetros científicos foram preservados: pH 7,4, caixa 24×24×24 Å, exhaustiveness 20, num_modes 10, receptores previamente preparados, cargas AM1 e minimização dos ligantes.

Foram verificados os hashes de **980 arquivos originais**, preservados integralmente. Nenhuma execução real foi cancelada/retomada por esta revisão e nenhum resultado do usuário foi excluído ou substituído. As execuções científicas de teste usaram apenas cópias temporárias. Não foram instaladas bibliotecas e o requirements.yml não precisou de alterações.

## Diagnóstico

A execução persistida estava interrompida no bloco Vina, sem job_id nem log de worker desse bloco. Os eventos científicos disponíveis não demonstram que o Vina tivesse sido iniciado nessa tentativa. A interrupção persistida, isoladamente, não identifica a pilha exata que estava ativa; os gargalos abaixo foram reproduzidos diretamente no caminho correspondente do código.

1. **Leituras auxiliares excessivas.** `docking_data.input_records` varria todas as tabelas disponíveis para procurar metadados de poses, inclusive quando a entrada escolhida já era um CSV molecular completo. Na reprodução de um único BBB+ com 3.362 linhas, leu 49 tabelas, somando 446.560 linhas. As combinações entre receptores e arquivos, além das sondagens de cache, repetiam esse trabalho.
2. **Normalização repetida.** `_cache_key` chama `_resolve(cache_only=True)`, que materializa as entradas. A mesma tabela era novamente lida, copiada e validada quimicamente em diferentes fases e para cada receptor, mesmo sem alterações no conjunto molecular.
3. **Bloqueio prolongado do coordenador.** `resume` mantinha o lock global durante a resolução de todas as combinações. Uma validação demorada também podia bloquear cancelamento e outras ações do serviço.
4. **Seletores grandes na interface.** O catálogo de identificação executava limpeza química completa para preencher um campo opcional de IDs; dropdowns construíam uma opção por composto. A seleção de vários arquivos ampliava simultaneamente o trabalho de leitura, a memória e o número de controles enviados ao frontend.

## Correções

- Tabelas auxiliares só são examinadas para entradas que realmente necessitam resolver poses avulsas. CSVs escolhidos são lidos diretamente. A autorização das poses referenciadas na própria tabela continua obrigatória, incluindo os aliases `vina_pose` e `dock6_pose`.
- A normalização de CSVs sem poses usa o **CalculationCache existente**, em namespace próprio. As chaves incluem arquivos selecionados, SHA-256 de seus conteúdos, seleção opcional de composto, versão da transformação e versão do RDKit. A tabela normalizada pode ser compartilhada entre receptores e entre blocos com a mesma necessidade de entrada.
- Entradas com poses seguem o caminho que preserva suas conformações e arquivos associados. As opções de preparo científico continuam sendo consideradas pelos caches científicos já existentes, depois da normalização da tabela.
- Os objetos do cache são verificados por hash. A publicação da tabela normalizada ocorre a partir de uma área temporária, com cópia atômica e índice transacional. Uma alteração da fonte durante a normalização impede publicar um checkpoint incorreto; corrupção provoca reconstrução/restauração validada.
- A validação onerosa da retomada ocorre fora do lock global. Antes de efetivar a retomada, o serviço confere novamente estado, permissões e snapshot. O cancelamento permanece disponível durante essa validação; entre combinações, o estado da execução é consultado para interromper trabalho desnecessário.
- `compound_identifiers` faz a leitura leve dos IDs para apresentação. A validação química dos dados efetivamente executados permanece no backend. A busca PubChem continua usando referências com SMILES, conforme seu contrato.
- O componente compartilhado `CompoundChoice` apresenta até 100 IDs por página, com pesquisa, navegação, opção de todos os compostos e preservação do ID previamente selecionado. Ambos os caminhos de configuração usam o mesmo componente. O preparo antecipado limita as referências de docking às origens conectadas.

## Medições executadas

| Cenário | Tempo | Tabelas moleculares lidas | Linhas lidas |
| --- | ---: | ---: | ---: |
| Leitura isolada anterior: um BBB+ com 3.362 linhas | 11,525 s | 49 | 446.560 |
| Mesma leitura após a correção | 0,088 s | 1 | 3.362 |
| 36 combinações, cache inicialmente vazio | 20,212 s | 12 | 47.536 |
| Mesmas 36 combinações, cache preenchido | 6,588 s | 0 | 0 |

As duas primeiras linhas medem apenas o carregador de registros. As duas últimas medem `_resolve`, incluindo materialização, normalização e validação de parâmetros; não incluem o cálculo de todas as cargas AM1 nem todos os dockings. O contador de leituras considera `docking_data.read_compounds`, não todas as consultas de filesystem ou cabeçalho.

A prévia real da configuração, com **12 arquivos BBB+ e três receptores selecionados**, mediu 10,603 s para preparar dados em thread e 0,201 s para construir os controles Python depois disso. Esse último número não é uma medição completa de renderização/rede. O maior dropdown teve 101 opções, incluindo Todos os compostos. As seleções permaneceram idênticas após redimensionamento.

Foram capturadas telas de 1440, 1024 e 360 pixels; o botão de configuração permaneceu visível. O único evento de console registrado foi Flutter app loaded, sem exceções JavaScript.

## Testes científicos com as ferramentas instaladas

| Origem | Composto | Preparo e Vina | Tempo | Scores das duas poses |
| --- | --- | --- | ---: | --- |
| ADMET BBB+ | CHEMBL58983 | 2 poses, zero falhas | 16,027 s | −8,334; −8,329 |
| Evaluate ADMET BBB+ | PUBCHEM72757 | 2 poses, zero falhas | 7,705 s | −6,590; −6,576 |

O Chimera aplicou o preparo configurado, incluindo AM1 e minimização. O receptor 6ZWI preparado foi reutilizado; os dois sítios existentes nos metadados desse receptor deram origem às duas poses por composto. Os valores são evidência de execução do motor, não validação experimental de afinidade.

Uma terceira amostra, **CHEMBL339995**, falhou no Chimera durante AM1: o log informa 153 elétrons e carga líquida zero, seguido de falha em SQM. A falha foi registrada por composto e não impediu o caso seguinte de Evaluate ADMET. Não foi removida validação, modificada a carga nem declarada inviabilidade intrínseca dessa molécula.

Uma tentativa ampliada com seis compostos foi interrompida manualmente somente no projeto temporário. Antes disso, CHEMBL377041 completou seu preparo em **233,797 s**; o processo SQM seguinte estava ativo a 100% de CPU. Esse teste ampliado não é declarado concluído e sua interrupção não foi classificada como falha científica. Isso demonstra a distinção entre cálculo AM1 caro e o gargalo de organização das entradas corrigido aqui.

**Não foram executados todos os dockings das 47.536 linhas.** O preparo/normalização de todos os arquivos e combinações foi validado, e o percurso científico completo foi exercitado nas duas amostras acima. AM1/minimização em grandes coleções continua tendo custo real significativo.

## Regressões

- `tests.log`: **133 testes direcionados aprovados**, cobrindo seleção/configuração, preservação de poses, cache incremental, persistência, invalidação seletiva, UI e os casos novos de segurança/desempenho.
- `atomic-tests.log`: **7 testes específicos aprovados** na versão final; sobrepõem a suíte anterior. Incluem alteração concorrente da fonte, corrupção do cache, reaproveitamento entre receptores, paginação de 50.000 IDs e cancelamento durante uma validação lenta.
- `vina-tests.log`: **12 testes específicos de Vina/reúso aprovados**, em processo separado para isolar mocks de importação.
- `broader-tests.log`: uma execução ampliada manteve uma falha em `test_docking_interoperability...test_both_engines_accept_actual_poses_and_keep_compound_identity`. O teste usa arquivos literais pose A/pose B e intercepta `docking_data.convert_structure`, enquanto o preparo atual chama o símbolo importado em `docking_preparation`. Esse teste não foi alterado nem declarado aprovado; sua separação consta explicitamente de `run-acceptance.py` e do log final. Os testes nativos acima usaram moléculas reais.

Os avisos RDKit de hidrogênios explícitos/depreciação foram mantidos. Uma primeira execução ampliada também identificou ausência do Open Babel no PATH do processo de teste; a repetição usou o PATH do ambiente BioMolExplorer já instalado. A primeira tentativa do smoke test apontou saídas fora do escopo do projeto temporário e foi corretamente rejeitada; os casos concluídos usam destinos dentro desse projeto.

## Arquivos modificados e continuação

Backend: `pipeline.py`, `docking_data.py` e `docking_inputs.py`. Interface: `ui/compound_choice.py`, `ui/compound_catalog.py`, `ui/file_selection.py`, `ui/input_editor.py` e os catálogos en/pt. Regressões: `tests/test_vina_multiple_inputs.py` e atualização do ponto de interceptação do catálogo em `test_editor_responsiveness.py`.

Carregue a versão atualizada e retome a execução com reutilização dos resultados válidos. A correção não exige apagar resultados anteriores. Os scripts anexos documentam a reprodução local e dependem dos snapshots/cópias temporárias usados nesta validação; os testes unitários são reproduzíveis sem os dados particulares do experimento.
