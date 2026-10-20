# Triagem de ligantes PDB para redocking — 18/10/2026

O processo anterior tratava quase todo resíduo heterogêneo não aquoso como
ligante. Isso incluía íons, solventes de cristalização, glicanos e cópias de
componentes. A recuperação agora produz uma seleção inicial reduzida e preserva
as evidências para revisão manual.

## Comportamento implementado

- Novo parâmetro `ligand_selection`: `recommended` por padrão; `all` mantém a
  seleção ampla. O campo aparece no formulário guiado e é validado na CLI/backend.
- Análise local do primeiro modelo, com um índice espacial Biopython por
  estrutura e os limites de download existentes. Não há requisições extras por
  componente, cálculo de afinidade nem execução de docking durante a triagem.
- Candidatos orgânicos precisam ter carbono, seis átomos pesados observados,
  três resíduos proteicos em contato e pelo menos 50% dos átomos pesados a até
  4 Å da proteína. Hidrogênios, águas, ocupação zero e coordenadas não finitas
  não contribuem para os contatos. Esses limites são heurísticos, não critérios
  universais de atividade biológica.
- Solventes/aditivos comuns têm baixa prioridade. Ligações depositadas com a
  proteína e seus ramos orgânicos, incluindo glicanos, impedem recomendação
  automática para o protocolo convencional. A análise não propaga a ligação
  através de íons, evitando confundir pontes de coordenação com ligações
  covalentes entre moléculas orgânicas. O campo `PROTEIN_LINKED` registra a
  ligação anotada; ele não classifica sozinho seu tipo químico.
- Cofatores ficam como alternativas quando outro candidato orgânico é adequado;
  podem ser recomendados quando são a única referência adequada. ATP não é
  excluído por uma regra genérica de cofatores.
- SITE/REMARK 800 dão evidência separada de sítios anotados. Descrições
  explicitamente ativas/catalíticas têm prioridade sobre outros bolsões. Um
  SITE genérico de ligação não se torna automaticamente evidência catalítica.
- Cópias só são agrupadas quando componente, sequência/numeração proteica e
  resíduos contatados são equivalentes. Outros componentes e bolsões não são
  reduzidos a um único vencedor arbitrário.

## Interface e persistência

`pdb_codes.csv` mantém o contrato anterior de cinco colunas e contém os
candidatos selecionados. É esse arquivo que alimenta preparação e redocking.
`ligand_candidates.csv` registra todos os componentes avaliados: nome depositado
em HETNAM, identidade, átomos observados, contatos, cadeias e resíduos contatados,
fração em contato, menor distância de contato, SITE, motivo da decisão e
representante das cópias equivalentes. `retrieval_report.json` inclui a versão
1 da política, limites, contagens e estruturas sem recomendação.

Os arquivos PDB preservam as coordenadas e registros LINK/SSBOND/SITE depositados.
Nenhuma molécula é apagada da estrutura para reduzir a lista. O relatório de
triagem documenta a seleção automática original; a curadoria posterior fica no
CSV de seleção e no histórico do projeto.

O diálogo **Ligantes** oferece filtro para os recomendados, acesso aos demais
resíduos, evidências do candidato selecionado e botão para repor a recomendação.
Salvar continua a exigir permissão de editor, revisão correspondente e ausência
de cálculo ativo, atualizando a integridade dos artefatos. Recuperações antigas
recebem sugestões ao abrir o diálogo sem alterar automaticamente sua curadoria.
A versão do código de triagem integra o digest de implementação usado pelo
reaproveitamento de resultados.

## Verificação com estruturas reais

Os PDBs públicos foram baixados da RCSB para `/tmp` e analisados localmente.
O [registro JSON](pdb_ligand_triage_2026-10-18.json) contém URLs, hashes SHA-256,
parâmetros, contagens, tempo de análise e evidências por componente.

| Estrutura | Componentes não aquosos avaliados | Recomendados | Componentes selecionados |
| --- | ---: | ---: | --- |
| [1ACJ](https://www.rcsb.org/structure/1ACJ) | 1 | 1 | THA, cadeia A |
| [1HVR](https://www.rcsb.org/structure/1HVR) | 1 | 1 | XK2, cadeia A |
| [4EY7](https://www.rcsb.org/structure/4EY7) | 14 | 2 | E20, cadeias A e B |

Em 4EY7, a redução é de aproximadamente 86%; os outros 12 componentes continuam
no relatório e na revisão manual. As cópias E20 foram mantidas porque as cadeias
não satisfizeram a equivalência estrita adotada. Essa amostra confirma os casos
exercitados, não mede sensibilidade/especificidade em todo o PDB.

## Limites científicos e do protocolo

Contato geométrico sustenta a hipótese de um bolso de ligação; não comprova sítio
ativo, afinidade ou seletividade. SITE pode descrever ambientes de ligantes,
sítios regulatórios e outros locais essenciais, conforme a
[especificação wwPDB](https://www.wwpdb.org/documentation/file-format-content/format33/sect7.html).
Sem descrição catalítica explícita, o sistema não anuncia confirmação de sítio
ativo. Aditivos, cofatores e fragmentos pequenos podem ser biologicamente
relevantes; a lista ampla permite rever as decisões.

A recuperação por alvo não garante que cada cadeia de um heterocomplexo seja
esse alvo. O sistema informa as cadeias e os resíduos contatados para revisão.
A análise não reconstrói a montagem biológica nem as cópias cristalográficas.
Em ligantes de interface, como o XK2 de 1HVR, é necessário conferir se o receptor
preparado contém todas as cadeias relevantes. O protocolo atual compartilha a
cadeia do receptor e do ligante: resíduos em outra cadeia, sem cadeia ou com
inserção não são recomendados automaticamente. A curadoria rejeita códigos de
inserção para impedir que sua omissão selecione um resíduo diferente.

Somente o primeiro modelo e as conformações alternativas selecionadas pelo
parser são analisados. Alvos exclusivamente de DNA/RNA não são cobertos pela
triagem proteica; nesses casos use a seleção ampla. Entradas somente mmCIF
continuam fora do contrato atual. Nenhum motor externo foi executado para
validar RMSD ou resultados de redocking nesta revisão.

## Testes e verificações

A suíte ampliada passou com **108 testes**. Após os ajustes finais de campos e
apresentação, passaram também as verificações focadas de recuperação/curadoria
(**39 testes**) e triagem/localização (**28 testes**). `git diff --check` passou.
O executável Open Babel disponível no ambiente científico foi incluído no PATH
para os testes de integração do pipeline; a primeira execução sem esse PATH
havia falhado em um teste de consenso, sem relação com a triagem.

```bash
PATH=/home/michel/progs/miniconda3/envs/BioMolExplorer/bin:$PATH \
PYTHONPATH=src:tests MPLCONFIGDIR=/tmp/biomol-mpl \
/home/michel/progs/miniconda3/envs/BioMolExplorer/bin/python -m unittest \
  test_pdb_ligand_triage test_flexible_retrieval test_retrieval_usability \
  test_redocking_configuration test_redocking_execution_review \
  test_pipeline_execution test_docking_preparation test_input_processing \
  test_localization test_pipeline_editor_regressions -q
```

As regressões cobrem limiares de contato, componentes não orgânicos, solventes,
cofatores, sítios anotados, cópias equivalentes e bolsões distintos, LINK/CONECT,
ramos de glicanos, primeiro modelo, ocupação zero, hidrogênios, inserções,
restrições de cadeia, ausência de recomendação, relatório de recuperação,
seletor da interface, alteração manual auditada e encaminhamento dos metadados.

Os HTMLs de recuperação em português/inglês e do manual em português foram
regenerados com o renderizador do projeto. O validador global de documentação
continua apontando duas âncoras ausentes nas páginas `frontend.html` em
português/inglês; ambas foram verificadas também na revisão HEAD e já existiam
antes desta alteração. As páginas do frontend não foram modificadas.
