# Validação de poses: critérios físicos e geométricos

[Portal de documentos](../index_doc.html) · [English](en/pose_validation.md) · [Parâmetros da API](backend_usage.md#api-pose-validation)

O bloco **Validar poses** avalia a plausibilidade física e geométrica das poses produzidas por AutoDock Vina, DOCK6 ou pelo consenso, antes que elas sigam para o Amber Score, para o consenso ou para um relatório. Ele **não realiza docking, não reposiciona o ligante, não minimiza a pose entregue e não altera scores**: cada pose é medida como foi recebida, e o resultado é um laudo por pose.

Os critérios seguem a formulação de PoseBusters (Buttenschoen, Morris e Deane, *Chemical Science* **15**, 3130–3139, 2024), implementada em `pose_geometry.py` com a identificação `posebusters_inspired_v1`. A implementação é própria, usa RDKit/NumPy/SciPy e não executa o pacote original; por isso o relatório registra a identificação do conjunto de critérios, não uma certificação PoseBusters.

**Leitura rápida:** uma pose aprovada satisfez todos os critérios avaliados. Uma pose reprovada traz os critérios que falharam e as medidas correspondentes. Aprovação significa plausibilidade física da geometria, **não** afinidade, atividade biológica ou acerto do modo de ligação.

## 1. Entradas e escopo

| Entrada | Conteúdo aceito | Uso |
| --- | --- | --- |
| `base_input_path` | Receptores preparados do redocking ou da preparação (`receptor_pdb`, `prepared_structures`, `structures`) | Contexto estérico do teste intermolecular |
| `base_selected_mols` | Poses e tabelas de Vina, DOCK6 ou consenso (`vina`, `dock6`, `molecular_poses`, `scores`) | Poses a validar |

As poses chegam pelos blocos conectados. O editor lista alvos pelo `receptor_id` dos receptores preparados e compostos pelas moléculas que efetivamente produziram resultados de docking, lidos do registro de artefatos e das tabelas de score já publicadas — o formulário não carrega o JSON de execução nem percorre a árvore de poses. Sem seleção explícita, todas as poses recebidas são validadas.

| Seleção | Efeito |
| --- | --- |
| Alvos selecionados (`selected_receptors`) | Restringe a validação aos `receptor_id` escolhidos |
| Ligantes selecionados (`selected_ligands`) | Restringe aos compostos escolhidos, por `compound_id` ou `ligand_id` |
| Complexos removidos (`excluded_pairs`) | Remove um complexo específico na forma `composto:receptor`, sem excluir o composto nos outros alvos |

Uma execução isolada pela CLI resolve as poses a partir da pasta informada em `base_selected_mols` ou `base_input_path`; as mesmas exclusões são aplicadas sobre as poses resolvidas.

## 2. Preparação da molécula medida

Cada arquivo de pose é lido com as ordens de ligação disponíveis:

| Formato | Leitura | Ordens de ligação |
| --- | --- | --- |
| `.sdf` | `SDMolSupplier` sem sanitização automática | Explícitas |
| `.mol2` | `MolFromMol2File` sem sanitização automática | Explícitas |
| `.pdbqt` | Conversão interna para PDB, modelo selecionável | Por proximidade |

A distinção importa: sem ordens de ligação explícitas, toda ligação seria lida como simples, e o campo de força e o teste de planaridade mediriam uma molécula diferente da real. Por isso os critérios de geometria interna, planaridade e energia **só são avaliados quando as ordens de ligação não vêm de proximidade**. O campo `bond_orders` do relatório registra qual leitura foi usada, e um critério sem medida fica registrado como indisponível, em vez de reprovar a pose silenciosamente.

## 3. Critérios avaliados

| Critério | Medida | Reprovação |
| --- | --- | --- |
| Química (`valid_chemistry`) | Sanitização RDKit da molécula lida | Valência ou aromaticidade inválida |
| Conectividade (`valid_connectivity`) | Número de fragmentos | Mais de um fragmento desconexo |
| Coordenadas | Finitude e unicidade das posições | Coordenadas ausentes, repetidas ou não finitas |
| Geometria interna (`valid_geometry`) | Maior desvio relativo de comprimento de ligação e maior desvio absoluto de ângulo contra o confôrmero relaxado | Desvio acima de `bond_tolerance` ou `angle_tolerance` |
| Planaridade (`valid_planarity`) | Maior distância ao plano ajustado em anéis aromáticos e ligações duplas fora de anel | Desvio acima de `planarity_tolerance` |
| Energia interna (`valid_energy`) | ΔE entre a pose ancorada e o mesmo ligante minimizado | ΔE acima de `delta_energy_limit` |
| Choque estérico (`valid_clash`) | Menor razão entre a distância observada e a soma dos raios de van der Waals, entre átomos pesados da pose e do receptor | Razão abaixo de `vdw_clash_scale` |

### 3.1 Referência interna relaxada

A referência geométrica e energética é uma **cópia da mesma conectividade**, minimizada com MMFF94 e, quando o campo não parametriza a molécula, com UFF. O método efetivo é registrado em `energy_method`. A pose original nunca é modificada: a minimização ocorre sobre a cópia.

ΔE = E(pose) − E(confôrmero relaxado), em kcal/mol, com o mesmo campo de força nas duas avaliações. ΔE é uma **medida de tensão interna da conformação entregue**, não energia de ligação, não energia livre e não termo de interação com o receptor.

### 3.2 Desvios internos

O desvio de ligação é relativo: |d_pose − d_relaxado| / d_relaxado, tomando o maior valor da molécula. O padrão 0,25 aceita desvios de até 25% por ligação. O desvio de ângulo é absoluto, em graus, com padrão 25°. Ambos comparam a pose com o confôrmero relaxado da mesma molécula, não com valores tabelados.

### 3.3 Planaridade

Para cada anel aromático e cada ligação dupla fora de anel, ajusta-se o plano do grupo e mede-se a maior distância dos átomos a esse plano. O número de grupos avaliados aparece em `planar_groups`.

### 3.4 Choque estérico intermolecular

O teste usa raios de van der Waals de Bondi (elementos não listados usam o raio do carbono, 1,70 Å). Solvente é descartado, e o ligante cristalográfico do próprio sítio é excluído por `exclude_residues`, para que a referência cristalográfica não seja contada como choque. A decisão usa a menor razão distância/soma de raios do conjunto; o relatório traz também a sobreposição absoluta do pior par (`min_vdw_overlap`), a menor distância entre átomos pesados (`min_heavy_distance`) e quantos pares ficaram abaixo da soma dos raios (`clashing_contacts`).

O padrão 0,75 preserva ligações de hidrogênio e contatos polares próximos, como na formulação de distância intermolecular de PoseBusters. **Sem receptor associado, o critério não é avaliado**: a pose não é reprovada por isso, e o relatório registra a ausência em `resolved_by` e no resumo (`unresolved_inputs`).

## 4. Tolerâncias e padrões

| Parâmetro | Padrão | Domínio aceito |
| --- | --- | --- |
| `strict_filter` | `true` | booleano |
| `vdw_clash_scale` | `0.75` | 0,3 a 1,0 |
| `delta_energy_limit` | `100.0` | 0 a 100000 kcal/mol |
| `bond_tolerance` | `0.25` | 0 a 2 (fração) |
| `angle_tolerance` | `25.0` | 0 a 180 graus |
| `planarity_tolerance` | `0.25` | 0 a 5 Å |
| `minimization_steps` | `400` | 1 a 20000 iterações |

Com `strict_filter` ativado, somente as poses aprovadas são entregues aos blocos seguintes. Desativado, todas as poses seguem disponíveis com os metadados da validação, e a decisão fica com quem interpreta o relatório. As tolerâncias aplicadas são gravadas no resumo da execução, de modo que um relatório antigo continua interpretável com os limites que o produziram.

## 5. Resultados publicados

| Arquivo | Conteúdo |
| --- | --- |
| `pose_validation.csv` | Uma linha por pose avaliada, com identificação, motor, score de origem, verdictos por critério, medidas, motivos de reprovação e arquivos de origem |
| `pose_validation.json` | Resumo: alvo, total de poses, aprovadas, reprovadas, publicadas, filtro estrito, contagem por motivo, tolerâncias aplicadas, entradas não resolvidas e identificação dos critérios |
| `poses/` | Cópias das poses entregues aos blocos seguintes |
| `receptors/` | Receptores usados no teste estérico, em PDB |

As portas de saída são `validated_poses` (`molecular_poses`, `vina`, `dock6`), `validation_report` (`scores`) e `receptors` (`receptor_pdb`). Cada pose publicada mantém um arquivo próprio: quando o motor nomeia a pose apenas pelo composto — como o DOCK6 faz — o nome da cópia é qualificado pelo receptor, para que três alvos do mesmo composto não colapsem em um único arquivo.

Ordenação do relatório: reprovadas primeiro, depois por alvo, composto e pose. A ordem existe para colocar na frente as poses que exigem decisão.

## 6. Colunas do relatório e sua origem no trabalho de referência

### 6.1 Testes definidos pelo trabalho

PoseBusters organiza a validade em três categorias. Os limiares abaixo são os do trabalho, não os desta implementação.

| Categoria | Teste | Definição e limiar no trabalho |
| --- | --- | --- |
| Química | Carregamento do arquivo | O RDKit consegue ler a molécula |
| Química | Sanitização | A molécula passa na sanitização do RDKit: valência, aromaticidade, radicais e hibridização |
| Química | Fórmula molecular | Igual à do ligante verdadeiro, comparada por InChI padrão, sem isótopos e com cargas neutralizadas quando possível |
| Química | Ligações | Iguais às do ligante verdadeiro, pela camada de conectividade do InChI |
| Química | Quiralidade tetraédrica | Igual à do ligante verdadeiro, pela camada estéreo do InChI |
| Química | Estereoquímica de ligações duplas | Igual à do ligante verdadeiro, pela camada de ligações duplas do InChI |
| Intramolecular | Comprimentos de ligação | Entre 0,75 do limite inferior e 1,25 do limite superior da geometria de distâncias do RDKit (tolerância de 25%) |
| Intramolecular | Ângulos de ligação | Mesma tolerância de 25% sobre os limites de geometria de distâncias |
| Intramolecular | Choque estérico interno | Pares não ligados covalentemente acima de 0,8 do limite inferior de distância |
| Intramolecular | Planaridade de anéis aromáticos | Anéis aromáticos de 5 e 6 membros com todos os átomos a até 0,25 Å do plano comum mais próximo |
| Intramolecular | Planaridade de ligações duplas | Ligações C=C não cíclicas e não aromáticas: os dois carbonos e seus quatro vizinhos a até 0,25 Å do plano |
| Intramolecular | Energia interna | Conjunto de 50 confôrmeros gerados por ETKDGv3 e relaxados com UFF (até 200 iterações); a razão entre a energia UFF da pose e a média do conjunto não pode passar de 100 — limite menos estrito que a razão 7 de Wills et al. |
| Intermolecular | Distância mínima proteína–ligante | Distância entre átomos pesados acima de 0,75 vez a soma dos raios de van der Waals |
| Intermolecular | Distância mínima a cofatores orgânicos | Acima de 0,75 vez a soma dos raios de van der Waals |
| Intermolecular | Distância mínima a cofatores inorgânicos | Acima de 0,75 vez a soma dos raios covalentes |
| Intermolecular | Sobreposição de volume | Abaixo de 7,5% com a proteína e com cofatores orgânicos (raios de van der Waals escalados por 0,8) e abaixo de 7,5% com cofatores inorgânicos (escala 0,5), pelo ShapeTverskyIndex do RDKit |

Uma pose que passa em todos os testes é **PB-valid**. A exatidão é medida separadamente: RMSD mínimo de átomos pesados, consciente de simetria (`GetBestRMS`), de até 2 Å ao ligante cristalográfico mais próximo. Validade e exatidão são critérios distintos — uma pose pode ser válida e estar errada, ou estar próxima da referência e ser fisicamente inválida.

### 6.2 Colunas de `pose_validation.csv`

| Coluna | O que registra aqui | Correspondência no trabalho |
| --- | --- | --- |
| `pose_id` | Identificador da pose, do produtor ou do nome do arquivo | Identificação; não é teste |
| `receptor_id` | Alvo da pose, pelo receptor preparado | Identificação |
| `ligand_id` | Ligante informado pela origem | Identificação |
| `compound_id` | Composto avaliado | Identificação |
| `engine` | Motor que produziu a pose: `vina`, `dock6` ou consenso | Identificação |
| `score` | Score bruto do motor de origem, preservado sem transformação | Nenhuma; o trabalho não usa o score para julgar validade |
| `pose_valid` | Aprovação em todos os critérios efetivamente avaliados | Equivale a **PB-valid** restrito ao subconjunto implementado |
| `valid_chemistry` | Sanitização RDKit da molécula lida | Teste de sanitização |
| `valid_connectivity` | Molécula em um único fragmento | Conectividade: a consistência de ligações do trabalho exige o ligante verdadeiro e não é reproduzida |
| `valid_geometry` | Comprimentos e ângulos dentro das tolerâncias configuradas | Comprimentos e ângulos de ligação, com referência diferente (ver 6.4) |
| `valid_planarity` | Planaridade de anéis aromáticos e ligações duplas | Planaridade de anéis aromáticos e de ligações duplas |
| `valid_energy` | Tensão interna dentro do limite configurado | Energia interna, com formulação diferente (ver 6.4) |
| `valid_clash` | Contato intermolecular acima da fração configurada dos raios | Distância mínima à proteína e a cofatores orgânicos |
| `min_vdw_ratio` | Menor razão distância/soma dos raios de van der Waals; é a métrica de decisão | Mesma grandeza do teste de distância mínima, cujo limite no trabalho é 0,75 |
| `min_vdw_overlap` | Sobreposição absoluta do pior par, em Å; negativo indica contato mais próximo que a soma dos raios | Medida auxiliar; não é um teste do trabalho |
| `min_heavy_distance` | Menor distância entre átomos pesados da pose e do receptor, em Å | Grandeza bruta do teste de distância mínima |
| `clashing_contacts` | Quantos pares ficaram abaixo da soma dos raios | Medida auxiliar; o trabalho decide pelo pior par e pela sobreposição de volume |
| `delta_energy` | ΔE entre a pose e o confôrmero relaxado, em kcal/mol | Teste de energia interna, em formulação diferente (ver 6.4) |
| `energy_method` | Campo de força efetivo: `MMFF94` ou `UFF` | O trabalho usa UFF |
| `max_bond_deviation` | Maior desvio relativo de comprimento de ligação | Teste de comprimentos de ligação |
| `max_angle_deviation` | Maior desvio absoluto de ângulo, em graus | Teste de ângulos de ligação |
| `max_planarity_deviation` | Maior distância ao plano ajustado, em Å | Testes de planaridade, com a mesma tolerância de 0,25 Å |
| `planar_groups` | Quantos grupos planares foram avaliados | Medida auxiliar de cobertura do teste |
| `atoms` | Átomos lidos na pose | Contexto da medida |
| `bond_orders` | Origem das ordens de ligação: `explicit` ou `proximity` | Condição de aplicabilidade dos testes intramoleculares |
| `resolved_by` | Como o receptor foi associado: `binding`, `identifier` ou `unresolved` | Condição de aplicabilidade do teste intermolecular |
| `rejection_reasons` | Chaves dos critérios reprovados, separadas por `;` | Conjunto de testes falhados |
| `rejection_details` | Descrição legível dos mesmos motivos | Apresentação |
| `error` | Diagnóstico de leitura ou de avaliação, quando houver | Equivale à falha do carregamento do arquivo |
| `pose_file` | Arquivo de pose avaliado | Proveniência |
| `receptor_file` | Receptor usado no teste intermolecular | Proveniência |
| `source_table` | Tabela de score que originou a pose | Proveniência |

### 6.3 Colunas da tabela exibida

| Coluna | Conteúdo | Origem |
| --- | --- | --- |
| Pose | `pose_id` | Identificação |
| Alvo | `receptor_id` | Identificação |
| Composto | `compound_id`, com `ligand_id` como alternativa | Identificação |
| Motor | `engine` em maiúsculas | Identificação |
| Score | `score` do motor de origem | Não participa da validade |
| Situação | `pose_valid` como aprovada ou reprovada | PB-valid do subconjunto implementado |
| Critérios reprovados | Rótulos dos critérios com verdicto falso | Testes falhados |
| Contato mín. · fração vdW | `min_vdw_ratio` | Distância mínima intermolecular; limite 0,75 no trabalho |
| ΔE interna | `delta_energy` | Energia interna, em formulação diferente |
| Ações | Laudo por pose e complexo receptor–pose em 3D | Apresentação |

O laudo de cada linha repete os seis critérios com seus verdictos, todas as medidas com unidade e as tolerâncias aplicadas naquela execução.

### 6.4 Divergências declaradas em relação ao trabalho

| Ponto | Trabalho | Esta implementação |
| --- | --- | --- |
| Referência de comprimentos e ângulos | Limites de geometria de distâncias do RDKit | Confôrmero da mesma molécula relaxado com MMFF94 ou UFF; a tolerância relativa de 25% coincide, a referência não |
| Energia interna | Razão entre a energia UFF da pose e a média de 50 confôrmeros ETKDGv3, limite 100 | Diferença ΔE contra um único confôrmero relaxado, limite padrão de 100 kcal/mol; é energia, não razão adimensional |
| Choque estérico interno | Pares não ligados acima de 0,8 do limite inferior | Não avaliado |
| Sobreposição de volume | Abaixo de 7,5%, com escalas 0,8 e 0,5 | Não avaliada; a decisão intermolecular usa apenas a distância |
| Consistência de identidade | Fórmula, ligações, quiralidade e estereoquímica contra o ligante verdadeiro, por InChI | Não avaliada: uma triagem prospectiva não tem ligante verdadeiro. O SMILES do composto, quando disponível, serve apenas de molde para atribuir ordens de ligação |
| Cofatores e águas | Testes separados para cofatores orgânicos, inorgânicos e águas, com raios covalentes para íons | Um único conjunto de átomos do receptor com raios de van der Waals, descartando solvente e o ligante cristalográfico do sítio |
| Planaridade de anéis | Anéis aromáticos de 5 e 6 membros | Todos os anéis aromáticos, sem restrição de tamanho |
| Exatidão | RMSD ≤ 2 Å ao ligante cristalográfico | Fora do escopo deste bloco; o redocking calcula o RMSD contra a referência cristalográfica |

Por essas diferenças, o relatório identifica os critérios como `posebusters_inspired_v1` e não afirma equivalência numérica com o PB-valid do trabalho. Um resultado desta implementação não deve ser publicado como taxa de PB-valid.

## 7. Apresentação dos resultados

Na aba **Execuções**, o bloco apresenta a tabela de validação como seu resultado:

- Cartões de resumo com poses avaliadas, aprovadas, reprovadas e entregues aos blocos seguintes, mais a contagem por motivo de reprovação e a identificação dos critérios.
- Tabela paginada com pose, alvo, composto, motor, score, situação, critérios reprovados, contato mínimo (fração dos raios de van der Waals) e ΔE interna, com filtro entre todas, aprovadas e reprovadas.
- Por linha, o laudo completo — os seis critérios, cada medida com sua unidade e as tolerâncias aplicadas — e a abertura do complexo receptor–pose em 3D, usando as cópias que o próprio bloco publicou.
- Abaixo da tabela, recolhido, o conjunto de arquivos do bloco para download, pré-visualização ou remoção.

## 8. Limites de interpretação

- Aprovação indica plausibilidade física da geometria entregue. Não indica afinidade, atividade, seletividade nem que o modo de ligação está correto.
- ΔE é tensão interna conformacional, no campo de força usado. Não é energia de ligação nem energia livre; valores de moléculas diferentes não são diretamente comparáveis.
- O choque estérico depende do receptor associado e dos resíduos excluídos. Sem receptor, esse critério fica sem medida.
- Poses lidas por proximidade (PDBQT) não recebem os critérios que dependem de ordens de ligação. O relatório registra essa condição em `bond_orders`.
- As tolerâncias são configuráveis: um conjunto de resultados só é comparável com outro quando as tolerâncias e os critérios registrados coincidem.
- O bloco não corrige poses. Reprovar uma pose não a reposiciona nem a minimiza para o próximo bloco.

## 9. Verificação

```bash
PYTHONPATH=src:tests python -m unittest test_pose_validation test_pose_validation_flow
```

`tests/test_pose_validation.py` cobre leitura por formato, cada critério com casos aprovados e reprovados, tolerâncias, filtro estrito, identidade das poses publicadas e os parâmetros declarados pelo bloco. `tests/test_pose_validation_flow.py` cobre o catálogo de alvos e compostos, o formulário, as conexões no canvas, a resolução das poses na execução, a recusa identificada de um bloco não configurado e a apresentação da tabela de resultados com seus arquivos.

## 10. Referências

- Buttenschoen, M., Morris, G. M., Deane, C. M. (2024). PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences. *Chemical Science*, 15, 3130–3139. DOI 10.1039/D3SC04185A. Os testes, categorias e limiares citados na seção 6 vêm deste trabalho.
- Bondi, A. (1964). van der Waals Volumes and Radii. *The Journal of Physical Chemistry*, 68(3), 441–451.
- Halgren, T. A. (1996). Merck molecular force field. I. Basis, form, scope, parameterization, and performance of MMFF94. *Journal of Computational Chemistry*, 17(5–6), 490–519.
- Rappé, A. K. et al. (1992). UFF, a full periodic table force field for molecular mechanics and molecular dynamics simulations. *Journal of the American Chemical Society*, 114(25), 10024–10035.
- Landrum, G. RDKit: Open-source cheminformatics. https://www.rdkit.org
