# Configuração do redocking

[Documentação](README.md) · Português · [English](en/redocking_configuration.md)

## Selecionar e preparar os pares

Após recuperar os PDBs, selecione os arquivos de entrada e abra a configuração do redocking. Em **Entradas**, escolha explicitamente cada par PDB / ligante / resíduo. A recuperação inicial não seleciona pares para execução sem uma escolha do usuário. É obrigatório selecionar pelo menos um antes de aplicar a configuração.

Ao revisar uma estrutura em **Resultados → Ligantes → Salvar ligantes**, a escolha explícita também atualiza os pares dos blocos de redocking conectados àquela recuperação. Isso inclui o primeiro ligante de um PDB que estava sem referência configurada. Se a entrada do redocking estava restrita a outros arquivos, a estrutura revisada é acrescentada à seleção. Os pares dos outros PDBs e suas opções de preparação são preservados; para um ligante novo da mesma cadeia, as opções compartilhadas do receptor são mantidas. Salvar novamente a mesma lista não duplica pares.

Depois, use **Run pipeline / Executar pipeline** e **Sim, reaproveitar**: os pares novos são processados e incorporados aos resultados acumulados. A curadoria não inicia cálculos por conta própria. Se a escolha foi salva numa versão anterior e não entrou no redocking, abra **Ligantes** e salve a lista novamente para reconciliar a configuração existente.

A tabela **Pares selecionados** exibe apenas PDB, Ligante, Número do resíduo, Cadeia e Ações. Use **Editar par** para abrir o formulário de uma linha; cancelar descarta a edição. **Remover par** retira a linha do conjunto, sem apagar o PDB de origem. A resolução acompanha os metadados e não precisa ser digitada. Para PDBs próprios sem metadados, use **Informar par manualmente**.

Expanda **Cofatores do receptor (opcional)** e ative **Considerar cofatores como parte do receptor** e informe os códigos, separados por vírgula (por exemplo, `FAD, MG`). Sem cofatores, a seleção do complexo usa `select #0:.{chain}`. Os cofatores declarados também são preservados na preparação do receptor, inclusive fora da cadeia selecionada. O ligante de redocking não pode ser declarado como cofator.

Configure remoção de solvente, remoção dos hidrogênios existentes, adição de hidrogênios, minimização e método de cargas para receptor e ligante. A seção **Preparação única do conjunto** oferece um painel comum para todos os receptores e outro para todos os ligantes. A preparação e a conformação do ligante usam essas mesmas opções. O pH da conversão PDBQT aparece nas opções Vina. Os painéis ficam lado a lado em telas largas e empilhados em telas estreitas; o campo de cofatores só aparece quando ativado.

Pares que usam o mesmo PDB e a mesma cadeia compartilham o arquivo de receptor preparado. A interface aplica uma preparação comum aos receptores e uma aos ligantes. Projetos antigos com opções diferentes por par continuam executáveis no backend; ao salvar pela interface, um aviso informa que a configuração comum substituirá essas diferenças. No processamento individual, cada execução recebe somente os pares e as configurações da estrutura correspondente.

Quando os pares já foram configurados, o pipeline usa essa seleção sem solicitar novamente. Antes dos cálculos, valida a presença do receptor, ligante, resíduo, cadeia e cofatores na estrutura. Complexos preparados devem conter centros válidos e os formatos do motor escolhido: PDBQT para Vina; MOL2 do receptor e ligante, mais receptor `.noH.pdb`, para DOCK6. A preparação utiliza Chimera clássico e Open Babel; as entradas já preparadas dispensam essas etapas. O motor selecionado deve estar instalado. DOCK6 exige também sphgen, sphere_selector, showbox, grid e seus arquivos de parâmetros na instalação 6.11. Verbosidade permanece desativada.

## Escolher o método

Em **Método de redocking**, escolha **AutoDock Vina** ou **DOCK6 6.11**. Apenas os campos do motor ativo são exibidos e validados. A troca conserva os valores do outro motor para uma próxima edição.

- **Vina:** caixa X/Y/Z em Å, pH da conversão, esforço de busca e quantidade de poses. O centro vem do ligante de referência. **Vina · opções do motor** reúne semente e função de score.
- **DOCK6:** instalação administrada, busca flexível/rígida, densidade SES, raio da sonda, distância de seleção das esferas, máximo de orientações, poses salvas e semente. As subseções expansíveis reúnem opções do motor, grade, caixa e esferas. O fluxo usa diretamente o ligante preparado, sem uma minimização adicional da referência cristalográfica.

```text
Selecione os pares e revise-os na tabela.
Defina uma preparação para receptores e outra para ligantes.
Escolha Vina ou DOCK6 6.11 e ajuste seus campos específicos.
Salve e execute o pipeline; confira RMSD e os diagnósticos.
```

Novas execuções com MOL2 disponível calculam RMSD de átomos pesados, com correção de simetria e sem sobreposição geométrica. A referência usa as coordenadas do ligante no PDB original; na ausência do PDB, o resultado registra explicitamente `RMSD_REFERENCE=prepared`. Entradas Vina legadas contendo apenas PDBQT mantêm o cálculo legado, identificado por `RMSD_METHOD`. Esses métodos não devem ser comparados como se fossem idênticos. `Engine`, `RMSD_method` e `RMSD_reference` acompanham a exportação CSV. RMSD finito não representa aprovação automática do protocolo.

O preparo usa o grafo CCD do componente antes de adicionar hidrogênios e cargas no Chimera. Tipos IDATM são fixados a partir desse grafo para evitar reinterpretação incorreta pela geometria. O MOL2 é escrito sem reabrir o PDB. Após a conversão PDBQT, nomes são preservados por correspondência de coordenadas antes do docking. O arquivo `.export.json` ou `.input.json` registra a definição e eventuais átomos ausentes; fragmentos observados não são completados com átomos pesados inventados. Falhas recuperáveis de um par no Vina/RMSD são registradas sem interromper os demais.

## Configuração pelo backend

Consulte o executável e a mensagem no log de execução. Scripts Chimera que falharam ficam preservados. Centro do ligante inválido ou arquivos preparados incompletos impedem o docking; corrija a configuração ou instalação e tente novamente.

Informe, por exemplo, `pdb_codes=[["4M0E", "1YL", 604, "A", 2.0]]`. `preparation_options` define a preparação comum (`receptor`, `ligand`, `cofactors`) e tem precedência. Para compatibilidade com clientes antigos, `preparation_pairs` usa a chave `4M0E|1YL|604|A`, contendo `cofactors`, `receptor` e `ligand`. As opções permitidas são `remove_solvent`, `remove_hydrogens`, `add_hydrogens`, `minimize` e `charge_type` (`gas` ou `am1`). Opções omitidas usam os padrões. O quinto elemento do registro é a resolução. Projetos legados com `ligand_chain` passam a usar a cadeia comum `CHAIN`.

Veja [validação do pipeline](pipeline_validation.md) para cobertura e limites da execução nativa.

A seleção de um subconjunto não remove PDBs nem arquivos preparados de outras estruturas da coleção. Quando o par não define seu método de cargas, a preparação do ligante conserva o método global do bloco.

## Falhas de preparação com Chimera clássico

Um processo Chimera pode terminar com código zero mesmo quando um comando `.com` falhou. O pipeline verifica essas mensagens e os arquivos esperados antes de avançar. Cada saída deve existir, estar preenchida e conter átomos; o mesmo vale para a conversão pelo Open Babel. No redocking, falhas de preparação de um par são registradas e esse par é ignorado; os demais continuam. O script Chimera é preservado para diagnóstico. Falhas globais, como instalação ausente ou erro de disco, continuam sendo erros da etapa.

Os caminhos de `open` e `write` são passados sem aspas de shell: no Chimera clássico, essas aspas fazem parte do nome do arquivo. Caminhos com espaços foram verificados com o executável real. A saída padrão e os erros das ferramentas são registrados nos logs científicos.

## Consultar resultados de redocking

Ao expandir uma etapa concluída em **Execuções**, a tabela mostra PDB, ligante, resíduo, cadeia e **RMSD (Å)**. **Ver simulação** abre um popup com os arquivos correspondentes àquela seleção. A lista distingue o receptor preparado, o ligante de referência, as poses e os metadados. Arquivos compartilhados pelo receptor ou pela coleção acompanham as simulações correspondentes.

Use o botão de download de cada linha para baixar um arquivo, ou **Baixar todos (ZIP)** para obter o conjunto da simulação. O ZIP conserva as pastas, inclusive quando o ligante de referência e a saída do Vina têm o mesmo nome. Arquivos de outras simulações ficam fora desse conjunto.

O botão **Visualizar estrutura 3D** abre o visualizador no navegador, como na consulta dos PDBs. Ele aceita PDB, PDBQT e MOL2 para examinar receptor, ligante e poses antes do download. Ao abrir um arquivo com várias poses, use **Modelo** no visualizador. A cena aberta pelo botão **3D** da tabela seleciona automaticamente a pose de menor score, ou a primeira sem score; essa cena não oferece troca de modelo. As legendas e ações acompanham o idioma selecionado. Leitores também podem consultar e baixar os resultados.

A tabela resume os valores finitos e não negativos da coluna `RMSD` de `pdb_codes.csv`, com três casas decimais na apresentação; o arquivo conserva a precisão original. O RMSD compara o ligante de referência com o resultado do redocking. A interface não aplica um limiar automático de aprovação científica. A identificação PDB / ligante / resíduo / cadeia permite distinguir cada simulação. Etapas com falha ou ainda em execução mantêm a lista geral de arquivos para diagnóstico.

[Uso do backend](backend_usage.md) inclui um exemplo completo do dicionário de parâmetros e a organização dos arquivos de saída.

## Reutilizar o receptor no docking de candidatos

Conecte o redocking a **Receptor PDB (retrieval ou redocking)** do bloco **Preparar para docking**. Escolha o receptor `.dockprep.pdbqt` (ou `.dockprep.mol2` quando só houver a preparação DOCK6) no campo **Receptor que deseja utilizar**; arquivos exclusivos do ligante não aparecem nessa seleção. O receptor, os formatos complementares e os centros são reutilizados, com as opções de preparo do receptor desabilitadas. Os candidatos vêm de uma ou mais fontes ChEMBL, PubChem, ZINC ou arquivos próprios e continuam sendo preparados. Escolha Vina, DOCK6 ou ambos como saída e conecte receptor e compostos a cada motor. Consulte o [manual do usuário](user_manual.md#preparacao-e-redocking) para processamento individual e merge.

## Estilos e interações em 3D

A cena **3D** compara pose, receptor e referência nas coordenadas originais.
Escolha estilos independentes para ligante e receptor e use **Hidrogênios do
ligante** para mostrar apenas os H existentes no arquivo. Interações classificadas
recebem linhas tracejadas coloridas, filtros e identificação pelo mouse. PDBQT
exige SMILES associado para recuperar topologia; resultados de redocking sem essa
informação continuam oferecendo contatos por distância, com a indisponibilidade
química indicada. Consulte o [guia do visualizador](molecular_viewer.md) para a
legenda, requisitos e limites da classificação.

## Reutilização e duração dos lotes posteriores

A preparação válida do redocking é reutilizada no docking dos candidatos; não é repetida por composto. O Vina mantém mapas por configuração dentro do lote quando a API Python suporta o template. O DOCK6 compartilha superfície, esferas e grid por receptor/sítio. A preparação de cada candidato continua sendo necessária quando não existe entrada preparada compatível.

Não há timeout científico implícito de 300 segundos. Cálculos continuam até conclusão, erro ou cancelamento; timeouts HTTP são separados. Candidatos com problemas individuais são registrados em **Compostos não utilizados**. Consulte [reutilização Vina](user_manual.md#reutilizacao-do-receptor-no-vina), [validação DOCK6](user_manual.md#validacao-e-exclusoes-no-dock6) e [duração no backend](backend_usage.md#duracao-das-execucoes).

## Redocking incremental e resultados acumulados

Ao executar o pipeline, escolha **Sim, reaproveitar** para conservar os cálculos compatíveis já concluídos. Para um ligante novo de um PDB/cadeia já processado, o receptor preparado é reutilizado: somente o ligante é preparado e seu redocking é executado. Um receptor novo recebe preparação completa. Pares já calculados com a mesma configuração não são executados novamente.

A tabela da nova execução reúne os pares novos e os resultados anteriores compatíveis da mesma origem de dados, inclusive receptores e ligantes anteriores que não tenham sido novamente marcados. A origem do PDB anterior precisa continuar disponível: arquivos excluídos não são recuperados do histórico automaticamente. Cada par aparece uma única vez. As poses e os arquivos preparados reutilizados são copiados para os resultados da nova execução, mantendo os downloads e a visualização 3D disponíveis.

A compatibilidade considera o conteúdo do PDB, a cadeia, cofatores, opções de preparação, pH, cargas, templates, parâmetros do Vina e implementação científica. Alterar a estrutura ou a preparação invalida o receptor; alterar os parâmetros do redocking exige novo cálculo da pose. Os arquivos são verificados por SHA-256 antes do reaproveitamento. Resultados nativos anteriores ao índice incremental também podem ser adotados quando sua configuração, arquivos preparados, pose, RMSD e manifesto de integridade estão disponíveis.

Somente execuções do pipeline concluídas com sucesso fornecem dados reutilizáveis. Execuções canceladas, interrompidas, com erro ou ainda aguardando entrada não fornecem resultados ao cache incremental, mesmo quando uma etapa individual terminou. Os arquivos parciais podem permanecer no histórico para diagnóstico, mas são descartados como fonte de novos cálculos. Falhas de preparação e resultados excluídos também não são reutilizados. **Não, executar novamente** ignora o reaproveitamento.

## Falhas parciais e exclusão do complexo de origem

A preparação do redocking isola cada par receptor/ligante. Um receptor preparado é reutilizado pelos outros ligantes da mesma cadeia. Se a preparação desse receptor falhar, os pares que dependem dele são ignorados; falhas de um ligante não excluem os outros ligantes de um receptor válido. Arquivos parciais do par são removidos e os demais complexos continuam. O comportamento estrito das outras funções de preparação permanece disponível.

Mesmo que nenhum par possa ser preparado, o redocking termina com aviso e não inicia o motor escolhido. `pdb_codes.csv` conserva todos os pares selecionados, incluindo os que falharam, com `STATUS`, `FAILED_STEP` e `ERROR`. O estado `preparation_failed` mantém RMSD vazio. `redocking_failures.csv` reúne os motivos. Esses pares não são oferecidos como complexos preparados válidos às etapas seguintes. Falhas globais de instalação, configuração, acesso ou armazenamento continuam interrompendo a operação.

Ao finalizar, a interface informa a quantidade de falhas. Em **Resultados**, os pares com erro aparecem em vermelho com **Falha na preparação**, RMSD **—** e **Ver erro**. As ações de pose 3D e contatos ficam desativadas porque não houve simulação concluída. Pares bem-sucedidos continuam mostrando o RMSD real.

**Excluir falha e complexo de origem** pede confirmação e remove a linha com erro, o PDB de entrada e sua cópia usada pelo redocking, os registros correspondentes nos metadados da origem e cópias verificadas no cache de entradas. A origem é determinada pelos arquivos de entrada e pela procedência do lote; arquivos de outra coleção com o mesmo nome são preservados. Se o PDB de origem foi alterado desde a execução, a exclusão é recusada até revisar os resultados.

Um arquivo PDB pode servir a vários ligantes. Excluir o arquivo remove todos os seus pares da próxima execução desse bloco, mas preserva os resultados históricos já concluídos e seus arquivos preparados. A confirmação explica esse efeito. Para recuperações PDB, o código entra em **PDBs excluídos de novas recuperações** (`excluded_pdb_ids`), impedindo que a mesma recuperação o baixe novamente. Essa lista pode ser editada no bloco de recuperação para permitir uma nova tentativa intencional. Para arquivos importados, o arquivo e sua referência no bloco de importação são removidos.

A ação exige proprietário/editor e execução encerrada. As alterações de arquivos, configuração e integridade dos resultados são auditadas no histórico; uma falha durante a operação restaura o estado anterior. As cópias históricas para restauração não são entradas ativas do redocking. Nenhuma exclusão é executada automaticamente ao ocorrer uma falha de preparação.


## Target e melhores RMSDs por alvo

A tabela de redocking inclui **Target**, obtido do bloco de recuperação que
forneceu o PDB, usando a proveniência salva dos arquivos. Em entradas mescladas,
o PDB é associado à sua origem específica. Metadados importados também podem
informar `Target` ou `TARGET`. Quando a origem não pode ser identificada, a célula
mostra **—**; o nome do bloco de redocking não é usado para inventar uma associação.

Use **Alvo terapêutico** para filtrar e a busca para localizar alvo, PDB ou ligante.
A ordem inicial agrupa por alvo e prioriza os menores RMSDs. Clique em **Target**
ou **RMSD (Å)** para ordenar. **Destacar melhores RMSDs por alvo** colore em verde
o mínimo de cada alvo e seus empates (tolerância numérica de 1e-9 Å); falhas continuam
vermelhas e não participam da comparação. O mínimo considera todos os resultados
do alvo, inclusive os que estão em outra página ou fora da busca textual.

**Baixar tabela RMSD** exporta os pares correspondentes aos filtros, com `Target`
e `Best_for_target`, revalidando o acesso. Esses rótulos não alteram poses, RMSDs,
metadados científicos originais ou cache. Resultados anteriores podem receber a
coluna ao serem abertos, quando sua proveniência está disponível.

```text
Abra os resultados e selecione ACHE no filtro de alvo.
Ative Destacar melhores RMSDs por alvo.
Compare as linhas verdes e abra a simulação desejada.
Baixe a tabela RMSD para guardar a seleção com sua origem.
```

## Integridade do preparo nativo

O preparo de ligantes usa a topologia CCD do componente observado antes de adicionar
hidrogênios no Chimera. Após a protonação, o grafo é validado e sua carga formal é
informada explicitamente à rotina nativa de AM1-BCC ou Gasteiger. A minimização usa
esse conjunto já preparado, sem repetir DockPrep nem substituir as cargas escolhidas.
O relatório `.lig.preparation.json` registra a carga, o método e as validações.

Se o FindHBond do Chimera apresentar o erro conhecido `accPhiPsi`, apenas a colocação
de hidrogênios é repetida com o método geométrico nativo; cargas e minimização não
são desativadas. Aminoácidos livres de cristalização registrados como HETATM não
são tratados como proteína, salvo se explicitamente incluídos como cofatores;
aminoácidos ligados por ligações peptídicas permanecem no receptor.

O RMSD usa somente átomos pesados na referência cristalográfica e na melhor pose,
sem alinhamento que esconda deslocamentos. RMSD alto é um resultado de avaliação,
não uma falha de execução. Erros de topologia, ferramentas ou configuração continuam
visíveis, com o passo e o diagnóstico, sem impedir os demais pares de executar.
