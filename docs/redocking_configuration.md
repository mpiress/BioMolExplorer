# Configuração do redocking

[Documentação](README.md) · Português · [English](en/redocking_configuration.md)

## Selecionar e preparar os pares

Após recuperar os PDBs, selecione os arquivos de entrada e abra a configuração do redocking. Em **Entradas**, escolha explicitamente cada par PDB / ligante / resíduo. A recuperação inicial não seleciona pares para execução sem uma escolha do usuário. É obrigatório selecionar pelo menos um antes de aplicar a configuração.

Ao revisar uma estrutura em **Resultados → Ligantes → Salvar ligantes**, a escolha explícita também atualiza os pares dos blocos de redocking conectados àquela recuperação. Isso inclui o primeiro ligante de um PDB que estava sem referência configurada. Se a entrada do redocking estava restrita a outros arquivos, a estrutura revisada é acrescentada à seleção. Os pares dos outros PDBs e suas opções de preparação são preservados; para um ligante novo da mesma cadeia, as opções compartilhadas do receptor são mantidas. Salvar novamente a mesma lista não duplica pares.

Depois, use **Run pipeline / Executar pipeline** e **Sim, reaproveitar**: os pares novos são processados e incorporados aos resultados acumulados. A curadoria não inicia cálculos por conta própria. Se a escolha foi salva numa versão anterior e não entrou no redocking, abra **Ligantes** e salve a lista novamente para reconciliar a configuração existente.

Cada seleção mostra a cadeia comum ao receptor e ao ligante, os cofatores e as opções de preparação. A resolução acompanha os metadados do par e não precisa ser digitada. Para PDBs próprios sem metadados, use **Informar par manualmente**.

Ative **Considerar cofatores como parte do receptor** e informe os códigos, separados por vírgula (por exemplo, `FAD, MG`). Sem cofatores, a seleção do complexo usa `select #0:.{chain}`. Os cofatores declarados também são preservados na preparação do receptor, inclusive fora da cadeia selecionada. O ligante de redocking não pode ser declarado como cofator.

Configure remoção de solvente, remoção dos hidrogênios existentes, adição de hidrogênios, minimização e método de cargas para receptor e ligante. A configuração do ligante é compartilhada entre preparação e conformação. O pH continua sendo definido uma vez para o bloco.

Pares que usam o mesmo PDB e a mesma cadeia compartilham o arquivo de receptor preparado; devem usar os mesmos cofatores e as mesmas opções do receptor. As opções de seus ligantes podem ser diferentes. No processamento individual, cada execução recebe somente os pares e as configurações da estrutura correspondente.

Quando os pares já foram configurados, o pipeline usa essa seleção sem solicitar novamente. Antes dos cálculos, valida a presença do receptor, ligante, resíduo, cadeia e cofatores na estrutura. Complexos preparados devem conter os arquivos PDBQT e centros válidos. Chimera, Open Babel e Vina precisam estar instalados no ambiente de execução; com entradas já preparadas, somente Vina é exigido. Verbosidade permanece desativada.

## Configuração pelo backend

Consulte o executável e a mensagem no log de execução. Scripts Chimera que falharam ficam preservados. Centro do ligante inválido ou arquivos preparados incompletos impedem o docking; corrija a configuração ou instalação e tente novamente.

Informe, por exemplo, `pdb_codes=[["4M0E", "1YL", 604, "A", 2.0]]`. `preparation_pairs` usa a chave `4M0E|1YL|604|A`, contendo `cofactors`, `receptor` e `ligand`. As opções permitidas são `remove_solvent`, `remove_hydrogens`, `add_hydrogens`, `minimize` e `charge_type` (`gas` ou `am1`). Opções omitidas usam os padrões. O quinto elemento do registro é a resolução. Projetos legados com `ligand_chain` passam a usar a cadeia comum `CHAIN`.

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

Conecte o redocking a **Receptor PDB (retrieval ou redocking)** do bloco **Preparar para docking**. Escolha o receptor `.dockprep.pdbqt` no campo **Receptor que deseja utilizar**; arquivos exclusivos do ligante não aparecem nessa seleção. O receptor, os formatos complementares e os centros são reutilizados, com as opções de preparo do receptor desabilitadas. Os candidatos vêm de uma ou mais fontes ChEMBL, PubChem, ZINC ou arquivos próprios e continuam sendo preparados. Escolha Vina, DOCK6 ou ambos como saída e conecte receptor e compostos a cada motor. Consulte o [manual do usuário](user_manual.md#preparacao-e-redocking) para processamento individual e merge.

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

Mesmo que nenhum par possa ser preparado, o redocking termina com aviso e não inicia o Vina. `pdb_codes.csv` conserva todos os pares selecionados, incluindo os que falharam, com `STATUS`, `FAILED_STEP` e `ERROR`. O estado `preparation_failed` mantém RMSD vazio. `redocking_failures.csv` reúne os motivos. Esses pares não são oferecidos como complexos preparados válidos às etapas seguintes. Falhas globais de instalação, configuração, acesso ou armazenamento continuam interrompendo a operação.

Ao finalizar, a interface informa a quantidade de falhas. Em **Resultados**, os pares com erro aparecem em vermelho com **Falha na preparação**, RMSD **—** e **Ver erro**. As ações de pose 3D e contatos ficam desativadas porque não houve simulação concluída. Pares bem-sucedidos continuam mostrando o RMSD real.

**Excluir falha e complexo de origem** pede confirmação e remove a linha com erro, o PDB de entrada e sua cópia usada pelo redocking, os registros correspondentes nos metadados da origem e cópias verificadas no cache de entradas. A origem é determinada pelos arquivos de entrada e pela procedência do lote; arquivos de outra coleção com o mesmo nome são preservados. Se o PDB de origem foi alterado desde a execução, a exclusão é recusada até revisar os resultados.

Um arquivo PDB pode servir a vários ligantes. Excluir o arquivo remove todos os seus pares da próxima execução desse bloco, mas preserva os resultados históricos já concluídos e seus arquivos preparados. A confirmação explica esse efeito. Para recuperações PDB, o código entra em **PDBs excluídos de novas recuperações** (`excluded_pdb_ids`), impedindo que a mesma recuperação o baixe novamente. Essa lista pode ser editada no bloco de recuperação para permitir uma nova tentativa intencional. Para arquivos importados, o arquivo e sua referência no bloco de importação são removidos.

A ação exige proprietário/editor e execução encerrada. As alterações de arquivos, configuração e integridade dos resultados são auditadas no histórico; uma falha durante a operação restaura o estado anterior. As cópias históricas para restauração não são entradas ativas do redocking. Nenhuma exclusão é executada automaticamente ao ocorrer uma falha de preparação.
