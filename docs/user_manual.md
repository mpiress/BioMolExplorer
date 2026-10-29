# Manual do usuário: do primeiro acesso aos resultados

[Documentação](README.md) · Português · [English](en/user_manual.md)

Este manual acompanha um estudo desde a escolha do idioma até a exportação do projeto. Os nomes dos botões abaixo correspondem à interface em português. Para conhecer parâmetros técnicos e APIs, consulte [Workspace Flet](frontend.md), [uso do backend](backend_usage.md), [projetos e versões](projects.md) e [validação do pipeline](pipeline_validation.md).

[Exemplos rápidos por funcionalidade](#exemplos-rapidos-por-funcionalidade)

## 1. Preparar o ambiente e iniciar

Baixe o [projeto no GitHub](https://github.com/mpiress/BioMolExplorer) e entre na pasta do código:

```bash
git clone https://github.com/mpiress/BioMolExplorer.git
cd BioMolExplorer
```

Sem Git, use **Code → Download ZIP**, extraia o arquivo e abra um terminal na pasta extraída. O [guia de instalação e configuração](installation.md) detalha o download, as ferramentas externas e os caminhos.

**Instale UCSF Chimera 1.17 e DOCK6 6.11 antes de iniciar a aplicação.** Eles precisam estar disponíveis no computador dos cálculos, inclusive em modo web. A instalação da interface e do ambiente Conda não instala essas ferramentas; sem elas, as etapas que dependem delas falharão. Conclua as verificações do guia de instalação antes de executar os comandos de inicialização.

1. Na raiz do código baixado, crie o ambiente científico.
2. Ative o ambiente e instale a interface.
3. Configure o `PATH` das ferramentas externas e substitua `/caminho/para/dock6-6.11` pela raiz real do DOCK6, que contém `bin/` e `parameters/`.
4. Inicie o modo navegador ou desktop após conferir os executáveis.

```bash
conda env create -f requirements.yml
conda activate BioMolExplorer
python -m pip install -e '.[ui]'
biomolexplorer-ui --web --language pt --dock6-path /caminho/para/dock6-6.11
```

No navegador, abra `http://127.0.0.1:8550`. Para desktop, omita `--web`. Sem `--language`, a aplicação inicia em inglês. Use `--language en` para escolher explicitamente inglês; essa opção também funciona no desktop.

Em um servidor, use `--web --no-browser --port 8550`. `--data-dir /caminho/workspace` escolhe o armazenamento das contas e dos índices de projetos; mantenha o mesmo caminho nas próximas inicializações. `--worker-python /caminho/python` seleciona outro Python para os cálculos. `--dock6-path /caminho/dock6` informa a instalação DOCK6.

Fingerprints, similaridade, grafos e ADMET usam o ambiente científico local. Recuperação consulta provedores externos. Preparação e docking precisam dos motores correspondentes: Chimera, Open Babel, Vina e, para DOCK6, os executáveis e recursos desse protocolo. Instalar a interface não instala esses motores. O guia do backend explica as dependências; valide o protocolo de docking com um caso de referência antes de usá-lo no estudo.

## 2. Escolher o idioma, criar uma conta e entrar

1. Na barra superior, clique na bandeira no canto direito para abrir **Selecionar idioma**.
2. Escolha a bandeira do Brasil para português ou dos Estados Unidos para inglês. A tela é atualizada imediatamente; campos já preenchidos são preservados.
3. Para o primeiro acesso, clique em **Ainda não tenho uma conta**.
4. Preencha nome, e-mail, senha e confirmação. A senha precisa ter pelo menos dez caracteres; a confirmação deve ser igual à senha.
5. Clique em **Criar conta**. Depois, use o mesmo e-mail e senha em **Entrar no workspace**.

O idioma escolhido vale para essa sessão da interface: menus, formulários, seleções, acompanhamento e mensagens conhecidas usam esse idioma. Nomes fornecidos pelo usuário, arquivos, códigos moleculares, SMILES, parâmetros científicos e chaves JSON não são traduzidos. Logs de motores externos conservam a linguagem produzida pelo motor. Para mudar de idioma durante o uso, utilize o mesmo menu de bandeiras na barra superior, sem sair da conta. A próxima abertura começa com o idioma definido no comando de inicialização.

As contas pertencem ao workspace utilizado na inicialização. Não há recuperação de senha por e-mail nesta versão. **Sair** encerra a sessão; uma sessão expirada retorna à autenticação.

## 3. Criar um projeto e escolher sua pasta

1. No workspace, clique em **Novo projeto**.
2. Informe um nome que identifique o estudo; acrescente descrição e escolha a cor do card na paleta visual. O campo Tags foi removido.
3. Em **Pasta do projeto**, clique no botão de pasta dentro do campo e selecione a pasta principal. Preencha o nome antes de abrir o seletor: o projeto ficará numa subpasta com esse nome. Uma subpasta já existente exige confirmação de substituição; os dados são removidos somente ao salvar. O campo não permite digitação manual. No desktop e no navegador, o botão abre a navegação de pastas do computador que executa o BioMolExplorer. Pastas ocultas não aparecem. **Criar pasta** solicita um nome e cria e seleciona a pasta dentro do local atual.
4. Confirme a criação e abra o projeto.

Exemplo Linux: `/home/pesquisador/estudos/enzima-a`. No navegador, esse caminho pertence ao computador do backend. Um caminho do seu notebook só funciona se o backend também estiver no notebook ou tiver acesso à pasta. A aplicação precisa ter permissão de escrita. Pastas de outros projetos, pastas sobrepostas e o diretório das contas não podem ser reutilizados como pasta de um novo projeto.

A pasta escolhida guarda `project.json`, arquivos em `assets/`, execuções em `runs/`, versões em `.history/` e as cópias de trabalho do projeto. O banco central conserva contas, sessões, permissões e índices; logs gerais de infraestrutura ficam no diretório de logs configurado. Resultados e logs associados às execuções podem ser consultados pelo projeto. Reabrir a aplicação com o mesmo workspace conserva a associação à pasta.

Use busca por nome ou descrição para localizar um projeto. Arquivar retira o projeto da lista ativa; restaurar permite voltar a usá-lo. Excluir exige confirmação e remove permanentemente a pasta do projeto, incluindo entradas, resultados e histórico. Exporte um backup antes; a exclusão não pode ser desfeita pelo histórico.

## 4. Conhecer as áreas do projeto

| Área | Quando usar |
| --- | --- |
| Pipeline | Adicionar, conectar, configurar e executar os blocos |
| Arquivos | Enviar dados próprios, consultar entradas e baixar arquivos |
| Execuções | Acompanhar estados, abrir resultados, examinar erros e logs |
| Compartilhar | Convidar colaboradores e revisar permissões |
| Histórico | Consultar alterações e, como proprietário, restaurar versões |

O pipeline define as dependências. A posição de um bloco na tela ajuda a organização visual; ela não determina a ordem dos cálculos. Duas ramificações independentes podem executar mesmo quando outra ramificação falha. Descendentes de uma etapa com falha ficam impedidos até que sua entrada esteja disponível.

## 5. Enviar e conferir arquivos próprios

1. Adicione **Importar meus arquivos** ao pipeline e abra sua configuração.
2. Escolha **Tipo dos arquivos a adicionar** e clique em **Selecionar arquivos no disco**. Selecione todos os arquivos necessários; cada arquivo é validado antes de entrar no projeto.
3. Para aproveitar arquivos existentes, use **Arquivo já enviado ao projeto** ou **Adicionar todos os arquivos do projeto**.
4. Confira a tabela **Tipo / Arquivo / Ações**. Você pode alterar o tipo por linha (com nova validação) ou usar o botão **Remover** na coluna **Ações**. A remoção dessa lista preserva o arquivo na biblioteca do projeto.
5. Para receptores preparados, inclua também os metadados e formatos complementares. Configure a pasta do alvo para estruturas e aplique a configuração. Um bloco pode publicar vários tipos, conectados às entradas correspondentes dos consumidores.

Um CSV mínimo de compostos é:

```csv
name,smiles
ETHANOL,CCO
BENZENE,c1ccccc1
PYRIDINE,c1ccncc1
```

Também são aceitos `molecule_chembl_id,canonical_smiles`. Códigos ausentes são gerados a partir da estrutura. Use códigos estáveis e seguros para nomes de arquivos, sem barras ou caminhos. Um mesmo código não deve representar moléculas diferentes. Use vírgula como separador de campos e ponto como separador decimal; salve em UTF-8. Valores que contêm vírgulas, como uma lista de fingerprint, precisam estar entre aspas no CSV.

| Tipo | Conteúdo necessário |
| --- | --- |
| Compostos | CSV com `canonical_smiles` ou `smiles`; código recomendado |
| Fingerprints | CSV com `molecule_chembl_id,fingerprint`; listas de bits 0/1 de igual tamanho |
| Similaridades | CSV com `source,target,value`; códigos e valores finitos entre 0 e 1 |
| Complexos PDB | `.pdb` com ATOM/HETATM e coordenadas válidas; registros de ligante/cadeia no formulário ou `pdb_codes.csv` |
| Receptores preparados | Receptores, ligantes e metadados completos compatíveis com o protocolo |
| Vina | `.pdbqt` com coordenadas e `REMARK VINA RESULT` |
| DOCK6 | `*_scored.mol2` com seções MOLECULE/ATOM e `Grid_Score` numérico |
| Scores | CSV com código (`molecule`, `molecule_chembl_id` ou `id`) e scores numéricos |
| Lista de downloads ZINC | Arquivo TXT/URI ou script exportado pelo navegador de tranches, com links SMI/MOL2, inclusive comprimidos |
| Outros | Arquivos de dados gerais; selecione um tipo específico para validar padrões científicos |

Para receptores preparados, preserve `<PDB>_<CHAIN>.dockprep.pdbqt`, `pdb_codes.csv` com `PDB_CODE,LIGAND,RESNUM,CHAIN` e `centers.csv` com três linhas de coordenadas por complexo. O protocolo DOCK6 também utiliza `.dockprep.mol2` e `.noH.pdb`; redocking utiliza os ligantes de referência. Os arquivos auxiliares produzidos pela preparação acompanham automaticamente os receptores selecionados nas conexões do pipeline. Na importação externa, envie o conjunto completo.

## 6. Montar e configurar o pipeline

1. Abra **Pipeline**. Use um modelo pronto ou clique no título de uma categoria da biblioteca para expandir seus blocos. As categorias começam recolhidas; a busca abre as categorias com resultados.
2. Arraste o bloco para a área quadriculada ou clique em **+**. Arraste seu cabeçalho para reposicioná-lo.
3. Ligue a saída da origem à entrada do consumidor. Também é possível clicar primeiro na saída e depois na entrada.
4. Use o nome da entrada para distinguir receptores, compostos, fingerprints e poses.
5. Abra a configuração com dois cliques no cabeçalho ou **Configurar seleção**.
6. No modo **Visual**, revise o nome, a ativação, os parâmetros, as entradas e os filtros. **Adicionar entrada** permite várias origens no mesmo campo.
7. Clique em **Aplicar configuração**. **Cancelar** conserva a configuração anterior.
8. Use **Salvar** para confirmar o rascunho quando necessário; alterações comuns também são salvas após uma breve pausa.

O ícone **ⓘ** de cada bloco, na biblioteca e no canvas, abre uma explicação breve com as entradas e saídas aceitas. Conexões incompatíveis exibem apenas “Dados incompatíveis. Escolha uma saída do tipo esperado por esta entrada.” e conservam o pipeline. Para similaridade, use **ChEMBL / PubChem / ZINC → Gerar fingerprints → Calcular similaridade**.

Entradas incompatíveis, ciclos e dependências de blocos desativados são recusados. **Avançado** permite editar JSON e templates; use-o quando precisar de opções que não estão no formulário. Preserve os nomes técnicos dos parâmetros e os placeholders exigidos pelos templates. Um erro de validação deve ser corrigido antes de aplicar.

**Ctrl+Z** desfaz, **Ctrl+Y** refaz e **Delete** remove a seleção quando nenhum formulário de configuração está aberto. **Escape** cancela uma conexão em andamento. Duplicar um bloco permite comparar parâmetros; renomeie as alternativas. Organização automática, navegação, zoom e recentralização ajudam a explorar fluxos extensos.

## 7. Entender as conexões disponíveis

| Consumidor | Entrada permitida no canvas |
| --- | --- |
| Expandir similares | Recuperação ChEMBL com os downloads originais; CSV genérico de compostos não substitui essa origem |
| Recuperar ZINC | Lista de downloads do tipo `zinc_urls` ou `other`, com links das tranches |
| Preparar para docking | PDB bruto ou receptor do redocking, e uma ou mais fontes de compostos |
| Avaliar ADMET / Gerar fingerprints | Compostos de recuperação, importação, ADMET ou grafos |
| Calcular similaridade | Fingerprints compatíveis; não misture tipos ou larguras diferentes |
| Filtrar por grafos | Saída de **Calcular similaridade**; CSV externo pode ser escolhido no próprio formulário |
| Executar redocking | PDBs brutos com preparação ativada; receptores preparados com preparação desativada |
| Docking com Vina | Receptores preparados e compostos selecionados |
| Docking com DOCK6 | Receptores preparados e compostos ou poses selecionados; poses Vina opcionais |
| Consenso de docking | Poses Vina e resultados DOCK6 com os mesmos identificadores |

Recuperar compostos e recuperar estruturas podem iniciar ramificações independentes. **Importar meus arquivos** publica os tipos dos arquivos selecionados. O modo **resultados prontos** de um bloco valida saídas já existentes e dispensa seu cálculo; isso não transforma qualquer CSV no resultado daquele bloco. Consulte a [matriz e os testes de validação](pipeline_validation.md) para as condições adicionais.

## 8. Executar, selecionar arquivos e escolher individual ou merge

1. Clique em executar o pipeline inteiro ou a etapa selecionada com suas dependências.
2. Se houver resultados concluídos salvos, responda à pergunta sobre reaproveitamento conforme a seção seguinte.
3. Acompanhe a janela de execução. Quando uma entrada conectada ainda exigir seleção automática, a execução pausa e abre **Selecionar arquivos**. Seleções explícitas salvas no bloco dispensam essa pergunta.
4. Marque somente os arquivos desejados em cada entrada. A lista mostra arquivos compatíveis das origens conectadas, sem selecionar automaticamente todos os resultados novos.
5. Escolha **Processar individualmente** ou **Mesclar arquivos (merge)**.
6. Para ajustar parâmetros, clique em **Configurar etapa**. As escolhas já feitas são incluídas no formulário, evitando outra seleção do mesmo conjunto.
7. Clique em **Continuar com os arquivos selecionados**. Repita a decisão somente nas próximas etapas que ainda precisarem de seleção.

| Escolha | Efeito | Exemplo |
| --- | --- | --- |
| Individual | Um processamento separado por arquivo, com saídas separadas | `alpha.csv` e `beta.csv` geram fingerprints próprios e continuam distinguíveis |
| Merge | Reúne os dados escolhidos para uma execução conjunta | Os compostos de `alpha.csv` e `beta.csv` participam da mesma análise |

Em operações com entradas distintas, como receptores e tabelas de compostos, individual forma combinações. DOCK6 limita essas combinações aos receptores/compostos com poses correspondentes. O consenso pareia Vina e DOCK6 pelo identificador da pose, evitando cruzamentos entre moléculas ou referências diferentes. Escolher merge não remove a necessidade de compatibilidade: tipos de fingerprint, códigos e metadados ainda precisam ser coerentes. Nomes iguais com conteúdos diferentes precisam ser desambiguados ou renomeados antes de combinar.

Metadados de estruturas acompanham as seleções em ambos os modos. **Selecionar depois** fecha a janela e mantém a pausa; **Configurar entradas e continuar** reabre a seleção. Fechar e reabrir a janela durante a sessão conserva arquivos marcados e modo. Escolhas confirmadas ficam na configuração do bloco e são reutilizadas ao abri-la posteriormente. Quando dois arquivos têm o mesmo nome, selecione o caminho que distingue cada resultado.

No docking Vina ou DOCK6, escolher um arquivo de compostos disponibiliza **Composto para docking (opcional)**. Selecione um identificador para usar somente essa molécula; **Todos os compostos** ou o campo vazio usa a tabela inteira. Cada arquivo tem sua própria seleção. Trocar o arquivo limpa a escolha anterior. As opções ficam disponíveis após a origem produzir seus dados; arquivos enviados ou importados podem ser configurados antes da execução.

Se receptor e compostos já estiverem definidos com arquivos explícitos na configuração de Vina/DOCK6, o pipeline executa sem repetir a seleção. Conexões com **Identificar automaticamente**, ainda sem arquivo escolhido, continuam pausando para configurar as entradas. Arquivos e identificadores são revalidados na execução; uma escolha que deixou de existir gera uma mensagem de validação.

A seleção ocorre antes dos blocos consumidores que precisam executar; blocos sem entrada não precisam dessa decisão. Etapas reaproveitadas não repetem a seleção. A última etapa disponibiliza resultados para consulta, sem exigir um próximo consumidor. Leitores acompanham; editores e proprietários confirmam as escolhas.

## 9. Reaproveitar dados depois de reiniciar

1. Encerre e reabra a aplicação usando o mesmo `--data-dir`.
2. Entre na conta e abra o projeto existente, sem criar outro projeto sobre sua pasta.
3. Clique em executar. Se houver resultados íntegros, a pergunta **Reaproveitar dados e resultados?** aparece antes da submissão.
4. Escolha **Sim, reaproveitar** para concluir imediatamente as etapas compatíveis. Isso inclui suas seleções anteriores e evita repetir os popups dessas etapas.
5. Escolha **Não, executar novamente** para recalcular as etapas solicitadas e confirmar as entradas novamente quando necessário.

O reaproveitamento verifica parâmetros, conexões, templates, entradas e integridade dos resultados. Mover ou renomear um bloco não exige recalcular. Alterar dados, modo individual/merge ou opções científicas invalida as etapas afetadas; as demais ainda podem ser reaproveitadas. Acrescentar um novo bloco permite usar as etapas anteriores concluídas. Artefatos reaproveitados permanecem na execução que os produziu.

Uma execução interrompida durante um cálculo não retoma a instrução interna do motor. Inicie uma nova execução e considere reaproveitar as etapas anteriores completas. Uma execução pausada aguardando seleção pode ser reaberta e continuada. Cancelar a pergunta de reaproveitamento não inicia uma nova execução.

## 10. Primeiro estudo completo com arquivos locais

Este roteiro exercita importação, fingerprints, similaridade, grafos e ADMET sem consultas a provedores ou docking.

1. Crie um projeto em uma pasta nova.
2. Salve o exemplo CSV da seção 5 como `molecules.csv` e envie-o como **Compostos**.
3. Adicione **Importar meus arquivos** e selecione esse CSV.
4. Adicione **Gerar fingerprints**; conecte a importação à entrada de compostos.
5. Escolha Morgan, raio 2 e 2.048 bits. Esses parâmetros são uma escolha de exemplo, não uma validação do estudo.
6. Adicione **Calcular similaridade** e conecte os fingerprints. Escolha Tanimoto e um limiar de 0% para observar as relações disponíveis nesse exemplo pequeno.
7. Adicione **Filtrar por grafos** e conecte a saída de similaridade.
8. Adicione **Avaliar ADMET** e conecte os compostos selecionados pelo grafo.
9. Execute. Em cada seleção, marque o arquivo produzido pelo bloco anterior e confirme o modo individual.
10. Em **Execuções**, abra fingerprints, similaridades, grafo e ADMET. Confira códigos e arquivos para verificar que o conjunto percorreu as etapas esperadas.
11. Execute novamente, escolha reaproveitar e verifique a indicação **Reaproveitadas**.
12. Para comparar conjuntos, envie outro CSV, acrescente-o à importação e repita escolhendo individual. Depois compare com merge, observando os resultados e a composição dos grafos.

Um conjunto pequeno pode produzir componentes isolados ou ser inteiramente removido pelos filtros ADMET. Confira os relatórios e não interprete a ausência de resultados como indicação de eficácia ou toxicidade experimental.

## 11. Configurar cada operação científica

### Recuperação e expansão

Em **Recuperar compostos**, escolha o modo de consulta: alvo por nome/IDs/UniProt/texto ou composto diretamente por nome/IDs/similaridade/subestrutura. Revise os filtros opcionais e os limites; escolha se deseja expandir com ChEMBL e/ou PubChem. Limiar e máximo de similares controlam essa ampliação. Após executar, examine `compounds.csv` e os downloads ChEMBL disponíveis; use a tabela adequada como entrada da análise seguinte.

**Expandir similares** utiliza os downloads ChEMBL de uma recuperação existente. Conecte essa origem e revise limiar/máximo; não use um CSV comum como substituto do contexto ChEMBL. **Recuperar ZINC** recebe o arquivo TXT/URI ou script de download exportado pelas tranches ZINC. Envie-o diretamente em **Lista de downloads das tranches ZINC** ou conecte uma importação com tipo **Lista de downloads ZINC**. As tranches SMI/MOL2, inclusive `.gz` e `.bz2`, geram `compounds.csv`; as conformações MOL2 acompanham a tabela. Conecte-a a ADMET/fingerprints ou ao preparo e docking. Consulte o [guia de recuperação](retrieval.md#tranches-zinc-2d-e-3d).

**Recuperar estruturas PDB** aceita texto, IDs PDB, acessos UniProt e códigos de ligantes; EC é opcional. Organismo, resolução, polímero e método refinam a busca. O nome da coleção é opcional. Consulte [recuperação flexível](retrieval.md) e `retrieval_report.json` para os critérios e resultados. Use o botão **Ligantes** de cada PDB para revisar, remover ou adicionar registros antes de preparar ou fazer redocking. O explorador 3D abre em uma aba do navegador, com rotação pelo mouse, zoom e opções de representação; o botão RCSB abre sua página oficial. `pdb_codes.csv` e `retrieval_report.json` são preservados no disco e ocultos nessa listagem. Uma recuperação sem estruturas compatíveis exige revisar os filtros; não avance com um conjunto vazio.

A recuperação PDB prioriza por padrão moléculas orgânicas em contato com a proteína e sítios anotados, reduzindo solventes, íons e cópias equivalentes na seleção para redocking. Consulte `ligand_candidates.csv` para os contatos e motivos da triagem. Em **Ligantes**, use **Mostrar também os demais resíduos** para incluir uma alternativa ou **Usar candidatos recomendados** para recuperar a sugestão; confirme com **Salvar ligantes**. Contato geométrico não comprova sítio ativo. A opção **Manter todos os ligantes** conserva a seleção ampla.

### Preparação e redocking

Em **Preparar para docking**, conecte **Receptor PDB (retrieval ou redocking)** e **Compostos selecionados** separadamente. Para um PDB bruto do retrieval, informe os registros `[PDB, ligante de referência, resíduo, cadeia]` ou forneça `pdb_codes.csv`; as opções de preparo do receptor permanecem disponíveis. Para receptores do redocking, escolha o arquivo `.dockprep.pdbqt` no campo **Receptor que deseja utilizar**: os arquivos exclusivos do ligante ficam ocultos e as opções de preparo do receptor são desabilitadas. Os arquivos complementares e os centros acompanham o receptor; ele é reutilizado sem novo preparo.

Adicione uma ou mais fontes de compostos de **ChEMBL**, **PubChem**, **ZINC** ou arquivos próprios em **Compostos selecionados**. O grupo **Preparação e conformação do ligante** configura esses candidatos externos e permanece disponível. O ligante de referência do PDB serve para definir o sítio de docking, não integra a lista de candidatos. Escolha **Preparar saída para → Vina, DOCK6 ou Vina e DOCK6**. O bloco prepara cada candidato e exporta os formatos selecionados a partir da mesma conformação. Conecte a saída deste bloco tanto à entrada de receptor quanto à entrada de compostos de cada ferramenta escolhida. Para reunir várias fontes em uma execução, selecione **Mesclar arquivos (merge)**; o modo individual prepara cada combinação separadamente. Confira os receptores, compostos e centros produzidos.

Use um bloco separado para cada modo de receptor: não combine PDBs brutos com receptores preparados no mesmo bloco. O pH continua disponível para o preparo dos candidatos quando o receptor é reutilizado. Antes de preparar os compostos, o bloco verifica os centros do sítio (três coordenadas finitas) e os arquivos do receptor exigidos pela saída escolhida. Vina requer `.dockprep.pdbqt`; DOCK6 requer `.dockprep.mol2` e `.noH.pdb`. A seleção aceita o receptor MOL2 quando não existe PDBQT; quando ambos estão presentes, mostra apenas o PDBQT. Se um consumidor Vina precisar do formato rígido, ele é exportado em uma cópia de trabalho sem refazer o preparo Chimera. Arquivos ausentes interrompem o processo com o nome do arquivo necessário, sem refazer o preparo do receptor.

Em **Executar redocking**, mantenha **Preparar complexos antes do redocking** ativado para PDBs brutos. Desative-o apenas ao fornecer o conjunto já preparado. Escolha Vina ou DOCK6 6.11 e confira os campos específicos de busca, sítio e poses. Examine os RMSDs e os logs; defina os critérios de aceitação do protocolo no estudo antes de usar seus receptores no docking de candidatos.

Em **Entradas**, selecione e valide ao menos um par PDB / ligante / resíduo / **Cadeia**; a cadeia é comum às duas estruturas. Revise a tabela e use as ações de editar/remover. Defina uma preparação comum para receptores e outra para ligantes, com cofatores opcionais recolhidos. A resolução vem dos metadados e as opções do ligante valem também para a conformação. Pares já configurados são reutilizados sem uma nova solicitação. Consulte [Configuração do redocking](redocking_configuration.md) para as opções, validações e diagnóstico.

### Consultar resultados de redocking

Falhas de preparação são isoladas por par; os outros complexos continuam e a conclusão informa os erros. A tabela mantém os pares que falharam em vermelho, com **Falha na preparação**, RMSD **—** e **Ver erro**. **Excluir falha e complexo de origem** remove o PDB da origem e os pares desse PDB da próxima execução, após confirmação, preservando resultados já concluídos. O código fica excluído da recuperação original até uma liberação intencional no campo **PDBs excluídos de novas recuperações**. Se todos os pares falharem, o processo conclui com aviso e não inicia o Vina.

Ao expandir uma etapa concluída em **Execuções**, a tabela mostra Target, PDB, ligante, resíduo, cadeia e **RMSD (Å)**. **Ver simulação** abre um popup com os arquivos correspondentes àquela seleção. A lista distingue o receptor preparado, o ligante de referência, as poses e os metadados. Arquivos compartilhados pelo receptor ou pela coleção acompanham as simulações correspondentes.

Use o botão de download de cada linha para baixar um arquivo, ou **Baixar todos (ZIP)** para obter o conjunto da simulação. O ZIP conserva as pastas, inclusive quando o ligante de referência e a saída do Vina têm o mesmo nome. Arquivos de outras simulações ficam fora desse conjunto.

O botão **Visualizar estrutura 3D** abre o visualizador no navegador, como na consulta dos PDBs. Ele aceita PDB, PDBQT e MOL2 para examinar receptor, ligante e poses antes do download. Ao abrir um arquivo com várias poses, use **Modelo** no visualizador. A cena aberta pelo botão **3D** da tabela seleciona automaticamente a pose de menor score, ou a primeira sem score; essa cena não oferece troca de modelo. As legendas e ações acompanham o idioma selecionado. Leitores também podem consultar e baixar os resultados.


O botão **3D** de cada resultado abre uma comparação com a melhor pose do Vina em ciano, o ligante cristalográfico em dourado e o receptor utilizado. O complexo cristalográfico completo pode ser ativado no painel de camadas. As coordenadas são mantidas, sem realinhamento que pudesse ocultar deslocamentos. Quando o PDB original não estiver disponível, uma referência preparada, se existente, aparece identificada dessa forma.

**Resíduos** apresenta as menores distâncias entre átomos pesados do receptor e do ligante, para a pose e a referência. O limite inicial é 4 Å, ajustável entre 2 e 8 Å. É possível baixar a tabela em CSV e selecionar um resíduo no visualizador para destacá-lo. Essa tabela indica proximidade geométrica. A classificação por tipo aparece separadamente nos traçados do visualizador, quando os dados químicos permitem, conforme a seção seguinte.

### Estilos, hidrogênios e interações da pose

No visualizador, **Estilo do receptor** e **Estilo do ligante** são independentes:
por exemplo, receptor em fitas e ligante em esferas. **Hidrogênios do ligante**
mostra ou oculta somente os H presentes no arquivo, sem gerar os ausentes no Vina
nem alterar a preparação ou os scores.

**Mostrar interações em 3D** desenha cada relação classificada com traçado
tracejado: verde para ligação de hidrogênio, rosa/coral para π–π paralelo/em T,
lilás para contato hidrofóbico e verde-lima para van der Waals geométrico. Passe o
mouse sobre o traçado para ver tipo, resíduo, cadeia e distância. A legenda e os
filtros usam as mesmas cores; a lista permite destacar o resíduo pelo clique.

A classificação depende de topologia química e resíduos válidos; ligações de
hidrogênio exigem H explícitos. Sem os dados necessários, a tela informa a
limitação e mantém os contatos por distância. A análise 3D é geométrica e difere
das energias por resíduo de **Footprint**. Consulte o
[guia do visualizador](molecular_viewer.md) para camadas, seleção de poses,
controles, interpretação e atualização.

### Fingerprints, similaridade e grafos

**Gerar fingerprints** escolhe um tipo por bloco: Morgan, MACCS ou farmacóforo. Raio e bits aplicam-se ao Morgan. Para comparar tipos, duplique o bloco e mantenha os resultados identificados.

**Calcular similaridade** lê fingerprints e identifica o tipo das saídas nativas. Para arquivos próprios, confira o seletor de tipo. Métrica e limiar são definidos nessa etapa. A opção aproximada altera a estratégia de busca; compare-a com a análise exata conforme o protocolo.

**Filtrar por grafos** preserva as relações recebidas e seleciona o MCC, o maior componente conectado. Configure o tempo máximo do MCS e as opções de anéis quando desejar examinar o fragmento comum. Os SMILES vêm da linhagem que produziu as similaridades; para dados externos, forneça a tabela molecular correspondente quando quiser visualizar estruturas.

### ADMET

**Avaliar ADMET** oferece perfis de pesquisa exploratória, oral, SNC, periférica, fragmentos, produtos naturais e critérios personalizados. Escolha entre **Anotar todos os compostos** e **Filtrar pela hipótese**; configure BBB/HIA, regras, alertas e intervalos em seções organizadas. RDKit fornece descritores e heurísticas locais. O método opcional ADMET-AI acrescenta 41 endpoints de modelos treinados e informa sua disponibilidade na configuração. Para instalar seu ambiente isolado, execute `python scripts/setup_admet_ai.py` na raiz do projeto.

O gráfico EGG adapta-se ao painel, identifica compostos ao passar o mouse e permite consultar estrutura e propriedades. As classes BBB/HIA locais não representam uma reprodução validada do BOILED-Egg publicado. Examine também o CSV `_evaluated` e o relatório de métodos e versões. Consulte a [análise científica ADMET](admet_research.md) para interpretar os métodos e seus limites.

Seleções de arquivos já salvas nos blocos dispensam nova escolha durante a execução. Em **Execuções**, cada etapa oferece **Remover todos os resultados**, com confirmação; os arquivos podem ser restaurados pelo histórico. A remoção preserva registro e log e fica bloqueada durante execuções ativas.

### Vina, DOCK6 e consenso

**Docking com Vina** e **Docking com DOCK6** recebem compostos de qualquer bloco que publique moléculas: ChEMBL, PubChem, ZINC, ADMET, fingerprints, grafos, consenso ou importação do usuário. Selecione a tabela na entrada de compostos e conecte separadamente o receptor preparado. Não é obrigatório passar por grafos, ADMET ou pelo outro motor.

Para reutilizar conformações, conecte **Vina → DOCK6** ou **DOCK6 → Vina** na entrada de compostos, escolhendo `docking_results.csv` ou a pose desejada. A identidade e o SMILES são preservados; a conversão utiliza a conformação selecionada, sem gerar outra a partir do SMILES. Sem uma pose de entrada, a preparação gera um conformero 3D. DOCK6 usa o centro do sítio do receptor para selecionar esferas e posicionar ligantes recém-gerados; a entrada opcional de poses Vina continua disponível para refinamento de candidatos correspondentes.

Confira a caixa, pH e esforço do Vina e as cargas, superfície, distância, raio, busca flexível/rígida e footprint do DOCK6. Ambos exigem receptores preparados, metadados e centros; DOCK6 também exige MOL2 do receptor e PDB sem hidrogênios. Seus resultados apresentam código molecular, receptor, score, SMILES, **3D** da pose calculada e **Remover**. O botão **3D** sobrepõe a melhor pose ao receptor utilizado; uma referência preparada, quando disponível, aparece em dourado. **Resíduos** oferece a mesma tabela de contatos do redocking. As novas execuções guardam cópias do receptor e da referência junto aos resultados para preservar essa associação.

O botão **Footprint** do DOCK6 abre um popup integrado com o gráfico das energias de Van der Waals e eletrostáticas por resíduo, comparando a referência com a pose final calculada. Use o zoom e arraste para explorar o gráfico; **Recentrar** restaura a visualização inicial e **Baixar PDF** salva o documento original. O mesmo popup está disponível no consenso quando há um footprint associado. A geração valida o conteúdo das energias e interrompe a etapa se o motor produzir um arquivo vazio. Resultados antigos podem ter usado o ligante minimizado antes do docking; a interface avisa sobre essa origem. Reexecute esses experimentos para obter o gráfico da pose final.

**Consenso de docking** reúne os resultados selecionados de cada motor, incluindo vários lotes e ramificações. Calcula somente a interseção por código molecular e receptor, usando o menor score quando há poses repetidas. Resultados do mesmo código com estruturas diferentes são recusados. Não havendo interseção, o bloco é marcado como não executado e explica que os motores não avaliaram compostos em comum para o mesmo receptor.

Após conectar os resultados Vina e DOCK6, o consenso segue automaticamente na execução do pipeline, sem confirmação de arquivos ou parâmetros adicionais. A tabela apresenta código do composto, score Vina, score DOCK6 e score normalizado (0 a 1). Clique no cabeçalho de uma coluna para alternar a ordenação crescente e decrescente; a ordenação considera toda a tabela antes da paginação. O score normalizado é a média das normalizações min-max dos dois motores, e valores maiores indicam posições melhores no ranking. O CSV também preserva z-score e min-max. Há botões **3D Vina**, **3D DOCK6** e **Remover** em cada linha. O consenso usa os scores originais dos dois motores, sem exigir receptor, ligante de referência ou pose. Arquivos estruturais disponíveis são anexos opcionais para visualização e sua ausência não bloqueia o cálculo. Para um composto ou scores constantes, as normalizações são zero. A remoção exige permissão de edição, fica registrada e afeta próximas entradas; não recalcula análises já concluídas. As tabelas principais de cada lote são usadas sem reintroduzir linhas removidas a partir das tabelas auxiliares.

## 12. Explorar resultados e navegar pelos grafos

1. Abra **Execuções** e expanda a etapa, ou use **Resultados** no bloco.
2. Se houver várias tabelas ou análises, selecione a desejada.
3. Use busca e paginação para localizar compostos; escolha 10, 25, 50 ou 100 itens por página quando disponível.
4. **2D** mostra a estrutura; **3D** gera um conformero local do SMILES, com rotação, zoom e recentralização. Essa visualização não é uma pose de docking.
5. Baixe o arquivo para conservar os dados completos da análise.

No grafo, escolha **Grafo completo** para ver todos os componentes. Cada componente recebe espaço próprio e os vértices são separados. Use o seletor de componentes para focar uma região; clique em um vértice para ver seus detalhes, vizinhos e pesos das relações. Buscar um código permite encontrar o nó e navegar até seu componente. Use arraste e zoom; **Ajustar à tela** e **Recentrar** ajudam a recuperar a visão geral.

**Restaurar grafo completo** limpa todos os filtros e a agregação e recupera a visão geral, inclusive após investigar uma comunidade. Em **Layout, centralidade e componentes**, escolha ForceAtlas2, LinLog, Fruchterman–Reingold ou outras organizações e configure cores e tamanhos. A legenda indica os valores ou os grupos de cada cor. **Laboratório de dados** permite pesquisar nós e arestas e exportar a tabela completa. PNG/SVG incluem legenda; GEXF/GraphML permitem continuar a análise no Gephi. Consulte [recursos, interpretação e limites](graph_explorer.md).

Ao escolher **MCC**, a exploração se restringe ao maior componente conectado. O painel mostra **SMILES do fragmento**, sua estrutura e o estado da busca do fragmento comum quando esses dados estão disponíveis. O SMILES é extraído da molécula de referência para o padrão encontrado por MCS; a conectividade comum não comprova identidade estereoquímica de todos os compostos. Busca interrompida pelo limite de tempo deve ser interpretada pelo estado informado. Sem SMILES correspondentes, a topologia pode ser explorada, mas não há estrutura molecular do fragmento a mostrar.

Resultados antigos podem ter apenas PNG. Para produzir os artefatos interativos, execute novamente escolhendo **Não, executar novamente**. O limite de abertura é de 128 MB para grafos e de 32 MB para outras visualizações `*.biomol-view.json`; arquivos maiores podem ser baixados. Em BOILED-Egg, pontos sobrepostos oferecem seleção dos compostos; busca também permite distingui-los.

## 13. Conferir exclusões, erros e cancelamentos

Durante o processamento, linhas inválidas são retiradas das cópias de trabalho, mantendo os originais. A etapa informa a quantidade excluída e oferece `molecule_exclusions.json`, com origem, linha, código e motivo. Confira esse relatório antes de interpretar o tamanho final do conjunto.

| Situação | Ação recomendada |
| --- | --- |
| CSV recusado | Conferir UTF-8, cabeçalho, separador, número de campos e formato mostrado no formulário |
| Nenhum arquivo compatível no popup | Conferir o tipo publicado pela origem e seus resultados; selecionar a origem correta |
| Nome de arquivo ambíguo | Selecionar o caminho completo que distingue a execução ou o lote |
| Fingerprints incompatíveis | Usar o mesmo tipo e a mesma largura; processar separadamente alternativas |
| Metadados/centros sem correspondência | Conferir PDB, ligante, resíduo e cadeia; reenviar o conjunto preparado completo |
| Poses não correspondem aos compostos/receptores | Conferir os códigos e selecionar Vina/candidatos da mesma referência |
| Falha de consulta externa | Ler o log com contexto da consulta; revisar filtros e repetir após verificar o provedor |
| Motor não encontrado ou comando falhou | Conferir ambiente do worker, executáveis, arquivos de entrada e log do comando |
| Etapa descendente impedida | Corrigir a origem que falhou e executar novamente |
| Reaproveitamento não ocorreu | Conferir se entradas, configuração, modo ou integridade das saídas mudaram |

**Minimizar** mantém o acompanhamento em uma faixa do projeto; **Acompanhar execução** reabre a janela. Para interromper, use **Cancelar execução** e aguarde o estado final. O pedido de cancelamento pode levar tempo enquanto o backend encerra o worker. Resultados completos de etapas anteriores podem ser reaproveitados em nova execução; uma etapa parcial não é tratada como concluída.

## 14. Curadoria, colaboração e histórico

**Remover** em uma tabela de compostos solicita confirmação e remove a linha apenas daquele CSV. Não sincroniza automaticamente outras tabelas. Escolha explicitamente qual tabela alimentará a próxima etapa. Remover um arquivo exige desvincular entradas que o utilizam. Edição/remoção requer editor ou proprietário, sem pipeline ativo; os originais e a autoria ficam no histórico. Dados alterados invalidam as análises dependentes.

Para compartilhar, o proprietário abre **Compartilhar**, informa o e-mail de uma conta já cadastrada e escolhe **Editor** ou **Leitor**. O destinatário precisa aceitar no workspace. O convite é interno, sem envio de e-mail. Leitores consultam e baixam; editores também configuram, selecionam arquivos e executam; o proprietário gerencia colaboradores e exclusão do projeto.

Alterações de colaboradores são atualizadas no projeto aberto. Formulários e rascunhos em edição são preservados; um conflito no mesmo campo precisa ser revisado antes de reaplicar. Use nomes de blocos e descrições para tornar as alternativas compreensíveis à equipe.

Em **Histórico**, confira data, usuário e alteração. Como proprietário, **Restaurar antes** retorna configurações, arquivos, resultados e permissões ao estado anterior ao registro escolhido; leia a confirmação e cancele qualquer execução ativa antes. A restauração também é registrada. Consulte [projetos e versões](projects.md) para os efeitos sobre colaboradores e versões importadas.

## 15. Exportar, fazer backup e migrar

1. Conclua ou cancele execuções ativas.
2. Use **Exportar** no projeto e conserve o pacote `.bme.zip`.
3. Para backup completo da instalação, encerre a aplicação e copie também o workspace das contas e as pastas dos projetos, incluindo arquivos ocultos.
4. Em outro computador, instale o ambiente, crie/entre em uma conta e use **Importar projeto**.
5. Escolha uma pasta nova ou vazia no computador do backend e importe o pacote.
6. Abra o projeto, confira entradas/resultados e instalações de motores. Ao executar, escolha se deseja reaproveitar os dados compatíveis.

`project.json` sozinho não contém os resultados nem as entradas. O pacote reúne dados, configurações e histórico; contas e senhas não são exportadas. O usuário que importa torna-se proprietário. Colaboradores precisam de novo aceite no destino. Pacotes grandes gerados no navegador ficam em `.exports/` para cópia direta; consulte os limites no guia de projetos.

## 16. Glossário e percurso de referência

| Termo | Significado nesta aplicação |
| --- | --- |
| Bloco / etapa | Uma operação configurada no pipeline |
| Entrada / binding | Arquivo ou resultado conectado a um campo de uma operação |
| Artefato | Arquivo produzido por uma etapa |
| Fingerprint | Representação molecular usada na comparação |
| Similaridade | Relação numérica entre dois códigos moleculares |
| Componente conectado | Grupo de vértices ligados por caminhos no grafo |
| MCC | Maior componente conectado |
| MCS | Busca de subestrutura comum; pode ter limite de tempo |
| SMILES | Representação textual de uma estrutura molecular |
| Redocking | Reposição de um ligante de referência para avaliar o protocolo |
| Reaproveitamento | Uso de resultados completos e compatíveis já persistidos |

Um percurso comum é recuperar/importar compostos → fingerprints → similaridade → grafos → ADMET, em paralelo a recuperar/importar PDBs → preparação/redocking. As ramificações convergem em Vina → DOCK6 → consenso. Escolha os arquivos e o modo de processamento em cada passagem; registre parâmetros, exclusões e critérios científicos junto ao estudo. O [relatório de validação](pipeline_validation.md) explica os testes e os limites de verificação dessa cadeia.

## Verificações antes da execução

A seleção de pares é obrigatória e usa uma única **Cadeia** para receptor e ligante. Um redocking já configurado reutiliza os pares escolhidos; o pipeline só pede a seleção quando ela ainda não existe. Resolução vem dos metadados, sem campo editável. A tabela lista os pares com ações de editar/remover. Cofatores e opções de solvente, hidrogênios, minimização e cargas são definidos uma vez para todos os receptores e uma vez para todos os ligantes. A preparação e a conformação do ligante usam as mesmas opções.

Antes de calcular, o sistema verifica os resíduos, a cadeia e os cofatores nos PDBs. Escolha **AutoDock Vina** ou **DOCK6 6.11** e preencha os campos específicos exibidos. Entradas já preparadas precisam de PDBQT para Vina ou MOL2 e receptor sem H para DOCK6, além de centros válidos. Chimera e Open Babel são necessários para a preparação; o motor escolhido e suas ferramentas devem estar instalados. O diretório do interpretador científico também integra o PATH dos processos.

Falhas de ferramentas mostram o executável, o código de saída e a mensagem do processo. Scripts Chimera que falharam são preservados para diagnóstico. Corrija a entrada ou a instalação indicada e tente novamente. Não há controle de verbose: a verbosidade do Vina permanece em zero.

Consulte [Configuração do redocking](redocking_configuration.md) para o fluxo completo.

## ADMET por objetivo de pesquisa

O bloco ADMET oferece perfis editáveis, modo descritivo ou filtragem, BBB/HIA positivos e negativos, regras físico-químicas, alertas estruturais e intervalos personalizados. Novos blocos começam em exploração ampla. Consulte a [revisão científica e o guia de configuração ADMET](admet_research.md) para os métodos, arquivos gerados, bibliotecas opcionais e limitações.

### Compostos não utilizados na preparação

A preparação de candidatos para Vina e DOCK6 trata cada composto separadamente. Quando a conformação, a conversão, a atribuição de cargas ou a minimização falha, esse composto é excluído das entradas do docking e os demais continuam. Arquivos parciais do composto são removidos. Se nenhum candidato puder ser preparado, o motor não é iniciado e a etapa informa o motivo.

Em **Execuções → Resultados por etapa**, use **Compostos não utilizados** para consultar a listagem paginada e **Baixar lista de compostos não utilizados** para baixar `preparation_failures.csv`. O relatório registra identificador, SMILES enviado, arquivo/linha de entrada, contexto da preparação e motivo. As falhas também entram no relatório geral de exclusões. As entradas originais são preservadas. Falhas de infraestrutura, como motor ausente ou disco cheio, continuam sendo erros da etapa.

## Reutilização do receptor no Vina

Ao conectar a saída preparada do redocking ao docking, o pipeline reutiliza o receptor e o centro do sítio em todos os lotes. Antes de iniciar, verifica os arquivos exigidos pelo motor e as coordenadas do centro. Um arquivo preparado ausente gera um erro de entrada; não dispara uma nova preparação silenciosa. Na preparação de estruturas brutas, a conversão do receptor ocorre uma vez por PDB/cadeia, mesmo quando existem vários ligantes de referência.

No fluxo padrão, a API Python do Vina mantém o receptor e os mapas de afinidade em memória durante o lote. Os mapas são calculados antes de carregar qualquer ligante, incluindo todos os tipos atômicos, conforme a [documentação oficial do Vina](https://autodock-vina.readthedocs.io/en/latest/docking_python.html). Mudanças no receptor, centro, tamanho da caixa ou parâmetros de inicialização criam uma nova sessão. As opções dos templates continuam sendo respeitadas; opções especiais não suportadas pela API ou a ausência dessa biblioteca usam o executável CLI e registram essa alternativa no log.

Os eventos `receptor.reused`, `ligand.preparation.started`, `vina.maps.started` e `vina.maps.reused` distinguem cada operação. `prepare.com` é a preparação de um candidato externo, que pode demorar devido a cargas AM1 ou minimização; não indica preparação do receptor. Candidatos com arquivos preparados são reutilizados. A semente configurada e os parâmetros de busca são preservados, com os compostos processados em ordem estável. A reutilização ocorre dentro do lote; os mapas em memória não são compartilhados entre workers.


## Validação e exclusões no DOCK6

O fluxo verifica instalação e parâmetros globais, receptor preparado, centro do sítio, esferas, caixa e arquivos de grid. Para cada ligante, verifica o MOL2 (contagens de átomos e ligações, identificadores, coordenadas e cargas finitas, topologia química e conectividade). Um ligante com fragmentos desconectados ou topologia inválida é excluído, sem reconstruir ligações por suposição. Isso é relevante na conversão de poses PDBQT do Vina, cujo formato não conserva toda a informação química de um MOL2. Uma pose Vina solicitada precisa corresponder ao composto e receptor selecionados; não é substituída silenciosamente por outra entrada.

Superfície, esferas e grid são compartilhados pelos compostos do mesmo receptor/sítio. Minimização, docking e footprint usam processos independentes por composto. Falhas nativas (incluindo SIGSEGV, código -11), arquivos ausentes, saídas vazias ou scores inválidos excluem aquele composto das etapas seguintes e dos resultados finais; os demais continuam. Falhas ao gerar o gráfico do footprint também são registradas individualmente. Configurações e logs das ferramentas ficam preservados nos resultados para diagnóstico. Identificadores com pontos são conservados integralmente.

A lista **Compostos não utilizados** inclui problemas de preparação e de execução, com etapa, identificador, SMILES, entrada e motivo. Ela mantém o arquivo `preparation_failures.csv` por compatibilidade com o relatório do Vina. Se todos forem excluídos, a etapa informa que nenhum resultado utilizável foi produzido. Erros globais — instalação ausente, parâmetros incompletos, receptor/sítio/grid indisponíveis ou falhas de disco — continuam interrompendo a etapa, pois não representam um problema de um composto específico.

A revisão seguiu o [manual oficial do DOCK6](https://dock.compbio.ucsf.edu/DOCK_6/dock6_manual.htm), o [tutorial de docking rígido e flexível](https://dock.compbio.ucsf.edu/DOCK_6/tutorials/ligand_sampling_dock/ligand_sampling_dock.html) e o [tutorial de geração de grids](https://dock.docking.org/DOCK_6/tutorials/grid_generation/generating_grid.html). As verificações estruturais não garantem ausência de defeitos no executável nativo; o isolamento permite continuar o lote e identificar os compostos afetados.

## Navegação e manutenção dos resultados

As ações do workspace ficam na lateral. No projeto, **Navegação do projeto** e **Pipeline e execução** são recolhíveis, com ajuda para cada ação. Use o canvas ampliado, zoom e **Organizar blocos** para pipelines grandes. O [guia da interface](frontend.md#navegacao-canvas-e-execucoes-concluidas) detalha os controles, seleção de arquivos salvos, cartões de execução e remoção de resultados.

Nos traçados 3D, **Mostrar nomes dos resíduos** permite uma exportação limpa sem remover as interações. Consulte o [visualizador](molecular_viewer.md#ocultar-os-nomes-dos-residuos).

[Veja um footprint DOCK6](../imgs/examples/dock6-footprint.png) com o [PDF original](../imgs/examples/dock6-footprint.pdf): o exemplo do tutorial 1ABE compara contribuições de energia da referência e da pose gerada por Vina → DOCK6. Esses perfis de energia são distintos dos contatos geométricos 3D.


## Exportação do pipeline em PNG

A câmera nos controles do canvas salva **todo o pipeline em PNG**, incluindo blocos fora da área visível. A imagem preserva posições, nomes, cores, portas e conexões, independentemente do zoom. O arquivo tem metadados de 300 dpi e resolução ampliada, com limite de 32 milhões de pixels para controlar o uso de memória. Leitores também podem exportar o diagrama para artigos e apresentações.

## Reconhecer as origens e acompanhar novos conjuntos

Os resultados exibem o bloco de origem e o identificador do alvo, como **ACHE · CHEMBL220**, **BCHE · CHEMBL1914** e **BACE1 · CHEMBL4822**. A origem acompanha os arquivos através de ADMET, PubChem, fingerprints, similaridade e grafos, mesmo quando arquivos diferentes se chamam `compounds.csv` ou `morgan_compounds.csv`.

Os seletores permitem pesquisar pelo alvo ou identificador ChEMBL. As tabelas podem ser filtradas por origem; os dashboards de ADMET e grafos mostram **Origem dos dados** e um painel **Arquivo e caminho no pipeline**. No grafo, a origem também acompanha a seleção da fonte molecular, a inspeção dos compostos e os relatórios exportados. Arquivos complementares pertencem ao lote selecionado, evitando misturar análises que usam o mesmo nome de arquivo.

Ao acrescentar as três saídas BBB+ do ADMET às três entradas vindas da PubChem, fingerprints deve produzir seis conjuntos separados. Uma seleção antiga de similaridade ou grafos que abrangia todas as saídas anteriores passa a acompanhar essa ampliação. A seleção efetiva é registrada no histórico e salva no bloco. Assim, resultam seis grafos: ADMET direto e PubChem para cada alvo.

A opção **Acompanhar novos conjuntos da origem**, nos blocos de similaridade e grafos, inclui as saídas compatíveis dos produtores conectados. Desative-a para manter somente os arquivos explicitamente selecionados. Seleções parciais antigas não são convertidas automaticamente em seleções completas. Para manter os seis conjuntos separados, use **Processar individualmente**; a opção merge reúne as entradas conforme a configuração escolhida.

A origem exibida é reconstruída a partir dos arquivos efetivamente usados e do histórico da execução. Não se atribui um alvo apenas porque o nome de um arquivo parece pertencer a ele. Resultados já calculados recebem os novos rótulos ao serem abertos.


## Exemplos rápidos por funcionalidade

Os exemplos abaixo são roteiros de quatro linhas para a interface atual. Os IDs entre parênteses correspondem ao catálogo; importação é um bloco do workspace. Use dados e critérios adequados ao seu estudo. Buscas externas exigem rede; preparo e docking exigem os motores indicados no guia de instalação.

### Importar meus arquivos (`import_results`)

```text
Adicione o bloco de importação e escolha o tipo de dado.
Envie compostos.csv com molecule_chembl_id e canonical_smiles.
Valide o arquivo e conecte a saída a Avaliar ADMET.
Execute o bloco consumidor e confira a tabela de resultados.
```

### Recuperar ChEMBL (`retrieve_compounds`)

```text
Adicione Recuperar ChEMBL e busque pelo alvo CHEMBL220.
Revise os filtros de bioatividade e a expansão por similares.
Salve a configuração e execute o bloco.
Selecione molecules ou similars na entrada do próximo bloco.
```

### Recuperar PubChem (`retrieve_pubchem`)

```text
Adicione Recuperar PubChem e escolha referência manual.
Selecione CID e informe 2244; use limiar 90 e máximo 10.
Salve e execute a busca, que requer acesso à rede.
Conecte a tabela de compostos a ADMET ou fingerprints.
```

### Expandir similares (`expand_similar_compounds`)

```text
Conecte uma recuperação ChEMBL ao bloco Expandir similares.
Escolha o conjunto de origem e informe o limiar desejado.
Salve e execute a expansão.
Confira novos compostos e relações antes de selecionar a saída.
```

### Recuperar estruturas PDB (`retrieve_structures`)

```text
Adicione Recuperar estruturas PDB e informe o ID 1P0I.
Defina o nome do alvo e revise filtros e triagem de ligantes.
Execute e abra os registros recuperados.
Confira cadeia, ligante e estrutura 3D antes do redocking.
```

### Recuperar ZINC (`retrieve_zinc`)

```text
Obtenha uma lista de URLs de tranches SMI ou MOL2 da ZINC.
Envie a lista e selecione-a em Recuperar ZINC.
Ajuste os downloads simultâneos e execute.
Confira o relatório e conecte a tabela de compostos ao próximo bloco.
```

### Preparar para docking (`prepare_structures`)

```text
Conecte o receptor e a tabela de compostos em entradas separadas.
Escolha Vina, DOCK6 ou ambos como formatos de saída.
Revise cofatores, hidrogênios, cargas gas/am1 e minimização.
Execute e conecte receptor e candidatos preparados ao motor escolhido.
```

### Avaliar ADMET (`admet`)

```text
Conecte uma tabela de compostos ao bloco Avaliar ADMET.
Escolha perfil oral, modo de anotação e método RDKit local.
Revise BBB/HIA, regras e alertas; salve e execute.
Abra o EGG e as tabelas; selecione a saída compatível para continuar.
```

### Gerar fingerprints (`fingerprints`)

```text
Conecte uma tabela com identificadores e SMILES válidos.
Escolha Morgan com raio 2 e 2048 bits.
Salve e execute; confira o conjunto gerado.
Conecte a saída ao bloco Calcular similaridade.
```

### Calcular similaridade (`similarity`)

```text
Conecte os fingerprints Morgan já calculados.
Escolha Tanimoto e limiar de 70%.
Revise o modo aproximado antes de executar.
Conecte os CSVs de relações ao bloco Filtrar por grafos.
```

### Filtrar por grafos (`graphs`)

```text
Conecte relações de similaridade prontas ao bloco de grafos.
Revise as opções de fragmento comum e execute.
Abra Full Graph ou MCC e inspecione componentes e moléculas.
Selecione a tabela molecular da análise para ADMET ou docking.
```

### Executar redocking (`redocking`)

```text
Conecte PDBs e selecione os pares receptor/ligante/resíduo/cadeia.
Valide a curadoria e confira cofatores e preparo com Chimera.
Escolha Vina ou DOCK6 6.11, revise suas opções e execute.
Avalie RMSD, logs, Resíduos e o botão 3D de cada par concluído.
```

### Docking com Vina (`docking_vina`)

```text
Conecte o receptor e os compostos nas respectivas entradas.
Selecione o receptor, a tabela molecular e o preparo necessário.
Confira caixa, pH, esforço de busca e número de poses.
Execute e abra scores e comparação 3D da pose com o receptor.
```

### Docking com DOCK6 (`docking_dock6`)

```text
Configure a instalação DOCK6 com bin/ e parameters/.
Conecte receptor compatível e compostos; poses Vina são opcionais.
Revise cargas, superfície, esferas e busca flexível ou rígida.
Execute e inspecione scores, poses 3D e footprint disponível.
```

### Consenso de docking (`consensus`)

```text
Conecte resultados Vina e DOCK6 às respectivas entradas.
Use compostos e receptores correspondentes nos dois resultados.
Revise o peso de repulsão e execute o consenso.
Ordene os scores normalizados e inspecione as poses disponíveis.
```

### Projetos, idioma e colaboração

```text
Escolha Português no login e entre com sua conta.
Crie um projeto e selecione a pasta de destino.
Abra Compartilhar e convide uma conta existente como leitor ou editor.
Consulte o histórico para identificar alterações e versões.
```

### Execução, cancelamento e reutilização

```text
Configure e salve os blocos; conecte suas dependências.
Execute o pipeline e acompanhe etapas e logs em Execuções.
Use Cancelar execução para interromper uma execução ativa.
Reabra o projeto e reutilize resultados compatíveis ao continuar.
```

### Seleção de entradas e merge

```text
Conecte duas fontes de compostos ao bloco consumidor.
Escolha apenas as saídas produzidas que deseja processar.
Use Processar individualmente para manter os conjuntos separados.
Escolha merge quando quiser reunir entradas compatíveis.
```

### Inspeção 2D, 3D e contatos

```text
Abra uma tabela de compostos e consulte a estrutura 2D.
Nos resultados de docking ou redocking, clique em 3D.
Na aba do visualizador, gire a cena e altere estilos e camadas.
Abra Resíduos para consultar os contatos geométricos.
```

### Resultados, limpeza e histórico

Ao lado de **Atualizar resultados**, há duas ações com ícones vermelhos:

| Ação | Remove | Preserva |
|---|---|---|
| **Remover todos / Remove All** | Resultados gerados e caches computacionais do projeto | Blocos, conexões e todas as configurações |
| **Limpar tudo / Clear All** | Os mesmos resultados, seleções de entrada, parâmetros científicos e templates personalizados | Blocos e conexões; os parâmetros voltam aos padrões |

Ambas preservam os arquivos originais da biblioteca, logs e histórico. Confira o
diálogo antes de confirmar; cancelar não apaga dados. Finalize ou cancele uma
execução ativa antes de limpar. Após Clear All, configure novamente os blocos.

Com a reutilização habilitada, lotes válidos permanecem disponíveis após reiniciar.
A equivalência entre configuração, entradas e resultados é verificada ao salvar e
antes de executar; não há sincronização manual. Adicionar entradas ou ações
independentes preserva resultados anteriores e processa somente o complemento.
Alterar parâmetros ou remover entradas/ações exige confirmação antes de salvar.
O modal lista os blocos com resultados afetados, incluindo dependentes. Confirmar
salva a alteração e remove somente resultados obsoletos; cancelar restaura a
configuração anterior e preserva os arquivos. Resultados de outras entradas e
arquivos compartilhados permanecem disponíveis. Uma execução que ainda utiliza
os dados afetados deve terminar ou ser cancelada antes dessa alteração.

Por exemplo, adicionar MACCS a um bloco com Morgan concluído executa apenas MACCS.
Alterar depois o raio de Morgan preserva MACCS. Remover um arquivo de um conjunto
individual remove apenas seu lote e os lotes dependentes. Um resultado calculado
por merge representa todo o conjunto: alterar seus integrantes exige reconstruir
essa combinação, sem substituir o cálculo por uma aproximação de lotes individuais.

```text
Finalize ou cancele a execução antes de remover seus dados.
Abra Resultados por etapa e escolha Remover todos os resultados.
Confira a confirmação e aplique a remoção desejada.
Use o histórico para restaurar arquivos de resultados se necessário.
```

### Exportar e migrar

```text
Use a câmera do canvas para exportar o pipeline em PNG.
Baixe tabelas, figuras ou o ZIP da simulação em Resultados.
Exporte o projeto pelo cartão para guardar configuração e dados.
Importe o pacote no workspace de destino e confira as ferramentas.
```


## Identificar os conjuntos nas saídas ADMET

As novas saídas ADMET incluem a origem do conjunto antes do tipo de resultado:

| Tag | Origem |
| --- | --- |
| `MOL` | Moléculas de referência da seleção ChEMBL ou arquivo identificado como moléculas |
| `SIM` | Moléculas similares da seleção ChEMBL ou arquivo identificado como similares |
| `PUBCHEM` | Compostos recuperados da PubChem |
| `ZINC` | Compostos recuperados da ZINC |
| `MIX` | Entradas mescladas de categorias diferentes ou coleção consolidada |
| `INPUT` | Arquivo cuja categoria não pôde ser determinada com segurança |

Exemplos: `ACHE_molecules_MOL_BBB+.csv`, `ACHE_similars_SIM_HIA-.csv` e
`ACHE_molecules_MOL_egg.png`. A mesma identificação acompanha a tabela avaliada,
selecionados/rejeitados, gráficos interativos, relatório de métodos e relatório
de exclusões, quando gerado. O nome do alvo vem da linhagem do pipeline; arquivos
avulsos usam uma identificação conservadora pelo nome. Caracteres impróprios e
nomes muito longos são normalizados com um identificador curto quando necessário.

Os finais `_evaluated`, `_BBB+`, `_HIA-`, `_egg` e `_admet_report` continuam sendo
reconhecidos pelos leitores e blocos seguintes. A mudança se aplica às novas
saídas; arquivos históricos não são renomeados nem precisam ser recalculados
para continuar sendo lidos.
