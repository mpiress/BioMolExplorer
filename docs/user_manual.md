# Manual do usuário: do primeiro acesso aos resultados

[Documentação](README.md) · Português · [English](en/user_manual.md)

Este manual acompanha um estudo desde a escolha do idioma até a exportação do projeto. Os nomes dos botões abaixo correspondem à interface em português. Para conhecer parâmetros técnicos e APIs, consulte [Workspace Flet](frontend.md), [uso do backend](backend_usage.md), [projetos e versões](projects.md) e [validação do pipeline](pipeline_validation.md).

## 1. Preparar o ambiente e iniciar

Baixe o [projeto no GitHub](https://github.com/mpiress/BioMolExplorer) e entre na pasta do código:

```bash
git clone https://github.com/mpiress/BioMolExplorer.git
cd BioMolExplorer
```

Sem Git, use **Code → Download ZIP**, extraia o arquivo e abra um terminal na pasta extraída. O [guia de instalação e configuração](installation.md) detalha o download, as ferramentas externas e os caminhos.

**Instale UCSF Chimera 1.17, DOCK6 6.11 e DMS antes de iniciar a aplicação.** Eles precisam estar disponíveis no computador dos cálculos, inclusive em modo web. A instalação da interface e do ambiente Conda não instala essas ferramentas; sem elas, as etapas que dependem delas falharão. Conclua as verificações do guia de instalação antes de executar os comandos de inicialização.

1. Na raiz do código baixado, crie o ambiente científico.
2. Ative o ambiente e instale a interface.
3. Configure o `PATH` das ferramentas externas e substitua `/caminho/para/dock6-6.11` pela raiz real do DOCK6, que contém `bin/` e `parameters/`.
4. Inicie o modo navegador ou desktop após conferir os executáveis.

```bash
conda env create -f environment.yml
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
2. Informe um nome que identifique o estudo; acrescente descrição, tags separadas por vírgulas e uma cor para o card.
3. Em **Pasta do projeto**, clique no botão de pasta dentro do campo e selecione uma pasta nova ou vazia. O campo não permite digitação manual. No desktop, o botão abre o seletor do sistema; no navegador, abre a navegação de pastas do computador que executa o BioMolExplorer e permite criar uma pasta.
4. Confirme a criação e abra o projeto.

Exemplo Linux: `/home/pesquisador/estudos/enzima-a`. No navegador, esse caminho pertence ao computador do backend. Um caminho do seu notebook só funciona se o backend também estiver no notebook ou tiver acesso à pasta. A aplicação precisa ter permissão de escrita. Pastas de outros projetos, pastas sobrepostas e o diretório das contas não podem ser reutilizados como pasta de um novo projeto.

A pasta escolhida guarda `project.json`, arquivos em `assets/`, execuções em `runs/`, versões em `.history/` e as cópias de trabalho do projeto. O banco central conserva contas, sessões, permissões e índices; logs gerais de infraestrutura ficam no diretório de logs configurado. Resultados e logs associados às execuções podem ser consultados pelo projeto. Reabrir a aplicação com o mesmo workspace conserva a associação à pasta.

Use busca por nome ou tags para localizar um projeto. Arquivar retira o projeto da lista ativa; restaurar permite voltar a usá-lo. Excluir remove o acesso normal ao projeto, preservando os arquivos no disco.

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

1. Abra **Arquivos**, escolha o tipo de dado e envie o arquivo. Cada arquivo pode ter até 200 MB pela interface.
2. Confirme que ele apareceu na biblioteca do projeto.
3. Adicione **Importar meus arquivos** ao pipeline ou use **Enviar meus arquivos** na configuração do bloco consumidor.
4. Escolha o tipo correto e os arquivos correspondentes. Aplicar a configuração confirma essa entrada; enviar um arquivo, sozinho, não o conecta a todas as etapas.
5. Confira o formato esperado mostrado pelo formulário. Para conjuntos preparados, selecione também os metadados necessários.

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
| Outros | Para ZINC, `.txt` com um endereço HTTPS autorizado por linha |

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

Entradas incompatíveis, ciclos e dependências de blocos desativados são recusados. **Avançado** permite editar JSON e templates; use-o quando precisar de opções que não estão no formulário. Preserve os nomes técnicos dos parâmetros e os placeholders exigidos pelos templates. Um erro de validação deve ser corrigido antes de aplicar.

**Ctrl+Z** desfaz, **Ctrl+Y** refaz e **Delete** remove a seleção quando nenhum formulário de configuração está aberto. **Escape** cancela uma conexão em andamento. Duplicar um bloco permite comparar parâmetros; renomeie as alternativas. Organização automática, navegação, zoom e recentralização ajudam a explorar fluxos extensos.

## 7. Entender as conexões disponíveis

| Consumidor | Entrada permitida no canvas |
| --- | --- |
| Expandir similares | Recuperação ChEMBL com os downloads originais; CSV genérico de compostos não substitui essa origem |
| Recuperar ZINC | Dados do tipo `other`, contendo a lista de endereços |
| Preparar meus complexos | Estruturas PDB brutas |
| Avaliar ADMET / Gerar fingerprints | Compostos de recuperação, importação, ADMET ou grafos |
| Calcular similaridade | Fingerprints compatíveis; não misture tipos ou larguras diferentes |
| Filtrar por grafos | Saída de **Calcular similaridade**; CSV externo pode ser escolhido no próprio formulário |
| Executar redocking | PDBs brutos com preparação ativada; receptores preparados com preparação desativada |
| Docking com Vina | Receptores preparados e compostos selecionados |
| Docking com DOCK6 | Receptores preparados, compostos selecionados e poses Vina correspondentes |
| Consenso de docking | Poses Vina e resultados DOCK6 com os mesmos identificadores |

Recuperar compostos e recuperar estruturas podem iniciar ramificações independentes. **Importar meus arquivos** publica o tipo escolhido. O modo **resultados prontos** de um bloco valida saídas já existentes e dispensa seu cálculo; isso não transforma qualquer CSV no resultado daquele bloco. Consulte a [matriz e os testes de validação](pipeline_validation.md) para as condições adicionais.

## 8. Executar, selecionar arquivos e escolher individual ou merge

1. Clique em executar o pipeline inteiro ou a etapa selecionada com suas dependências.
2. Se houver resultados concluídos salvos, responda à pergunta sobre reaproveitamento conforme a seção seguinte.
3. Acompanhe a janela de execução. Quando uma etapa precisar de entradas de seus blocos de origem, a execução pausa e abre **Selecionar arquivos**.
4. Marque somente os arquivos desejados em cada entrada. A lista mostra arquivos compatíveis das origens conectadas, sem selecionar automaticamente todos os resultados novos.
5. Escolha **Processar individualmente** ou **Mesclar arquivos (merge)**.
6. Para ajustar parâmetros, clique em **Configurar etapa**. As escolhas já feitas são incluídas no formulário, evitando outra seleção do mesmo conjunto.
7. Clique em **Continuar com os arquivos selecionados**. Repita a decisão nas próximas etapas que precisam de processamento.

| Escolha | Efeito | Exemplo |
| --- | --- | --- |
| Individual | Um processamento separado por arquivo, com saídas separadas | `alpha.csv` e `beta.csv` geram fingerprints próprios e continuam distinguíveis |
| Merge | Reúne os dados escolhidos para uma execução conjunta | Os compostos de `alpha.csv` e `beta.csv` participam da mesma análise |

Em operações com entradas distintas, como receptores e tabelas de compostos, individual forma combinações. DOCK6 limita essas combinações aos receptores/compostos com poses correspondentes. O consenso pareia Vina e DOCK6 pelo identificador da pose, evitando cruzamentos entre moléculas ou referências diferentes. Escolher merge não remove a necessidade de compatibilidade: tipos de fingerprint, códigos e metadados ainda precisam ser coerentes. Nomes iguais com conteúdos diferentes precisam ser desambiguados ou renomeados antes de combinar.

Metadados de estruturas acompanham as seleções em ambos os modos. **Selecionar depois** fecha a janela e mantém a pausa; **Configurar entradas e continuar** reabre a seleção. Fechar e reabrir a janela durante a sessão conserva arquivos marcados e modo. Escolhas confirmadas ficam na configuração do bloco e são reutilizadas ao abri-la posteriormente. Quando dois arquivos têm o mesmo nome, selecione o caminho que distingue cada resultado.

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

Em **Recuperar compostos**, informe nome ou ID ChEMBL do alvo, revise filtros de alvo/bioatividade/moléculas e escolha se deseja ampliar com PubChem. Limiar e máximo de similares controlam essa ampliação. Após executar, examine `compounds.csv` e os downloads ChEMBL disponíveis; use a tabela adequada como entrada da análise seguinte.

**Expandir similares** utiliza os downloads ChEMBL de uma recuperação existente. Conecte essa origem e revise limiar/máximo; não use um CSV comum como substituto do contexto ChEMBL. **Recuperar ZINC** exige a lista de endereços autorizados; selecione o `.txt` e examine os compostos produzidos.

**Recuperar estruturas PDB** oferece critérios como EC, organismo, resolução, ligante, tipo de polímero e método experimental. Confira `pdb_codes.csv`, ligantes e cadeias antes de preparar ou fazer redocking. Uma recuperação sem estruturas compatíveis exige revisar os filtros; não avance com um conjunto vazio.

### Preparação e redocking

Em **Preparar meus complexos**, selecione PDBs brutos e informe registros `[PDB, ligante, resíduo, cadeia]` ou forneça `pdb_codes.csv`. Configure pH e método de cargas. Confira os receptores, ligantes e centros produzidos.

Em **Executar redocking**, mantenha **Preparar complexos antes do redocking** ativado para PDBs brutos. Desative-o apenas ao fornecer o conjunto já preparado. Confira registros, dimensões da caixa, esforço de busca e número de poses. Examine os RMSDs e os logs; defina os critérios de aceitação do protocolo no estudo antes de usar seus receptores no docking de candidatos.

### Fingerprints, similaridade e grafos

**Gerar fingerprints** escolhe um tipo por bloco: Morgan, MACCS ou farmacóforo. Raio e bits aplicam-se ao Morgan. Para comparar tipos, duplique o bloco e mantenha os resultados identificados.

**Calcular similaridade** lê fingerprints e identifica o tipo das saídas nativas. Para arquivos próprios, confira o seletor de tipo. Métrica e limiar são definidos nessa etapa. A opção aproximada altera a estratégia de busca; compare-a com a análise exata conforme o protocolo.

**Filtrar por grafos** preserva as relações recebidas e seleciona o MCC, o maior componente conectado. Configure o tempo máximo do MCS e as opções de anéis quando desejar examinar o fragmento comum. Os SMILES vêm da linhagem que produziu as similaridades; para dados externos, forneça a tabela molecular correspondente quando quiser visualizar estruturas.

### ADMET

**Avaliar ADMET** calcula propriedades e aplica os filtros moleculares disponíveis. Examine tabelas, exclusões e o gráfico BOILED-Egg. BBB/HIA e os demais indicadores são estimativas e regras computacionais; não substituem validação experimental. Consulte o backend para conhecer os cálculos efetivamente implementados.

### Vina, DOCK6 e consenso

**Docking com Vina** recebe receptores preparados e a tabela de candidatos. Confira alvo/pasta, complexo de referência, caixa, pH, esforço e poses. Em individual, cada conjunto recebe saídas próprias. Nomes nativos das poses incluem PDB, ligante, resíduo, cadeia e composto para distinguir referências.

**Docking com DOCK6** exige receptores e auxiliares preparados, candidatos selecionados, poses Vina e uma referência explícita. Confira instalação, cargas, superfície, distância, raio, busca flexível/rígida e parâmetros de footprint. A seleção de candidatos restringe as poses que passam ao refinamento. Uma tabela de candidatos diferente da que produziu as poses precisa conter códigos correspondentes.

**Consenso de docking** recebe resultados Vina e DOCK6 correspondentes. Revise o peso de repulsão. São gerados scores e normalizações z-score/min-max; para um único composto ou scores constantes, a contribuição normalizada é zero, pois não existe dispersão para comparar. Preserve os scores brutos na interpretação. Arquivos sem score finito ou sem a pose correspondente são recusados.

## 12. Explorar resultados e navegar pelos grafos

1. Abra **Execuções** e expanda a etapa, ou use **Resultados** no bloco.
2. Se houver várias tabelas ou análises, selecione a desejada.
3. Use busca e paginação para localizar compostos; escolha 10, 25, 50 ou 100 itens por página quando disponível.
4. **2D** mostra a estrutura; **3D** gera um conformero local do SMILES, com rotação, zoom e recentralização. Essa visualização não é uma pose de docking.
5. Baixe o arquivo para conservar os dados completos da análise.

No grafo, escolha **Grafo completo** para ver todos os componentes. Cada componente recebe espaço próprio e os vértices são separados. Use o seletor de componentes para focar uma região; clique em um vértice para ver seus detalhes, vizinhos e pesos das relações. Buscar um código permite encontrar o nó e navegar até seu componente. Use arraste e zoom; **Ajustar à tela** e **Recentrar** ajudam a recuperar a visão geral.

Ao escolher **MCC**, a exploração se restringe ao maior componente conectado. O painel mostra **SMILES do fragmento**, sua estrutura e o estado da busca do fragmento comum quando esses dados estão disponíveis. O SMILES é extraído da molécula de referência para o padrão encontrado por MCS; a conectividade comum não comprova identidade estereoquímica de todos os compostos. Busca interrompida pelo limite de tempo deve ser interpretada pelo estado informado. Sem SMILES correspondentes, a topologia pode ser explorada, mas não há estrutura molecular do fragmento a mostrar.

Resultados antigos podem ter apenas PNG. Para produzir os artefatos interativos, execute novamente escolhendo **Não, executar novamente**. Visualizações `*.biomol-view.json` acima do limite de abertura de 32 MB podem ser baixadas. Em BOILED-Egg, pontos sobrepostos oferecem seleção dos compostos; busca também permite distingui-los.

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
