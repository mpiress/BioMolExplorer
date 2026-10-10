# Workspace Flet

[Documentação](README.md) · Português · [English](en/frontend.md)

A interface reúne autenticação, projetos privados, convites, biblioteca de arquivos,
editor de pipeline e acompanhamento das execuções. Ela usa os serviços Python do
backend; os cálculos científicos continuam em processos separados.

## Projetos portáteis e colaborativos

A versão atual solicita a pasta do projeto, oferece **Importar projeto**, **Exportar**
e **Histórico**, com restauração de versões completas. Configurações são salvas
automaticamente e atualizações de colaboradores aparecem no projeto aberto. Cada
bloco pode receber várias entradas e usar resultados prontos validados. Consulte
[o guia de projetos e versões](projects.md) para o passo a passo e os formatos.

Ao criar um projeto, use o botão de pasta dentro de **Pasta do projeto** para
selecionar a pasta principal; o nome deve estar preenchido e será usado para criar a subpasta do projeto. O campo não permite digitação manual.
No desktop e na versão web, o botão abre a navegação de pastas visíveis
no computador que executa o BioMolExplorer. Criar pasta solicita um nome e seleciona automaticamente a nova pasta no local atual. A subpasta com o nome do projeto recebe entradas, configurações,
resultados e histórico. Destinos existentes exigem confirmação para substituição ao salvar. O caminho continua associado ao projeto ao reabrir a aplicação.

## Iniciar

Antes de iniciar, baixe o código do [GitHub](https://github.com/mpiress/BioMolExplorer) e instale **UCSF Chimera 1.17 e DOCK6 6.11** no computador dos cálculos. Siga o [guia de instalação e configuração](installation.md) para preparar o ambiente, configurar o `PATH` e conferir os executáveis. Instalar a interface não instala essas ferramentas; etapas que dependem delas falharão se estiverem ausentes.

Na raiz do código baixado, com o ambiente científico preparado:

```bash
conda activate BioMolExplorer
python -m pip install -e '.[ui]'
biomolexplorer-ui --web --language pt --dock6-path /caminho/para/dock6-6.11
```

Acesse `http://127.0.0.1:8550`. Para servir sem abrir automaticamente o navegador:

```bash
biomolexplorer-ui --web --no-browser --data-dir ./workspace-data --port 8550 --dock6-path /caminho/para/dock6-6.11
```

`--worker-python /caminho/do/python` seleciona o ambiente dos cálculos.
`--dock6-path /caminho/do/dock6` define a instalação DOCK6 disponível às etapas.
Chimera, Open Babel, Vina, DOCK6 e suas dependências continuam necessários para os
respectivos protocolos. A instalação da interface não instala esses motores.

O diretório de dados padrão é `~/.local/share/biomolexplorer`. Ele contém o banco
SQLite de contas/projetos e a área temporária de uploads. Entradas, configuração,
versões e resultados ficam na pasta escolhida para cada projeto. Use um diretório persistente
e faça backup dele com a aplicação encerrada. Uma instância supervisiona cada
diretório; dois processos não podem executar o mesmo workspace simultaneamente.

## Escolher o idioma

Use `--language pt` ou `--language en`; sem essa opção a interface inicia em inglês. Na barra superior, o menu de bandeiras no canto direito seleciona português (Brasil) ou inglês (Estados Unidos) para a sessão, preservando os campos preenchidos. Nomes do usuário, arquivos, dados moleculares e chaves de parâmetros conservam seus valores originais. Para trocar durante o uso, abra o mesmo menu de bandeiras, sem sair da conta. Consulte o [manual detalhado](user_manual.md).

## Contas e colaboração

Crie uma conta na tela inicial. A senha tem pelo menos 10 caracteres e é armazenada
com salt individual e scrypt. A sessão dura oito horas e é invalidada ao sair.

Cada projeto começa privado. Na aba **Compartilhar**, o proprietário informa o
e-mail de uma conta já cadastrada e escolhe a permissão. O convite aparece no
workspace do destinatário, que precisa aceitá-lo. O proprietário pode revogá-lo.
O convite informa se o papel é **Editor** ou **Leitor**. Convites recusados ou
revogados não podem ser aceitos depois. Mudar o papel de um colaborador cria um
novo convite e suspende o acesso anterior até a aceitação. Repetir um convite com
a mesma permissão de um colaborador ativo não remove seu acesso.
Se o papel mudar enquanto o convite estiver na tela, a aplicação mostra o convite
atualizado antes de permitir sua aceitação. A revogação é conferida nas consultas:
ao detectar perda de acesso, a interface fecha as visualizações e retorna ao
workspace. Uma sessão expirada retorna à autenticação. Respostas de consultas
iniciadas numa sessão anterior não reabrem projetos ou gráficos após sair da conta.

| Permissão | Capacidades |
| --- | --- |
| Leitor | Consultar pipeline, arquivos, resultados e logs; baixar artefatos |
| Editor | Além das consultas, configurar etapas, enviar arquivos e executar/cancelar |
| Proprietário | Além da edição, convidar/revogar colaboradores e excluir o projeto |

Os convites são internos à aplicação. Não há envio de e-mail, validação de posse
do endereço, recuperação de senha ou integração com um provedor institucional de
identidade nesta versão. Para acesso institucional pela internet, configure HTTPS,
identidade verificada e controles de infraestrutura antes de publicar o serviço.
Os workers locais compartilham o usuário do sistema operacional: o isolamento
implementado é de autorização e caminhos na aplicação, não de containers.

## Projetos e pipelines

Crie projetos com nome, descrição e cor escolhida na paleta. Busque pelo nome ou descrição.
Projetos podem ser arquivados, restaurados ou excluídos. A exclusão exige confirmação mostrando o caminho e remove permanentemente o projeto, seus registros e a pasta completa no disco. Conclua ou cancele execuções ativas antes de excluir.

Na aba **Pipeline**, clique no título de uma categoria da biblioteca para expandir seus blocos; as categorias começam recolhidas. A busca expande as categorias com resultados. Arraste os blocos para a área
quadriculada, ou use **+**. Arraste o cabeçalho para reposicionar um bloco. Conecte
com o mouse a saída de um bloco à entrada de outro; também é possível clicar na
saída e depois na entrada. Cores identificam recuperação, análise, docking e
importação. Entradas incompatíveis e ciclos são rejeitados sem alterar o fluxo.

Dê dois cliques no cabeçalho ou use **Configurar seleção** para abrir a configuração.
O modo **Visual** oferece formulários, listas, seletores e controles para filtros
ChEMBL, registros PDB, caixa de docking e opções dos motores. Os scripts são
gerados a partir dessas escolhas. O modo **Avançado** permite editar parâmetros e
conexões em JSON e os templates científicos completos. **Aplicar configuração**
valida e confirma a edição; **Cancelar** preserva o bloco anterior.

Exceto pelo redocking com pares previamente configurados e pelo docking Vina/DOCK6 com arquivos de entrada já selecionados, ao concluir uma etapa o pipeline pausa antes do próximo bloco conectado e abre
automaticamente **Selecionar arquivos**. Marque os arquivos desejados em cada
entrada e clique em **Continuar com os arquivos selecionados**. O popup lista
somente resultados compatíveis das origens conectadas nesta execução. Em **Como
processar os arquivos selecionados?**, escolha **Processar individualmente**
(padrão do popup) ou **Mesclar arquivos (merge)**. Individual executa cada arquivo
em separado, preserva seu nome e disponibiliza resultados separados para a próxima
seleção. Merge reúne os arquivos selecionados em uma entrada. Com vários tipos de
entrada (por exemplo, receptores e compostos), individual executa cada combinação
de arquivos. DOCK6 restringe as combinações aos receptores/compostos/poses correspondentes; o consenso sempre reúne as seleções e calcula a interseção por código e receptor. Metadados necessários acompanham as estruturas escolhidas em ambos
os modos. A escolha pode mudar a cada etapa e fica registrada na execução.

**Selecionar depois** fecha o popup e mantém a execução pausada. Para reabrir,
use **Configurar entradas e continuar** no acompanhamento ou **Configurar [nome do
bloco]** na aba **Execuções**. **Configurar etapa** permite ajustar parâmetros e
selecionar outras entradas antes de continuar, já com as escolhas feitas no popup
preenchidas. Fechar e reabrir o popup conserva suas seleções e o modo de processamento.
As escolhas confirmadas também ficam salvas no bloco. A retomada preserva os resultados
concluídos e a ordem definida pelas conexões. Após a última etapa, os resultados
ficam disponíveis na aba **Execuções**.

A seleção é solicitada para etapas que precisam ser processadas, inclusive em projetos
antigos. Etapas reaproveitadas são concluídas sem repetir a seleção. Uma pausa continua disponível ao reabrir a aplicação e
pode ser cancelada. Leitores acompanham o fluxo; editores selecionam e retomam.

Durante a execução, as tabelas moleculares passam por uma limpeza comum às etapas.
Registros com SMILES ausente ou inválido, códigos ambíguos, fingerprints malformados
ou valores de similaridade inválidos são excluídos das cópias usadas no processamento.
Os arquivos originais são preservados. Nos grafos, os SMILES vêm dos fingerprints
que produziram a similaridade; relações sem um composto correspondente são removidas.
Arquivos externos de similaridade sem tabela de compostos continuam permitindo a
análise da topologia, sem estruturas moleculares.

Quando houver exclusões, **Execuções** mostra a contagem do bloco e disponibiliza
`molecule_exclusions.json`, com arquivo de origem, linha, códigos e motivo de cada
registro excluído. A limpeza é reaplicada nas reexecuções. Cabeçalhos incorretos,
arquivos indisponíveis e falhas de serviços ou motores continuam sendo informados.

Use zoom, navegação da área, organização automática, duplicação, desfazer e refazer.
O pipeline, os gráficos e as imagens permitem zoom de 0,1% a 10.000%, inclusive
quando o conteúdo fica menor que a área visível. Use a recentralização para
restaurar a visão inicial. O controle de zoom das moléculas 3D usa a mesma faixa,
com escala logarítmica para facilitar ajustes pequenos e grandes.
**Ctrl+Z**, **Ctrl+Y** e **Delete** atuam no editor quando a configuração está fechada;
**Escape** cancela a conexão em andamento. O histórico de desfazer é da sessão de
edição; posições e conexões são salvas automaticamente após uma breve pausa.
**Salvar** também confirma as mudanças. O menu de seleção e os
botões oferecem alternativas ao arraste. A posição visual não define a ordem de
execução: as conexões definem as dependências. Uma etapa ativa não pode depender
de outra desativada.

Você pode executar o pipeline inteiro ou a etapa selecionada com suas dependências.
Cada execução registra uma cópia das configurações. Se houver resultados salvos,
a aplicação pergunta **Reaproveitar dados e resultados?** antes de iniciar. **Sim,
reaproveitar** considera concluídas as etapas compatíveis, inclusive depois de
reiniciar a aplicação ou acrescentar novos blocos. O sistema
compara configurações, templates, entradas e integridade dos arquivos. Mover ou
renomear blocos não repete os cálculos; mudar parâmetros ou entradas recalcula as
etapas afetadas. **Não, executar novamente** inicia os cálculos sem usar o cache.
Os resultados reaproveitados continuam no diretório da execução que os produziu.

Ao iniciar, uma janela mostra a etapa atual, o tempo decorrido, os estados de todos
os blocos e quantas etapas foram concluídas. Na recuperação ChEMBL, ela informa
também a fase da consulta e a contagem dos registros moleculares processados.
O indicador animado representa trabalho em andamento, sem estimar uma porcentagem
para cálculos cuja duração não é conhecida. **Minimizar** mantém o acompanhamento
em uma faixa do projeto; **Acompanhar execução** reabre a janela. Ao terminar, a
janela apresenta conclusão, cancelamento ou falha, com acesso ao log e aos resultados.
Editores e proprietários podem solicitar **Cancelar execução**. O pedido permanece
visível enquanto o backend interrompe o processo; leitores apenas acompanham.

A aba **Execuções** mostra estados, erros, logs e arquivos para download. O
cancelamento interrompe o worker; após uma interrupção da aplicação, a execução
anterior fica registrada como interrompida. Uma nova execução reaproveita as etapas
concluídas compatíveis e executa as demais. O popup e o histórico identificam as
etapas **Reaproveitadas**.
A cada dois segundos, quando não há formulário aberto nem rascunho local, a interface carrega também as últimas
alterações dos colaboradores. Se a permissão passar para leitor, rascunhos de
edição e janelas privadas são descartados quando a mudança é detectada.

Os modelos incluem recuperação → ADMET, seleção por grafos → ADMET, compostos
próprios → ADMET, PDBs próprios → preparação e um percurso completo até consenso.
Eles são pontos de partida: revise alvos, entradas e escolhas científicas antes de
executar. DOCK6 exige a instalação configurada e um complexo de referência.

## Explorar gráficos e compostos

Na aba **Execuções**, expanda **Recuperar compostos** ou **Expandir similares**:
a tabela aparece diretamente, sem uma janela intermediária. O campo **Tabela de compostos** lista os CSVs disponíveis:
`<alvo>_FULL`, `<alvo>_MOLS`, `<alvo>_SIMS` e `compounds.csv`, o conjunto integrado
ChEMBL/PubChem. A tabela mostra código e SMILES, com busca, paginação de 25 linhas
e download do arquivo selecionado. **Itens por página** permite escolher 10, 25,
50 ou 100 linhas; o padrão é 25. A lista contém somente arquivos produzidos;
etapas sem similares podem não gerar todos os CSVs.

**2D** abre a estrutura molecular. **3D** gera um conformero local a partir do
SMILES, com rotação por arraste, zoom e recentralização; não representa uma pose
de docking. Leitores também podem consultar essas visualizações.

**Remover** solicita confirmação e exclui apenas a linha do arquivo selecionado,
preservando suas demais colunas. A ação exige permissão de editor/proprietário e
nenhum pipeline ativo no projeto. A cópia original e o registro de autoria ficam
no histórico privado de curadoria (`.curation/` e tabela SQLite `compound_edits`).
A coleta permanece reaproveitável após a curadoria; etapas que consomem os dados
alterados serão recalculadas na próxima execução. A remoção não sincroniza outros
CSVs: escolha na conexão do próximo bloco qual tabela deverá ser utilizada.

Após executar uma etapa, clique no ícone **Resultados** do bloco ou abra
**Execuções** e expanda a etapa. Os resultados de fingerprints, similaridade e
demais listas de arquivos aparecem em tabelas paginadas, com ações de download,
visualização quando disponível e remoção à direita. A aba **Arquivos** segue o
mesmo padrão. A remoção exige permissão de edição e nenhum pipeline ativo;
entradas associadas a um bloco precisam ser desvinculadas antes da exclusão.
O histórico permite restaurar os arquivos removidos. Uma etapa com saídas removidas
não é reaproveitada como um resultado completo na próxima execução.

No bloco **Filtrar por grafos**, a única entrada do canvas é **Similaridades**.
Conecte um ou mais blocos **Calcular similaridade** e selecione os arquivos
produzidos por cada um. Cada CSV gera uma análise independente. Não é possível
conectar fingerprints, recuperação de compostos ou blocos genéricos de importação.
Métrica, tipo de fingerprint e limiar são configurados em **Calcular similaridade**;
o bloco de grafos preserva as relações e os pesos recebidos.

Opcionalmente, use **Enviar meus arquivos** para fornecer CSVs externos UTF-8 com
as colunas `source,target,value`, códigos válidos e pesos numéricos finitos entre
0 e 1. Exemplo: `MOL1,MOL2,0.85`. Cabeçalhos incorretos, códigos vazios e pesos
inválidos são rejeitados com uma mensagem explicando o padrão esperado.

Os SMILES das entradas do pipeline são identificados automaticamente. Para
entradas externas, **Arquivo externo de compostos e SMILES · opcional** aceita
uma tabela com `molecule_chembl_id,canonical_smiles`, usando os mesmos códigos
das arestas. Essa tabela permite a visualização 2D, a busca do fragmento comum e
a inclusão de compostos isolados. Sem ela, o grafo apresenta a topologia e informa
que as estruturas estão indisponíveis; seus CSVs não podem alimentar etapas que
exijam SMILES até que a tabela correspondente seja fornecida.

Projetos antigos preservam as conexões de similaridade e os resultados salvos.
Conexões diretas de fingerprints e parâmetros de cálculo do bloco de grafos são
removidos ao abrir o projeto. Se não houver uma conexão de similaridade,
conecte **Calcular similaridade** antes de executar. Exportações antigas recebem
o mesmo ajuste ao serem importadas.

Expanda a etapa em **Execuções** e escolha **Análise de grafos · entrada e origem**.
O grafo aparece na própria tela. Selecione **Grafo completo**, incluindo compostos
isolados, ou **Máximo componente conectado (MCC)**. O contorno azul evidencia os
nós do MCC no grafo completo. Os vértices têm tamanho pequeno e uniforme,
inclusive na apresentação MCC exportada. As cores Viridis indicam o grau (número
de conexões); o seletor também oferece conectividade normalizada, definida como
grau / (n − 1) para o conjunto exibido. A escala explica o intervalo. Essa medida
por nó é diferente da densidade global de arestas, registrada no resultado.
Os componentes aparecem separados, com distância mínima entre os vértices.
Redes maiores expandem a área navegável. Em **Navegar entre componentes**, escolha
um componente para explorar seus compostos e relações; **Todos os componentes**
retorna à rede completa. Arraste o fundo, amplie, recentre, use **Ajustar à tela**
para ver o conjunto ou escolha a organização em anéis. Ative os códigos junto
aos nós quando desejar. Passe o mouse para identificar o composto; o clique
abre sua molécula 2D e suas propriedades, destaca as conexões e oferece links
para navegar até cada vizinho com seu valor de similaridade. A busca por código
localiza o composto mesmo quando outro componente está selecionado.

Ao selecionar o MCC, o painel à direita mostra o fragmento molecular comum,
sua imagem, o **SMILES do fragmento**, os números de átomos e ligações e o estado
da busca. O SMILES e a imagem representam o fragmento encontrado em uma molécula
de referência; o SMARTS preserva as condições da correspondência no arquivo de
resultado. Os resultados salvos em versões anteriores recebem o novo layout e
exibem o SMILES ao serem abertos, preservando os arquivos originais.
O tempo máximo (30 segundos
por padrão) e as opções de comparação de anéis são
configuráveis no bloco. A busca exige correspondência em todos os compostos do
MCC. Se atingir o limite, a interface e a figura identificam o melhor fragmento
encontrado como parcial; o tamanho máximo não fica confirmado. Essa distinção
segue a [documentação do RDKit FindMCS](https://www.rdkit.org/docs/source/rdkit.Chem.rdFMCS.html).

**Baixar apresentação MCC** exporta um PNG com o grafo por grau, a imagem do
fragmento, o ranking, o histograma e a distribuição de graus. **Baixar compostos
do MCC** exporta o CSV da análise selecionada, preservando seus códigos e
metadados. **Arquivos desta análise** reúne os arquivos associados ao resultado escolhido
em tabela paginada,
com 25 itens por página inicialmente. Cada análise tem seu próprio CSV, figura,
modelo interativo e arestas. `Molecules/molecules.csv` mantém a união dos MCCs
para pipelines existentes; se códigos iguais identificarem estruturas diferentes
entre análises, essa união recebe códigos com sufixo e uma coluna `original_code`.
Os resultados individuais preservam os códigos originais. Escolha um resultado
específico ao conectar a próxima etapa para trabalhar com apenas uma análise.

Em empates de tamanho, o MCC é escolhido de forma determinística pelo código
dos compostos. Sem arestas, um composto isolado forma um MCC de tamanho um.
Redes com mais de 1.000 nós usam layout circular para reduzir o custo de
organização, preservando todos os nós e relações.

Nos resultados de **Avaliar ADMET**, escolha **Arquivo de resultados ADMET** e
**Compostos exibidos**: todos os compostos avaliados, BBB+, BBB− ou HIA+.
O EGG aparece na própria etapa, com TPSA × WLOGP e regiões HIA/BBB.
Os limites incluem pontos fora da faixa padrão. O clique em um ponto abre um
popup com a molécula 2D e suas propriedades. Os botões à direita baixam o gráfico
EGG ou o CSV correspondente ao filtro selecionado; o CSV filtrado preserva as
propriedades dos compostos. Os arquivos JSON continuam armazenados internamente
para a exploração interativa, sem aparecer na lista de arquivos do ADMET.

**Gerar fingerprints** oferece um único tipo por bloco: Morgan, MACCS ou
farmacóforo. Raio e número de bits aparecem somente ao escolher Morgan.
**Calcular similaridade** identifica e bloqueia o tipo de fingerprint a partir
do arquivo gerado selecionado na entrada. Para arquivos próprios, o seletor
permanece disponível. Entradas de tipos diferentes não podem ser combinadas.
Em projetos antigos que geravam vários tipos no mesmo bloco, selecione um arquivo
específico antes de configurar a similaridade.

Em ambos os visualizadores:

1. Passe o mouse sobre um ponto para ver o código do composto.
2. Clique para abrir sua estrutura 2D, SMILES e propriedades disponíveis.
   Quando vários compostos ocupam o mesmo ponto, uma lista permite escolher qual
   deles consultar; o indicador `(+N)` informa os compostos sobrepostos.
3. Busque um código para selecionar o composto pelo teclado, inclusive quando
   pontos coincidem; no MCC, a busca se restringe aos compostos desse componente.
4. Arraste o fundo para navegar, use o zoom e **Recentrar** para voltar ao início.

Leitores também podem explorar e baixar gráficos; editar e executar continuam
restritos a editores e proprietários. A autorização é verificada ao abrir o
arquivo e ao consultar detalhes. As visualizações são arquivos versionados
`*.biomol-view.json`, gerados junto com os PNGs, sem serviços externos ou execução
de HTML/scripts. O limite de abertura é 32 MB; arquivos maiores podem ser baixados.
Resultados antigos mantêm os PNGs: escolha **Não, executar novamente** na pergunta de reaproveitamento e
execute a etapa para gerar os artefatos interativos.

## Usar arquivos próprios

Os formulários de cada bloco oferecem **Enviar meus arquivos**, **Adicionar entrada**
e um modo separado para **resultados prontos**. A mensagem de validação informa o
formato esperado. O bloco de importação continua disponível como alternativa:

1. Adicione **Importar meus arquivos** ao pipeline e abra sua configuração.
2. Escolha **Tipo dos arquivos a adicionar** e clique em **Selecionar arquivos no disco**. Selecione todos os arquivos necessários; cada arquivo é validado antes de entrar no projeto.
3. Para aproveitar arquivos existentes, use **Arquivo já enviado ao projeto** ou **Adicionar todos os arquivos do projeto**.
4. Confira a tabela **Tipo / Arquivo / Ações**. Você pode alterar o tipo por linha (com nova validação) ou usar o botão **Remover** na coluna **Ações**. A remoção dessa lista preserva o arquivo na biblioteca do projeto.
5. Para receptores preparados, inclua também os metadados e formatos complementares. Configure a pasta do alvo para estruturas e aplique a configuração. Um bloco pode publicar vários tipos, conectados às entradas correspondentes dos consumidores.

Cada arquivo pode ter até 200 MB. No navegador, o upload usa uma URL assinada e
temporária; o arquivo entra na biblioteca apenas após nova verificação da sessão,
permissão e tamanho. O servidor limita a transferência pelo tamanho configurado,
com nova validação antes de incorporar o arquivo ao projeto. Os arquivos não são
servidos por uma pasta pública de assets; downloads passam pela autorização.

Compostos aceitam CSV com `canonical_smiles` e `molecule_chembl_id`, ou `smiles` e
`name`. A importação normaliza os nomes das colunas e gera identificadores quando
faltam. Identificadores devem ser seguros para nomes de arquivos. O original
enviado é preservado. Para filtrar por um resultado específico, informe o nome do
CSV no seletor **Resultado usado nesta entrada**. Ele lista arquivos importados e
resultados das execuções anteriores; durante a execução, confirme os arquivos no popup. O modo
avançado permite especificar outros seletores (`selector`).

Para PDBs próprios, use importação **Complexos PDB**, seguida de **Preparar para docking**, e adicione os complexos no formulário (PDB, ligante, resíduo e cadeia). Você também
pode enviar `pdb_codes.csv` com `PDB_CODE,LIGAND,RESNUM,CHAIN`. Assim, não é necessário
incluir recuperação PDB ou redocking.

Para receptores já preparados, escolha **Receptores preparados** e envie o conjunto
esperado pelo protocolo existente: arquivos `<PDB>_<CHAIN>.dockprep.pdbqt`,
`centers.csv` e `pdb_codes.csv`. A importação cria `<alvo>/Prepared` e coloca os
metadados em `<alvo>/pdb_codes.csv`. Conecte essa importação diretamente ao Vina e
conecte os compostos à outra entrada. Um PDB bruto sozinho ainda precisa passar
pela preparação molecular.

Fingerprints, relações de similaridade, resultados Vina e DOCK6 também podem ser
importados. Eles devem respeitar os formatos do backend: fingerprints possuem
`molecule_chembl_id,fingerprint`, relações possuem `source,target,value`, resultados
Vina incluem `.lig.pdbqt` e DOCK6 inclui `_scored.mol2`. O consenso aceita entradas
separadas dos dois motores, inclusive vindas de importações.

## Diagnóstico e recuperação ChEMBL

A tela inicial usa os logos originais de `imgs`, empacotados com a aplicação, e
adapta a apresentação a telas menores. A configuração de cada bloco é organizada
em cartões de identificação, entradas e parâmetros. Filtros e opções dos motores
abrem em seções com espaçamento próprio e campos responsivos.

Os erros são gravados em `logs` na raiz do repositório. Na versão instalada fora
do repositório, o padrão é `logs` no diretório de execução. Para escolher outro
local, use `--log-dir /caminho/logs` na interface ou `BIOMOL_LOG_DIR`.

| Arquivo | Conteúdo |
| --- | --- |
| `logs/frontend.log` | Falhas de ações, autenticação e eventos da interface |
| `logs/backend.log` | Falhas do supervisor e dos pipelines, com IDs de execução |
| `logs/errors.log` | Erros agregados do processo principal, com traceback |
| `logs/jobs/<job-id>/` | Logs científicos, erros e cópia de `execution.log` de cada worker |
| `logs/chembl-check.log` | Diagnóstico do teste de consulta ao vivo |
| `logs/chembl-check-report.json` | Resultado mais recente do teste, com status, alvo e erro/contagens |

Os arquivos de diagnóstico rotacionam a cada 5 MB, mantendo três cópias anteriores.
A cópia de `execution.log` é feita quando o worker termina; o log privado original
continua disponível na aba **Execuções** durante a execução. Não publique `logs`
como assets: os diagnósticos contêm caminhos e identificadores de projetos.

A recuperação ChEMBL usa endpoints REST paginados, com timeout e tentativas
limitadas, sem depender da inicialização do endpoint `/spore`. Uma indisponibilidade
da API continua sendo informada como erro; ela não é convertida em um conjunto vazio.
Os limites de atividade usam `standard_value` e `standard_units`; a seleção de
produtos naturais respeita tanto “Sim” quanto “Não”.

O bloco **Recuperar compostos** chama `wrappers.crawlers.retrieve_compounds` em um
worker Python separado. O arquivo `workflow/1-InformationRetrieval/retrieve_compounds.py`
é um exemplo manual que chama a mesma função; suas opções e pasta de saída podem
ser diferentes das configurações do projeto. O início do worker registra o
interpretador utilizado em `logs/jobs/<job-id>/backend.log`. `request.json`, no
diretório privado do job, registra a operação e os parâmetros enviados. Os filtros
do bloco também podem vir dos templates salvos com a etapa.

Um HTTP 500 recebido do endpoint ChEMBL ou um timeout durante a coleta é propagado
ao frontend com o contexto da consulta, após tentativas automáticas. Uma falha
depois de baixar alguns arquivos não significa que a recuperação foi concluída:
os arquivos parciais permanecem no diretório privado daquela execução, e as etapas
dependentes não são executadas. O popup e **Ver log da etapa** ajudam a distinguir
essa situação de uma falha de inicialização do worker.

Para verificar CHEMBL220 com apenas IC50, sem expansão por similaridade:

```bash
PYTHONPATH=src python scripts/check_chembl.py --target CHEMBL220 --limit 5
```

Esse teste baixa uma amostra de até cinco atividades IC50 em nM, até 5000 nM e
com pChEMBL, e seus compostos. Exporta as tabelas em `datasets/chembl-check/CHEMBL220`
quando a API responde. HTTP 500, timeout e ausência de resultados são registrados
no relatório e fazem o comando terminar com código 1. O comando não depende dos
testes offline, que usam respostas controladas para verificar filtros e paginação.

## Configuração avançada e extensão

Cada estágio apresenta os argumentos do wrapper correspondente em formulários; o
editor JSON e os scripts ficam no modo avançado. Entradas e saídas são administradas pelo projeto para preservar seu escopo.
Filtros ChEMBL e templates de Chimera, Vina e DOCK6 podem ser editados por estágio;
o worker recebe uma cópia própria. Preserve os marcadores de formatação. A edição
aceita configurações científicas e rejeita comandos externos e caminhos livres.

O catálogo reside em `biomolexplorer/catalog.py`; contratos e execução científica
em `operations.py`; autenticação/arquivos em `workspace.py`; ordenação e execução
em `pipeline.py`; interface em `ui/app.py`, canvas em `ui/flow_canvas.py`, formulários em
`ui/guided.py` e contratos visuais em `flow.py`. Uma nova operação registrada pode usar
o mesmo editor de parâmetros, com título e ajuda adicionados ao catálogo.

O executor atual usa SQLite e uma fila local limitada, com duas execuções
simultâneas e no máximo 20 execuções pendentes/ativas. Para múltiplos servidores,
a interface pode conservar esses serviços como contrato, substituindo o
armazenamento e o supervisor por serviços distribuídos. Esta versão não oferece
uma API HTTP pública dos serviços científicos.

Para consultar os contratos e a evidência de testes, veja [validação do pipeline](pipeline_validation.md).


## Paleta, pastas e estruturas recuperadas

Escolha a cor por amostra visual; o formulário não apresenta Tags. O navegador de pastas é o mesmo no desktop e na versão web e oculta diretórios cujo nome começa com ponto e os sinalizados como ocultos pelo sistema. Em **Criar pasta**, informe somente o nome: a pasta é criada dentro do local atual e selecionada automaticamente. Por exemplo, em `/home/michel/Downloads`, criar `teste` seleciona `/home/michel/Downloads/teste`. Para novos projetos, selecione a pasta principal após preencher o nome: o destino será uma subpasta com o nome do projeto. Uma subpasta existente exige confirmação para substituição ao salvar.

A listagem de PDBs oferece edição de ligantes, visualização local e acesso ao RCSB. Os metadados são ocultados na interface, preservando seu uso no pipeline. Consulte [recuperação](retrieval.md) para os filtros ChEMBL e regras de curadoria.


## Explorar um PDB em 3D

O botão **Visualizar estrutura 3D** abre uma aba com o explorador WebGL local. Arraste para rotacionar, use a roda para zoom e o botão direito para deslocamento. O painel oferece representações, cores, cadeias, modelos, ligantes e água; a barra superior oferece tela inteira e exportação PNG. Consulte [o guia de recuperação](retrieval.md) para controles, permissões e requisitos do navegador. Nenhum PDB é enviado a serviços externos.

## Verificações antes da execução

A seleção de pares é obrigatória e usa uma única **Cadeia** para receptor e ligante. Um redocking já configurado reutiliza os pares escolhidos; o pipeline só pede a seleção quando ela ainda não existe. Resolução vem dos metadados, sem campo editável. Cofatores e opções de solvente, hidrogênios, minimização e cargas são configurados por par. A preparação e a conformação do ligante usam as mesmas opções.

Antes de calcular, o sistema verifica os resíduos, a cadeia e os cofatores nos PDBs. Entradas já preparadas precisam dos PDBQT e de três coordenadas finitas em `Prepared/centers.csv`. O ambiente científico deve disponibilizar Chimera, Open Babel (`obabel`) e Vina; sem preparação, somente Vina é exigido. O diretório do interpretador científico também integra o PATH dos processos.

Falhas de ferramentas mostram o executável, o código de saída e a mensagem do processo. Scripts Chimera que falharam são preservados para diagnóstico. Corrija a entrada ou a instalação indicada e tente novamente. Não há controle de verbose: a verbosidade do Vina permanece em zero.

Consulte [Configuração do redocking](redocking_configuration.md) para o fluxo completo.

## Consultar resultados de redocking

Ao expandir uma etapa concluída em **Execuções**, a tabela mostra PDB, ligante, resíduo, cadeia e **RMSD (Å)**. **Ver simulação** abre um popup com os arquivos correspondentes àquela seleção. A lista distingue o receptor preparado, o ligante de referência, as poses e os metadados. Arquivos compartilhados pelo receptor ou pela coleção acompanham as simulações correspondentes.

Use o botão de download de cada linha para baixar um arquivo, ou **Baixar todos (ZIP)** para obter o conjunto da simulação. O ZIP conserva as pastas, inclusive quando o ligante de referência e a saída do Vina têm o mesmo nome. Arquivos de outras simulações ficam fora desse conjunto.

O botão **Visualizar estrutura 3D** abre o visualizador no navegador, como na consulta dos PDBs. Ele aceita PDB, PDBQT e MOL2 para examinar receptor, ligante e poses antes do download. Ao abrir um arquivo com várias poses, use **Modelo** no visualizador. A cena aberta pelo botão **3D** da tabela seleciona automaticamente a pose de menor score, ou a primeira sem score; essa cena não oferece troca de modelo. As legendas e ações acompanham o idioma selecionado. Leitores também podem consultar e baixar os resultados.

Consulte [Logs e diagnóstico](logging.md) para o formato comum, contexto por execução, códigos de falha, resumo do job e o comando `python -m biomolexplorer.log_report`.


## Entradas e resultados de docking

**Docking com Vina** e **Docking com DOCK6** recebem compostos de qualquer bloco que publique moléculas: ChEMBL, PubChem, ZINC, ADMET, fingerprints, grafos, consenso ou importação do usuário. Selecione a tabela na entrada de compostos e conecte separadamente o receptor preparado. Não é obrigatório passar por grafos, ADMET ou pelo outro motor. PubMed fornece referências bibliográficas; moléculas obtidas dessas referências devem ser importadas com código e SMILES ou estrutura válida.

Para reutilizar conformações, conecte **Vina → DOCK6** ou **DOCK6 → Vina** na entrada de compostos, escolhendo `docking_results.csv` ou a pose desejada. A identidade e o SMILES são preservados; a conversão utiliza a conformação selecionada, sem gerar outra a partir do SMILES. Sem uma pose de entrada, a preparação gera um conformero 3D. DOCK6 usa o centro do sítio do receptor para selecionar esferas e posicionar ligantes recém-gerados; a entrada opcional de poses Vina continua disponível para refinamento de candidatos correspondentes.

Confira a caixa, pH e esforço do Vina e as cargas, superfície, distância, raio, busca flexível/rígida e footprint do DOCK6. Ambos exigem receptores preparados, metadados e centros; DOCK6 também exige MOL2 do receptor e PDB sem hidrogênios. Seus resultados apresentam código molecular, receptor, score, SMILES, **3D** da pose calculada e **Remover**. O botão abre a estrutura obtida no docking.

**Consenso de docking** reúne os resultados selecionados de cada motor, incluindo vários lotes e ramificações. Calcula somente a interseção por código molecular e receptor, usando o menor score quando há poses repetidas. Resultados do mesmo código com estruturas diferentes são recusados. Não havendo interseção, o bloco é marcado como não executado e explica que os motores não avaliaram compostos em comum para o mesmo receptor.

A tabela de consenso apresenta os scores Vina e DOCK6, SMILES, receptor e normalizações z-score/min-max, com botões **3D Vina**, **3D DOCK6** e **Remover** em cada linha. O score DOCK6 usado no consenso é `min(0, Grid_Score + repulsion_weight × Internal_energy_repulsive)`; a repulsão ausente vale zero. Para um composto ou scores constantes, as normalizações são zero. A remoção exige permissão de edição, fica registrada e afeta próximas entradas; não recalcula análises já concluídas. As tabelas principais de cada lote são usadas sem reintroduzir linhas removidas a partir das tabelas auxiliares.


Em **Preparar para docking**, conecte **Receptor PDB (retrieval ou redocking)** e **Compostos selecionados** separadamente. Para um PDB bruto do retrieval, informe os registros `[PDB, ligante de referência, resíduo, cadeia]` ou forneça `pdb_codes.csv`; as opções de preparo do receptor permanecem disponíveis. Para receptores do redocking, escolha o arquivo `.dockprep.pdbqt` no campo **Receptor que deseja utilizar**: os arquivos exclusivos do ligante ficam ocultos e as opções de preparo do receptor são desabilitadas. Os arquivos complementares e os centros acompanham o receptor; ele é reutilizado sem novo preparo.

Adicione uma ou mais fontes de compostos de **ChEMBL**, **PubChem**, **ZINC** ou arquivos próprios em **Compostos selecionados**. O grupo **Preparação e conformação do ligante** configura esses candidatos externos e permanece disponível. O ligante de referência do PDB serve para definir o sítio de docking, não integra a lista de candidatos. Escolha **Preparar saída para → Vina, DOCK6 ou Vina e DOCK6**. O bloco prepara cada candidato e exporta os formatos selecionados a partir da mesma conformação. Conecte a saída deste bloco tanto à entrada de receptor quanto à entrada de compostos de cada ferramenta escolhida. Para reunir várias fontes em uma execução, selecione **Mesclar arquivos (merge)**; o modo individual prepara cada combinação separadamente. Confira os receptores, compostos e centros produzidos.

Use um bloco separado para cada modo de receptor: não combine PDBs brutos com receptores preparados no mesmo bloco. O pH continua disponível para o preparo dos candidatos quando o receptor é reutilizado. Antes de preparar os compostos, o bloco verifica os centros do sítio (três coordenadas finitas) e os arquivos do receptor exigidos pela saída escolhida. Vina requer `.dockprep.pdbqt`; DOCK6 requer também `.dockprep.mol2` e `.noH.pdb`. O PDBQT permanece como arquivo de seleção do receptor nos dois casos. Arquivos ausentes interrompem o processo com o nome do arquivo necessário, sem refazer o preparo do receptor.

## Ajuda dos blocos e seleção de composto

O ícone **ⓘ** aparece na biblioteca e em cada bloco do canvas, inclusive para leitores. O popup explica a função do bloco, suas entradas e suas saídas. Similaridade recebe fingerprints; conectar diretamente ChEMBL ou PubChem apresenta um aviso padrão de incompatibilidade, sem traceback. A recusa não altera conexões nem o histórico de desfazer.

Na configuração do docking e no popup de arquivos, cada tabela de compostos oferece **Composto para docking (opcional)**. O seletor permite busca por identificador e usa todos quando vazio. A escolha é conservada ao reabrir a configuração; trocar o arquivo limpa o identificador. Com receptor e arquivos de compostos explícitos, Vina e DOCK6 não repetem o popup durante a execução. Origens ainda automáticas solicitam a escolha quando seus dados estiverem disponíveis.

## Atualizações pela branch master

O botão **Atualizações**, no cabeçalho, consulta a branch `master` de
[mpiress/BioMolExplorer](https://github.com/mpiress/BioMolExplorer) ao abrir a
interface e a cada hora. Quando há uma revisão nova, passa a mostrar
**Atualização disponível**. Clique para consultar o commit e escolher
**Baixar atualização** ou **Lembrar em 24 horas**. O adiamento é salvo por conta
no workspace; o botão continua permitindo uma consulta manual.

O download é um ZIP do commit mostrado no aviso e vai para o computador do
usuário, inclusive quando a interface está em um servidor web. A instalação
continua sendo uma etapa manual: encerre as análises, extraia o código e siga o
[guia de instalação](installation.md), depois reinicie o aplicativo. O download
não altera o código em execução, as contas ou os projetos. Falhas de conexão não
impedem o uso da interface; a consulta manual permite tentar novamente.

Checkouts Git usam o commit instalado como referência; pacotes wheel construídos
pelo projeto incluem essa referência. ZIPs baixados pelo botão preservam o commit
no nome da pasta `BioMolExplorer-<commit>`: mantenha esse nome. ZIPs sem referência
(como `BioMolExplorer-master`) registram o primeiro commit consultado e avisam
sobre mudanças posteriores. As preferências ficam em `updates.json` no diretório
de dados da aplicação.

Seleções de resultados guardam o lote e o caminho lógico do arquivo, sem o
identificador temporário do worker. Assim, uma nova execução do PDB ou do docking
mantém a seleção do mesmo resultado. Configurações anteriores também são aceitas;
resultados ausentes ou ambíguos continuam exigindo uma nova seleção.

## Visualizador molecular e interações 3D

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
