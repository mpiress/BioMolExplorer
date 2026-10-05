# Projetos, arquivos e versões

[Documentação](README.md) · Português · [English](en/projects.md)

## Escolher a pasta do projeto

Em **Novo projeto**, informe nome, descrição, tags, cor e **Pasta do projeto**.
Escolha uma pasta nova ou vazia pelo botão dentro do campo, que não permite digitação manual. A interface desktop oferece um seletor do sistema; na versão web, o botão abre a navegação de pastas e permite criar uma pasta;
no navegador, o caminho pertence ao computador que executa o backend, e não ao
computador do navegador. A aplicação precisa poder escrever nesse caminho.
Pastas sobrepostas a outros projetos ou ao armazenamento de contas são rejeitadas.

Essa pasta guarda `project.json` (configuração portátil), `assets/` (entradas
originais), `runs/` (experimentos e resultados) e `.history/` (versões).
`project.json` acompanha as edições; referências internas são expressas com
`@project`, permitindo mudar a pasta durante a importação. Use **Exportar** para
migrar o projeto completo: o JSON sozinho não contém os dados científicos.
Contas, sessões, permissões correntes e índices continuam no SQLite do workspace.
Os logs continuam centralizados na pasta de logs configurada pela aplicação.

## Migrar para outro computador

1. Conclua ou cancele a execução ativa.
2. No card do projeto, use **Exportar** e salve o pacote `.bme.zip`.
3. No novo computador, instale o ambiente científico e entre na aplicação.
4. Use **Importar projeto**, selecione o pacote e informe uma pasta nova ou vazia.
5. Abra o projeto importado e execute o pipeline normalmente.

O pacote inclui configurações, arquivos enviados, resultados, metadados das
execuções e histórico. A importação verifica hashes e caminhos antes de criar o
projeto, ajusta as referências e evita colisões de identificadores. Etapas
concluídas com dados e configurações compatíveis são reaproveitadas; arquivos
faltantes, alterações de entrada ou parâmetros exigem execução das etapas afetadas.
Uma versão diferente das ferramentas científicas deve ser avaliada pelo pesquisador;
escolha **Não, executar novamente** na pergunta de reaproveitamento quando desejar recalcular os experimentos.

Os caminhos de instalação são configurações do computador. Novas execuções DOCK6
usam a instalação configurada no backend; reutilizar um resultado concluído não
exige o caminho da instalação anterior.

O usuário que importa torna-se proprietário. Contas, senhas e sessões não são
exportadas. Colaboradores presentes no workspace de destino recebem acesso apenas
após novo aceite; o histórico conserva os nomes e e-mails dos autores anteriores.
Restaurações de versões importadas também preservam essa exigência de aceite.

No navegador, envio e download pela interface são limitados a 200 MB. Um pacote
de exportação maior fica disponível em `.exports/` para cópia direta; use a
aplicação desktop para importá-lo. A importação aceita até 4 GB descompactados e
100 mil arquivos. Exportações exigem permissão de leitura; importar cria um
projeto privado do usuário autenticado.

## Compartilhar e trabalhar em equipe

Dentro do projeto, o proprietário usa **Compartilhar** para convidar contas já
cadastradas como **Editor** ou **Leitor**. O destinatário aceita no workspace.
Todos consultam o mesmo projeto persistido; não são criadas cópias do pipeline.

Aplicar configurações salva imediatamente. Inclusão, remoção, conexões e movimentos
dos blocos são salvos após uma breve pausa de edição. O botão **Salvar** continua
disponível. A interface verifica atualizações de colaboradores a cada dois segundos
enquanto o projeto está aberto. Formulários abertos e rascunhos locais não são
substituídos por essa atualização.

Edições independentes são reunidas ao salvar. Se duas pessoas alterarem o mesmo
campo, a interface informa o conflito e conserva o rascunho local. Confira a versão
compartilhada antes de reaplicar esse campo. Revogações e mudanças de papel são
verificadas também no backend; leitores não podem alterar nem executar o projeto.

## Consultar o histórico e restaurar uma versão

O botão **Histórico** no card ou dentro do projeto abre uma tabela da alteração
mais recente para a mais antiga, com data, hora, usuário e descrição. O registro
inclui configurações, metadados, convites, uploads, curadoria e conclusão de execuções.
Datas são exibidas no fuso do computador que executa a interface Python.

O proprietário pode usar **Restaurar antes**. Após a confirmação, o projeto volta
ao estado anterior à alteração escolhida: configurações, arquivos, resultados e
permissões são restaurados. A revisão interna aumenta para invalidar rascunhos
antigos; o proprietário e o histórico de auditoria são mantidos. A restauração
também aparece no histórico. Versões posteriores permanecem preservadas e podem
ser recuperadas através das versões anteriores às restaurações subsequentes.
Uma execução ativa impede a restauração.

As versões usam cópias por conteúdo em `.history/blobs/`: arquivos iguais são
armazenados uma vez; arquivos modificados geram novos conteúdos. Preserve essa
pasta ao fazer backup. Não há retenção automática nem expurgo das versões nesta
edição. A exclusão do card continua sendo lógica; projetos excluídos ficam fora
do acesso normal e não são restaurados pelo botão de histórico.

## Escolher várias entradas de um bloco

Abra a configuração do bloco e escolha a origem em **Entradas**. Para uma etapa
já concluída, **Resultado usado nesta entrada** lista arquivos da última execução
bem-sucedida, com seu caminho relativo para distinguir nomes iguais.
Na recuperação CHEMBL220, por exemplo, selecione `CHEMBL220_FULL.csv`,
`CHEMBL220_MOLS.csv`, `CHEMBL220_SIMS.csv` ou o `compounds.csv` integrado.

**Adicionar entrada** inclui arquivos da mesma etapa ou de outras origens. Também é possível misturar saídas e arquivos enviados; uma nova conexão acrescenta uma origem. Na execução, o popup permite **Processar individualmente** ou **Mesclar arquivos (merge)**. Individual mantém resultados por arquivo; merge forma uma união normalizada de CSVs. Códigos repetidos da mesma molécula são deduplicados; códigos conflitantes são registrados e excluídos pela limpeza molecular. As demais colunas são preservadas, e a primeira linha escolhida prevalece nos duplicados equivalentes. Originais são conservados. Arquivos estruturais com nomes iguais e conteúdos diferentes precisam ser renomeados antes de combinar.

Blocos mostram suas entradas e execuções registram os arquivos resolvidos. Individual forma combinações entre portas distintas; DOCK6 mantém apenas combinações correspondentes e o consenso pareia poses/scores pelo identificador. O popup conserva escolhas ao abrir a configuração, e seleções confirmadas ficam salvas no bloco. Projetos continuam limitados a 100 blocos. Ao executar com resultados anteriores íntegros, escolha se deseja reaproveitar; etapas compatíveis não repetem a seleção. Consulte o [manual](user_manual.md).

## Entradas próprias e resultados prontos

Em cada entrada, **Enviar meus arquivos** disponibiliza dados para o cálculo do
bloco. O envio valida colunas, SMILES, fingerprints, scores ou coordenadas conforme
o tipo escolhido. Se houver erro, a mensagem informa o formato esperado e o arquivo
não é incorporado à biblioteca. O original válido é preservado; a normalização
ocorre na cópia usada pela execução.

Se o cálculo já foi realizado fora da aplicação, escolha **Como usar este bloco? →
Usar meus resultados prontos · etapa concluída**, selecione ou envie os resultados
e aplique a configuração. O bloco passa a indicar **Resultados fornecidos**.
A execução apenas organiza os arquivos em seu diretório de resultados; o cálculo
científico desse bloco não é chamado. As etapas seguintes consomem esses resultados.
Entradas e dependências antigas do bloco ficam inativas nesse modo.

| Resultado | Formato esperado |
| --- | --- |
| Compostos | CSV UTF-8 com `molecule_chembl_id,canonical_smiles`; aliases `name,smiles` são aceitos e códigos ausentes são gerados |
| ADMET pronto | CSV com código, SMILES e valores finitos de `TPSA,WLOGP`; outras propriedades são preservadas |
| Fingerprints | CSV com código e `fingerprint`, lista de bits 0/1 de tamanho consistente |
| Similaridade | CSV `source,target,value`, com score finito entre 0 e 1 |
| PDB | Registros ATOM/HETATM com coordenadas válidas; metadados opcionais em `pdb_codes.csv` |
| Receptores preparados | PDBQT, `pdb_codes.csv` e `centers.csv` com as três coordenadas por complexo |
| Vina / DOCK6 | PDBQT com `REMARK VINA RESULT` ou `*_scored.mol2` com seções de molécula, átomos e Grid_Score |

ADMET pronto gera um EGG interativo usando os descritores fornecidos. Nos resultados
do bloco, expanda a etapa e selecione o arquivo; passe o mouse para identificar o composto e clique
para ver a estrutura 2D e suas propriedades. O gráfico não recalcula os descritores.

O bloco de grafos recebe relações prontas dos blocos de similaridade ou arquivos externos pelo formulário; fingerprints não são uma entrada direta. Os SMILES vêm da linhagem ou de tabelas externas correspondentes. Individual mantém análises independentes; merge combina relações selecionadas em uma análise. Selecione um MCC específico para alimentar a etapa seguinte. Modelos interativos, imagens de fragmentos e apresentações de graus acompanham a exportação. O fragmento é apresentado como SMILES extraído da referência.
