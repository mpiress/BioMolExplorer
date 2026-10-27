# Organizar as estruturas selecionadas para redocking

A seleção de arquivos PDB apresenta uma tabela com arquivo, alvo/origem,
ligantes selecionados, revisão e ações. O identificador do bloco de origem
distingue estruturas com o mesmo código provenientes de recuperações diferentes.

1. Abra **Adicionar PDBs à seleção** para incluir um arquivo, todos os arquivos
   da origem ou todas as origens. Após incluir, esse painel se recolhe para
   liberar espaço para a tabela.
2. Use **Alvo dos arquivos selecionados** para trabalhar com um alvo de cada vez.
   O seletor informa quantas estruturas foram revisadas em cada grupo.
3. Pesquise pelo PDB, alvo ou ligante e filtre por **A revisar**, **Revisados**
   ou **Indisponíveis**. Os filtros preservam os demais arquivos selecionados.
4. Em cada linha, use as ações para editar os ligantes, visualizar a estrutura
   em 3D, abrir sua página no RCSB PDB ou remover o arquivo da seleção.
5. Ao salvar os ligantes, a linha fica verde e recebe o indicador **Revisado**.
   Fechar o editor sem salvar mantém o estado anterior.

**Remover desta visualização** remove da seleção todos os itens que correspondem
ao alvo, estado e busca atuais, inclusive em outras páginas. **Limpar toda a
seleção** atua sobre todas as origens. Essas ações mantêm os arquivos de origem.

A tabela renderiza até 25 linhas por página e cria as ações das estruturas sob
demanda. Ela ocupa a largura disponível do painel. Em telas estreitas, agrupa
PDB e alvo na primeira coluna, reúne ligantes e revisão na segunda e distribui
as ações em duas linhas. Os campos se empilham com espaçamento entre eles.
Somente em larguras muito reduzidas é necessária a rolagem horizontal.
Os filtros e a seleção são preservados ao selecionar depois e reabrir a janela.

A revisão é registrada por estrutura junto à edição auditada dos ligantes.
Continua visível ao reabrir a seleção ou reutilizar seus resultados. Alterações
na estrutura ou nos ligantes tornam a revisão anterior pendente. A cor verde
indica uma revisão salva dos ligantes; a cadeia, cofatores e demais opções de
preparação continuam sendo configurados e validados na etapa de redocking.
