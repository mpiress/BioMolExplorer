# Redocking: acesso às entradas e limites das tabelas

## Diagnóstico

A montagem de `RedockingFilePicker` dependia de posições incorretas nos controles do seletor base: incluía a mesma coluna de arquivos selecionados duas vezes e omitia a navegação. Só os primeiros 25 itens eram visíveis. O cabeçalho também apontava para o controle de feedback errado.

Na execução atual, auditada com SQLite em modo somente leitura, existem 106 PDBs disponíveis nas três origens conectadas: AChE (1), BChE (5) e BACE1 (100). A configuração continha 104 seleções, das quais 55 não estão nos resultados atuais dessas origens. Essas referências não foram substituídas, descartadas ou consideradas válidas automaticamente.

## Correções

- Cabeçalho e paginação restaurados, com uma única instância da tabela.
- Pesquisa, filtros, seleção em lote, revisão de ligantes e navegação entre páginas preservados.
- Aviso explícito para referências indisponíveis: filtrar por Indisponíveis, revisar/remover as referências e adicionar as estruturas atuais desejadas. Alterações só são efetivadas pelo fluxo normal de confirmação.
- Agrupamento usa o nome da origem quando o metadado de alvo é apenas o nome genérico Estruturas/Structures; evita duplicar grupos de uma mesma origem entre arquivos atuais e referências antigas.
- O diálogo inicial usa o mesmo dimensionamento responsivo empregado no redimensionamento da janela.
- Tabela de pares ocupa a largura delimitada da seção; rolagem horizontal interna permanece disponível em telas estreitas. Ações de configuração permanecem no rodapé do diálogo.
- Traduções das mensagens do seletor e seus contadores em português/inglês.

## Validação

O log `tests.log` contém a execução final dos testes de seleção, curadoria, handoff, configuração e layout. As regressões novas verificam acesso às 61 entradas em três páginas, ausência de controles duplicados, remoção explícita de referências indisponíveis preservando as entradas atuais e agrupamento/tradução dos contadores.

O navegador local capturou oito telas em larguras de 1440, 1024 e 360 pixels. Os botões de paginação e do rodapé ficaram acessíveis; em telas estreitas, a tabela mantém rolagem interna. O único evento de console registrado foi `Flutter app loaded`, sem exceções JavaScript. Capturas precedem apenas o ajuste final dos rótulos de agrupamento/contagem.

Na prévia em memória do pipeline real, a revisão explícita das referências indisponíveis e a seleção das três origens permitiram construir 175 pares a partir dos 106 PDBs. Isso valida o handoff dos metadados, não a viabilidade científica da execução de todos os pares.

Nenhuma configuração do experimento, resultado, arquivo científico ou estado de execução foi alterado. Não foram instaladas dependências nem reexecutados motores de docking nesta revisão.

## Continuação

Na interface atualizada, reabra a configuração da execução aguardando entrada. Use o filtro Indisponíveis para revisar as seleções antigas; se desejar removê-las, use Remover desta visualização enquanto esse filtro estiver ativo. Em Adicionar PDBs à seleção, escolha as estruturas das origens atuais. Confirme a curadoria/configuração para continuar a mesma execução. Não é necessário limpar resultados ou criar uma nova execução.
