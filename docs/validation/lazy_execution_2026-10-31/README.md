# Espera lazy e paginação do redocking — 31/10/2026

## Comportamento

O monitor de execução abria `configure_pending` automaticamente ao detectar `awaiting_input`. Agora ele apenas prepara os dados e mantém o estado de espera. Os botões existentes em Runs e Follow run continuam sendo as entradas explícitas para configurar o próximo bloco. Pausa manual e permissões de visualizador continuam respeitadas.

A preparação usa uma tarefa assíncrona com trabalho de leitura em thread, sem criar controles, abrir diálogos, mostrar splash ou avançar o pipeline. Mantém um único cache de sessão, identificado pela execução, bloco, configuração e saídas das dependências. Reúne catálogo/identidades, proveniência, saídas públicas, cabeçalhos CSV, metadados de pares, revisões PDB e referências de compostos. No clique, confere assets e os atributos de tamanho/mtime/ctime dos arquivos; mudanças provocam nova preparação. A validação científica do backend não foi substituída por esse cache de apresentação.

Falhas de preparação antecipada são registradas sem interromper o acompanhamento ou navegar para outro workspace. Ao clicar, uma preparação que falhou pode ser tentada novamente. Fechar o projeto cancela a tarefa de interface; a leitura em thread não publica controles nem altera arquivos científicos.

## Redocking

O contador usava o total filtrado como quantidade mostrada, embora a tabela exibisse no máximo 25 linhas. A navegação estava depois da tabela larga, cujo conteúdo podia ocupar a largura do componente pai. Agora:

- O contador distingue intervalo visível, total filtrado, total geral e página atual/total.
- A navegação está antes da tabela, limitada à largura disponível.
- A tabela tem rolagem horizontal interna e os filtros usam a grade responsiva compartilhada.
- O deslocamento é limitado à última página válida; filtros e ordenação continuam locais, sem novas consultas de dados.

## Evidências

`real-results-audit.json`: leitura somente do banco e arquivos do projeto atual. Foram encontradas 50 simulações de AChE PDB, BChE PDB e BACE1 PDB. As duas páginas de 25 foram percorridas; o botão Next está habilitado na primeira e desabilitado na última.

`browser-report.json`: seis capturas e cliques reais em Next nas larguras de 1440, 1024 e 360 pixels. O botão ficou dentro da área visível nas três larguras; o estado avançou de 0 para 25. O console teve apenas o registro normal Flutter app loaded, sem exceções JavaScript. A prévia usa cópia dos resultados em memória e não altera o experimento.

`tests.log`: suíte final de seleção/espera, preparação em thread, cache de sessão, paginação, filtros, ordenação, i18n, pausa/retomada, cancelamento corrente, curadoria e navegação 3D.

Uma execução exploratória mais ampla (`broader-regression.log`) passou 99 de 101 testes. Os dois restantes em `test_execution_cancellation` esperam que fechar/cancelar uma configuração cancele a execução inteira. Isso diverge do fluxo atual preservado, no qual cancelar uma edição mantém a execução aguardando, e Cancel Run é uma ação distinta. Não foram alterados os botões de cancelamento para satisfazer essas expectativas antigas. Esses dois testes não são declarados aprovados.

A primeira tentativa de navegador precisou ajustar o adaptador de prévia (faltava docking_view_url); outra precisou reiniciar apenas a página da prévia entre tamanhos. Ambas foram corrigidas no harness temporário antes da captura final. A tentativa de testes dentro do sandbox ficou bloqueada no suporte a threads e foi encerrada; a validação final foi executada com suporte local a threads.

## Limites e preservação

Nenhuma execução real foi retomada/cancelada, nenhum resultado foi apagado e nenhum parâmetro científico foi modificado. Motores de docking não foram reexecutados: esta revisão altera a apresentação e a preparação de dados para configuração. Não foram instaladas dependências. O cache é restrito à sessão e reconstruído após reiniciar a interface.
