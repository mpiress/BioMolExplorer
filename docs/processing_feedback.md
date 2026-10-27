# Feedback visual de processamento

As operações da interface utilizam o splash padrão do BioMolExplorer: logo,
anel animado e mensagens que identificam o trabalho em andamento. Ele aparece
acima da janela de edição e permanece até concluir o processamento e a
atualização da interface. Falhas e cancelamentos liberam o splash antes de
apresentar o diagnóstico ou permitir nova tentativa.

## Ações abrangidas

| Ação | Feedback |
| --- | --- |
| Identificar ligantes PDB | Analisando ligantes PDB |
| Salvar ligantes, inclusive escolhidos manualmente | Salvando ligantes PDB, antes da validação e gravação |
| Configurar um bloco | Carregando configurações / Aplicando configuração |
| Selecionar entradas do próximo bloco | Preparando seleção de arquivos |
| Salvar pipeline | Salvando pipeline |
| Criar, abrir ou atualizar projeto | Mensagem correspondente ao projeto |
| Importar, exportar ou restaurar projeto | Mensagem correspondente à operação |
| Ler arquivos, resultados, histórico ou logs | Mensagem correspondente ao conteúdo |
| Excluir arquivos, compostos ou resultados | Mensagem correspondente à exclusão |
| Visualizar moléculas, estruturas, ADMET ou footprint | Mensagem correspondente à análise ou visualização |
| Calcular métricas, comunidades, disposição do grafo ou MCS | Mensagem correspondente ao cálculo |
| Iniciar, continuar ou cancelar pipeline | Mensagem correspondente à execução |

No salvamento de ligantes, os botões ficam desabilitados durante a operação;
cliques repetidos não duplicam a gravação. O splash permanece durante a
atualização dos pares de redocking e da tabela. Em caso de erro, os campos
preenchidos são mantidos e o usuário pode tentar novamente.

O acompanhamento das execuções usa o mesmo emblema animado, incluindo o resumo
minimizado. Os controles de cancelar, minimizar e configurar continuam
disponíveis. A animação para quando a execução termina ou aguarda uma decisão.
O cálculo interativo de MCS também mantém a opção de interrupção no splash.

Consultas automáticas de acompanhamento e salvamentos automáticos do rascunho
não abrem janelas sobre o trabalho do usuário. Durante seletores nativos de
arquivo, o splash é suspenso enquanto o usuário escolhe o destino ou a entrada;
ele volta se houver outro processamento depois da escolha.

## Padrão para novas ações

- Operações de backend passam por `WorkspaceUI.call()`, que mostra o splash
  antes de iniciar o trabalho em outra thread e escolhe a mensagem pelo serviço.
- `guard()` agrupa chamadas de uma mesma ação. `@with_processing()` define
  mensagens específicas e mantém um único splash durante suas várias etapas.
- Para botões que precisam publicar o feedback antes de qualquer validação,
  use `async with processing(ui, título, descrição)` no início da operação.
- Forneça `on_cancel` para cálculos que suportam interrupção.
- Use `save_file()` / `pick_files()` do módulo de processamento para os
  seletores nativos e `@background` para tarefas automáticas.

Requisições simultâneas compartilham o componente. Uma ação que termina não
remove o indicador de outra que ainda está trabalhando. As mensagens e os
rótulos de acessibilidade possuem traduções em português e inglês.

A consulta periódica de atualizações também usa `@background`: mesmo com uma janela de configuração, remoção ou cancelamento aberta, essa tarefa não cria splash. A consulta manual de atualizações continua com feedback quando solicitada.
