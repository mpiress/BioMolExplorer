# Ordem de execução, pausa, resultados e screenshot — 2026-10-30

## Comportamento implementado

- A primeira execução pela interface solicita confirmação da ordem. As setas movem o bloco e todos os seus dependentes transitivos; movimentos que atravessariam uma dependência são bloqueados. O agendador respeita essa prioridade dentro do DAG, inclusive ao executar somente a cadeia até o bloco selecionado.
- **Pipeline e execução → Ordem de execução** permite rever a última ordem. A preferência é gravada atomicamente no SQLite, separada da configuração científica. Novos blocos são incorporados respeitando dependências. Projetos que já possuem execuções aproveitam a ordem anterior sem uma nova solicitação automática. Alterar somente a ordem não invalida resultados.
- **Runs → Pausar após operação atual** conclui a operação científica em andamento e pausa antes da próxima operação/lote. **Continuar execução** retoma o mesmo run e preserva etapas e lotes já concluídos, mesmo após reiniciar o serviço. O estado usa `awaiting_input` com `pause_reason=user`, mantendo compatibilidade com os controles existentes. A pausa manual não abre automaticamente a curadoria.
- A pausa permite minimizar o acompanhamento, consultar resultados e configurar entradas antes da continuação. **Cancelar execução** continua disponível para encerramento definitivo. A pausa não suspende Vina/DOCK6 no meio de um cálculo nem promete um checkpoint interno desses executáveis.
- Uma solicitação tardia de pausa não pode reabrir uma execução terminada. Leitores não podem pausar, continuar ou salvar uma ordem.

## Causas dos erros de interface

O log da interface registrava `NameError: unsupported is not defined` ao criar a tabela de redocking. A revisão também identificou chamadas a `view_button` e `view_url` inexistentes e ausência do import de `save_file` no download ZIP. Esses caminhos foram corrigidos, preservando navegação 3D direta no navegador, autorização, ordenação e diagnósticos de falhas.

O screenshot falhava com `KeyError: Dinâmica Molecular`: o exportador tinha uma cópia incompleta da paleta do canvas. Ambos agora compartilham a mesma paleta. A exportação continua independente do viewport e limitada a 32 milhões de pixels.

## Auditoria do projeto atual

O banco original foi aberto exclusivamente com SQLite `mode=ro`; não houve atualização de credenciais, migração do banco real, reexecução científica ou escrita em arquivos do projeto nesta auditoria.

A execução `384ffe70b0484994ab021a823aa78efb` contém 161 entradas de redocking de AChE, BChE e BACE1 e 1.542 artefatos presentes. **As 161 entradas estão registradas como `preparation_failed`**: a tabela corrigida exibe os diagnósticos; não existem RMSDs válidos nessa execução para essas entradas. Em 82 registros, o erro persistido é `O preparo do ligante alterou o grafo molecular validado`; outros registros incluem falhas de kekulização. Corrigir a apresentação não converte esses registros em simulações científicas aprovadas.

A tabela foi construída com esses dados reais sem erro, exibindo a primeira página de 25 entradas. Também foi exportado um PNG de 323.781 bytes com os 17 blocos atuais. [Auditoria em JSON](project-audit.json).

## Testes executados

Ambiente existente: `/home/michel/progs/miniconda3/envs/BioMolExplorer/bin/python`; `PYTHONPATH=src:tests`, `MPLCONFIGDIR=/tmp/biomol-mpl`. Nenhum pacote instalado.

- **79 testes passaram** na suíte focada: `test_execution_order_pause`, `test_redocking_results`, `test_redocking_viewer_navigation`, `test_pipeline_png`, `test_run_progress.RunFeedbackTests`, `test_runs_pipeline_controls`, `test_runtime_messages`, `test_result_sorting`, `test_published_results` e `test_md_downstream_persistence`.
- **52 testes passaram** na suíte adicional do agendador, mensagens, ordenação/pausa e exclusões de falhas de redocking: `test_pipeline_execution`, `test_runtime_messages`, `test_execution_order_pause` e `test_redocking_partial_failures.FailedResultDeletionTests`. As suítes se sobrepõem; os totais não devem ser somados.
- O teste adicional de reordenação após uma execução confirmou que as saídas mantêm os mesmos caminhos e que nenhuma nova chamada ao supervisor ocorre quando o reúso está habilitado.
- Pausa/retomada foram verificadas com supervisor controlado, inclusive reinicialização do serviço, pausa repetida entre lotes com reúso global desabilitado, cancelamento pausado e corrida com conclusão. Esses testes não executam novas simulações científicas.
- Revisão renderizada em Chrome local nas larguras **1440, 1024 e 360 px**, em projeto temporário: diálogo legível, rolagem interna e ações acessíveis; nenhuma exceção de navegador. O único registro de console foi a mensagem informativa `Flutter app loaded`. [Desktop](order-desktop.png), [celular](order-mobile.png).
- `git diff --check` sem erros.

## Limites da validação

A verificação mais ampla com `test_localization` encontrou 11 falhas de subcasos por mensagens sem tradução nos formulários existentes e na validação de exportação ADMET. Uma execução exploratória de `test_run_progress.WorkerRetrievalTests` também encontrou duas expectativas incompatíveis com o cache isolado quando o reúso está desabilitado e com os nomes publicados de arquivos ChEMBL. Esses resultados não são apresentados como aprovação global. As mensagens novas de ordem/pausa constam dos catálogos pt/en, e os testes de mensagens de sessão passaram.

Não foram refeitos cálculos de redocking, docking ou Amber Score, nem alterados parâmetros químicos para contornar as falhas científicas persistidas.
