# Redocking incremental e exportação do pipeline — 2026-10-19

## Comportamento implementado

- O reaproveitamento ocorre por receptor e por par, além do cache existente de etapas completas.
- Um ligante novo reutiliza o receptor compatível; apenas o ligante é preparado e seu redocking é executado. Um receptor novo recebe preparação completa.
- A nova tabela consolida os resultados anteriores compatíveis da mesma origem e os pares novos, mesmo quando apenas as novas entradas são selecionadas. Lotes individuais distribuem os resultados anteriores para evitar duplicação.
- Os resultados reutilizados têm arquivos próprios na nova execução: estrutura, receptor preparado, ligante preparado, centro e pose. O cache não depende de poses soltas ou de um RMSD sem saída correspondente.
- SHA-256 verifica os arquivos. A assinatura científica inclui estrutura PDB, cadeia, preparação, cofatores, pH, cargas, templates, parâmetros de redocking e implementação. Receptores com arquivos preparados diferentes não são misturados sob o mesmo nome.
- Somente execuções do pipeline com estado `succeeded` fornecem dados. Uma etapa bem-sucedida de uma execução interrompida também é recusada pelo cache de etapas completas. Execuções parciais permanecem apenas como histórico diagnóstico.
- Resultados nativos antigos sem índice incremental podem ser adotados a partir da configuração, metadados, arquivos preparados, poses e manifesto original verificado. Resultados fornecidos/importados não são adotados implicitamente.
- Entradas apagadas não são recuperadas automaticamente do histórico. Mudanças na origem de dados não misturam coleções. A escolha de executar novamente ignora o cache.
- Os arquivos temporários que transportam o cache são removidos após uso. Arquivos compartilhados do receptor são verificados e copiados uma vez por conjunto, evitando trabalho repetido entre ligantes.
- A câmera dos controles do canvas exporta o diagrama completo como `pipeline.png`, preservando posições, nomes, categorias, portas e conexões, independentemente do zoom ou da área visível. A renderização não depende de capturar widgets fora da tela.
- PNG RGB, escala de até 2×, metadados de 300 dpi e limite de 32 milhões de pixels. Leitores também podem exportar. Trocas de projeto/sessão durante a geração impedem o download do projeto anterior.

## Validação

Ambiente: Python 3.12 do ambiente Conda BioMolExplorer, Flet 1.0.3, Pillow e dependências científicas instaladas.

Resultado final: **103 testes aprovados em 46,277 segundos**. `git diff --check` sem problemas de espaçamento. A cobertura inclui retornar a uma seleção antiga sem perder resultados acrescentados depois e descartar um redocking aparentemente completo quando a execução do pipeline foi interrompida.

```bash
PATH=/home/michel/progs/miniconda3/envs/BioMolExplorer/bin:$PATH \
PYTHONPATH=src:tests MPLCONFIGDIR=/tmp/biomol-mpl \
python -m unittest \
  test_redocking_incremental test_pipeline_png test_stage_cache \
  test_restart_reuse test_redocking_execution_review \
  test_redocking_configuration test_redocking_results \
  test_redocking_partial_failures test_canvas_drag_regression \
  test_pipeline_editor_regressions test_run_reuse_ui test_localization -q
```

A suíte cobre receptor compartilhado, novo ligante, novo receptor, acumulação sem repetir seleções, origens removidas, mudança de coleção, separação de lotes, recusa de execuções incompletas, pose corrompida, mudanças de estrutura/preparação/parâmetros, execução forçada e adoção de resultados antigos. A integração usa o armazenamento e o agendador reais, executa o adaptador de operação e consulta o serviço da tabela final.

A exportação foi verificada por decodificação do PNG, pixels dos blocos fora da tela, conexão entre blocos, dimensões, metadados de resolução, limite de memória, controle de câmera para leitores e proteção durante mudança de projeto. Um exemplo com dois blocos e conexão foi inspecionado visualmente. Os Markdown e HTML afetados foram atualizados em português e inglês.

Limites: Chimera, Open Babel e Vina foram substituídos por respostas determinísticas nos testes. Esta revisão valida orquestração, reaproveitamento e apresentação; não representa uma nova validação científica dos motores externos nem uma sessão manual do cliente Flet com uma execução molecular real.
