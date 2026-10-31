> Registro histórico. As instruções de ambiente/protocolos desta validação foram substituídas pela [revisão de Amber Score nativo de 2026-10-29](native_amber_score_2026-10-29/README.md). Não representam as opções atuais da aplicação.

# Validação do bloco de dinâmica — 2026-10-23

## Causas e correções

- A ausência de `selected_pairs` significava executar todas as combinações. Agora a seleção é explícita; adicionar entradas não seleciona novos pares automaticamente. Há seleção múltipla, exclusão e pré-visualização individual.
- O campo de compostos aceitava apenas poses MOL2. A resolução agora aceita etapas compatíveis, CSV com referências a poses, SDF e PDBQT, utilizando o conversor existente. Tabelas com o mesmo composto em receptores distintos não perdem poses por deduplicação. O índice do conformero é preservado.
- A identidade de uma molécula CSV dependia de dados globais da tabela. O cache por par agora considera o conteúdo da pose selecionada e seu índice, sem depender de caminhos de snapshots ou linhas não utilizadas. Alterar uma pose referenciada invalida seu resultado.
- Os checkpoints históricos não eram consultados pelo índice compartilhado. Resultados completos e preparações de receptor, ligante e complexo agora utilizam `CalculationCache`, hashes SHA-256 e o índice SQLite já adotado no docking. O controle global de reuso é respeitado.
- A preparação era repetida ao mudar parâmetros de simulação. Os três checkpoints de preparação não incluem temperatura/duração. Um teste nativo identificou e corrigiu a perda da permissão executável dos helpers locais usados para pular preparações já disponíveis.
- A análise antiga ajustava o complexo e misturava componentes. Foi substituída por ajuste ao backbone proteico, RMSD separado do ligante sem reajuste, RMSF C-alfa em torno da média alinhada, Rg ponderado por massa, ligações H apenas da interface e ocupâncias. Foram acrescentados SASA LCPO, distância COM ao sítio inicial e ocupância de contatos por resíduo.
- A leitura energética agora exclui as seções estatísticas de médias/RMS do sander, que não são quadros adicionais, e calcula a média sobre amostras efetivamente gravadas.
- Validação de topologias, ordenação ligante–receptor, massas, cargas e coordenadas finitas antecede os cálculos. Cargas preservadas devem somar aproximadamente um inteiro; isso não substitui revisão química de protonação.
- Falhas moleculares são isoladas por par, com fase registrada e logs AMBER preservados. Cofatores/metais/resíduos desconhecidos recebem orientação de diagnóstico; não são silenciosamente removidos.

## Resultados e rastreabilidade

`dynamics_results.json` indexa os pares; `analysis.json` consolida séries, máscaras, definições, resíduos e ocupâncias. `protocol.json` registra parâmetros, hashes de ferramentas/entradas, origem, campos de força, solvente e restrições. O painel responsivo apresenta resumo, métricas, interações persistentes, gráficos, energias e arquivos. A comparação descritiva agrupa pelo conteúdo do receptor, parâmetros e ferramentas; não gera ranking de afinidade.

Após sucesso, somente cópias locais dos helpers e o arquivo transitório de cargas são removidos. `BIOMOL_MD_DEBUG=1` os conserva. Topologias, coordenadas, configurações científicas, logs e trajetórias permanecem disponíveis. O cache deduplica objetos por conteúdo e utiliza cópia por reflink quando suportada, sem hardlinks mutáveis.

## Testes executados

Ambiente principal: Python 3.12 do BioMolExplorer. DOCK6 nativo: `/home/michel/progs/dock6`. AmberTools 26 foi instalado exclusivamente em `/tmp/biomol-md-validation`, a partir de `requirements-ambertools.yml`, sem modificar o ambiente principal.

- Regressão: **133 testes aprovados em 78,755 s**, cobrindo dinâmica, métricas, painel, localização, execução do pipeline, canvas/editor, preparação, handoffs, Vina, redocking e caches incrementais.
- Rodada final após ajustes de leitura energética/entradas: **31 testes aprovados em 14,902 s**, nos módulos `test_md_analysis`, `test_md_dashboard`, `test_molecular_dynamics` e `test_localization`. Há sobreposição com a regressão acima.
- Execução nativa completa do exemplo DOCK6 `amber_score_1lgu`: preparação AMBER com cargas preservadas, AMBER Score, minimização sander, trajetória NetCDF de 20 passos com stride 5, quatro quadros, análises e PNGs finitos. Teste adicional executa AM1-BCC nativo e verifica carga total da topologia.
- Execução nativa repetida a 301 K: os logs dos três helpers são vazios na segunda execução, confirmando que receptor, ligante e complexo foram restaurados sem repetir sua parametrização.
- Casos sintéticos controlados: transformação rígida resulta em RMSD/RMSF aproximadamente zero; deslocamento do ligante de 10 Å permanece detectável; contato presente em um de dois quadros tem ocupância 50%; Rg conhecido, tempo/stride e rejeição de NaN/quadros insuficientes.
- Lotes com múltiplos pares e falha parcial: testes com motor controlado verificam preservação dos sucessos, retomada apenas da falha, inclusão/remoção, alteração de parâmetros, cache corrompido, reuso sem seed e desativação do reuso.
- Painel construído com largura de 360 px e comparação agrupada testados sem navegador. Não foi realizada inspeção visual interativa em dispositivo físico.
- Documentação HTML regenerada; 32 páginas com arquivos locais, âncoras e IDs válidos. `git diff --check` sem erros.

## Limites científicos e operacionais

O AMBER Score nativo continua distinto da trajetória sander adicional. Os helpers DOCK6 testados adotam **ff94/GAFF**; não houve troca silenciosa por ff14SB/GAFF2. A trajetória é refinamento curto em GB implícito, sem uma etapa de equilibração validada ou controle de pressão. Não é dinâmica de produção, MM/GBSA nem energia livre de ligação.

A máscara sander padrão `:1` foi preservada: somente o primeiro resíduo se move. O usuário deve revisar ligantes multirresíduo e não interpretar uma proteína congelada como evidência de baixa flexibilidade. Um composto sem pose tridimensional associada exige docking ou pose fornecida; não é colocado arbitrariamente no receptor.

As validações nativas são testes de integração curtos, não validação científica de todos os sistemas, metais ou cofatores. Interações persistentes são contatos geométricos e ligações H; não foram classificados automaticamente contatos hidrofóbicos ou pares iônicos. A suíte completa do repositório não foi reexecutada nesta rodada; os testes acima são o escopo efetivamente validado.
