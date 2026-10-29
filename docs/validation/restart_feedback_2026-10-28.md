# Feedback durante o reinício — 28/10/2026

O fluxo de reinício não acionava o splash durante salvamento, hashing e validação
dos resultados. Agora utiliza o splash padrão imediatamente, com a mensagem
“Preparando a execução…” e a descrição da validação, no idioma da sessão.
Após confirmar uma alteração, o splash permanece visível até surgir o menu de
execução com ou sem reúso. Sem esse menu, cobre também a submissão até a transição
para o progresso da execução. Erros e cancelamentos sempre liberam o indicador.
O bloqueio existente de submissões repetidas foi preservado.

Validação: 42 testes aprovados, incluindo os dois idiomas, presença do splash antes
dos acessos, continuidade durante a validação, transição para o menu, cancelamento,
falhas e cliques repetidos. No Chrome real em 360 px, um projeto temporário validou
confirmação → splash → menu de reúso, sem erros JavaScript. O backend foi atrasado
de forma controlada nesse teste visual. Nenhuma biblioteca foi instalada.

[Logs e capturas](restart_feedback_2026-10-28/).
