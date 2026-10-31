# Resultados preservados e saídas de redocking — 2026-10-29

## Caso real: 4M0E / CHEMBL1179697

A execução recente `04047e5c` foi cancelada após iniciar novamente a recuperação de compostos. A última execução concluída `15577f2d` permanece no histórico, com **15 etapas e 6.409 arquivos distintos**. Todos os arquivos catalogados existem e seus conteúdos correspondem aos manifestos SHA-256 persistidos. As configurações solicitadas dessas 15 etapas continuam equivalentes às atuais.

Portanto, neste incidente não houve exclusão física desses resultados. Houve recusa do reúso pelo cache de implementação, seguida da apresentação da execução cancelada/vazia como foco de Runs. Os resultados anteriores ficavam escondidos no histórico, dando a impressão de perda, enquanto a execução tentava iniciar os cálculos novamente.

A validação final usou o caminho real de `_reuse_stage` com materialização de entradas proibida durante a auditoria: **15 etapas reaproveitáveis, zero jobs científicos iniciados e zero arquivos do projeto regravados**. Somente a MD nova aparece como trabalho pendente. A seleção persistida foi resolvida como **4M0E_A.noH.pdb + CHEMBL1179697 · 4M0E_A**, sem exigir nova curadoria dos blocos anteriores. [Auditoria](project-audit.json).

O banco original foi aberto exclusivamente com SQLite `mode=ro`. Não foi necessário recriar, copiar, alterar status de execuções ou inventar manifestos: os originais estavam íntegros. A correção restaura sua utilização e acesso a partir do histórico durável.

## Correções

### Experimentos concluídos versus cache de cálculos novos

`published_results.py` valida resultados concluídos independentemente da assinatura global da versão atual do aplicativo. Exige requisição científica equivalente, escolhas explícitas de protocolo compatíveis, manifesto completo das saídas, integridade dos arquivos, entradas originais comprovadas e equivalência das dependências atuais. Resultados removidos, incompletos, obsoletos, com falhas ou manifestos inválidos não são certificados como uma etapa inteira concluída.

Essa retenção preserva a **versão do produtor original**, o run de origem, a assinatura da requisição e a política de verificação. Não rotula um experimento antigo como se tivesse sido executado pelo algoritmo atual. O cache de cálculos novos continua usando sua verificação de implementação/runtime. A opção de executar sem reúso é respeitada, produzindo uma nova execução sem remover os arquivos anteriores.

O mesmo mecanismo opera por lote para reaproveitar unidades válidas quando chegam novas entradas. Uma escolha automática só pode herdar seletores previamente confirmados se o conjunto completo da origem continuar igual; entradas novas não desaparecem por restauração da seleção antiga. Configurações alteradas ou entradas/saídas adulteradas rejeitam o reúso. Uma etapa já concluída também pode ser compartilhada com outra equivalente dentro da execução ativa, sem aceitar etapas pendentes ou em processamento.

As leituras de histórico são compartilhadas dentro da execução. Hashes são cacheados com limite de memória e assinatura de stat incluindo inode, tamanho, mtime e ctime; o escopo do projeto é resolvido uma vez por verificação, mantendo a proteção contra caminhos externos e evitando consultas de diretório por arquivo.

### Acesso aos resultados

Runs mostra uma área **Resultados disponíveis do pipeline**, com acesso direto a cada resultado preservado. O status da execução cancelada continua verdadeiro e separado do estado dos dados disponíveis. Mesmo que a origem tenha saído da janela das 100 execuções recentes, o botão busca aquela execução e abre seus resultados. Novos blocos inválidos/incompletos não ocultam os resultados válidos anteriores. Registros legados inválidos não impedem encontrar uma versão íntegra mais antiga.

Seletores e rótulos automáticos deixam de considerar uma etapa vazia/removida como substituta de um resultado anterior disponível. Novas execuções registram o evento `stage.results_preserved` com proveniência, facilitando futuras auditorias.

### Redocking somente como receptor

Removido o papel de moléculas/poses das saídas públicas de redocking. Conexões para portas de compostos, ligantes ou escores são rejeitadas no canvas e novamente na validação do backend. As escolhas de arquivos expõem somente receptores preparados/PDBs macromoleculares compatíveis; ligantes cristalográficos não são oferecidos como entrada molecular aos demais blocos.

Os arquivos próprios do redocking e seu contexto de referência não foram apagados: permanecem necessários ao acervo, RMSD, visualização e rastreabilidade da execução original. O filtro de handoff não remove arquivos científicos ou metadados internos.

## Testes executados

| Suíte | Resultado | Evidência |
|---|---:|---|
| Pipeline, persistência, invalidação seletiva, execução incremental e controles | 85 aprovados | [pipeline-regression.log](pipeline-regression.log) |
| Reúso, redocking/docking incremental, seletores, limpeza e feedback | 79 aprovados | [redocking-ui-regression.log](redocking-ui-regression.log) |
| Integridade, escopo dos arquivos, reinício e cache | 56 aprovados | [integrity-regression.log](integrity-regression.log) |
| Reabertura após mais de 100 cancelamentos, dados legados e navegação | 17 aprovados | [display-regression.log](display-regression.log) |

Há sobreposição entre as suítes; os números não representam testes únicos nem aprovação da suíte completa. Os testes usam projetos temporários e supervisores simulados; a auditoria do projeto real foi somente leitura. O fixture genérico de ADMET foi ajustado para fornecer colunas reais de saída, e o teste de curadoria de pares usa recuperação de PDB como origem, conforme o novo contrato de redocking.

Chrome local: visualização da área de resultados preservados junto à execução cancelada em 360, 1024 e 1440 px, sem exceções JavaScript de runtime. [360 px](runs-360.png), [1024 px](runs-1024.png), [1440 px](runs-1440.png), [metadados](browser.json).

Verificados a sintaxe dos 207 arquivos Python de produção e `git diff --check`. Nenhuma biblioteca instalada. Nenhuma execução científica do projeto original foi iniciada. Reinicie o BioMolExplorer para carregar a correção e execute com reúso para continuar pela nova MD.
