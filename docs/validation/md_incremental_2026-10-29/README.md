# MD: persistência incremental, seletores e tabelas — 2026-10-29

## Correções

- A reconciliação compara a edição com o pipeline salvo antes de reconciliar históricos. Adicionar/reconfigurar MD não autoriza limpeza de históricos upstream inalterados. Alterações reais de arquivos/ativos continuam detectadas; exclusões seguem com confirmação e transação existentes.
- Assinaturas de outros blocos não percorrem os módulos específicos de MD. A migração de fingerprints/similaridade foi auditada: trocar somente o validador atual por `58b8707^` reproduz exatamente os hashes do certificado anterior. A diferença é a introdução de ramos MD/CIF que esses algoritmos não utilizam. Aliases valem somente para os hashes exatos auditados, preservando versões dos produtores, validação de runtime, entradas e manifestos; mudanças científicas posteriores não herdam aprovação.
- Corrigidos import de `file_digest` e contador `docking_failures` ausentes no pipeline. O contador é persistido/restaurado por lote, evitando falhas em etapas que não executam docking.
- A projeção dos resultados respeita a invalidação seletiva já aplicada transacionalmente: a remoção de uma fonte não oculta o arquivo sobrevivente. Uma etapa historicamente vazia não é marcada obsoleta só por estar vazia.
- MD sem pares confirmados aguarda configuração. O handoff abre o formulário científico para escolher pares antes de executar. Corrigida também uma variável ausente que impedia a abertura de seletores de etapas.
- Receptores são agrupados por origem e cadeia, com preferência por arquivo preparado, sem descartar receptores de outros alvos. No projeto persistido examinado em modo somente leitura: 34 receptores existentes, provenientes de AChE, BChE e BACE1; plano de exclusão vazio.
- Tabela própria para receptores com identificador, visualizar e remover. Tabela de ligantes/pares sem número de resíduo ou visualização redundante de receptor; apenas visualizar ligante/remover. Ações fixas à direita, rolagem horizontal interna, alturas de linha iguais, elipse/tooltip, busca e filtro de alvo. IDs dos pares preservados.
- Remover a última entrada é permitido durante edição e não restaura silenciosamente uma seleção automática. A execução mantém validação estrita. Traduções de seletores e mensagens dinâmicas completadas em português/inglês.
- Fontes nativas Vina/DOCK6/consenso oferecem tabelas consolidadas de compostos/escores na seleção primária. Contextos, logs e poses intermediárias ficam ocultos; seleções explícitas antigas de poses continuam válidas e removíveis. Consenso com poses dos dois motores produz escolhas distintas, sem sobrescrever uma com a outra.

## Testes executados

Todos usam arquivos/projetos temporários; nenhum job foi enviado ao pipeline científico do usuário.

| Suíte | Resultado | Evidência |
|---|---:|---|
| Persistência, restart, crescimento, remoção seletiva, invalidação e migração | 42 aprovados | [cache.log](cache.log) |
| Seleção/handoff, preservação de entradas, splash e tabelas | 49 aprovados | [handoff.log](handoff.log) |
| Preparação, docking/redocking, login e feedback | 82 aprovados | [preparation-login.log](preparation-login.log) |
| MD, configuração, análise, fontes e tabelas | 48 aprovados | [md-ui.log](md-ui.log) |
| Fontes consolidadas e poses de consenso por motor | 3 aprovados | [consensus.log](consensus.log) |

Há sobreposição entre as suítes; esses números não representam testes únicos. Os fixtures de dois testes antigos de handoff foram alinhados aos nomes de arquivo que o gerenciador simulado realmente fornece. Sintaxe dos 205 arquivos Python de produção e `git diff --check` também verificados.

## Execução científica e navegador

Depois da restauração das alterações da sessão, executado novamente o exemplo 1LGU: Chimera → DOCK6 6.11 nativo, com carga do ligante preservada, em diretório temporário. Amber Score −2,798773 kcal/mol; RMSD final do ligante 0,491373 Å; preparação 1,07 s; cálculo nativo 7,21 s. [Log](native-1lgu.log). O refinamento nativo é estocástico: não se exige igualdade numérica bit a bit nem se interpreta o score como energia livre. Não há trajetória/RMSF nesta modalidade.

Chrome local: inspeção das tabelas com dados fictícios em 360, 1024 e 1440 px, sem exceções JavaScript de runtime. Capturas `pairs-*.png` e `table-*.png`; ações permanecem visíveis em telas estreitas. Isto não substitui uma auditoria de todas as telas da aplicação.

## Integridade e limites

Nenhuma biblioteca instalada; nenhuma dependência adicional necessária. Banco original aberto apenas com SQLite `mode=ro`; resultados originais não foram removidos, regravados nem recalculados. Durante a auditoria houve extração equivocada de uma revisão Python antiga sobre a cópia de trabalho; a ocorrência foi informada, os arquivos foram restaurados e as regressões acima foram executadas depois da restauração. A operação atingiu código Python, não banco nem resultados.

Reuso exige configuração, conteúdo e contrato científico compatíveis. A mera existência de um arquivo não autoriza reuso de saída corrompida, incompleta ou produzida por um protocolo incompatível. Nenhuma execução científica real do projeto atual do usuário foi utilizada como teste.
