# Revisão de resultados, controle de execução e regressão do redocking

Validação de 30/10/2026. Nenhum pacote foi instalado; `requirements.yml` não precisou de alterações.

## Causas e correções

### Runs, persistência e cache

A interface construía cartões de cada execução e um painel separado de resultados publicados. Isso fragmentava a apresentação: uma execução recente cancelada podia esconder resultados válidos de execuções anteriores.

Runs agora apresenta somente **Latest Run**, com uma projeção dos resultados atuais do pipeline. A seção **Available Pipeline Results** e os cartões de execuções anteriores foram removidos dessa tela. O histórico durável continua disponível pelo acesso de histórico do projeto.

A projeção usa o catálogo existente de resultados publicados, seus manifestos e as assinaturas de requisição/entrada. Cada resultado mantém o ID de sua execução produtora para autorização, downloads, inspeção molecular e rastreabilidade. Não é criado um run sintético no banco, nem são copiados arquivos para agregar resultados.

Quando há somente checkpoints, o catálogo verifica cada unidade concluída. Unidades complementares de diferentes execuções, inclusive da execução pausada atual, são reunidas sob o mesmo bloco. Unidades equivalentes são deduplicadas por assinatura e hashes de entradas; sua versão original de processamento é preservada. A política `verified_completed_units_v1` distingue esses checkpoints de uma etapa integralmente concluída e registra a assinatura, versão de processamento e chave de cada unidade validada. Em agregações parciais com múltiplas origens, os arquivos são apresentados por conjunto usando o componente existente `FileTable`, sem inventar uma tabela científica final. A conclusão da etapa mantém suas visualizações especializadas.

A integridade também é conferida antes de exibir os resultados da execução mais recente. Arquivos alterados externamente ou ausentes não são apresentados como resultados válidos; versões compatíveis anteriores podem continuar disponíveis. Essa verificação não apaga arquivos nem modifica registros históricos.

Os mecanismos existentes de invalidação seletiva, confirmação e cache por unidade foram mantidos. O histórico usa os blobs endereçados por conteúdo já existentes. A ordem de execução passou a integrar os snapshots e sua alteração gera um evento reversível, sem duplicar blobs de resultados. Snapshots antigos continuam legíveis.

Durante a regressão foi encontrada uma segunda implementação de `remove_project_results`, que sobrescrevia a anterior. A implementação efetiva protegia os bytes de fontes legadas não registradas, mas apagava suas referências. O caminho unificado agora preserva ambos em **Remove All**. A definição morta foi removida. **Clear All** mantém sua política de restauração de padrões e limpeza, com confirmação, feedback e transação existentes.

### Pause, Resume e Cancel

Foram reutilizados os mecanismos existentes de pausa segura e retomada. **Pause Run** e **Resume Run** estão ligados aos handlers em Runs e no acompanhamento. A pausa termina a operação atual e para entre operações/lotes. **Ver resultados** também fica disponível durante a pausa manual. O cancelamento encerra o avanço e conserva unidades concluídas.

Os testes confirmam: retomada do mesmo run pausado após reiniciar o serviço; recuperação de checkpoints de execuções interrompidas/canceladas; ausência de recálculo de unidades válidas; execução nova quando o reúso está desabilitado; bloqueio de leitores e de ações incompatíveis. Uma interrupção abrupta continua sendo registrada como interrupção, conforme a política existente; uma nova execução com reúso recupera o trabalho concluído.

### Redocking

Três causas foram comprovadas:

1. **Vina apagava uma pose que acabara de gerar.** `VinaSessions.run()` acessava `fields` sem defini-lo ao registrar conclusão ou erro. O supervisor tratava o `NameError` como falha e removia a pose. O arquivo real `vina_execution.json` de 4M0E/1YL registrava esse erro. O contexto de logging agora é definido antes do cálculo; erros reais do motor continuam sendo propagados.
2. **A opção `prep false` era emitida duas vezes.** A configuração já a incluía, e a fase de carga a acrescentava novamente. A reprodução de 6ZWI/QRH com `minimize spec #0 prep false prep false`, seguida de exportação PDB/MOL2, alterou a protonação: a carga formal observada passou de +1 para +2. A validação corretamente rejeitou a mudança. A geração agora mantém exatamente uma opção `prep false`, conservando os demais parâmetros de minimização. Com a mesma sequência de exportação corrigida, o grafo e a carga +1 foram preservados.
3. **A validação química era antecipada no intermediário errado.** `validate_mol2(..., sanitize=False)` ainda chamava `Chem.SanitizeMol()` incondicionalmente. O intermediário MOL2 anterior ao AddH, embora derivado de um grafo CCD já validado, podia não representar corretamente o H de nitrogênios aromáticos. O parâmetro agora funciona conforme seu contrato: nesse ponto são verificadas estrutura do arquivo, contagens, coordenadas, cargas numéricas e conectividade. Após adicionar H, carregar e minimizar, a sanitização e as validações de grafo/carga continuam estritas e habilitadas por padrão. O caso real 6QAA/HUN concluiu o preparo com 66 átomos, carga formal +1 e soma de cargas parciais 1,0012.

A identidade do protocolo de ligantes inclui os helpers modificados. Isso impede reutilizar um intermediário de preparo incompatível; a identidade dos receptores permanece independente desses helpers e sua reutilização foi testada.

A tabela diferencia falhas operacionais de falhas de preparação, incluindo registros antigos identificados pelo campo `FAILED_STEP`. RMSDs altos continuam sendo resultados numéricos, sem serem mascarados ou convertidos em falhas operacionais.

### Configurações e ADMET

`GuidedForm.layout()` envolvia novamente os conteúdos de `ExpansionTile` em cada chamada, acumulando margens. O layout agora é idempotente. Outra causa era a ausência de registro da janela em `settings_dialog` e no rastreador de configuração existente: janelas podiam se empilhar e callbacks antigos continuavam habilitados. O fluxo agora registra, invalida e fecha os diálogos corretos, inclusive após confirmação; falhas de confirmação permitem tentar novamente.

O componente compartilhado mantém ações separadas da rolagem do conteúdo, com margens e espaçamento uniformes. Rótulos de switches/checkboxes podem quebrar linha, preservando acessibilidade. Não foram impostas dimensões fixas novas aos formulários.

As listas de arquivos de origem e a curadoria dos pares de redocking foram mantidas porque representam entidades diferentes. A preparação continua única para receptores e única para ligantes. Não foi removida funcionalidade de seleção, edição ou inspeção.

Os cinco ícones de **Alerts for Chemical Review** tinham somente tooltip, sem evento de clique. Cada um agora abre a descrição de seu próprio catálogo pelo componente de ajuda existente. As descrições são traduzidas para a sessão pt/en. Nenhum critério científico do ADMET foi modificado.

Também foi corrigido o import ausente de `save_file` usado para baixar relatórios de exclusões e outros arquivos da interface.

## Testes e evidências

Ambiente existente: `/home/michel/progs/miniconda3/envs/BioMolExplorer/bin/python`, `PYTHONPATH=src:tests`, `MPLCONFIGDIR=/tmp/biomol-mpl`, sem novos pacotes.

- **181 testes passaram** na suíte de aceitação de resultados, cache, histórico, invalidação, pausa/retomada/cancelamento, limpeza, configuração e visualização. [Log](acceptance-tests.log).
- **12 testes dedicados ao Vina passaram**, incluindo propagação de erro, validação de geometria, reutilização de mapas/receptor e supervisor nativo isolado. [Log](vina-tests.log).
- **72 testes adicionais de integridade da visualização passaram** após acrescentar a detecção de arquivos alterados externamente. [Log](view-integrity-tests.log). As suítes se sobrepõem; seus totais não devem ser somados.
- O teste de checkpoints complementares cria duas execuções interrompidas com resultados diferentes, confirma um único bloco com ambos os conjuntos, verifica a visibilidade de unidades da execução pausada e retoma sem novas chamadas científicas para os lotes já concluídos.
- O teste de histórico confirma restauração da ordem e igualdade do conjunto de blobs antes/depois de mudanças de ordem. A regressão de limpeza mantém as verificações de originais, referências, restauração, permissões e indicadores.
- Os fixtures de testes antigos foram alinhados aos contratos atuais: o pipeline de apresentação contém os blocos que geraram seus resultados; contagens de limpeza usam o retorno estruturado atual; os indicadores vêm do catálogo de idioma; o objeto de teste do Vina inicializa o mesmo atributo do construtor. Nenhum teste foi ignorado para aprovar as alterações.
- **48 telas renderizadas:** 16 blocos em 1440, 1024 e 360 px; os dois botões de ação ficaram dentro do viewport. [Verificação](all-forms-checks.json).
- **Reaberturas repetidas:** posições de Cancel/Apply permaneceram estáveis nos quatro formulários representativos, nas três larguras. [Verificação](repeated-layout-checks.json).
- **15 popups ADMET renderizados:** cinco catálogos nas três larguras. Os testes de eventos também cobrem pt/en.
- Latest Run após cancelamento e a composição de checkpoints foram renderizados nas três larguras. Um primeiro teste de acessibilidade do navegador excedeu o timeout de 15 s; a repetição com tolerância de inicialização maior concluiu. As execuções finais não registraram exceções JavaScript; o registro de console foi a mensagem informativa `Flutter app loaded`.
- `git diff --check` sem erros.

Imagens: [configuração redocking no celular](redocking-settings-mobile.png), [alerta PAINS](admet-info-mobile.png), [Latest Run](latest-run-desktop.png), [checkpoints complementares](complementary-checkpoints-mobile.png).

## Execuções científicas reais

Todas ocorreram em cópias temporárias. Os parâmetros de busca Vina foram caixa 24 Å, exhaustiveness 20, 10 modos, seed 73 e um núcleo. Para os ligantes corrigidos, o centro foi calculado pelo mecanismo de centro de massa do próprio pipeline. O preparo conservou AM1 e a minimização configurados. DOCK6 usou o protocolo rígido e a versão 6.11 comprovada no arquivo nativo.

| Par | Motor | RMSD (Å) | Score nativo | Resultado computacional |
|---|---|---:|---:|---|
| 4M0E / 1YL / 604 / A | Vina 1.2.5 | 0,194442 | -11,998 | Concluído |
| 6ZWI / QRH / 604 / A | Vina 1.2.5 | 4,625933 | -8,513 | Concluído; reprodução cristalográfica desfavorável |
| 6QAA / HUN / 613 / A | Vina 1.2.5 | 6,965525 | -9,356 | Concluído; reprodução cristalográfica desfavorável |
| 4M0E / 1YL / 604 / A | DOCK6 6.11, rígido | 0,277706 | -47,111782 | Concluído |

Scores de motores diferentes não são tratados como diretamente comparáveis. Nenhum limiar de RMSD, caixa, regra de seleção, validação de carga ou validação final de topologia foi relaxado para obter aprovação.

[Casos científicos](scientific-cases.json), [Vina 4M0E](vina-4M0E.json), [DOCK6](dock6-redocking.json), [Chimera antes](chimera-before.log), [Chimera depois](chimera-after.log), [preparo aromático](aromatic-preparation.log).

Os scripts de reprodução usados estão nesta pasta. Contêm caminhos deste checkout e escrevem somente em `/tmp`; não devem ser interpretados como execução do pipeline inteiro nem como atualização dos resultados do projeto real.

## Integridade dos dados originais e limites

O banco real foi aberto somente com SQLite `mode=ro`. A execução original continua com seus 161 diagnósticos de falha e 1.542 artefatos presentes. Os SHA-256 de receptores/ligantes/complexo selecionados permanecem iguais às cópias feitas antes das validações. Nenhuma limpeza, reescrita de manifesto ou substituição de resultados foi realizada no projeto real. [Auditoria](original-integrity.jsonl).

Os resultados corrigidos dos testes não foram injetados no experimento do usuário. Para produzir novos resultados na plataforma, é necessário reiniciar o aplicativo atualizado e executar com reúso; etapas válidas são preservadas e o trabalho inválido/pendente é processado.

Não foi refeita a coleção completa de 161 pares nem o pipeline científico integral. Não há aprovação global de todas as químicas, configurações ou motores. Uma varredura exploratória mais ampla encontrou expectativas legadas fora deste recorte em `test_project_portability` (nome de saída ADMET antigo) e em `test_redocking_chimera_repair` (texto de erro e fallback químico antigo). Esses testes não são declarados aprovados. As validações científicas não foram desabilitadas para satisfazê-los.

## Arquivos alterados nesta revisão

- Resultados/persistência: `published_results.py`, `result_visibility.py`, `workspace.py`, `project_state.py`, `execution_order.py`, `clear_results.py`, `result_files.py`.
- Preparo/execução: `chimera_preparation.py`, `dock6_validation.py`, `structure_cache.py`, `vina_session.py`, `redocking_results.py`.
- Interface: `ui/app.py`, `ui/stage_results.py`, `ui/run_progress.py`, `ui/guided.py`, `ui/settings_layout.py`, `ui/admet_settings.py`, `ui/redocking_results.py`.
- Traduções: catálogos `en.json`, `pt.json`, `runtime_en.json`, `runtime_pt.json`.
- Regressões: novo `tests/test_unified_results_review.py` e ajustes em `test_published_results.py`, `test_execution_order_pause.py`, `test_vina_reuse.py` e `test_results_display.py`.
