# Controles de Runs, persistência, portas e seleção MD — 2026-10-29

## Diagnóstico e implementação

1. **Clear All / Remove All:** `ui/removal.py` chamava `asyncio.sleep(0)` sem importar `asyncio`. O helper importava o contexto de processamento, mas não o utilizava. Agora importa o módulo explicitamente e inicia o splash padrão imediatamente após a confirmação, antes de qualquer remoção. Desabilita os botões da confirmação, bloqueia confirmações concorrentes e restaura os controles em falhas/cancelamentos. Uma confirmação aberta não é duplicada por cliques repetidos na barra de Runs; o estado é descartado ao sair do projeto. Exceções seguem pelo guard da interface, com aviso e possibilidade de tentar novamente.

2. **Save:** o binding existia, porém `save_draft` podia retornar silenciosamente se `dirty=False`, inclusive depois do autosave. O botão agora usa `save_pipeline_action`: mostra splash antes de esperar o lock/armazenamento, compara o pipeline com a base salva, cancela apenas o timer pendente de debounce, impede chamadas concorrentes e confirma o sucesso mesmo quando o autosave já persistiu o estado. Preserva o lock de escrita, revisão e merge existentes. Erros não deixam a ação bloqueada; uma resposta de outro contexto não gera confirmação de sucesso no projeto atual.

3. **Templates:** validação do DAG antes de aplicar, conexão automática pelo mesmo contrato tipado usado pelo canvas/backend, aplicação com splash e bloqueio de repetição, fechamento do diálogo específico e persistência imediata. Os sete modelos foram aplicados e relidos do banco. As entradas externas de grafos não recebem mais uma conexão automática indevida. Durante a regressão também foi corrigido o salvamento de um rascunho de redocking vazio: a seleção é obrigatória ao executar/retomar, mas não ao montar o pipeline.

4. **Run through selected stage:** o clique real transmite somente o ID selecionado ao serviço. O serviço existente já calcula o fechamento transitivo de ancestrais e os ordena pelo DAG. Testado com um ramo independente e um descendente: ambos ficam fora da execução. A ausência de seleção agora produz aviso explícito antes da submissão.

5. **Portas:** receptores MD usam o papel `receptor_pdb`; ligantes usam `compounds` / `molecular_poses`. Aliases estreitos mantêm estruturas macromoleculares e poses Vina/DOCK6 compatíveis, inclusive uploads válidos. Arquivos genéricos, arquivos genéricos de dinâmica e tabelas de escores sem papel molecular não são aceitos como substitutos desses dados. O consenso nativo continua fornecendo seu conjunto tipado de compostos. Não há wildcard `other` nas entradas do grafo. Portas ocultas e PubChem em modo manual não recebem novas arestas. Rejeições são transacionais no canvas e verificadas novamente no backend.

6. **Conteúdo molecular:** na seleção/validação MD, PDBs receptores exigem evidência de resíduos de proteína/ácido nucleico; MOL2 de ligantes deve declarar `SMALL`, evitando oferecer um `PROTEIN` como molécula pequena. A classificação é cacheada por assinatura de arquivo; mudança de conteúdo invalida a classificação. Se uma preparação candidata estiver inválida, um receptor bruto válido continua disponível. Essa identificação de papel não substitui a validação científica de coordenadas, topologia, cargas e parâmetros na execução.

7. **Available Files → Selected Files:** o vínculo automático mantém os arquivos disponíveis, mas não seleciona todos. A tabela de receptores e a matriz de pares começam vazias. Somente escolhas explícitas produzem candidatos/pares. Seleções explícitas já salvas são preservadas; projetos antigos com vínculo automático e IDs de pares salvos recuperam somente os receptores efetivamente associados àqueles IDs. Sem evidência de escolha salva, nada é selecionado implicitamente. Rótulos e avisos do seletor também foram completados em inglês/português.

## Validação executada

| Execução | Resultado | Evidência |
|---|---:|---|
| Regressão de controles, editor/autosave, limpezas, MD, cache, arquivos e canvas | 100 aprovados | [regression.log](regression.log) |
| Callbacks reais, tipagem, rascunhos, splash, login e poses consolidadas | 75 aprovados | [feedback-handoff.log](feedback-handoff.log) |
| Verificação final dos uploads tipados, seletores, controles e remoções | 58 aprovados | [selection-controls.log](selection-controls.log) |

As suítes se sobrepõem; os números não representam testes únicos nem aprovação da suíte completa do repositório. Falhas de disco/salvamento são injetadas intencionalmente para verificar captura de exceções e retomada; seus logs de erro são esperados. Os testes de limpeza utilizam o armazenamento real em projetos temporários e verificam confirmação, cancelamento, execução única, preservação dos originais e os diferentes efeitos de Remove All/Clear All. A execução dirigida usa o scheduler real com supervisor científico simulado, sem consultar provedores ou executar cálculos do usuário.

O pipeline atual, com 16 blocos, passou na validação das novas portas usando SQLite em **modo somente leitura**. Não foi alterado, executado ou limpo.

Chrome local em 360, 1024 e 1440 px: verificado o estado inicial com **34 arquivos disponíveis, 0 selecionados e 0 pares**, sem exceções JavaScript de runtime. Capturas [360 px](md-empty-360.png), [1024 px](md-empty-1024.png) e [1440 px](md-empty-1440.png); [metadados](browser.json). A auditoria visual usou dados fictícios.

Também verificados sintaxe dos 206 arquivos Python de produção e `git diff --check`. Nenhuma biblioteca foi instalada e nenhuma dependência adicional foi necessária. Banco, credenciais e resultados científicos originais permaneceram intocados.
