# Dinâmica DOCK6: configuração e validação — 12/11/2026

Escopo: apenas configuração e execução de dinâmica (`dynamics_dock6`), com os contratos científicos e resultados dos outros blocos preservados. O serviço geral de salvamento ativa a otimização de auditoria somente quando todos os blocos alterados são de dinâmica e a revisão é compatível.

## Configuração e idiomas

- Categoria em maiúsculas e mensagens de progresso traduzidas no inglês; opções da região móvel recebem rótulos legíveis em português/inglês, preservando os valores de API `ligand`, `everything`, `nothing`, `distance` e `nab_atom_expression`.
- Histórico consultado pelo modelo de leitura compacto. Construção e leitura do formulário saem do loop da interface, com verificação de sessão/requisição após operações assíncronas.
- O seletor físico oculto da dinâmica deixa de construir menus para todas as poses. Tabelas consolidadas são priorizadas antes de inspecionar moléculas intermediárias.
- Catálogos de entradas são reutilizados enquanto os arquivos/referências não mudarem. A rastreabilidade é calculada para receptores, sem repetir a descrição do mesmo CSV para cada ligante.
- Pares e índice textual de busca são reutilizados. A tabela continua paginada em 50 linhas, sem recalcular entradas ao selecionar pares.
- Aplicação sem mudanças fecha o formulário sem novo salvamento. Edições exclusivas de dinâmica usam o arquivo imutável dos JSON nativos de execução na auditoria, preservando snapshots/rollback. Uma prévia compacta só evita carregar JSON completo se nenhuma remoção for necessária e as revisões dos runs permanecerem idênticas. Invalidação ou concorrência retornam ao plano nativo completo; provas compactas nunca substituem a proveniência autoritativa.

## Medição com o conjunto real

Medição isolada do formulário com os arquivos reais, em inglês, sem profiler e sem reconstruir resultados. Mesmos 5.187 pares candidatos e três pares selecionados; não representa tempo de rede, navegador ou todo o salvamento do projeto.

| Ação | Antes | Depois |
|---|---:|---:|
| Construir formulário | 6,510 s | 3,088 s |
| Atualizar pares, mesmas entradas | 2,347 s | 0,009 s |
| Buscar composto | 0,059 s | 0,002 s |
| Selecionar pares filtrados | 0,110 s | 0,001 s |
| Ler configuração final | 0,0012 s | 0,0012 s |

O benchmark anterior usa as classes de configuração do HEAD anterior, no mesmo processo de teste com entradas equivalentes; as fontes atuais são medidas em outro processo. Arquivos estão aquecidos pelo sistema operacional. Não é um benchmark de renderização física no navegador. O fechamento sem alterações e os contratos de salvamento/auditoria são cobertos por testes automatizados.

## 4M0E–CHEMBL92663

A pose representativa salva no consenso pertence a **2QMG_A**. A dinâmica agora recupera a pose **DOCK6 de 4M0E_A** pela tabela de origem autorizada e seus hashes, mantendo o ID do par e calculando sua identidade científica a partir da pose efetivamente usada. A resolução é feita somente para pares selecionados e também no preview. Ausência de uma pose correspondente falha explicitamente por par, sem executar a pose de outro alvo.

O formulário real foi parametrizado com exatamente um par; `GuidedForm.read` → `MDInputs.resolve` → wrapper nativo foi exercitado em uma área isolada. Os arquivos e a configuração persistidos do Alzheimer não foram modificados.

Parâmetros: 100 ciclos iniciais, 1.000 passos nativos (1 ps), 10.000 ciclos finais, 300 K, GB 5, região móvel somente do ligante, cálculo AM1-BCC do DOCK6. Foram mantidos os parâmetros do bloco salvo, sem reduzir o teste para facilitar aprovação.

A execução inicial identificou ligações SG–SG em CONECT com resíduos CYS e duas interrupções C–N de 11,35 Å e 4,52 Å. Na cópia de execução, resíduos de dissulfetos explícitos são identificados como CYX e interrupções de backbone (>2,0 Å) recebem TER. Hidrogênio de tiol incompatível ou coordenação ambígua são rejeitados. Não se movem átomos, removem cofatores ou inventam resíduos; a normalização fica em protocol.json.

O retorno total do DOCK6 6.11 usa float32, enquanto as contribuições são impressas em precisão dupla (`src/dock/score_amber.cpp`). A conferência da soma considera esse limite de arredondamento e a precisão decimal impressa, mantendo rejeição de valores não finitos/incompatíveis.

- Execução nativa: **25,35 s**, um par concluído, nenhuma falha.
- Repetição: **0,134 s**, um resultado reutilizado; nenhum novo refinamento.
- Validação adicional da seleção pelo formulário: **3,47 s** para construir/configurar; execução nativa em **25,93 s**, sem falhas, e repetição com cache em **0,152 s**.
- Amber Score: **+2391,533691 kcal/mol**; energia do complexo −17393,140872; contribuições do receptor +19769,696901 e do ligante +14,977599.
- RMSD final do ligante: **0,5840 Å** após ajuste do backbone; Rg inicial/final: **22,49696 / 22,49667 Å**.

A energia positiva é preservada exatamente como produzida pelo motor; conclusão da execução não significa interação favorável. É uma validação de um único sistema e de refinamento curto em solvente implícito, sem demonstração de equilíbrio, dinâmica convencional de produção ou energia livre de ligação. Não foi feita aceitação visual manual em um navegador real nesta rodada.

## Verificações automatizadas

- Suíte ampliada: 157 testes aprovados, incluindo interface, configuração, resultados, portabilidade e histórico.
- Regressão final focada: 75 testes aprovados; mensagens de progresso verificadas separadamente após a atualização dos catálogos.
- Pacote instalado fora do checkout: 25 testes de configuração, seleção de pares, resolução de consenso, precisão numérica, auditoria e invalidação aprovados.
- Auditoria de traduções dos módulos MD sem pendências. A identidade científica dos demais blocos permaneceu inalterada; a identidade do bloco MD mudou por causa das correções de preparação, seleção de pose e precisão do score.

## Revisão da porta semântica e dos rótulos Targets

O run real que falhou tinha 409 descritores do Redocking. Descritores dos targets 4M0E, 2WSL e 2QMG estavam presentes e disponíveis, com IDs corretos, mas sem `output_port`; outros descritores já tinham portas. A compatibilidade global só classificava manifests inteiramente legados e, portanto, recusava esse manifest misto.

A entrada do Amber Score agora normaliza apenas o contexto transitório dos produtores conectados em modo corrente. A porta ausente é recuperada pelo contrato tipado do produtor e pela validação estrutural PDB já usada na seleção MD. IDs, paths e proveniência são preservados. Não se altera um port declarado, não se inventam arquivos, não se consulta outro run como substituto e referências pinned continuam estritas. Os descritores científicos persistidos não são regravados. A correção não mudou nenhuma identidade científica e preserva caches existentes.

Os rótulos de canvas, formulário, cabeçalho de pares, busca e ajuda NAB usam Targets/target em português e inglês. Nomes de arquivos e chaves internas `receptor` permanecem intactos para compatibilidade.

A validação anterior começava após a resolução de portas. Esta revisão reproduziu o erro com descritores mistos e exerceu `PipelineService._resolve` sobre o snapshot real que falhou, antes da resolução MD e do motor nativo.

- As três referências passaram pela resolução semântica corrigida.
- 4M0E–CHEMBL92663 e 2QMG–CHEMBL92663 concluíram o Amber Score.
- 2WSL contém o resíduo modificado CSS, sem parâmetros na biblioteca DOCK6 instalada. A preparação desse par falhou de forma explícita; a estrutura não foi convertida artificialmente em CYS nem teve átomos removidos. Esse caso requer parametrização validada específica para CSS.
- A repetição dos dois pares concluídos reutilizou ambos os resultados em 0,256 s, sem novo refinamento.
- Testes exercitam manifests mistos, todos os produtores com contrato de target PDB, IDs/portas incorretos, arquivo ausente, porta declarada conflitante e referências históricas pinned.

A validação real foi feita em área isolada com as entradas originais somente para leitura; não substituiu a configuração ou os resultados do projeto do usuário.

## Instalação concluída e integridade após a retomada

A instalação anteriormente interrompida foi concluída no ambiente BioMolExplorer. A revisão final limitou a compatibilidade legada aos produtores ligados à entrada de targets (`base_input_path`), evitando inspecionar produtores usados somente para ligantes.

- 71 testes de regressão no pacote instalado e 32 testes finais de configuração/semântica/seleção passaram.
- As 15 identidades científicas foram preservadas.
- 185 arquivos de resultados, entradas e objetos de cache mantiveram os mesmos SHA-256 após a instalação.
- 4M0E e 2QMG foram recuperados do cache em 0,159 s. O teste bloqueou `md_engine.run_pair` com uma exceção para comprovar que nenhum motor nativo era invocado na recuperação.
- A execução de 2WSL foi conferida separadamente na versão instalada: a falha CSS continua explícita, sem resultado de energia fabricado e sem afetar os pares concluídos.

As verificações usam dados reais em área isolada, além de testes automatizados; não constituem garantia universal para qualquer estrutura, parametrização ou ambiente. A sessão já aberta da aplicação precisa ser reiniciada para carregar os módulos atualizados.
