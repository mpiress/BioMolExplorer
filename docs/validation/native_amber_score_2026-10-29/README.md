# Amber Score nativo, preparação Chimera e acesso ao workspace — 2026-10-29

## Entrega e escopo

A categoria **Dinâmica Molecular / Molecular Dynamics** já existia. Seu bloco foi restringido a **Amber Score · DOCK6**, exclusivamente com o Amber Score nativo do **DOCK6 6.11**. O identificador `dynamics_dock6`, os vínculos e a leitura dos resultados históricos foram preservados. Não foram instalados pacotes.

- Removidos o protocolo executável `dock6_sander`, as chamadas de trajetória `sander/cpptraj`, seus controles na interface e a criação de um ambiente AmberTools pelo instalador.
- Removidos campos invisíveis/obsoletos de configuração de preparação e ramos de cache de um backend standalone que os blocos ativos já normalizavam para Chimera.
- Cargas AM1-BCC/Gasteiger do adaptador molecular usam utilitários distribuídos com **UCSF Chimera**, com executáveis absolutos e bibliotecas apenas no ambiente do processo filho. Não se procura Conda, `BIOMOL_AMBERTOOLS_*` ou executáveis AmberTools no PATH.
- Minimização desse adaptador utiliza Chimera/MMTK, tipagem GAFF, 100 passos SD e 10 CG, `prep false`, mantendo cargas, conectividade e mapeamento. Seu checkpoint é durável, separado do diretório temporário; entradas, coordenadas completas, cargas e implementação integram sua assinatura.
- O adaptador de topologia GAFF usa apenas as ferramentas distribuídas com DOCK6 6.11 e verifica as cargas fornecidas. O formato de cargas em grupos de oito valores foi corrigido para o leitor antigo incluído no DOCK6. Não se anuncia suporte GAFF2.
- O antigo protocolo experimental de receptor ff14SB/sander não foi falsamente relabelado como ff94 ou Chimera. Seus pontos de compatibilidade recusam execução com instruções claras. A entrada pública antiga de preparo de estruturas direciona ao preparo Chimera já utilizado pelos blocos atuais.
- A escolha histórica `dock6_sander` é recusada até seleção explícita do protocolo nativo. Ao aplicar essa escolha na interface, campos aposentados são descartados da nova configuração; resultados antigos continuam disponíveis. Um certificado congelado de cache legado não autoriza reutilização de uma versão científica alterada: arquivos anteriores são preservados e o cálculo é refeito conservadoramente.

## Limites científicos e parâmetros reais

O executável selecionado precisa informar `DOCK v6.11`. São usados `prepare_amber.pl`, `amberize_*`, `teLeap`, `antechamber` e `parmchk` **incluídos no DOCK6**, com cópias locais dos helpers relocados e bibliotecas da própria instalação. Não se modifica a instalação global. `AMBERHOME` é retirado do ambiente filho.

O Amber Score usa ff94 para o receptor e GAFF para o ligante. A opção de conservar cargas MOL2 refere-se ao ligante; a proteína é parametrizada pelo DOCK6. A interface informa essa diferença. Histidinas genéricas `HIS` exportadas pelo Chimera são identificadas como HID/HIE/HIP quando há evidência explícita de HD1/HE2, sem remover átomos ou mover coordenadas; a decisão é registrada no protocolo.

Configurações expostas: ciclos de minimização antes/depois, convergência, passos do refinamento interno (inclusive zero), temperatura, GB 1/2/5, corte não ligado, região móvel `ligand/everything/nothing/distance/nab_atom_expression`, distância da região móvel, expressões NAB do receptor/complexo, aborto de pares não preparados e conservação das cargas do ligante. Prefixos de arquivos e esferas do sítio são gerados a partir do par selecionado. Valores inválidos/não finitos são recusados.

O passo de integração nativo é fixo em 0.001 ps, Langevin gamma 2, RATTLE desativado. O cálculo final inclui GB/SA. **Não há caixa de água explícita, barostato, aquecimento/equilíbrio de produção, sander ou substituição por outra função de pontuação.** O refinamento curto interno faz parte do próprio Amber Score; não constitui simulação convencional validada ou energia livre de ligação.

As métricas disponíveis são:

| Métrica | Disponibilidade e definição |
|---|---|
| Amber Score e contribuições | Quatro valores finitos; validação da soma das contribuições assinadas produzidas pelo DOCK6 |
| RMSD de backbone e ligante | Comparação do complexo preparado inicial com o final, após ajuste rígido do backbone proteico; ligante sem reajuste ou correção de simetria |
| Raio de giro | Ponderado por massa, estruturas inicial/final do complexo; somente inicial quando não existe estrutura final |
| Hidrogênios | Contagem de átomos/isótopos pela massa na topologia; não é contagem de ligações H |
| RMSF, ocupância e séries temporais | Indisponíveis: o caminho nativo verificado não publica trajetória utilizável |

A ausência de estrutura final, como no modo `nothing`, é informada; não se fabricam RMSD/RMSF. Receptores sem backbone proteico reconhecido mantêm suas energias e métricas estáticas possíveis. Metais/cofatores/resíduos sem parâmetros compatíveis permanecem falhas explícitas, sem exclusão silenciosa ou scores aproximados. Os analisadores Python de trajetórias históricas foram preservados para compatibilidade de leitura.

Verificação primária: fonte **instalada do DOCK6 6.11**, `src/dock/score_amber.cpp`, `src/dock/nab/sff.c` e `src/dock/version.h`. O caminho de coordenadas em `mme_init` não abre um arquivo de trajetória; os reinícios finais são usados como estruturas, não como amostras temporais. Referências oficiais: [tutorial Amber Score](https://dock.compbio.ucsf.edu/DOCK_6/tutorials/amber_score/amber_score.htm) e [FAQ do DOCK6](https://dock.compbio.ucsf.edu/DOCK_6/faq.htm). O manual online atual pode corresponder a outra versão; não foi usado para inventar parâmetros da 6.11.

## Validação científica executada

O exemplo **1LGU**, distribuído com DOCK6 6.11, foi preparado por Chimera (receptor PDB e ligante MOL2 com AM1-BCC) e executado com os **3.000 passos nativos padrão**. A fase nativa recebeu PATH `/usr/local/bin:/usr/bin:/bin`, sem o ambiente Conda externo, e `AMBERHOME` intencionalmente inválido no pai; a preparação do DOCK6 o removeu corretamente.

| Resultado de uma execução | Valor |
|---|---:|
| Amber Score | −2.761784 kcal/mol |
| Complexo | −5082.720620 kcal/mol |
| Contribuição do receptor | +5000.012698 kcal/mol |
| Contribuição do ligante | +79.946138 kcal/mol |
| RMSD final do ligante | 0.487560 Å |
| Rg inicial/final | 16.465162 / 16.465705 Å |
| Átomos de hidrogênio | 1315 |
| Preparo Chimera | 1.03 s |
| Amber Score + preparo/análise | 6.73 s |

Também executados os modos `nothing`, `distance`, `nab_atom_expression` e `everything`, com refinamentos curtos e energias finitas. Foram testados conservação de cargas, cargas AM1-BCC, reutilização de preparação após alteração de temperatura, inclusão/remoção de pares, corrupção/invalidação e isolamento de falhas.

Esses valores são uma observação em uma máquina e um sistema, não benchmark estatístico, afinidade experimental ou aceitação científica para todos os complexos. Etapas de preparação podem gerar pequenas diferenças geométricas entre execuções; conservação dos inputs, hashes, topologias, cargas e protocolo é necessária para rastreabilidade.

Dados: [scientific-validation.json](scientific-validation.json). Exemplos de gráficos reais: [RMSD estrutural](structural_rmsd.png) e [raio de giro](structural_rg.png).

## Correção de acesso à conta

Os logs locais registraram autenticação aceita seguida de falha **durante a abertura do workspace**: `RuntimeError: Dialog is already opened`. O Flet conserva o diálogo fechado no registro até receber o evento de dismiss do cliente. A reutilização do controle (inclusive com igualdade estrutural antes da montagem) causava a falha.

O splash agora usa uma rota nova com chave única a cada ciclo, compartilha o diálogo ainda aberto entre requisições simultâneas e preserva títulos/idioma configurados. O feedback de salvar foi unificado para evitar sobreposição de duas rotas. Erros posteriores à autenticação recebem uma mensagem de abertura de workspace, não de credenciais incorretas. Senhas e contas reais não foram alteradas.

Testes usam banco exclusivamente temporário: cadastro, fechamento/reabertura do store, autenticações consecutivas, rejeição de senha errada e clique real no botão de login até renderizar **Your workspace**, antes de qualquer dismiss do cliente.

## Regressão e dependências

- **39 testes científicos/contratos** aprovados: preparação, cargas, GAFF, minimização Chimera com reuso, pares nativos, cache e configurações/resultados de Amber Score. [Log](scientific-tests.log).
- **57 testes de regressão** aprovados: Chimera, preparação, docking/redocking, cache incremental e configurações. O teste de certificado legado foi atualizado para exigir recálculo de versão científica aposentada e verificar a preservação dos arquivos antigos. [Log](preparation-regression.log).
- **26 testes de login, splash e confirmação/reinício** aprovados, incluindo conta persistida e abertura do workspace pelo formulário real. [Log](login-and-splash.log).
- **9 testes de análise/interface** aprovados após complementar as métricas iniciais no modo sem estrutura final. [Log](analysis-ui-tests.log).
- Auditoria AST de **356 arquivos Python**: nenhum import de produção sem declaração nem pacote órfão comprovado no ambiente principal. `poppler` e dependências de build/runtime foram mantidos com justificativa. Playwright é apenas ferramenta opcional de auditoria visual e não está instalado; não foi incorporado como dependência de runtime. [Auditoria](dependency-audit.json).
- `requirements-ambertools.yml` já estava ausente e não era rastreado neste checkout; não havia arquivo para `git rm`. `requirements.yml` não declarava AmberTools; não foram removidos pacotes utilizados. O instalador e a documentação/index foram corrigidos para não solicitar esse ambiente inexistente.
- Reconstruídas 30 páginas HTML e dois portais; `bash -n install.sh` e `git diff --check` concluíram sem erros.

Os testes também revelaram código incompleto já presente: referências de splash/imports ausentes, esquema `block_ids`, campos de cache do AppConfig, `_prepare_one`, mapa de ligantes do Chimera e publicação de checkpoints de redocking. Foram corrigidos os pontos necessários à execução e regressão desses fluxos, mantendo confirmações de exclusão e os dados persistidos.

Não foi declarada aprovação da suíte completa do repositório ou de uma auditoria visual em navegador. Uma execução exploratória do módulo global de i18n encontrou traduções ausentes em outras telas; o novo formulário Amber Score foi validado isoladamente em inglês e a troca de idioma preserva suas configurações. Nenhum projeto científico real do usuário foi executado, apagado ou regravado nesta validação. Reinicie o processo BioMolExplorer para carregar a correção do login.
