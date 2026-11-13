# Auditoria documental — 13/11/2026

## Escopo e proteção do código

Esta auditoria usa os arquivos presentes no checkout, incluindo implementações não commitadas que já existiam antes da solicitação. **Somente documentação foi escrita nesta tarefa.** As alterações de código, configuração e testes que já apareciam no workspace foram preservadas; não foram revertidas, instaladas nem ampliadas.

A conferência inicial/final por SHA-256 protegeu **613 arquivos de código/configuração**, sem diferenças durante a auditoria. `pyproject.toml`, `requirements.yml`, arquivos Python/C/C++/JavaScript, scripts de build e testes não foram editados. CSS/JavaScript externos dos documentos também foram preservados. Não houve instalação de pacotes, execução de motores científicos ou edição de dados de projetos.

Foram inventariados e lidos automaticamente **139 documentos já existentes**, entre fontes Markdown, HTML, textos de instrução e relatórios históricos. A revisão semântica priorizou os guias ativos, confrontados com o registro de operações, parsers, assinaturas, validadores, resolvedores, bibliotecas de recursos e fluxos de persistência/UI. Os HTML derivados foram reconstruídos; evidências históricas foram mantidas como registros datados, sem simular novas validações.

## Fonte da verdade inspecionada

- `pyproject.toml`: versão 2.2.0, Python >=3.12, comandos `biomolexplorer`/`biomolexplorer-ui`, extra `ui` e ausência deliberada de dependências científicas de base.
- `requirements.yml`: ambiente científico/UI e ferramentas de compilação para suporte AMBER local.
- `__main__.py`, `ui/app.py`, `log_report.py`: comandos, flags, obrigatoriedade, padrões e idioma.
- `operations.py`, assinaturas dos wrappers e `catalog.py`: **15 operações científicas**, mais o bloco de importação exclusivo do workspace; chaves aceitas, campos internos e padrões.
- `dock6_protocol.py`, `consensus_model.py`, módulos `md_*` e `structure_sanitation.py`: modos de busca/avaliação, ECR, suporte estrutural e limites de parametrização.
- `updates.py` e `ui/updates.py`: instalação automática do pacote no host, SHA do commit, logs, ausência de atualização de dependências e reinício manual.
- `workspace.py`, `result_files.py`, `clear_results.py`, `pipeline.py`, `result_catalog.py`, `result_index.py`, `calculation_cache.py` e `ui/removal.py`: autorização, transações, caches, leitura e exclusão.
- Parsers dos scripts em `scripts/`: inventário de automações, sem executá-las como protocolos de validação científica.
- Guias de recuperação, ADMET, grafos, visualização, localização e exemplos/fixtures: coerência de caminhos, escopo e separação entre instruções atuais e evidências datadas.

## O que foi removido da orientação ativa

1. A descrição antiga do consenso como z-score/min-max, normalização constante para zero e fórmula que truncava energias positivas. Esses comportamentos não são apresentados como ativos; a compatibilidade de leitores históricos é explicitada.
2. A orientação de que o botão de atualização apenas baixa ZIP e exige extração/substituição manual de código no cliente web.
3. A alegação de uma limitação TOML atual que impediria o empacotamento, e comentários de inicialização que apontavam para essa limitação inexistente no código atual.
4. A afirmação de que toda proteína no Amber Score necessariamente usa ff94, independentemente de PTMs/perfil selecionado.
5. A promessa de reconstrução geral de cadeias laterais pelo DockPrep e de retry genérico automático após falha de carga quando o caminho ativo não os executa.
6. A apresentação do cache científico como invariavelmente exclusivo de um projeto.
7. A interpretação de status/testes antigos ou de uma instalação passada como certificação do checkout atual.

Relatórios históricos e referências bibliográficas não foram apagados. Foram identificados como históricos, porque remover a evidência destruiria rastreabilidade. Identificadores legados aceitos pelo código permanecem documentados como compatibilidade, não como novos recursos.

## O que foi adicionado

- Referências PT/EN derivadas do AST com os comandos públicos, todos os parâmetros registrados/padrões, opções da API `AppConfig`, variáveis lidas pelo código e flags declaradas nas automações.
- `dock6_mode`: distinção entre `auto`, `guided`, `denovo`, `score_only`, `restrained` e o contexto de `redock`; poses vinculadas a alvo versus conformeros genéricos; separação de grid scoring e Amber Score.
- Consenso ECR por composto: empates, singleton, agregação geométrica multialvo, cobertura, linha única por composto, precisão interna e três casas na apresentação; combobox de alvo antes de ações estruturais.
- Perfis `residue_parameter_policy`, saneamento genérico, 32 modelos PTM, parametrização de componentes orgânicos isolados e kernel de torções compatível em cópia local.
- Requisitos de fontes licenciados DOCK6/compiladores e escopo **CC BY-NC 4.0** dos parâmetros PTM, separado do código MIT.
- Limites explícitos: MSE sem modelo validado, metais/crosslinks especiais, inferência geral de protonação/tautômeros e reconstrução de loops não são garantidos.
- Exclusão MD na transação existente, fechamento do splash antes do diagnóstico e geração interna `result_cache_epoch`/`md_cache_epoch`; preparos e outros consumidores preservados.
- Compatibilidade MD para portas semânticas ausentes em manifests legados mistos, sem inventar arquivos ou reescrever artefatos persistidos.
- Índice de relatórios históricos e novos links nos manuais, índice documental e apresentação pública.

## O que foi corrigido

Comandos de instalação passaram a registrar efetivamente o aplicativo antes do uso dos console scripts. As tabelas required/optional foram regeneradas a partir do registro aceito, em vez de inferidas de rótulos. As versões PT/EN foram alinhadas; árvores e diagramas refletem os módulos e fluxos atuais. A documentação separa controles científicos de metadados internos e parâmetros Python de flags CLI. Os textos do atualizador distinguem host, cliente web, pacote, dependências e reinício. HTML, navegação e avisos de licença foram sincronizados com as fontes.

## Verificações realizadas

| Verificação | Resultado |
| --- | --- |
| Proteção de arquivos de código/configuração | 613/613 hashes preservados |
| Tabelas de parâmetros versus `OPERATIONS` | 30 tabelas válidas: 15 operações em PT e EN |
| Sintaxe de exemplos Bash/sh | 71 blocos aprovados por `bash -n`, sem executar comandos |
| Validador documental existente | 40 HTML com arquivos locais, IDs e âncoras válidos; 15 operações e importação presentes nos dois manuais e na apresentação |
| HTML reconstruído | 36 páginas de guias; portais mantidos compatíveis |
| `git diff --check` | Sem erros de whitespace |
| Instalação limpa / resolução Conda | Não executada nesta auditoria |
| Motores científicos / interface desktop interativa | Não executados nesta auditoria |
| Verificação de disponibilidade de provedores/DOIs externos | Não executada; observações anteriores permanecem identificadas como históricas |

O gerador existente regrava também CSS. Para respeitar a restrição, foram utilizadas suas funções de renderização em memória, gravando **somente HTML**, sem alterar nem executar o trecho de build que escreve CSS. As três referências complementares receberam páginas PT/EN sem modificar o registro executável de guias. O portal mantém seu registro de 15 guias principais; as referências complementares são alcançáveis pelos manuais e pela navegação das páginas geradas. Uma futura inclusão permanente desses guias no registro do gerador exigiria alteração de código e ficou fora desta tarefa. As fontes Markdown permanecem a referência de manutenção.

## Relação completa de arquivos alterados

Foram criados/atualizados **71 arquivos de documentação**: 33 Markdown e 38 HTML. A tabela inclui cada arquivo, inclusive os que mudaram apenas por regeneração da navegação/rodapé.

| Arquivo | Alteração |
| --- | --- |
| `README.md` | Consenso por composto; instalação e atualização reais; licença MIT do código separada dos dados PTM; links de referência; observação bibliográfica corretamente histórica. |
| `docs/ARTIFACT_REFACTOR_STATUS.md` | Status/testes de 06/11 identificados como registro histórico, não certificação atual. |
| `docs/README.md` | Índice bilíngue de referências atuais, relatório e escopo dos documentos históricos. |
| `docs/amber_score_structural_support.md` | Validação/instalação de 13/11 explicitamente histórica; não atesta instalação de mudanças posteriores. |
| `docs/architecture.md` | Árvore/fluxo atualizados; retirada da falsa limitação TOML; cache compartilhável, módulos novos e portas semânticas MD legadas. |
| `docs/backend_usage.md` | Tabelas sincronizadas com o registro real; parâmetros novos, consenso ECR, modos DOCK6, suporte AMBER e exclusão/cache. |
| `docs/chimera_migration.md` | Suporte estrutural AMBER integrado e referência ao caminho atual, conservando compatibilidade dos nomes. |
| `docs/current_reference.md` | Novo inventário derivado do AST: CLIs, 15 operações, parâmetros/padrões, AppConfig, variáveis de ambiente e flags dos scripts; legado versus ativo. |
| `docs/documentation_audit_2026-11-13.md` | Escopo, inventário, verificações, limites e relação completa dos arquivos desta auditoria. |
| `docs/en/admet_research.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/amber_score_structural_support.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/amber_score_structural_support.md` | Nova referência em inglês de saneamento genérico, perfis PTM/ff03, componentes isolados, compatibilidade local, cache, licenças e limites. |
| `docs/en/architecture.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/architecture.md` | Árvore/fluxo atualizados; retirada da falsa limitação TOML; cache compartilhável, módulos novos e portas semânticas MD legadas. |
| `docs/en/backend_usage.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/backend_usage.md` | Tabelas sincronizadas com o registro real; parâmetros novos, consenso ECR, modos DOCK6, suporte AMBER e exclusão/cache. |
| `docs/en/chimera_migration.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/chimera_migration.md` | Suporte estrutural AMBER integrado e referência ao caminho atual, conservando compatibilidade dos nomes. |
| `docs/en/current_reference.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/current_reference.md` | Novo inventário derivado do AST: CLIs, 15 operações, parâmetros/padrões, AppConfig, variáveis de ambiente e flags dos scripts; legado versus ativo. |
| `docs/en/dms_migration.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/docking_consensus.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/docking_consensus.md` | Nova referência em inglês do ECR multialvo, cobertura, estereoquímica, linha única por composto e ações por alvo. |
| `docs/en/frontend.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/frontend.md` | Retirada de consenso antigo e instalação manual via botão; comportamento real do atualizador, seletores por alvo e ciclo de exclusão/splash. |
| `docs/en/installation.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/installation.md` | Requisitos de compilação e perfis AMBER; parâmetros externos não instalados pelo pip; referência do código atual. |
| `docs/en/logging.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/molecular_preparation.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/molecular_preparation.md` | Modos de operação DOCK6 e suporte AMBER; distinção entre DockPrep e preenchimento pelo template AMBER; correção de retry histórico em inglês. |
| `docs/en/molecular_viewer.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/pipeline_validation.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/pipeline_validation.md` | Referência explícita ao estado atual do código e separação de verificação documental/científica. |
| `docs/en/projects.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/projects.md` | Transação única na exclusão; nova geração de resultados MD; preservação de preparos e distinção da exclusão de pares. |
| `docs/en/redocking_configuration.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/redocking_configuration.md` | Busca versus score-only/minimização restrita e exigência de coordenadas; separação de grid e Amber Score. |
| `docs/en/retrieval.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/user_manual.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/en/user_manual.md` | Perfil ff94/ff03 condicional, consenso atual, suporte AMBER e links detalhados; retirada de comentário sobre limitação de empacotamento. |
| `docs/frontend.md` | Retirada de consenso antigo e instalação manual via botão; comportamento real do atualizador, seletores por alvo e ciclo de exclusão/splash. |
| `docs/installation.md` | Requisitos de compilação e perfis AMBER; parâmetros externos não instalados pelo pip; referência do código atual. |
| `docs/molecular_preparation.md` | Modos de operação DOCK6 e suporte AMBER; distinção entre DockPrep e preenchimento pelo template AMBER; correção de retry histórico em inglês. |
| `docs/output_datasets.md` | Diferença entre remover resultados AMBER e excluir pares da seleção; geração de cache e preparos preservados. |
| `docs/pipeline_validation.md` | Referência explícita ao estado atual do código e separação de verificação documental/científica. |
| `docs/processing_feedback.md` | Fechamento do splash antes de apresentar erro e restauração dos controles de exclusão. |
| `docs/projects.md` | Transação única na exclusão; nova geração de resultados MD; preservação de preparos e distinção da exclusão de pares. |
| `docs/pt/admet_research.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/amber_score_structural_support.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/architecture.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/backend_usage.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/chimera_migration.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/current_reference.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/dms_migration.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/docking_consensus.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/frontend.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/installation.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/logging.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/molecular_preparation.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/molecular_viewer.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/pipeline_validation.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/projects.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/redocking_configuration.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/retrieval.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/pt/user_manual.html` | Página HTML sincronizada com a fonte Markdown; navegação para referências complementares e aviso de licença do código/dados PTM. |
| `docs/redocking_configuration.md` | Busca versus score-only/minimização restrita e exigência de coordenadas; separação de grid e Amber Score. |
| `docs/result_indexing.md` | Separação entre leitura derivada e exclusão transacional; diretório real e possível compartilhamento do cache. |
| `docs/user_manual.md` | Perfil ff94/ff03 condicional, consenso atual, suporte AMBER e links detalhados; retirada de comentário sobre limitação de empacotamento. |
| `docs/validation/README.md` | Índice de evidências datadas; orientação para não usar protocolos antigos como instrução atual. |
| `index.html` | Descrição pública do consenso por composto e suporte AMBER; limites científicos e aviso de licença PTM. |
| `index_doc.html` | Links para referências atuais e auditoria; scripts e estilos originais preservados. |

## Limites de conclusão

A conclusão é documental: contratos descritos e exemplos estaticamente verificáveis foram sincronizados com o checkout, preservando todo o código. As métricas e a aceitação científica de relatórios anteriores não foram reexecutadas. A auditoria não declara a plataforma inteira pronta para produção e não afirma que alterações pendentes de implementação tenham sido instaladas.
