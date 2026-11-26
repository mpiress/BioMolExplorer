# Auditoria da documentação — 26/11/2026

## Resultado

Esta revisão parte do inventário de funcionalidades derivado do código presente neste checkout, incluindo implementações ainda não commitadas, e confronta a documentação ativa com ele. Nenhum módulo científico foi alterado por esta auditoria. Os relatórios datados em [validation](validation/README.md) e as auditorias de [13/11](documentation_audit_2026-11-13.md) e [15/11](documentation_audit_2026-11-15.md) descrevem o estado daquelas datas e foram preservados como evidência histórica.

## Método do inventário

| Fonte consultada | Inventário extraído |
| --- | --- |
| `operations.py`, `catalog.py` e assinaturas dos wrappers | 16 operações científicas registradas e o bloco `import_results`, com obrigatórios, opcionais e padrões das funções |
| `artifact_contracts.py` | Portas de entrada/saída e tipos aceitos por bloco |
| `ui/app.py`, `ui/stage_results.py`, `ui/pose_results.py`, `ui/run_progress.py` | Abas, ações da barra, decisões, apresentação de resultados por bloco e acompanhamento de execuções |
| `resources/viewer/*` e `resources/i18n/viewer_*.json` | Estilos, controles e rótulos dos visualizadores |
| `pipeline.py`, `result_catalog.py`, `workspace.py` | Persistência do estado da execução, projeções de leitura e journal do banco |
| `scripts/*.py` | Scripts de automação e suas flags declaradas |

O validador documental do próprio projeto (`scripts/validate_docs.py`) foi usado como teste de cobertura: ele exige que cada operação ativa apareça exatamente uma vez nos cartões de `index.html` e nos exemplos dos dois manuais.

## Funcionalidade não documentada

O bloco **Validar poses** (`pose_validation`) estava implementado, registrado no catálogo, coberto por testes e sem nenhuma menção na documentação, na apresentação da plataforma ou nos portais. Antes desta auditoria, `scripts/validate_docs.py` acusava três ausências por causa disso.

Documentação criada:

- `docs/pose_validation.md` e `docs/en/pose_validation.md`: guia científico com entradas e escopo, leitura por formato e o papel das ordens de ligação, os sete critérios com suas medidas, a referência interna relaxada (MMFF94/UFF), o teste estérico com raios de Bondi, tolerâncias e domínios, arquivos publicados e portas de saída, apresentação dos resultados, limites de interpretação, verificação e referências bibliográficas.
- Registro no gerador (`scripts/build_docs.py`), no catálogo do portal (`docs/assets/portal-catalog.json`) e no índice `docs/README.md`, com páginas HTML geradas nos dois idiomas.

Também não havia documentação do estilo **Sítio de ligação (superfície)** do visualizador 3D, nem do contrato de persistência do estado da execução e das decisões durante o processamento. Esses pontos foram documentados nos guias correspondentes.

## Divergências corrigidas

| Divergência | Correção |
| --- | --- |
| `README.md` anunciava 16 tipos de bloco e 15 operações científicas | 17 tipos de bloco e 16 operações; o ramo estrutural passa a citar a validação de poses |
| `current_reference.md` e sua versão inglesa declaravam 15 operações registradas e não traziam a tabela do bloco novo | Contagem corrigida e tabela `pose_validation` acrescentada com obrigatórios, opcionais e padrões das funções |
| `backend_usage.md` e sua versão inglesa declaravam 15 operações no CLI, com `import_results` como 16º bloco, e não listavam o bloco novo | Contagens corrigidas, linha na tabela de operações e seção `API pose_validation` |
| `artifact_naming.md` declarava 16 tipos de bloco e não citava o produtor público do bloco novo | 17 tipos, `PoseValidation` incluído entre produtores e tipos de conteúdo |
| `output_datasets.md` não identificava os conjuntos do bloco novo | Linha com relatório, poses aprovadas e receptores usados |
| A matriz de contratos de `pipeline_validation.md` omitia **Amber Score · DOCK6** e o bloco novo | Duas linhas acrescentadas, com entradas e saídas para outras etapas |
| Os manuais não tinham a tabela de conexões nem o exemplo do bloco novo, exigidos pelo validador | Linha na tabela de conexões, seção narrativa e exemplo em PT e EN |
| `index.html` não tinha o cartão da funcionalidade | Cartão acrescentado, com link para o exemplo do manual |
| `molecular_viewer.md` listava apenas fitas, bastões, esferas e linhas como estilo de receptor | Estilo de superfície do sítio documentado: escala Kyte–Doolittle, malhas sólida e translúcida, controle de opacidade, identificação do sítio por contexto, indisponibilidade explicada e limites |
| `result_indexing.md` atribuía a revisão de leitura apenas aos triggers da coluna autoritativa | Registro do anúncio direto feito pela execução quando só o progresso muda, e das leituras por `review_run`/`list_review_runs` |
| `architecture.md` e sua versão inglesa não descreviam a persistência do estado da execução | Nova seção: o que vai para a linha autoritativa, o que vai para a sobreposição de leitura, auditoria após o commit e journal WAL |
| `processing_feedback.md` afirmava que rascunhos automáticos não abrem janelas, sem o mecanismo que garante isso | Regra explicitada e nova seção sobre decisões do usuário durante o processamento |
| `projects.md` e sua versão inglesa descreviam a confirmação destrutiva sem dizer quem a pede e o efeito de descartá-la | Texto completado: confirmação pedida pelas ações explícitas, descarte equivale a cancelar, demais ações continuam respondendo |
| `scientific_documentation.md` não cobria a validação física de poses nem a escala da superfície do sítio | Seção de critérios implementados, nota sobre a escala de hidrofobicidade, quatro referências e dois caminhos de rastreabilidade |
| A tabela de scripts de `current_reference.md` não listava `scripts/check_ui_installation.py` | Script acrescentado, registrado como utilitário sem flags argparse |

## Descrições obsoletas

A varredura dos identificadores citados nos guias contra `operations.py`, `catalog.py` e o código fonte não encontrou operações, parâmetros ou caminhos de implementação inexistentes: nenhuma descrição precisou ser removida por obsolescência. Os identificadores verificados incluem os caminhos de rastreabilidade científica (`src/caad/admet.py`, `src/kernel/descriptors.py`, `consensus_model.py`, `footprints.py`, `docking_interactions.py`, `network_analysis.py`, `network_centrality.py`, `preparation_tools.py`), as opções de preparação e os campos de artefato citados nos textos.

Afirmações de contagem que ficaram desatualizadas pela entrada do bloco novo foram corrigidas, não removidas: elas descrevem a plataforma atual. Os relatórios datados e as auditorias anteriores foram mantidos sem edição, conforme a convenção do repositório de que evidência datada não é instrução atual.

## Arquivos alterados

| Arquivo | Alteração |
| --- | --- |
| `docs/pose_validation.md` | Novo guia científico do bloco de validação de poses |
| `docs/en/pose_validation.md` | Tradução integral do guia |
| `docs/README.md` | Guia novo nos índices PT/EN; resumo do visualizador atualizado |
| `docs/user_manual.md`, `docs/en/user_manual.md` | Tabela de conexões, seção narrativa, exemplo prático e guias de referência |
| `docs/frontend.md`, `docs/en/frontend.md` | Entradas, configuração e apresentação dos resultados do bloco novo |
| `docs/backend_usage.md`, `docs/en/backend_usage.md` | Tabela de operações, contagens e seção de API |
| `docs/current_reference.md`, `docs/en/current_reference.md` | Contagem, tabela de parâmetros, script novo e data de revisão |
| `docs/pipeline_validation.md`, `docs/en/pipeline_validation.md` | Matriz de contratos completada com dois blocos |
| `docs/molecular_viewer.md`, `docs/en/molecular_viewer.md` | Estilo de superfície do sítio de ligação |
| `docs/architecture.md`, `docs/en/architecture.md` | Seção de persistência do estado da execução |
| `docs/scientific_documentation.md`, `docs/en/scientific_documentation.md` | Critérios de validação de poses, escala de hidrofobicidade, referências e rastreabilidade |
| `docs/projects.md`, `docs/en/projects.md` | Contrato da confirmação destrutiva |
| `docs/result_indexing.md` | Revisão de leitura e projeções usadas pela interface |
| `docs/processing_feedback.md` | Decisões do usuário durante o processamento |
| `docs/output_datasets.md` | Conjuntos de saída do bloco novo |
| `docs/artifact_naming.md` | Contagem de blocos e produtor público novo |
| `README.md` | Contagem de blocos e ramo estrutural |
| `index.html` | Cartão da funcionalidade |
| `scripts/build_docs.py` | Registro do guia novo com títulos e resumos bilíngues |
| `docs/assets/portal-catalog.json` | Categoria, palavras-chave e ícone do guia novo |
| `docs/en/*.html`, `docs/pt/*.html`, `docs/index.html`, `index_doc.html` | Páginas e portais regerados pelo gerador |

## Complemento do mesmo dia: colunas da tabela e o trabalho de referência

O guia de validação de poses passou a registrar, coluna por coluna, o que cada valor mede e qual teste do trabalho de referência lhe corresponde.

- Nova seção 6 em `docs/pose_validation.md` e `docs/en/pose_validation.md`: os testes do trabalho nas três categorias que ele define, com os limiares publicados; as 32 colunas de `pose_validation.csv`; as 10 colunas da tabela exibida; e uma tabela de divergências declaradas entre a implementação e o trabalho. As seções seguintes foram renumeradas e a referência bibliográfica passou a trazer o DOI.
- As divergências declaradas são: referência de comprimentos e ângulos, formulação da energia interna, ausência do choque estérico interno e da sobreposição de volume, ausência da consistência de identidade por InChI, tratamento único de cofatores e águas, planaridade sem restrição de tamanho de anel e o RMSD de exatidão fora do escopo do bloco.
- A interface passou a explicar cada coluna da tabela de resultados por dica de contexto, assim como cada critério e cada medida do laudo por pose. As 27 frases novas receberam tradução nos catálogos; a auditoria de localização permanece nos 76 erros pré-existentes.
- `tests/test_pose_validation_flow.py` cobre a presença da explicação em todas as colunas, nos seis critérios e nas medidas do laudo.

Os limiares citados foram conferidos no texto integral do trabalho: tolerância de 25% sobre os limites de geometria de distâncias do RDKit para comprimentos e ângulos, 0,8 do limite inferior para o choque interno, 0,25 Å de planaridade, razão de energia 100 contra 50 confôrmeros ETKDGv3 relaxados com UFF, 0,75 da soma dos raios para a distância mínima intermolecular, 7,5% de sobreposição de volume com escalas 0,8 e 0,5, e RMSD de até 2 Å consciente de simetria para exatidão.

## Verificações e limites

- `scripts/build_docs.py`: 38 páginas de guia e três entradas de portal geradas.
- `scripts/validate_docs.py`: 44 HTML com arquivos, IDs e âncoras locais válidos; 16 operações científicas mais importação documentadas nos dois manuais e na apresentação. Antes desta auditoria, o mesmo validador acusava três ausências.
- Os três portais publicam 19 cartões pesquisáveis, agora coerentes com o registro do gerador.
- O complemento acima foi validado pelo mesmo gerador e validador documental, sem novos erros.
- Esta é uma auditoria documental. Não substitui testes científicos, execução de motores externos nem aceitação visual interativa, e não declara prontidão de produção.
