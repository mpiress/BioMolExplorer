# Auditoria integrada da documentação — 15/11/2026

## Resultado

Esta revisão parte do inventário fonte-documentação de 13/11/2026 e reabre os contratos científicos e do portal solicitados. O código científico, o registro de operações, a CLI, a inicialização Flet, autenticação/sessão, artefatos, templates DOCK 6, visualizadores e documentação ativa foram confrontados novamente. Nenhum módulo científico foi alterado.

## Divergências corrigidas

- Não havia uma referência científica central, bilíngue e rastreável. Foram criados `docs/scientific_documentation.md` e `docs/en/scientific_documentation.md`, publicados também em HTML.
- Morgan, MACCS e Gobbi Pharm2D não estavam reunidos com seus parâmetros e limitações exatos. A referência agora registra `radius=2`, 2048 bits, quiralidade/features, `useBondTypes=False`, 166 chaves públicas e a codificação topológica RDKit.
- O consenso foi fixado como `target-balanced-ecr-v1`, com quantis de empate, `tau=0.1`, cobertura, penalização de ausências e `repulsion_weight` sem efeito. Z-score/min-max não são atribuídos ao caminho ativo.
- Footprint DOCK 6 foi separado de detecção geométrica: a plataforma lê energias por resíduo do motor, mas não inventa uma série independente de H-bond.
- Os thresholds 3D foram documentados a partir de `docking_interactions.py`; ponte salina, cátion-pi, halogênio e coordenação metálica foram marcados como não implementados.
- PCA, t-SNE, UMAP, AUC-ROC e Enrichment Factor foram auditados e registrados como ausentes do código ativo.
- O portal linear misturava engenharia e ciência. `index_doc.html`, `docs/index.html` e `docs/pt/index.html` agora são gerados por uma única fonte com dois hubs, filtros independentes e busca global agrupada.

## Estrutura publicada

O Hub Técnico cobre primeiros passos, workspace/fluxos e arquitetura/CLI. O Hub Científico cobre propriedades/ADMET, fingerprints/similaridade, docking/consenso e energia/interações 3D. Durante busca, cada contador é atualizado separadamente, hubs sem resultados são ocultados e o estado vazio aparece apenas quando ambos zeram.

Os guias de consenso e AMBER, antes complementares, passaram ao registro publicado do gerador. O validador documental exige que cada guia apareça exatamente uma vez no conjunto dos dois hubs, sem impor a ordem linear antiga.

## Verificações e limites

- `scripts/build_docs.py`: 36 páginas de guia e três entradas de portal geradas.
- `scripts/validate_docs.py`: 42 HTML com arquivos, IDs e âncoras locais válidos; 15 operações mais importação documentadas.
- JSON dos catálogos/traduções, compilação Python e `git diff --check`: válidos.
- Testes científicos, motores externos e aceitação visual interativa não são substituídos por esta auditoria documental. Não há declaração de prontidão de produção.
