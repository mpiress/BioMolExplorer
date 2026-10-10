# Documentação / Documentation

[Apresentação da plataforma](../index.html) · [Portal no navegador / Browser portal](index.html)

## Português

| Guia | Conteúdo | Ler no navegador |
| --- | --- | --- |
| [Instalação e configuração](installation.md) | Download do GitHub, Chimera 1.17, DOCK6 6.11 e verificações antes de iniciar | [HTML](pt/installation.html) |
| [Manual do usuário](user_manual.md) | Passo a passo completo: idioma, projeto, arquivos, operações, resultados e migração | [HTML](pt/user_manual.html) |
| [Workspace Flet](frontend.md) | Instalação da interface, contas, compartilhamento, projetos, pipelines e arquivos próprios | [HTML](pt/frontend.html) |
| [Uso do backend](backend_usage.md) | Ambiente científico, CLI, serviços de tarefas, parâmetros e operações | [HTML](pt/backend_usage.html) |
| [Arquitetura e revisão técnica](architecture.md) | Organização do código, contratos, correções, execução, verificações e limites | [HTML](pt/architecture.html) |
| [Projetos e versões](projects.md) | Pastas, migração, colaboração, histórico, rollback e múltiplas entradas | [HTML](pt/projects.html) |
| [Configuração do redocking](redocking_configuration.md) | Pares, cadeia, cofatores, preparação, RMSD, downloads e visualização 3D | [HTML](pt/redocking_configuration.html) |
| [Recuperação flexível da informação](retrieval.md) | Critérios PDB, modos ChEMBL, filtros e rastreabilidade | [HTML](pt/retrieval.html) |
| [Visualizador molecular e interações 3D](molecular_viewer.md) | Poses Vina/DOCK6, receptor e referência, estilos, hidrogênios, cores e identificação pelo mouse | [HTML](pt/molecular_viewer.html) |
| [Logs e diagnóstico](logging.md) | Contexto, causas, códigos de erro, resumos e filtros | [HTML](pt/logging.html) |
| [Validação do pipeline](pipeline_validation.md) | Matriz de conexões, regressões e limites da verificação | [HTML](pt/pipeline_validation.html) |

[Avaliação da substituição do Chimera](chimera_migration.md) · [HTML](pt/chimera_migration.html): inventário, alternativas e decisão de preservar o protocolo.

[Porte nativo DMS e validação](dms_migration.md) · [HTML](pt/dms_migration.html): superfície SES em Python, contrato sphgen e referências do algoritmo C.

## English

| Guide | Contents | Read in your browser |
| --- | --- | --- |
| [Installation and configuration](en/installation.md) | GitHub download, Chimera 1.17, DOCK6 6.11 and checks before starting | [HTML](en/installation.html) |
| [User manual](en/user_manual.md) | Complete walkthrough: language, projects, files, operations, results and migration | [HTML](en/user_manual.html) |
| [Flet workspace](en/frontend.md) | UI installation, accounts, sharing, projects, pipelines and user-supplied files | [HTML](en/frontend.html) |
| [Backend usage](en/backend_usage.md) | Scientific environment, CLI, job services, parameters and operations | [HTML](en/backend_usage.html) |
| [Architecture and technical review](en/architecture.md) | Code organization, contracts, fixes, execution, verification and limitations | [HTML](en/architecture.html) |
| [Projects and versions](en/projects.md) | Folders, migration, collaboration, history, rollback and multiple inputs | [HTML](en/projects.html) |
| [Redocking configuration](en/redocking_configuration.md) | Pairs, chain, cofactors, preparation, RMSD, downloads and 3D viewing | [HTML](en/redocking_configuration.html) |
| [Flexible information retrieval](en/retrieval.md) | PDB criteria, ChEMBL modes, filters and traceability | [HTML](en/retrieval.html) |
| [Molecular viewer and 3D interactions](en/molecular_viewer.md) | Vina/DOCK6 poses, receptor and reference, styles, hydrogens, colors and hover identification | [HTML](en/molecular_viewer.html) |
| [Logs and diagnostics](en/logging.md) | Context, causes, error codes, summaries and filters | [HTML](en/logging.html) |
| [Pipeline validation](en/pipeline_validation.md) | Connection matrix, regressions and verification limits | [HTML](en/pipeline_validation.html) |

[Chimera replacement assessment](en/chimera_migration.md) · [HTML](en/chimera_migration.html): inventory, alternatives and decision to preserve the protocol.

[Native DMS port and validation](en/dms_migration.md) · [HTML](en/dms_migration.html): native Python SES surfaces, sphgen contract and C reference comparisons.

## Manutenção / Maintenance

Os documentos Markdown são a fonte. Após editar uma versão, atualize sua tradução
e gere novamente as páginas HTML na raiz do projeto:

The Markdown documents are the source. After editing a guide, update its
translation and regenerate the HTML pages from the project root:

```bash
python scripts/build_docs.py
```

O gerador usa apenas a biblioteca padrão Python. As páginas geradas são estáticas,
não dependem de serviços externos e podem ser abertas diretamente pelo navegador
ou publicadas junto com `index.html`.

The generator uses only the Python standard library. Generated pages are static,
require no external services and can be opened directly in your browser or
published alongside `index.html`.
