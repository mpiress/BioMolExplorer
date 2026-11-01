# Revisão exclusivamente documental — 31/10/2026

## Escopo e fonte primária

Foram lidos o registro `operations.py`, o catálogo, os contratos tipados em `flow.py`, os serviços de ordem/execução, persistência/artefatos, os módulos e validadores Amber Score, os caminhos ativos de preparação, manifests e argumentos CLI. A auditoria não modifica a aplicação nem executa cálculos científicos.

O workspace já continha alterações de código/testes de tarefas anteriores. A restrição documental foi verificada pela comparação SHA-256 dos arquivos protegidos antes/depois desta tarefa, e não pela exigência de um `git status` limpo. O arquivo não versionado `index (Cópia).html` já existia e não foi usado nem alterado.

## Atualizações

- Site institucional: CSS/JavaScript locais, sem Tailwind/fontes/ícones carregados por CDN, navegação com teclado, menu móvel, zoom de figuras, suporte a movimento reduzido e impressão. As âncoras principais e URLs dos guias foram preservadas.
- Catálogo público: 16 tipos de bloco, incluindo importação e Amber Score nativo; distinção entre 15 operações CLI e 16 blocos do workspace. Links de cada capacidade apontam para exemplos da implementação atual.
- Fluxograma: `imgs/BioMolExplorer.png` foi redesenhado como figura científica; `imgs/BioMolExplorer.svg` fornece uma representação vetorial acessível. Ramos são ilustrativos, não uma cadeia obrigatória. Redocking fornece receptores aos consumidores; Amber Score requer estruturas de receptor e poses de ligante.
- README: apresentação concisa, requisitos por protocolo, exemplos CLI, execução pelo checkout, persistência, documentação e licença. O DOI do dataset associado não é apresentado como DOI do software.
- Guias PT/EN: instalação, preparação, referência de todos os parâmetros registrados, Amber Score (requisitos/defaults/unidades/limitações), seleção tardia, prioridade topológica, pausa/cancelamento, operações de limpeza e nomenclatura com IDs/linhagem.
- Documentos históricos de validação permanecem evidências datadas; não foram convertidos em promessas de capacidades atuais.

## Divergências da aplicação registradas, sem correção

1. `pyproject.toml:10` contém `dependencies = [dock6, chimera, miniconda]`: a leitura com `tomllib` produz `TOMLDecodeError` (linha 10, coluna 17). Instruções que presumiam instalação editável funcional foram substituídas por execução a partir do checkout, para um ambiente existente. Não foi alterado o empacotamento.
2. `requirements.yml` inclui ADMET-AI 2 no ambiente principal enquanto o guia/instalador opcional prevê isolamento pela compatibilidade RDKit. A resolução completa dos pins não foi testada; os manifests foram mantidos.
3. O gerador documental existente ainda contém, em seus metadados internos, um resumo de preparação citando backends alternativos. O script foi preservado por restrição de escopo. Os dois portais gerados receberam somente a correção documental desse resumo; o procedimento e o risco de regeneração estão registrados em `docs/README.md`.
4. Os sete casos de regressão preexistentes do relatório de nomenclatura continuam como pendências registradas naquela tarefa. Esta revisão documental não os corrigiu nem reexecutou e não certifica uma release científica.

## Verificações e evidências

As evidências finais deste diretório registram cobertura do catálogo/parâmetros, links locais/externos, ajuda CLI, navegador, responsividade e integridade dos arquivos protegidos. Não houve instalação de pacotes, criação de novo ambiente Conda, lançamento de execução científica ou alteração de dados de pesquisa. Verificação de link externo demonstra somente acessibilidade no momento da consulta; não valida o conteúdo científico da página.

### Resultado das verificações

- `python scripts/build_docs.py`: 30 páginas de guias e dois portais regenerados; o resumo legado dos portais foi ajustado documentalmente após a geração.
- `python scripts/validate_docs.py`: **34 páginas HTML**, caminhos, âncoras e IDs válidos; 15 operações científicas e o bloco de importação cobertos exatamente uma vez pelos exemplos dos dois manuais e cartões do site.
- Análise HTML estrita: nenhuma tag de fechamento inesperada ou aberta ao fim. Blocos de código Markdown equilibrados nos documentos correntes.
- Chrome real: **cinco larguras** (320, 390, 768, 1024 e 1440 px) para o site, mais **16 visualizações** dos guias PT/EN. Sem overflow horizontal da página e sem exceções JavaScript. Menu móvel, âncoras, ampliação da figura e fechamento por Escape exercitados; o conteúdo principal permanece disponível sem JavaScript. Capturas foram inspecionadas visualmente.
- Referências externas finais: **52 verificadas, 46 acessíveis, seis não confirmadas por timeout** (páginas DOCK e wiki ZINC). Timeouts não foram classificados como links definitivamente inexistentes nem como verificação aprovada. `external-links.json` contém os destinos e resultados.
- Dois destinos que retornavam 404 foram substituídos: instalação Chimera passou à página oficial de download; catálogo Gephi passou ao caminho atual da documentação Desktop. O DOI de dataset informado no projeto também retornou 404: o identificador e a referência foram preservados como texto com ressalva, sem outro DOI inventado. As consultas não validam a bibliografia completa.
- `--help` do CLI científico e do workspace executados no ambiente existente com `PYTHONPATH=src`; logs isolados em `/tmp`. Isso verifica os comandos e argumentos documentados, não a instalação de um ambiente limpo.
- **506 arquivos protegidos** (código, testes, scripts e arquivos operacionais) comparados por SHA-256: **zero alterações e zero novos arquivos nessas áreas**. As alterações de implementação já presentes no início da tarefa permaneceram byte a byte idênticas.
- `git diff --check`: sem erros de espaços. Nenhum pacote instalado; nenhum teste, algoritmo, dependência ou configuração de execução modificado. Não houve teste científico nem execução do pipeline.

### Limitações da publicação e da manutenção

Não foi publicado/deployado um site externo: a validação ocorreu na cópia estática local. O domínio e o repositório públicos não foram modificados. O documento não garante disponibilidade futura das referências externas ou compatibilidade de plataformas não testadas. A checagem de Markdown cobre a geração do portal e a inspeção visual das páginas geradas; não certifica todos os renderizadores de terceiros.

O fluxograma atual é ilustrativo e não substitui as regras de conexão/seleção. Versões antigas em arquivos de imagem históricos não foram apagadas. A aparência da aplicação não foi modificada; as capturas de workspace existentes no site são identificadas como ilustrações.

### Reprodução das verificações documentais

```bash
python scripts/build_docs.py
# Ajuste o resumo dos portais conforme docs/README.md; o gerador não foi alterado.
python scripts/validate_docs.py
PYTHONPATH=src python -m biomolexplorer --help
PYTHONPATH=src python -m biomolexplorer.ui.app --help
```

Os dados do navegador, inventário extraído, saídas de ajuda e resultados de links estão neste diretório. Nenhum desses resultados deve ser interpretado como aprovação científica da aplicação.

## Documentos e recursos atualizados

- `README.md`
- `docs/README.md`
- `docs/architecture.md`
- `docs/backend_usage.md`
- `docs/chimera_migration.md`
- `docs/en/architecture.html`
- `docs/en/architecture.md`
- `docs/en/backend_usage.html`
- `docs/en/backend_usage.md`
- `docs/en/chimera_migration.html`
- `docs/en/chimera_migration.md`
- `docs/en/frontend.html`
- `docs/en/frontend.md`
- `docs/en/installation.html`
- `docs/en/installation.md`
- `docs/en/logging.html`
- `docs/en/logging.md`
- `docs/en/pipeline_validation.html`
- `docs/en/pipeline_validation.md`
- `docs/en/projects.html`
- `docs/en/projects.md`
- `docs/en/user_manual.html`
- `docs/en/user_manual.md`
- `docs/frontend.md`
- `docs/graph_explorer.md`
- `docs/index.html`
- `docs/installation.md`
- `docs/logging.md`
- `docs/output_datasets.md`
- `docs/pipeline_validation.md`
- `docs/projects.md`
- `docs/pt/architecture.html`
- `docs/pt/backend_usage.html`
- `docs/pt/chimera_migration.html`
- `docs/pt/frontend.html`
- `docs/pt/index.html`
- `docs/pt/installation.html`
- `docs/pt/logging.html`
- `docs/pt/pipeline_validation.html`
- `docs/pt/projects.html`
- `docs/pt/user_manual.html`
- `docs/user_manual.md`
- `imgs/BioMolExplorer.png`
- `index.html`
- `docs/assets/site.css` e `docs/assets/site.js` (novos, exclusivos da vitrine).
- `imgs/BioMolExplorer.svg` (nova representação vetorial).
- Este relatório e suas evidências documentais.
