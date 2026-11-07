# Relatório e localização — 7 de novembro de 2026

## Diagnóstico e correções

A restauração anterior reativou o gerador completo existente. O clique carregava
até 50 execuções, 250 registros de artefatos, tabelas, imagens e análises antes de
mostrar uma prévia que usava apenas uma pequena parte desses dados. O orçamento
de oito gráficos era aplicado após a renderização; gráficos adicionais eram
produzidos e descartados.

A prévia agora consulta o projeto e duas contagens no banco, sem abrir resultados,
desserializar execuções ou desenhar figuras. A geração completa acontece no
primeiro download/impressão e é reutilizada pelos outros formatos naquele diálogo.
A prévia informa esse comportamento. Não há cache global que possa misturar
projetos; troca de projeto/sessão invalida os callbacks. Uma nova abertura consulta
o estado atual. Limites de tabelas, linhas e figuras, métricas e contratos de
artefatos foram preservados. O orçamento impede a renderização de gráficos que
seriam descartados, sem descartar suas métricas.

Os splashes tinham muitas mensagens novas fora dos dicionários. Overlays não
participavam das atualizações da localização; tooltips estruturados e alguns
outros campos de apresentação também não eram percorridos. Mensagens compostas
traduziam o prefixo e deixavam descrições internas em português. Foram completados
dicionários de mensagens, rótulos, erros e descrições da interface e dos serviços,
e corrigidos esses caminhos de apresentação. Os textos do visualizador 3D, cenas
de docking e métodos dos relatórios passaram a dicionários. O seletor e o argumento
`--language` leem `locales.json`; o teste de um terceiro idioma comprova a extensão
sem alterar a interface. Veja [LOCALIZATION.md](../../LOCALIZATION.md).

## Medições

Cenário sintético: quatro blocos ADMET, 32 CSVs de 1.000 linhas, em diretório
temporário. Nenhum projeto do usuário foi usado.

- Antes, o caminho completo executado na abertura levou 5,47 s sob cProfile.
- Após limitar os desenhos efetivos, a geração completa levou 2,06 s sob cProfile.
- Benchmark sem profiler, registrado em `benchmark.json`: geração completa 1,21 s;
  consulta da prévia com mediana de 0,81 ms e máximo de 1,56 ms em dez repetições.
- Foram mantidos oito gráficos e métricas para os 32 conjuntos.
- No Chrome, do clique até a prévia disponível na árvore de acessibilidade:
  aproximadamente 0,11–0,25 s nas quatro larguras verificadas. Essa medida inclui
  comunicação com o cliente, mas não a duração total da animação de abertura.

Reprodução do benchmark, no ambiente da aplicação:

```bash
PYTHONPATH=src MPLCONFIGDIR=/tmp/biomol-mpl python scripts/benchmark_report.py
PYTHONPATH=src python scripts/audit_localization.py
```

## Validação e limites

- 148 testes integrados passaram: relatórios, desempenho, idiomas, mensagens dos
  workers, artefatos semânticos, responsividade, navegação, visualizadores 3D,
  grafos, layout, atualizações, feedback de processamento e login persistido.
  O registro está em `tests.log`; avisos de moléculas inválidas são esperados nos
  casos de teste. Após os últimos ajustes de textos/eixo dos gráficos, as suítes
  focadas de relatório e localização foram executadas novamente.
- A auditoria estática percorreu 2.006 expressões de texto reconhecíveis no módulo
  BioMolExplorer, sem pendências nos dicionários ou parâmetros de interpolação.
  As contagens dos seis grupos de dicionários estão em `catalog-audit.json`.
- Todos os títulos e descrições cadastrados de processamento foram verificados
  em inglês, incluindo reutilização de diálogo/overlay e retorno ao português.
  Também foram testados tooltips, isolamento entre sessões, dados protegidos,
  idioma adicional por arquivos e reordenação de argumentos de templates.
- A inspeção no Chrome usou controles reais em uma página de auditoria com dados
  sintéticos, em 1440×1000, 1024×800, 560×760 e 360×760. Verificou acesso ao botão,
  prévia, splash modal em inglês e ausência de erros JavaScript. As capturas de
  360 px foram preservadas. O caminho de overlay tem teste de controles; sua
  renderização no navegador não foi confirmada nesta auditoria.
- A construção do wheel confirmou os 13 arquivos de localização, incluindo os
  namespaces novos. Não houve mudança de dependências do BioMolExplorer.
- Dois simuladores de página de testes foram atualizados para suportar os diálogos
  já usados pelo processamento; uma expectativa antiga de quatro categorias foi
  corrigida para cinco. A categoria de dinâmica molecular foi preservada.

Os tempos dependem de hardware, quantidade de arquivos e tamanho do projeto. A
análise estática usa reconhecimento lexical e não prova cobertura de toda frase
possível; novos textos compostos ainda precisam de revisão/teste. Logs nativos,
dados moleculares, nomes e valores do usuário permanecem sem reescrita. Esta
validação não constitui uma nova aceitação científica de pipelines reais.
