> Registro histórico. As instruções de ambiente/protocolos desta validação foram substituídas pela [revisão de Amber Score nativo de 2026-10-29](native_amber_score_2026-10-29/README.md). Não representam as opções atuais da aplicação.

# Retorno à preparação exclusiva com UCSF Chimera

Validação local em 26/10/2026. Docking e redocking normalizam configurações
antigas `auto`/`python` para `legacy_chimera`; a interface oferece apenas Chimera.
AmberTools continua disponível para os fluxos de dinâmica molecular.

## Correções verificadas

- Preparação, pré-validação, configuração, reaproveitamento de estruturas e
  preparação de candidatos usam Chimera. Resultados históricos identificados
  como nativos não são adotados como preparação Chimera.
- A minimização usa `prep false`: o DockPrep implícito do comando `minimize`
  substituía cargas `gas` por AM1-BCC. A protonação e as cargas explicitamente
  configuradas agora são preservadas, inclusive na exportação para PDBQT.
- A validação de estruturas preparadas consulta o diretório de saída, mesmo
  quando ele é diferente do diretório dos PDBs de entrada.
- A consulta periódica de atualizações é uma tarefa de fundo. Três consultas
  automáticas com um diálogo de configuração aberto não criaram splash.

## Execuções científicas reais

Ferramentas: UCSF Chimera 1.17.3, Open Babel e Vina do ambiente
`BioMolExplorer`, DOCK6 instalado em `/home/michel/progs/dock6`.
Os ensaios usaram cópias temporárias; os dados originais não foram alterados.

- Etanol, preparação com minimização e exportação MOL2/PDBQT: `gas` e `am1`
  concluíram. Carga do oxigênio: -0.3967/-0.397 para gas e -0.6028/-0.603
  para am1, respectivamente; a diferença corresponde à precisão dos formatos.
- Redocking Vina de 4M0E/1YL, resíduo 604, cadeia A: preparação a partir do
  PDB, receptor sem minimização, ligante com minimização, caixa 24 Å,
  exhaustiveness 1, um modo. Ambos concluíram sem falhas; RMSD gas
  5.8665 Å e am1 0.1988 Å. Esses parâmetros validam execução e artefatos,
  não constituem um benchmark de acurácia nem aceitação científica geral.
- Docking Vina de etanol em receptor sintético: ambos os métodos concluíram
  e produziram poses (scores -1.260 e -1.261).
- Docking rígido DOCK6 de etanol no receptor 4M0E preparado com Chimera:
  grade, docking, pose e footprint concluíram para gas e am1, com scores
  -17.92396 e -18.78022. Os scores entre motores não são comparáveis.
- A preparação de 4M0E em pasta separada produziu centros e metadados do
  par sem relatório de falha.

Evidências temporárias: `/tmp/biomol-chimera-real-charge-fixed.log`,
`/tmp/biomol-chimera-redock-target.log`, `/tmp/biomol-chimera-vina-final.log`,
`/tmp/biomol-chimera-dock6-final.log` e
`/tmp/biomol-chimera-prepared-folder-fixed.log`.

## Interface e regressão

A interface de configuração foi renderizada em larguras 560, 1024 e 1440 px,
incluindo a tabela de pares e as configurações gerais. Seis capturas foram
inspecionadas sem erros do navegador. A consulta automática sem splash foi
validada por teste determinístico, sem aguardar o temporizador real.

Na execução abrangente de 146 testes, 145 passaram e um fixture do teste novo
não criava o CSV de entrada exigido pela limpeza. Após corrigir esse fixture,
os 16 testes de preparação Chimera e isolamento de falhas foram repetidos e
passaram. A suíte abrangente não foi repetida após essa alteração exclusiva
do fixture. Logs: `/tmp/biomol-chimera-complete-tests-corrected.log` e
`/tmp/biomol-chimera-subset-corrected.log`.

A limpeza foi também verificada com uma pose antiga de um composto cuja
conversão falha: ela é removida e os demais compostos continuam disponíveis.

A documentação HTML foi regenerada; 34 páginas passaram na validação de
arquivos locais, âncoras e IDs. A validação é amostral: não executa todos os
pares dos projetos do usuário ou todos os modos flexíveis do DOCK6.
