# Dependências e idioma de execução — 2026-10-26

## Inventário e manifestos

A varredura AST percorreu 268 arquivos Python do projeto, incluindo workflows,
scripts e testes. Ambientes virtuais, caches, diretórios de build e metadados Git
foram excluídos. Foram considerados imports locais, biblioteca padrão Python 3.12,
imports de terceiros e chamadas dinâmicas literais de importação.
O inventário com ocorrências está em `dependency_imports_2026-10-26.json`.

- Removido `chembl_webresource_client`: não há imports; o cliente ChEMBL atual usa
  `requests` e `urllib3`.
- Declarados Flet e sua implementação desktop, Flet Web, FastAPI, Uvicorn e
  websockets, utilizados pela aplicação e scripts de validação.
- Preservados os pacotes científicos importados e a infraestrutura Python, pip e
  setuptools. Open Babel e Poppler são usados por `obabel` e `pdftoppm`;
  ausência de import Python não os torna órfãos.
- Versões diretas adicionais foram fixadas a partir do ambiente inspecionado.
  Os pins científicos existentes foram mantidos. `nodefaults` limita a resolução
  ao canal conda-forge declarado.
- O import opcional `admet_ai` está declarado em `requirements-admet.yml` e no
  extra existente `admet-ml`. Os metadados de ADMET-AI 2.0.1 exigem RDKit >= 2025.9.5,
  incompatível com o pin 2023.09.6 do protocolo principal. O ambiente opcional
  usa PyTorch CPU; o instalador CPU existente continua disponível.

O auditor termina com erro para imports não declarados ou dependências diretas
órfãs. Na execução final: zero ausências e zero órfãos.

## Resolução e reprodutibilidade

A resolução Conda simulada de `requirements.yml` passou, sem criar ambiente.
A resolução pip simulada dos requisitos da interface e recuperação também passou.
O arquivo `requirements-linux-64.lock.yml` registra 219 pacotes Conda com versões,
builds e subdiretórios, mais 39 distribuições pip com versões e hashes dos wheels
resolvidos para Linux x86_64 / Python 3.12. O próprio lock passou na resolução
Conda simulada. Seu cabeçalho registra o SHA-256 do manifesto de origem.

O ambiente atual passou em `python -m pip check`. Não houve instalação, remoção
ou atualização de pacotes do ambiente utilizado pelo usuário. As resoluções e o
lock não constituem uma instalação limpa nem validação científica dos motores.
A parte opcional ADMET continua isolada; não foi criada nem executada nesta revisão.

## Mensagens da sessão

O português possui catálogo explícito `pt.json` e padrões próprios, em vez de
inversão do inglês durante a renderização. `runtime_pt.json` e `runtime_en.json`
compartilham chaves e contratos de interpolação. O acompanhamento mantém o idioma
capturado da sessão e obtém títulos, ações, estados, contadores, tempos e dicas
por chave. Mensagens de reutilização, seleção pendente e resumo da execução usam
o mesmo catálogo. A troca de idioma é bloqueada durante execuções ativas.

Os workers publicam `message_key` e `message_args`, mantendo `message` para leitores
anteriores. Mensagens históricas continuam sendo traduzidas por frase/padrão.
Os controles de progresso já traduzidos não sofrem uma segunda tradução que
possa alterar nomes escolhidos pelo usuário, arquivos e identificadores.
Metadados de acompanhamento inválidos não interrompem o cálculo científico, e
mensagens grandes respeitam o limite de leitura do arquivo de progresso.

O log exibido traduz eventos e mensagens internas reconhecidas na sessão leitora.
Saída nativa de motores e tracebacks são dados técnicos preservados. A nova ação
**Baixar log original** entrega exatamente os bytes armazenados, sem tradução.

## Validação

- 123 testes focados passaram: catálogos, idiomas simultâneos, polling repetido,
  estados do pipeline, acessibilidade, entradas históricas, download do log,
  seleção/reutilização, resultados, navegação, configuração, recuperação e Vina.
- No Chrome local, o popup real foi inspecionado em português (1200 px) e inglês
  (390 px). Cinco observações sucessivas por sessão conservaram o título e o idioma;
  não houve exceções JavaScript. As capturas foram inspecionadas visualmente.
- A suíte ampliada anterior executou 100 testes, com uma falha em
  `CandidatePreparationTests.test_exported_candidates_are_reused_without_conversion`:
  o teste espera rejeição de um formato ausente, mas o código permite continuar.
  O mesmo caso falhou ao carregar `biomolexplorer.docking_data` diretamente do HEAD.
  Esse contrato científico preexistente não foi alterado por esta revisão.

Para repetir a auditoria de imports:

```bash
python scripts/audit_dependencies.py --output /tmp/biomol-imports.json
python scripts/build_docs.py
python scripts/validate_docs.py
```
