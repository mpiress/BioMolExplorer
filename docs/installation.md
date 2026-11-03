# Instalação e configuração do BioMolExplorer

[Documentação](README.md) · Português · [English](en/installation.md)

Siga esta ordem: baixar o projeto → instalar as ferramentas externas → preparar o ambiente Conda → configurar e verificar os executáveis → iniciar o BioMolExplorer. Os exemplos de terminal deste guia usam Linux e Bash, ambiente de execução do backend atual. Instale previamente Git e Anaconda ou Miniconda.

**Antes de iniciar a aplicação, instale e configure UCSF Chimera 1.17 e DOCK6 6.11 no computador que executará os cálculos.** Instalar o pacote Python ou a interface não instala essas duas ferramentas. A interface pode abrir mesmo sem elas, mas etapas de preparação, redocking e docking que dependem delas falharão durante a execução.

## 1. Baixar o projeto do GitHub

Abra o [repositório oficial do BioMolExplorer](https://github.com/mpiress/BioMolExplorer). Com Git instalado, escolha a pasta onde deseja manter o código e execute:

```bash
git clone https://github.com/mpiress/BioMolExplorer.git
cd BioMolExplorer
```

Se preferir baixar sem Git, no repositório clique em **Code → Download ZIP**, extraia o arquivo e abra um terminal dentro da pasta extraída, normalmente `BioMolExplorer-master`. Execute os próximos comandos nessa pasta, que deve conter `requirements.yml` e `pyproject.toml`; não use a pasta `src` nem execute comandos dentro do ZIP.

A pasta do código é diferente da pasta de cada projeto de pesquisa. Você escolherá a pasta que recebe entradas e resultados ao criar um projeto na interface.

## 2. Instalar as ferramentas externas obrigatórias

| Ferramenta | Versão exigida neste guia | Uso no BioMolExplorer | Obtenção e instalação |
| --- | --- | --- | --- |
| UCSF Chimera | 1.17 | Preparação ativa de receptores e ligantes | [Downloads de versões anteriores da UCSF](https://www.cgl.ucsf.edu/chimera/olddownload.html) e [instruções de instalação](https://www.cgl.ucsf.edu/chimera/download.html) |
| DOCK6 | 6.11 | Docking, refinamento e cálculo de scores, com seus programas auxiliares | [Página oficial DOCK 6](https://dock.compbio.ucsf.edu/DOCK_6/index.htm), [notas da versão 6.11](https://dock.compbio.ucsf.edu/DOCK_6/new_in_6.11.txt) e manual incluído na distribuição |

Instale UCSF Chimera 1.17, Open Babel e Vina para a preparação/redocking; instale DOCK6 6.11 para seus blocos de docking. Não é necessário ambiente AmberTools separado. ChimeraX não interpreta os scripts Chimera.

Instale o Chimera com o instalador adequado ao sistema. Para DOCK6, siga as instruções da distribuição 6.11, compile também os programas auxiliares e mantenha a instalação completa, incluindo `bin/` e `parameters/`. Execute os testes fornecidos pelo DOCK6 conforme seu manual.

A superfície molecular é calculada nativamente em Python com NumPy e SciPy. Não é necessário instalar DMS; consulte o [porte e sua validação](dms_migration.md).

As ferramentas exigidas pelo protocolo devem estar disponíveis no computador dos workers, inclusive em implantação web. Chimera é exigido pelos protocolos ativos de preparação; recuperação, descritores e grafos não o exigem.

## 3. Criar o ambiente científico e instalar a interface

Na raiz do código baixado:

```bash
conda env create -f requirements.yml
conda activate BioMolExplorer
python -m pip install -e '.[ui]'
python -m pip check
command -v biomolexplorer-ui
biomolexplorer-ui --help
```

Execute uma linha por vez e aguarde terminar sem erro antes de continuar. O comando `pip install` é obrigatório: instalar as bibliotecas pelo Conda não registra os comandos do BioMolExplorer. A opção `-e` mantém a instalação ligada à pasta do código; não a remova nem a mova depois. `.[ui]` inclui a interface web e desktop. O comando `--help` deve mostrar as opções sem abrir a aplicação.

Se o ambiente já existe, atualize-o com `conda env update -f requirements.yml --prune` e depois ative-o. Repita `python -m pip install -e '.[ui]'` e as verificações acima. O ambiente inclui Open Babel, Vina e a interface; ele não substitui Chimera e DOCK6, que continuam instalações externas. O arquivo mantém faixas compatíveis — em especial NumPy 1.x e pandas 2.x — para impedir que extras de aprendizado de máquina alterem silenciosamente o worker científico. O pacote requer Python 3.12.

## 4. Configurar caminhos e conferir a instalação

Com o ambiente Conda ativo, adicione ao `PATH` os diretórios que contêm os executáveis. Substitua todos os caminhos de exemplo pelos caminhos reais da sua instalação:

```bash
export BIOMOL_DOCK6_ROOT="/caminho/para/dock6-6.11"
export PATH="/caminho/para/chimera-1.17/bin:$BIOMOL_DOCK6_ROOT/bin:$PATH"
```

`BIOMOL_DOCK6_ROOT` é uma variável de conveniência deste guia. A aplicação recebe esse caminho pela opção `--dock6-path`, mostrada abaixo. Informe a raiz do DOCK6, que contém `bin/` e `parameters/`, e não apenas o executável ou a pasta `bin/`. O mesmo terminal deve iniciar a aplicação para que os workers herdem o `PATH`. Para manter a configuração entre sessões, inclua os exports com os caminhos reais no arquivo de inicialização do seu shell e abra um novo terminal.

Confira os executáveis antes de iniciar:

```bash
for ferramenta in chimera dock6 sphgen showbox grid obabel vina; do
    if command -v "$ferramenta"; then
        printf 'OK: %s disponível\n' "$ferramenta"
    else
        printf 'FALTA: %s — instale ou corrija o PATH antes de iniciar\n' "$ferramenta"
    fi
done
if test -d "$BIOMOL_DOCK6_ROOT/parameters"; then
    echo 'OK: diretório parameters do DOCK6 encontrado'
else
    echo 'FALTA: diretório parameters — confira a raiz da instalação DOCK6'
fi
```

**Não execute um protocolo enquanto faltar alguma ferramenta ou pasta exigida por ele.** Confirme os caminhos das ferramentas realmente usadas; Chimera 1.17 é exigido para preparação e DOCK6 6.11 apenas por seus protocolos, especialmente se houver várias versões instaladas. A localização de um executável confirma sua presença no `PATH`, mas não valida sua versão nem o funcionamento do protocolo. Confira as versões nas ferramentas/distribuições e faça uma execução com um caso de referência antes de iniciar seu estudo.

## 5. Iniciar a aplicação configurada

Ainda no mesmo terminal, com o ambiente ativo e as verificações concluídas:

```bash
biomolexplorer-ui --web --language pt --dock6-path "$BIOMOL_DOCK6_ROOT"
```

Abra `http://127.0.0.1:8550`. Para desktop, omita `--web`. Para servir sem abrir um navegador automaticamente, use `--web --no-browser`. Se selecionar outro interpretador com `--worker-python`, prepare também esse ambiente científico e garanta que os workers tenham acesso às ferramentas externas.

Continue no [manual do usuário](user_manual.md) para criar uma conta, selecionar a pasta de um projeto e configurar o pipeline. Para execução por CLI e parâmetros técnicos, consulte [uso do backend](backend_usage.md).

## Verificações antes da execução

O redocking com preparação verifica Chimera, Open Babel e Vina. Com entradas já preparadas, somente Vina é exigido. Configurações antigas automáticas/Python passam a utilizar Chimera.

**Preparar para docking** usa Chimera para receptores, ligantes, cargas `gas`/`am1` e minimização. Open Babel converte os formatos de Vina e DOCK6. Não é necessário configurar um ambiente AmberTools nesses blocos. Consulte [preparação molecular](molecular_preparation.md).

**Amber Score · DOCK6** requer instalação completa do DOCK6 **6.11**, incluindo seus helpers AMBER e biblioteca de parâmetros. Não usa ambiente AmberTools separado, sander ou cpptraj. Consulte a seção Amber Score nativo no [manual](user_manual.md).

Consulte [Configuração do redocking](redocking_configuration.md) para seleção de pares, cofatores, entradas preparadas, diagnóstico e resultados.

## Preditor ADMET-AI opcional

O backend **RDKit local** funciona com o ambiente científico principal. Para escolher **ADMET-AI**, instale também o preditor no computador do worker. Ele permanece deliberadamente em um ambiente isolado para não substituir pandas, RDKit, NumPy ou outras dependências do docking.

Na raiz do repositório, usando Python 3.12 compatível:

```bash
python scripts/setup_admet_ai.py
```

O script cria `.admet-env`, instala PyTorch CPU e `admet-ai==2.0.1`, realiza uma inferência e grava `.ready` após sucesso. A instalação requer acesso às distribuições e modelos. O worker descobre esse ambiente local automaticamente. Para instalar em outra pasta:

```bash
python scripts/setup_admet_ai.py --directory /caminho/para/admet-env
export BIOMOL_ADMET_PYTHON="/caminho/para/admet-env/bin/python"
```

Inicie a aplicação no terminal que contém esse export, para o worker herdar a configuração. A variável identifica o Python do preditor, independentemente de `--worker-python`. Na migração para outro computador, recrie o ambiente e atualize caminhos; exportar um projeto não instala bibliotecas científicas.

Se aparecer **ADMET-AI is not installed in the worker environment**, confira o estado mostrado na configuração ADMET, o caminho, a instalação e o acesso do worker aos modelos. Reinicie a aplicação após configurar o ambiente. Como alternativa, selecione explicitamente **RDKit local**; isso fornece descritores, estimativas e regras locais, sem os endpoints `ML_`. Não há substituição silenciosa de modelos. Use o ambiente isolado descrito acima.

Consulte [ADMET por hipótese](admet_research.md) para unidades, perfis, alertas, arquivos e limites dos modelos.


## Reabrir a aplicação e resolver erros comuns

Em um novo terminal, execute `conda activate BioMolExplorer`, reaplique os caminhos da etapa 4 se não os salvou no shell e execute o comando da etapa 5. Não é necessário reinstalar a cada abertura.

| Mensagem ou situação | Como resolver |
| --- | --- |
| `conda: command not found` | Instale Miniconda/Anaconda e abra um novo terminal Bash. Se necessário, inicialize Bash conforme as instruções do instalador. |
| `CondaError: Run conda init before conda activate` | Execute `conda init bash`, feche o terminal e abra outro. |
| `biomolexplorer-ui: command not found` | Ative `BioMolExplorer`, volte à pasta com `pyproject.toml` e repita `python -m pip install -e '.[ui]'`. Confira `python -m pip --version` e `command -v biomolexplorer-ui`: ambos devem apontar para o ambiente ativo. |
| `pyproject.toml` não encontrado | Entre na pasta extraída/clonada do código, não em `src`. |
| Erro na instalação ou em `pip check` | Pare e registre a mensagem completa. Não use `sudo pip` nem instale no Python do sistema. Em ambientes antigos, prefira criar um novo com `conda env create -n BioMolExplorer-clean -f requirements.yml` e ativá-lo antes da instalação editável. |
| Interface desktop não abre | Tente o modo `--web --no-browser` e abra o endereço mostrado no terminal. |

## Escopo e ambiente científico

O caminho recomendado usa `requirements.yml` para as dependências científicas e `python -m pip install -e '.[ui]'` para registrar o aplicativo. Chimera, DOCK6 e Conda são ferramentas externas, não dependências pip. O extra `ui` sozinho não instala o ambiente científico completo. ADMET-AI permanece opcional e isolado; não instale `.[admet-ml]` no ambiente de docking.

Este guia usa Linux/Bash e Python 3.12. Não há lock Linux distribuído neste checkout; não use um arquivo de lock inexistente. Uma instalação limpa requer rede e disponibilidade dos pacotes declarados. `pip check` verifica metadados Python, não executa os protocolos científicos. Valide as ferramentas externas e um caso de referência antes de usar os resultados em pesquisa.
