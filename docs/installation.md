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
| UCSF Chimera | 1.17, opcional | Fallback legado escolhido explicitamente | [Downloads de versões anteriores da UCSF](https://www.cgl.ucsf.edu/chimera/olddownload.html) e [instruções de instalação](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/installation.html) |
| DOCK6 | 6.11 | Docking, refinamento e cálculo de scores, com seus programas auxiliares | [Página oficial DOCK 6](https://dock.compbio.ucsf.edu/DOCK_6/index.htm), [notas da versão 6.11](https://dock.compbio.ucsf.edu/DOCK_6/new_in_6.11.txt) e manual incluído na distribuição |

Instale UCSF Chimera 1.17, Open Babel e Vina para a preparação/redocking; instale DOCK6 6.11 para seus blocos de docking. AmberTools separado é opcional para protocolos de dinâmica, não para preparar docking. ChimeraX não interpreta os scripts Chimera.

Instale o Chimera com o instalador adequado ao sistema. Para DOCK6, siga as instruções da distribuição 6.11, compile também os programas auxiliares e mantenha a instalação completa, incluindo `bin/` e `parameters/`. Execute os testes fornecidos pelo DOCK6 conforme seu manual.

A superfície molecular é calculada nativamente em Python com NumPy e SciPy. Não é necessário instalar DMS; consulte o [porte e sua validação](dms_migration.md).

As ferramentas exigidas pelo protocolo devem estar disponíveis no computador dos workers, inclusive em implantação web. Não instale Chimera apenas para satisfazer o caminho padrão: ele é opcional.

## 3. Criar o ambiente científico e instalar a interface

Na raiz do código baixado:

```bash
conda env create -f requirements.yml
conda activate BioMolExplorer
python -m pip install -e '.[ui]'
```

Se o ambiente já existe, atualize-o com `conda env update -f requirements.yml --prune` e depois ative-o. O ambiente inclui Open Babel, Vina e a interface; ele não substitui Chimera e DOCK6, que continuam instalações externas. O arquivo mantém faixas compatíveis — em especial NumPy 1.x e pandas 2.x — para impedir que extras de aprendizado de máquina alterem silenciosamente o worker científico. O pacote requer Python 3.12.

## 4. Configurar caminhos e conferir a instalação

Com o ambiente Conda ativo, adicione ao `PATH` os diretórios que contêm os executáveis. Substitua todos os caminhos de exemplo pelos caminhos reais da sua instalação:

```bash
export BIOMOL_DOCK6_ROOT="/caminho/para/dock6-6.11"
export PATH="/caminho/para/chimera-1.17/bin:$BIOMOL_DOCK6_ROOT/bin:$PATH"
```

`BIOMOL_DOCK6_ROOT` é uma variável de conveniência deste guia. A aplicação recebe esse caminho pela opção `--dock6-path`, mostrada abaixo. Informe a raiz do DOCK6, que contém `bin/` e `parameters/`, e não apenas o executável ou a pasta `bin/`. O mesmo terminal deve iniciar a aplicação para que os workers herdem o `PATH`. Para manter a configuração entre sessões, inclua os exports com os caminhos reais no arquivo de inicialização do seu shell e abra um novo terminal.

Confira os executáveis antes de iniciar:

```bash
for ferramenta in dock6 sphgen showbox grid obabel vina; do
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

**Não inicie a aplicação enquanto houver uma mensagem FALTA ou o diretório `parameters/` estiver ausente.** Confirme os caminhos das ferramentas realmente usadas; Chimera 1.17 é opcional e DOCK6 6.11 é exigido apenas por seus protocolos, especialmente se houver várias versões instaladas. A localização de um executável confirma sua presença no `PATH`, mas não valida sua versão nem o funcionamento do protocolo. Confira as versões nas ferramentas/distribuições e faça uma execução com um caso de referência antes de iniciar seu estudo.

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

**Dinâmica com DOCK6** usa o AMBER Score e os utilitários AMBER incluídos na instalação completa do DOCK6. O protocolo `DOCK6 + sander` requer também `sander`; a análise opcional da trajetória requer `cpptraj`. Instale-os em ambiente isolado, pois o AmberTools atual conflita com as bibliotecas nativas de Open Babel/PyMOL do worker principal:

```bash
conda env create -f requirements-ambertools.yml
```

Nenhuma ativação adicional é necessária. O bloco usa **Ambiente Conda AmberTools** = `BioMolExplorer-AmberTools` por padrão. Informe **Instalação AmberTools** apenas para substituir o ambiente por uma raiz explícita.

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

Se aparecer **ADMET-AI is not installed in the worker environment**, confira o estado mostrado na configuração ADMET, o caminho, a instalação e o acesso do worker aos modelos. Reinicie a aplicação após configurar o ambiente. Como alternativa, selecione explicitamente **RDKit local**; isso fornece descritores, estimativas e regras locais, sem os endpoints `ML_`. Não há substituição silenciosa de modelos. Outra opção é `python -m pip install -e '.[admet-ml]'` em um worker dedicado, respeitando as dependências desse ambiente.

Consulte [ADMET por hipótese](admet_research.md) para unidades, perfis, alertas, arquivos e limites dos modelos.

### AmberTools para protocolos de dinâmica

Crie também o ambiente isolado, sem ativá-lo como worker principal:

    conda env create -f requirements-ambertools.yml
    export BIOMOL_AMBERTOOLS_ENV=BioMolExplorer-AmberTools

Confira nessa raiz bin/antechamber, bin/parmchk2, bin/tleap e bin/sander. cpptraj é necessário somente para análise de trajetória. O backend inicia cada executável com ambiente próprio; não exporte LD_LIBRARY_PATH do AmberTools para o Python principal.
