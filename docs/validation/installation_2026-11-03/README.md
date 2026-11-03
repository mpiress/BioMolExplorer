# Revisão da instalação — 2026-11-03

## Diagnóstico e modificações

O ambiente Conda instalava bibliotecas, mas não registrava os entry points do projeto. O manual chamava `biomolexplorer-ui` sem executar a instalação editável. Além disso, `pyproject.toml` continha TOML inválido: `dependencies = [dock6, chimera, miniconda]`.

- Corrigido o TOML: ferramentas externas não são dependências pip. O ambiente científico continua declarado em `requirements.yml`.
- O extra `ui` inclui Flet desktop e web.
- Restaurado `python -m pip install -e '.[ui]'` no README e nos dois idiomas do guia de instalação, com verificações de registro do comando e ajuda.
- Acrescentadas instruções para iniciantes: aguardar cada comando, preservar a pasta da instalação editável, atualizar ambientes existentes, reabrir a aplicação e diagnosticar problemas de Conda, PATH e diretório incorreto.
- Acrescentado Chimera à verificação de executáveis.
- Removidos ADMET-AI, PyTorch CPU e seu índice do ambiente principal: o guia já exige preditor isolado devido à incompatibilidade com o RDKit científico. Ambientes existentes não perdem pacotes com esta edição; para evitar herdar conflitos, use um ambiente novo.
- Removida a indicação de um lock Linux inexistente no checkout e o aviso de empacotamento agora corrigido.
- Regeneradas as páginas HTML de instalação nos dois idiomas.

## Testes e resultados

Comando reproduzível, usando o Python do ambiente científico:

```bash
python -m unittest discover -s tests -p test_installation_guide.py -v
python scripts/build_docs.py
python scripts/validate_docs.py
```

Os três testes passaram: comandos bilíngues e sintaxe Bash dos blocos; contrato dos metadados e separação de ADMET; instalação editável real em venv temporário e execução dos dois entry points com `--help`, fora do checkout e sem PYTHONPATH. As opções usadas no guia foram confirmadas na ajuda da interface.

A instalação de teste usa `--no-deps --no-build-isolation` e herda as bibliotecas do ambiente científico existente. Isso verifica o empacotamento e o registro dos comandos sem baixar novamente o ambiente científico. Não reproduz a resolução de dependências feita pelo comando completo do manual.

Também passou `python -m pip check` no venv de teste: `No broken requirements found`. Documentação: 34 páginas HTML válidas, incluindo arquivos locais, âncoras, IDs e catálogo das operações. `git diff --check` passou.

## Limites da validação

Não foi criado um ambiente Conda científico do zero; disponibilidade dos pins e resolução completa em uma máquina nova permanecem sem validação. A tentativa de smoke test HTTP foi impedida pelo sandbox ao criar um socket (`PermissionError: Operation not permitted`), antes de iniciar o servidor. Não houve validação visual de desktop/navegador nem execução de protocolos Chimera/DOCK6 nesta revisão. Ajuda e metadados não comprovam aceitação científica.
