# Sincronização da documentação — 2026-10-26

## Escopo

Revisão da tarefa 1: documentação, exemplos e apresentação do site. O inventário
foi cruzado com `src/biomolexplorer/operations.py`, `catalog.py`, os wrappers,
serviços do workspace, configuração de preparação e pontos de entrada da CLI.
O catálogo contém 14 operações científicas e um bloco de importação.

## Alterações

- README concentrado em recursos ativos, instalação, CLI e navegação para guias;
  tutoriais duplicados e caminhos antigos foram retirados.
- Manuais português/inglês com um exemplo de quatro linhas para cada um dos 15
  blocos e seis roteiros adicionais para tarefas do workspace.
- Guia de preparação descreve Chimera, formatos e pontos de entrada ativos;
  propostas de backends alternativos foram retiradas. O endereço anterior foi
  preservado para manter links existentes.
- Referências a PubMed como orientação de entrada molecular foram retiradas.
  Citações bibliográficas de métodos científicos foram preservadas.
- Manual corrigido para informar que a exclusão do projeto remove a pasta inteira.
- A tabela técnica de operações agora inclui a recuperação PubChem.
- Site com 15 cartões vinculados aos exemplos, grade de três/duas/uma colunas,
  espaçamento consistente e diagrama das camadas atuais da aplicação.
- Nenhuma biblioteca adicionada; dependências e protocolos científicos preservados.

## Verificações

`python scripts/build_docs.py` gerou 28 guias HTML e dois portais de idioma.
`python scripts/validate_docs.py` validou 32 HTMLs, arquivos e âncoras locais,
links Markdown e cobertura do catálogo. A checagem compara os 15 cartões e os
exemplos dos dois manuais com o registro e confirma seus pontos de entrada em AST,
sem importar motores científicos. Os exemplos por bloco têm de três a cinco linhas.
O JSON de exemplo da CLI no README passou na validação de parâmetros, sem consulta
externa. Os blocos Python dos guias foram analisados sintaticamente e seus imports
internos apontam para módulos existentes.

No Chrome local, as larguras 320, 390, 768, 1024 e 1440 px apresentaram a grade
esperada sem rolagem horizontal da página. O menu móvel e a troca para ADMET no
exemplo interativo funcionaram; não houve exceções JavaScript. Capturas em 320 e
1440 px foram inspecionadas visualmente. Os dois manuais HTML também foram
verificados em 320 px. Imagens com carregamento adiado não foram consideradas
falhas enquanto fora da região carregada; seus caminhos foram validados.

Essas verificações cobrem documentação e interface estática. Não representam uma
nova execução ou validação científica de todos os protocolos descritos.
