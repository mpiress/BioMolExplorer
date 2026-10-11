# Correção de conexões após a execução de blocos

## Origem da falha

Os tracebacks de `logs/frontend.log` e `logs/errors.log`, registrados em
11/10/2026 às 23:28 e 23:31 UTC, apontam para o callback `accept` de
`ui/flow_canvas.py`: `self.link(e.src.data, ...)`.

O Flet resolve `DragTargetEvent.src` consultando o controle de `src_id` em seu
registro da página. O canvas reconstruía todos os nós e seus `Draggable`s em
`refresh()`. Essa atualização também ocorre ao selecionar nós e ao terminar o
salvamento automático, que pode substituir o objeto do projeto. Um controle de
origem substituído durante o arraste deixa de estar disponível para o evento de
aceitação. O callback acessava `data` sem resolver essa alteração de identidade.

A exceção ocorre antes de `flow.connect` e da descoberta de arquivos. O traceback
não aponta para metadados de artefatos nem para a invalidação de resultados.

## Correção

- Controles de saída são preservados por UUID do bloco ao reconstruir o canvas.
  Conteúdo e callbacks são atualizados; controles de blocos removidos são descartados.
- O início do arraste registra ID do controle, grupo e identidade da origem.
  Se o registro vivo do Flet já não localizar a origem, o evento só é resolvido
  quando seu `src_id` e grupo correspondem exatamente ao arraste registrado.
- A aceitação consome esse registro. Escape cancela também o arraste registrado.
- Uma origem não resolvível apresenta feedback para repetir o arraste ou conectar
  por clique. Blocos removidos continuam sendo rejeitados pela validação
  transacional existente. Não se escolhe automaticamente outra origem.
- A mesma resolução protege a inclusão de blocos arrastados da biblioteca.

Não foram alterados armazenamento, cache, algoritmos científicos ou descoberta
de artefatos. Resultados históricos continuam sendo recuperados de SQLite ao
abrir a configuração, inclusive em uma nova instância da interface.

## Verificação

`tests/test_canvas_drag_regression.py` contém seis testes de regressão:

1. Identidade do controle após redesenho e substituição do estado pelo salvamento.
2. Aceitação com `src=None` e identificação exata do arraste iniciado.
3. Origem desconhecida ou removida sem mutação indevida do pipeline.
4. Inclusão pela biblioteca mesmo sem o controle vivo de origem.
5. Conexão a um produtor ainda não executado, com seleção automática posterior.
6. Execução inicial, inclusão e conexão de novo consumidor, recuperação dos
   arquivos persistidos pela janela real de configuração, seleção por identidade
   e execução incremental com reutilização do produtor. Os metadados da execução
   anterior são comparados integralmente, e a desconexão remove a seleção disponível.

Os seis testes passaram. A rodada com dez módulos de canvas, navegação,
configuração guiada, artefatos, seleção, tradução e layout passou com 76 testes;
o sexto teste específico foi acrescentado e executado depois dessa rodada.
`compileall` e `git diff --check` também passaram.

O teste de integração usa workers controlados e armazenamento real. Não realiza
cálculos científicos externos nem simula gestos no navegador; exercita os callbacks
Flet com o mesmo estado do evento identificado no traceback.
