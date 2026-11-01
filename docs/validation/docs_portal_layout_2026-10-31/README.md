# Integração do novo portal de documentação

O layout fornecido em `docs/index.html` foi preservado como base visual. A revisão corrigiu os destinos em português que apontavam para inglês, os atalhos de manual que voltavam à página inicial e a declaração DOCTYPE duplicada. O menu agora também funciona entre 821 e 1120 pixels, faixa em que antes a navegação ficava oculta sem botão de acesso.

Os 15 guias são gerados a partir de `GUIDES`, com títulos/resumos PT/EN, grupos, palavras-chave e os ícones do layout. Cartões, acessos rápidos e percursos abrem os documentos existentes no idioma selecionado. Busca, filtros, contador, teclado e links sem JavaScript foram mantidos. Os rótulos acessíveis também acompanham o idioma.

Fontes de manutenção:

- `docs/assets/portal.template.html`: layout fornecido e interações locais.
- `docs/assets/portal-catalog.json`: categorias, ícones e palavras-chave.
- `docs/assets/portal-translations.json`: textos e rótulos da interface.
- `scripts/build_docs.py`: integra esses recursos, publica os guias e ambos os portais sem apagar o novo design.
- `scripts/validate_docs.py`: valida o catálogo e os destinos dinâmicos PT/EN, além dos links estáticos.

`docs/index.html` mantém a URL canônica inglesa; `docs/pt/index.html` mantém a portuguesa. `?lang=pt` e `?lang=en` permitem compartilhar a escolha. Os dois portais usam a mesma proposta visual.

Verificações executadas:

- Geração de 30 páginas de guias e dois portais; segunda geração com conteúdo idêntico.
- Validação existente: 34 páginas HTML, caminhos/âncoras/IDs e catálogo válido nos dois idiomas.
- Chrome local: dez combinações (dois portais × 320, 390, 768, 1024 e 1440 px), sem overflow horizontal ou exceções JavaScript.
- Troca PT/EN, busca, estado vazio, filtros, limpeza, menu móvel/tablet e Escape.
- Destinos de 15 cartões e seis acessos rápidos/percursos conferidos nos dois idiomas; abertura real dos manuais com conteúdo visível.
- Acesso aos 15 guias com JavaScript desativado e idioma por parâmetro de URL.
- Integridade dos arquivos protegidos: aplicação científica e `index_doc.html` fornecido pelo usuário permaneceram idênticos. A exclusão prévia de `index original.html` não foi revertida.

Correções adicionais de links documentais detectadas pelo validador: a âncora Amber Score do site principal foi alinhada ao título real do manual; a referência a um SVG inexistente foi retirada do índice Markdown, mantendo o PNG disponível. Nenhum algoritmo, resultado de pesquisa, dependência ou configuração do pipeline foi alterado. Não houve instalação de pacotes ou execução científica.

As evidências do navegador e da geração estão neste diretório. A validação ocorreu localmente; não houve publicação em servidor externo.
