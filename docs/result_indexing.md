# Indexação dos resultados

Os arquivos científicos e os registros originais em `runs.stages` continuam sendo
a fonte de verdade. A camada de leitura não executa preparação molecular,
docking, alteração de poses ou conversão química.

## Catálogo de execuções

O banco do workspace mantém três tabelas derivadas:

- `result_read_versions`: revisão por execução, atualizada por triggers quando
  `runs.stages` muda, inclusive em exclusões que preservam o horário da execução.
- `result_read_catalog`: representação compacta e comprimida para apresentação.
- `result_read_stages`: consulta por execução e bloco, com descritores de exibição
  separados das provas de integridade e compatibilidade.

Entradas de origem e listas de identificadores repetidas são compartilhadas em
grupos. As provas conservam configurações, manifests e entradas necessárias às
validações. Históricos são percorridos sob demanda; descritores são carregados
apenas para os resultados efetivamente encontrados. O JSON legado é convertido
em fluxo, por descritor, sem materializar todo o histórico em memória.

A API `get_run` mantém seu contrato completo. A interface de consulta usa
`get_result_run`, `list_result_runs` e `result_revision`. Operações de escrita,
histórico, exportação e restauração mantêm os dados originais. O cache de leitura
nunca espera por uma transação de exclusão para preencher seus próprios dados.

## Índice de tabelas

CSVs a partir de 512 KiB recebem um índice SQLite em
`.biomolexplorer/cache/result-index/`. Tabelas pequenas conservam o leitor em
memória. O índice guarda o ordinal original, os dados de cada linha, textos de
busca e, sob demanda, os valores usados na ordenação. Há no máximo oito visões
de busca/ordenação por tabela; páginas e prévias seguintes buscam apenas as
linhas solicitadas.

O índice é invalidado por dispositivo, inode, tamanho, mtime e ctime. Na
construção, calcula-se o SHA-256 usado no contrato de versão da tabela. Uma
alteração concorrente impede a publicação do índice. Arquivos derivados
corrompidos são reconstruídos a partir do CSV.

Os mesmos ordinais são usados na paginação, visualização molecular e curadoria.
Valores ausentes permanecem no fim em ambas as direções de ordenação; empates
preservam a ordem original. Buscas continuam literais e respeitam os campos de
cada visualizador. Autorização, associação ao projeto e possibilidade de edição
são verificadas novamente em cada chamada. Exclusões continuam conferindo o
hash do CSV original e usando as transações auditadas existentes.

Esses índices são descartáveis. Ficam fora dos artefatos publicados e das
capturas de histórico; sua ausência não elimina resultados científicos.

## Verificação

`scripts/benchmark_result_index.py --output /tmp/result-index.json` compara
paginação indexada com releitura/ordenação do CSV. As opções `--database` e
`--run-id` avaliam também um histórico real, com acesso somente de leitura.

`scripts/prepare_result_indexes.py` prepara antecipadamente o catálogo e até
oito tabelas de uma execução concluída, usando `--database`, `--project-id` e
`--report`. Ele verifica que o JSON científico da execução permaneceu idêntico
e não executa motores. Os demais arquivos são indexados quando consultados.

O custo de construção e a latência posterior são medidos separadamente. Os
tempos representam o backend de leitura, não FPS da interface nem execução
científica. `tests/test_result_index.py` cobre ordenação, buscas, versões,
reutilização entre instâncias, recuperação de corrupção e autorização.

## Leitura, escrita e exclusão

Os índices de consulta são derivados. A remoção usa `ResultFiles._transaction_stage` na transação de escrita já aberta; não chama o leitor de projeção para preencher caches por outra conexão durante a exclusão. `result_read_versions` continua sendo atualizado pelos triggers quando o JSON autoritativo muda, e a próxima leitura reconstrói uma projeção obsoleta. Consulte [exclusão e nova execução AMBER](projects.md#excluir-resultados-amber-e-solicitar-novo-calculo).

`CalculationCache` não é necessariamente exclusivo de um projeto: com `BIOMOL_PREPARATION_CACHE`, usa `Path(valor).parent / calculations`; sem essa variável, usa `.calculation-cache` no destino. A configuração do worker pode compartilhar caches entre jobs. Não apague todo o índice compartilhado para solicitar novo cálculo de um único bloco; a geração de resultados MD trata esse caso sem invalidar preparos ou outros consumidores.
