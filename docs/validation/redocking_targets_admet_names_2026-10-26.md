# Target no redocking e identidade das saídas ADMET — 2026-10-26

## Redocking

`RedockingResults.simulations` acrescenta `target`, `targets` e `best_rmsd_for` a
cada par, preservando sua chave. O alvo vem da proveniência salva do bloco que
recuperou o PDB. A comparação pelo código PDB separa as origens quando um único
metadado reúne vários alvos. Arquivos importados podem informar `Target`/`TARGET`.
Sem identificação segura, o alvo fica vazio; não se infere pelo diretório nem pelo
nome do bloco consumidor.

A interface inclui a coluna Target, filtro de alvo, busca, ordenação por alvo ou
RMSD e checkbox de destaque verde dos menores valores por alvo. Empates usam
tolerância absoluta 1e-9 Å. O mínimo é global para o alvo, independentemente da
paginação ou busca textual; falhas permanecem vermelhas. A exportação CSV dos pares
filtrados inclui Target e Best_for_target e revalida acesso. A atribuição de alvos
não regrava metadados científicos, poses, resultados calculados ou índices de cache.

Uma leitura sem escrita dos resultados existentes identificou 35 pares:
19 BChE PDB, 6 AChE PDB e 10 BACE1 PDB, com um mínimo por alvo. Não houve nova
execução dos motores nem alteração dos arquivos ou banco desse projeto.

## ADMET

O pipeline fornece a categoria e o rótulo de origem ao worker usando a linhagem
exata das entradas. As categorias são MOL, SIM, PUBCHEM, ZINC, MIX e INPUT. MIX
identifica categorias mescladas/coleções consolidadas; INPUT preserva a incerteza
quando não é possível determinar uma categoria. O parâmetro técnico dataset_origin
é validado pela API e não aparece como campo extra no formulário científico.

Tabelas, subconjuntos BBB/HIA, avaliações, gráficos PNG e interativos e relatórios
novos recebem a tag, por exemplo ACHE_molecules_MOL_BBB+.csv e
ACHE_similars_SIM_egg.png. Os relatórios de exclusões também recebem identificação;
seus leitores continuam aceitando o nome antigo. Normalização de nomes inseguros
ou longos usa um digest curto para evitar colisões decorrentes da normalização.

Os sufixos de finalidade permanecem reconhecíveis: _evaluated, _egg,
_admet_report, _BBB+ e similares. Resultados históricos continuam legíveis;
seus arquivos não foram renomeados. A renomeação das novas saídas não cria cópias
adicionais nem altera classificações, propriedades ou critérios científicos.

## Validação

- 103 testes passaram na suíte ampliada: nova identidade dos conjuntos, RMSD,
  exportação, filtros, destaque de empates, ADMET, leitores históricos, qualidade
  molecular, proveniência, processamento individual/merge, i18n e mensagens.
- Após os últimos ajustes de classificação por caminho e cabeçalho, passaram
  29 testes focados; os seis testes de identidade também passaram na conferência final.
- Um pipeline real em projeto temporário executou importação → ADMET → fingerprints
  com entradas molecules e similars. Produziu dois lotes ADMET com nomes distintos,
  transferiu os dois conjuntos aos fingerprints e preservou seus códigos moleculares.
- O Chrome local renderizou a tabela em português a 1200 px e inglês a 390 px.
  O checkbox foi alternado por cliques reais e o estado foi confirmado pelo navegador.
  As capturas foram inspecionadas; não houve exceções JavaScript. Em tela estreita,
  a tabela conserva suas colunas com rolagem horizontal interna.
- A auditoria final percorreu 270 Python, sem imports ausentes ou dependências
  órfãs. Nenhuma nova biblioteca foi necessária para estas duas funcionalidades.
- Os 32 documentos HTML passaram nas verificações de arquivos, links e âncoras.

A validação cobre associação de resultados, interface e processamento ADMET.
Não representa nova validação científica das poses ou dos motores de docking.
