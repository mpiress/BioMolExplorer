# Convenção global de artefatos

A publicação dos 16 tipos de bloco utiliza uma única implementação em
`biomolexplorer.artifact_naming`, sem modificar algoritmos ou formatos científicos.

`TG[ALVO]-SRC[ORIGEM]-OP[OPERAÇÃO]-TIPO-ID[identificador].ext`

- **TG**: alvo biológico identificado na recuperação/linhagem. Docking e Amber Score usam o alvo do receptor, independentemente da coleção dos compostos. Coleções independentes sem alvo usam `NA`. Resultados agregados podem representar vários alvos; listas extensas usam `MULTI-hash`, com a lista completa nos metadados.
- **SRC**: categoria da entrada imediata: `MOL`, `SIM`, `BBB+`, `Fingerprints-Morgan`, `REC`, `PDB` etc. Misturas extensas usam `MIX-hash`; os tipos e IDs de todas as entradas permanecem registrados. Arquivos auxiliares de um conjunto estrutural não inventam novas categorias biológicas.
- **OP**: produtor, como `ChEMBL`, `PubChem`, `ADMET`, `Fingerprints`, `Similarity`, `Graphs`, `Redocking`, `Preparation`, `Docking-Vina`, `Docking-DOCK6`, `Consensus`, `AmberScore-DOCK6` e `Import`.
- **TIPO**: conteúdo: `MOL`, `SIM`, `BBB+`, `BBB-`, `HIA+`, `HIA-`, `Rules-Lipinski`, `Alerts`, `Evaluated`, `Fingerprints-Morgan`, `Similarity-Tanimoto`, `Molecules-C2`, `Edges-LCC`, `Scores`, `RMSD`, `Receptor-4M0E_A`, `Pose-CHEMBL1179697`, `AmberScore`, `Protocol`, `Topology`, `Coordinates`, `Exclusions`, `Methods`, `Figure` etc. Uma trajetória só recebe esse tipo se existir no formato correspondente; a convenção não cria métricas ou arquivos científicos ausentes.
- **ID**: prefixo determinístico da ocorrência do artefato (normalmente oito dígitos hexadecimais). Está sempre presente para distinguir também resultados históricos. Colisões de prefixo detectadas no catálogo expandem esse sufixo. Não substitui o identificador interno completo nem a validação de conteúdo.
- **ext**: formato original, inclusive extensões compostas, como `.csv.gz` e `.biomol-view.json`. Os sinais `+` e `-` são preservados.

## Publicação, arquivos físicos e compatibilidade

O nome público é um alias persistido no descritor do artefato e, para novas saídas próprias, em `.artifact-publication.json`, gravado atomicamente. Listas de arquivos, seletores, tabelas e downloads usam esse alias; continuam selecionando/validando pelo ID e pelo caminho autorizados, nunca somente pelo texto exibido. O índice é metadado, não uma entrada científica do cache.

Nomes físicos operacionais como `centers.csv`, `pdb_codes.csv`, `*.dockprep.pdbqt` e referências relativas dentro das tabelas permanecem no protocolo nativo. Os nomes analíticos físicos já usados pela publicação antiga também continuam compatíveis. Isso evita romper ferramentas, bundles e resultados previamente validados, sem duplicar arquivos ou renomear dados existentes.

Downloads individuais recebem o nome público. ZIPs científicos preservam os caminhos internos necessários à execução e acrescentam `artifact-names.json`. Backups `.bme.zip` preservam a estrutura e contêm `artifact_names` no manifesto; a restauração conserva os descritores, inclusive quando precisa remapear IDs de execução.

O parser compartilhado reconhece a convenção nova e as publicações `B01__...`. Nomes legados arbitrários continuam aceitos pelos leitores existentes. Estruturas exportadas com nomes canônicos têm uma camada de adaptação para os nomes operacionais; a validação de átomos, coordenadas e cargas continua obrigatória. Arquivos isolados não recriam por nome uma linhagem que não foi exportada: para rastreabilidade completa use o pacote do projeto.

## Metadados e validade

O registro mantém arquivo original, caminho, ID do artefato e da ocorrência, parâmetros do produtor, versão de processamento, SHA-256 disponível, entradas e suas dependências. Nomes não decidem validade, identidade, permissões ou invalidação seletiva.

A alteração de nomenclatura tem compatibilidade de versões exata e auditada em `processing_compatibility.json`, sem ignorar mudanças futuras de algoritmos. Apenas versões específicas de alterações de publicação/roteamento recebem equivalência. Dados científicos alterados, ausentes ou inválidos continuam sujeitos às verificações existentes.

Consulte o [relatório de validação](validation/global_artifact_naming_2026-10-31/README.md) para exemplos reais e limites dos testes.
