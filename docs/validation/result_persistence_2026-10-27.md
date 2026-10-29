# Persistência e reutilização incremental — 27/10/2026

## Contrato implementado

Os resultados de cada execução permanecem no projeto e podem ser consultados em
execuções posteriores. Uma nova execução não apaga os artefatos anteriores.
Reutilização exige hashes de conteúdo íntegros, requisição científica compatível,
versão de processamento e dependências verificáveis. Falha ou cancelamento de uma
execução não descarta os lotes completos de outras etapas. Lotes incompletos,
artefatos alterados e execuções ainda ativas não são usados como cache.

A identidade de processamento distingue o código científico de alterações de
interface e infraestrutura. Inclui parâmetros, templates, entradas selecionadas,
bibliotecas do interpretador científico e executáveis relevantes. Alterações em
opções de um motor inativo não invalidam resultados do motor selecionado.
Os registros dos artefatos incluem hash, versão científica, runtime, produtor e
manifesto de dependências. Proveniência desconhecida continua explicitamente
desconhecida.

## Granularidade e compartilhamento

- Etapas: reutilização entre execuções, com validação do manifesto completo.
- Lotes: persistência imediata de cada lote concluído e execução somente dos novos
  ou invalidados; falhas posteriores não descartam os anteriores.
- Redocking: reaproveitamento por par e migração conservadora de pares anteriores.
- Receptores: cache de preparação compartilhado entre docking e redocking quando
  estrutura, cadeia, parâmetros, cofatores, pH, templates e protocolo equivalem.
  Opções exclusivas de ligantes não invalidam o receptor. As conversões necessárias
  ao motor continuam obrigatórias; não se reutiliza uma extração de complexo como
  substituta da extração local do ligante.
- Fingerprints: cache SQLite por molécula canônica, método e parâmetros, preservando
  os identificadores e a ordem de cada conjunto. Moléculas equivalentes em arquivos
  diferentes compartilham o cálculo. Desabilitar reutilização ignora esses caches.

## Exclusão e sincronização

| Operação | Remove | Preserva |
|---|---|---|
| Remove All | Resultados gerados e caches computacionais | Todas as configurações, blocos, conexões, originais, logs e histórico |
| Clear All | Resultados, caches, seleções, templates e parâmetros personalizados | Blocos, conexões, originais, logs e histórico; parâmetros voltam aos padrões |
| Remover fonte/bloco | Resultados exclusivamente associados, após confirmação | Resultados independentes e arquivos ainda usados por outros blocos |

Ambos os botões globais têm ícones vermelhos e confirmação descritiva. Cancelar
não modifica arquivos nem configurações. A limpeza rejeita execução ativa.
Salvar um pipeline via API sem `prune_results=True` nunca apaga resultados.
Se a sincronização for adiada, a interface identifica a pendência após reiniciar e
oferece a ação de sincronização com confirmação. Resultados compostos históricos
que deixaram de representar as entradas atuais são preservados e marcados como
históricos; não são tratados como cache válido do novo conjunto.

Arquivos compartilhados são preservados e referências do produtor removido são
destacadas. Exclusões usam journal durável coordenado com o commit SQLite.
Interrupção antes do commit restaura arquivos e metadados; interrupção depois do
commit conclui a exclusão ao abrir o workspace. O espelho `project.json` é
ressincronizado com o banco. Conflitos não sobrescrevem dados existentes.
O histórico deliberadamente preservado permite recuperação auditável; essas
operações não pretendem apagar definitivamente todas as versões históricas.

## Validação

As evidências ficam em [result_persistence_2026-10-27](result_persistence_2026-10-27/).
Todos os cenários destrutivos foram executados em projetos temporários.
Nenhuma biblioteca foi instalada e os resultados do projeto original não foram
limpos para testar a implementação.

A suíte ampla executou 276 testes: 275 passaram e um identificou a tradução
ausente de `Arquivo fora do projeto.` no novo journal. A tradução foi corrigida.
A execução final dos módulos afetados aprovou **117 testes em 115,059 s**,
incluindo internacionalização, persistência, sincronização, limpeza, workspace e
portabilidade. O cenário adicional de remoção de uma fonte importada e os dois
cenários de encerramento abrupto também passaram isoladamente. Os logs
`initial-regression.log` e `final-affected-regression.log` preservam esses resultados
sem ocultar a falha inicial. A sintaxe dos 26 arquivos Python alterados/adicionados,
os quatro JSONs de tradução e `git diff --check` foram verificados. A documentação
estática foi regenerada: 28 páginas e dois portais de idiomas.

| Cenário | Evidência observada |
|---|---|
| Reiniciar serviço e workspace | Resultados íntegros reutilizados sem novos jobs |
| Adicionar fonte | Somente o novo lote processado |
| Alterar fonte em pipeline real de fingerprints → similaridade → grafos | Primeira execução: 6 jobs; após alterar uma das duas fontes e reiniciar: somente 3 jobs adicionais; cadeia da outra fonte reutilizada |
| Alterar parâmetro científico | Somente artefatos atingidos invalidados |
| Alterar código de RMSD | Não invalida fingerprints; alteração do worker de fingerprints invalida |
| Falha em lote posterior | Lote anterior concluído persiste e é reutilizado após reinício |
| Remover arquivo de um lote | Outro lote continua reutilizável |
| Remover fonte/bloco | Cascata segue linhagem real; originais e arquivos compartilhados preservados |
| Remover uma fonte da lista de importação | Cópia normalizada e lote dependente removidos; outra fonte reutilizada sem novo job; ambos os originais mantidos |
| Falha de auditoria durante exclusão | Arquivos, configurações e metadados revertidos |
| Encerrar subprocesso antes/depois do commit | Recuperação idempotente de arquivos, banco e `project.json` |
| Cache versão 1 verificável | Migração integral, crescimento imediato de fontes e novos pares sem repetir os anteriores |
| Interface real | Chrome, inglês 1200 px e português 360 px; confirmações e cancelamentos corretos; sem erros JavaScript |
| Preparação nativa compartilhada | 1P0I/BUA: receptor preparado uma vez, reutilizado por segunda preparação Vina e preparação DOCK6, mesmo SHA-256 |

Na preparação nativa, os tempos observados foram aproximadamente 4,00 s, 2,09 s
e 1,61 s. A minimização foi explicitamente desabilitada nessa configuração de
teste. Foram validados estrutura e formatos de preparação; esses números não são
um benchmark geral nem uma nova validação científica de poses Vina/DOCK6.

## Compatibilidade e limites da validação

A migração automática de cache versão 1 exige build e ambiente explicitamente
auditados no recurso `legacy_cache_compatibility.json`, requisição original e
hashes válidos. Builds desconhecidos e pesos antigos não verificáveis do ADMET-AI
continuam legíveis, mas podem exigir recálculo. Não se inventa versão para um
produtor antigo nem se aceita resultado alterado para evitar processamento.

Os testes cobrem contratos, persistência, invalidação, caminhos nativos de
preparação e a interface descrita. Não constituem aceitação científica completa
de todos os provedores remotos, inferência ADMET-AI ou novas poses de todos os
pares do usuário.
