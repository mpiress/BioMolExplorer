# Equivalência automática de estado e resultados — 28/10/2026

## Comportamento

O controle manual de sincronização foi removido de Runs. A aplicação verifica
equivalência no salvamento e nos protocolos `submit`/`resume`, sem varreduras de
sincronização durante o polling ou enquanto a tela está parada.

| Alteração | Confirmação | Efeito |
|---|---|---|
| Nome/posição do bloco ou opção de motor inativo | Não | Resultados científicos preservados |
| Inclusão de entradas ou ações independentes | Não | Resultados anteriores preservados; somente novas unidades executadas |
| Remoção de uma entrada/ação | Sim | Exclusão da unidade correspondente e dependências afetadas |
| Alteração científica de parâmetros | Sim | Exclusão das unidades cujo protocolo mudou e de seus dependentes |
| Conteúdo de entrada alterado | Sim | Comparação dos hashes persistidos identifica as unidades afetadas |
| Cancelamento da confirmação | — | Configuração anterior restaurada na interface e no banco; arquivos preservados |

O modal lista os blocos cujos resultados serão invalidados, com os nomes padrão
traduzidos e os nomes personalizados preservados. A efetivação da mudança e a
limpeza ocorrem na mesma transação auditada. O hash de confirmação é revalidado
sob lock para impedir que uma confirmação antiga autorize outra exclusão.
O backend recusa alterações destrutivas sem autorização explícita. Uma execução
que ainda utiliza os resultados afetados deve terminar ou ser cancelada; a mudança
não é salva silenciosamente para uma sincronização manual futura. Configurar um
bloco pendente confirma e salva antes de retomar o processamento.

## Estado, granularidade e desempenho

`result_state.py` normaliza parâmetros, templates, entradas, modo de processamento
e eixos comprovadamente independentes. Assinaturas SHA-256 excluem apresentação e
opções irrelevantes de motores/algoritmos inativos. Não se interpreta uma lista
arbitrária de parâmetros como trabalho independente. As chaves de execução
continuam incorporando conteúdo, versão científica, runtime e dependências.

Normalizações possuem cache limitado; checksums de entrada são reutilizados por
inode/dispositivo/tamanho/mtime/ctime. Mudanças desses atributos refazem o hash.
Os valores omitidos são resolvidos pelos padrões reais do worker, preservando a
distinção em relação aos presets de interface. Opções padrão explícitas e rótulos
de perfil ADMET semanticamente equivalentes não causam exclusões desnecessárias.
Importações persistem hashes dos arquivos originais, permitindo detectar mudanças
de conteúdo sem confundi-las com simples adição de arquivos. Escolhas de curadoria
passam a ser configuração durável do projeto, além de estado da execução.

A cascata utiliza uma fila e um índice de consumidores por arquivo/produtor.
Linhagem real identifica os lotes afetados, evitando invalidar outro arquivo do
mesmo bloco. Para registros antigos sem linhagem suficiente, dependências do DAG
são tratadas conservadoramente. Arquivos compartilhados ainda utilizados por outro
bloco permanecem no disco.

Fingerprints são unidades por entrada e algoritmo. A preparação conjunta para
Vina/DOCK6 mantém ambos os formatos na mesma conformação, evitando duas minimizações
independentes. A escolha dos formatos de preparação continua um parâmetro do protocolo.
Morgan e MACCS podem ser reutilizados de um job antigo que os produziu juntos,
desde que os manifests, entradas, protocolo e runtime sejam verificáveis. Contagens
de exclusão de moléculas não são somadas novamente para cada algoritmo da mesma
entrada. Remover um par de redocking preserva arquivos exclusivos de pares válidos
e receptores compartilhados. Tabelas agregadas obsoletas são descartadas; o índice
original é mantido como proveniência e filtrado ao preparar os seeds, impedindo
que pares removidos reapareçam.

`processing_compatibility.json` declara equivalência científica apenas para hashes
exatos auditados da mudança de roteamento de seeds. Preparação, scoring, RMSD e
chaves científicas dos pares não foram alterados. Alterações posteriores no código
produzem outros hashes e não herdam essa declaração automaticamente. Os registros
mantêm a proveniência original do produtor.

Resultados de operações de merge ou métricas globais representam um conjunto
completo: quando seu escopo muda, uma nova combinação precisa ser calculada.
Resultados antigos continuam preservados na adição, mas não substituem um cálculo
global diferente por uma aproximação cientificamente incompatível. A granularidade
por arquivo preserva as unidades independentes no modo individual.

## Evidências

Os logs e a captura ficam em [automatic_result_sync_2026-10-28](automatic_result_sync_2026-10-28/).
Os testes destrutivos e de navegador utilizaram somente projetos temporários.
Nenhuma biblioteca adicional foi instalada e os resultados originais do usuário
não foram removidos para testar a implementação.

A suíte geral final aprovou **282 testes em 249,661 s**. Após preservar a preparação
conjunta em uma chamada e acrescentar sua regressão específica, a suíte final de
sincronização, preparação, handoffs, configuração e traduções aprovou **52 testes
em 30,069 s**. As suítes têm sobreposição; esses números não representam 334 testes
distintos. Os logs são `acceptance.log` e `preparation-final.log`. A documentação
foi regenerada em 28 páginas e dois portais. Sintaxe Python, JSONs de tradução e
`git diff --check` foram verificados.

- Cancelamento: configuração restaurada em memória e no banco, resultados intactos.
- Confirmação: raiz e descendentes listados e removidos, ramo independente mantido.
- Confirmação antiga: recusada quando a requisição efetiva muda.
- Algoritmos: adicionar MACCS após Morgan executa apenas MACCS; mudar o raio de
  Morgan depois preserva MACCS. Um job antigo combinado fornece seus componentes
  sem novos jobs para esses componentes.
- Redocking: remoção de um par preserva outra pose e receptor comum. Um seed parcial
  foi reutilizado sem novo docking e sem ressuscitar o par removido. O backend de
  docking é controlado nesse teste; não é uma nova validação de poses científicas.
- Pipeline científico real de duas fontes: fingerprints → similaridade → grafos;
  seis jobs iniciais e somente três adicionais após alterar uma fonte e confirmar
  a mudança, com reutilização integral do ramo da outra fonte após reinício.
- Navegador real: Chrome, inglês, largura 360 px, ausência do botão manual,
  confirmação antes da exclusão, cancelamento preservando configuração e arquivos,
  confirmação removendo somente os resultados previstos; sem erros JavaScript.
- Benchmark de metadados: 100 blocos, dez execuções, 4.000 unidades de entrada;
  previsão da cascata completa em aproximadamente 0,11–0,12 s em três medições.
  Não mede cálculos químicos, exclusão física nem desempenho geral de pipelines.

Os testes de persistência, recuperação, permissões, interface e execução continuam
separados da aceitação científica de provedores remotos ou de novas poses nativas.
