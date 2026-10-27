# Preparação de ligantes antes do Vina — 20/10/2026

O gargalo observado ocorreu antes da execução do Vina. No job
`94d08f66298d4b3889400c7d5d8a80e4`, os logs registram 112 preparações iniciadas,
101 concluídas, 10 exceções e nenhuma chamada de docking/mapas do Vina. A última
preparação estava em andamento quando o job foi cancelado, após aproximadamente
8 h 51 min. A mediana das preparações concluídas foi 262,541 s e o percentil 95,
664,570 s. Os tempos foram calculados pelos pares de eventos iniciada/concluída
do processamento sequencial antigo.

A configuração selecionada usava cargas AM1 e minimização dos ligantes. O
processamento individual também combina receptores e conjuntos de moléculas;
anteriormente, cada combinação preparava novamente os mesmos compostos. A tabela
de preparados só era publicada ao concluir o conjunto inteiro.

Agora, os ligantes são preparados em paralelo com o orçamento de CPUs do worker
(duas tarefas por padrão no aplicativo). Cada Chimera recebe uma pasta temporária
própria. O coordenador mantém a ordem dos registros e reúne as exclusões por
composto, preservando a identificação e a linha original. Erros de infraestrutura
continuam interrompendo a etapa.

O cache compartilhado do projeto fica em
`.biomolexplorer/cache/ligand-preparation`. Sua chave considera código molecular,
SMILES ou conteúdo/index da pose selecionada, pH, opções efetivas do ligante,
formatos de saída e identidade do código e das ferramentas de preparação. Cada
saída recebe SHA-256; arquivos ausentes ou alterados provocam novo processamento.
Locks por entrada impedem gravações concorrentes da mesma preparação. O manifesto
é publicado após os arquivos. O CSV de preparados é atualizado atomicamente a
cada resultado, e os sucessos do cache ficam disponíveis mesmo após interrupções.

Novos receptores e conjuntos parcialmente sobrepostos reutilizam os compostos
compatíveis, mantendo apenas as moléculas solicitadas no resultado atual. A opção
de executar sem reaproveitamento desabilita a leitura desse cache. Arquivos antigos
sem manifesto não são incorporados automaticamente; o reaproveitamento incremental
vale para preparações produzidas com esta atualização.

Os parâmetros científicos selecionados foram preservados, incluindo cargas,
minimização, protonação, exhaustiveness e número de poses. A interface explica o
custo dos métodos de cargas e distingue a preparação do início do docking. Os logs
passam a registrar identificador molecular e duração da geração da conformação,
do Chimera e da preparação total.

Foi executado um teste real com Chimera 1.17.3, Open Babel e dois compostos pequenos
(etanol e propilamina), mantendo AM1 e minimização:

| Cenário | Workers | Tempo total de preparação |
| --- | ---: | ---: |
| Série, cache vazio | 1 | 3,1490 s |
| Paralelo, cache vazio | 2 | 2,1755 s |
| Paralelo, cache preenchido | 2 | 0,5063 s |

Os hashes dos dois PDBQT foram idênticos nos três cenários. Esta medição valida o
adaptador com ferramentas reais; não estima o ganho das moléculas maiores do
projeto nem mede o docking completo. O benchmark inclui a criação do cache e a
identificação das ferramentas. Os dados estão no JSON ao lado deste relatório.

Dez testes novos verificam orçamento de paralelismo, ordem das saídas, subconjuntos
selecionados, formatos previamente preparados, invalidação por parâmetros,
estrutura/pose e ferramenta, corrupção de arquivos, execução sem reaproveitamento,
exclusões e retomada após interrupção. Uma suíte ampliada executou 130 testes: 129
passaram; permanece a falha anterior
`test_preparation_settings.CandidatePreparationTests.test_exported_candidates_are_reused_without_conversion`,
já reproduzida no HEAD original: o teste exige uma exceção que o tratamento de
falhas por composto registra em relatório.

Na verificação final, os 49 testes específicos de desempenho/preparação,
acompanhamento da execução, diagnósticos e localização passaram, incluindo a
propagação do cache do projeto e da opção de reaproveitamento até o processo filho.

Para usar as mudanças, reinicie o aplicativo e execute o bloco novamente. A
execução longa anterior não foi reiniciada nem alterada durante esta investigação.
