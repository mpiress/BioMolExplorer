# Cancelamento de execuções recuperadas

Uma execução aguardando configuração continua disponível depois de reiniciar
o BioMolExplorer. Ela pode ser retomada ou cancelada sem depender de um worker
ativo em memória. Cancelar preserva os blocos concluídos, marca o bloco pendente
como cancelado e ignora os blocos que ainda estavam na fila.

A transição para o estado terminal e seu registro no histórico usam a mesma
transação SQLite. A leitura da versão anterior do histórico recebe a conexão
já aberta, evitando uma segunda conexão durante a gravação. Essa conexão extra
podia bloquear a própria transação quando o SQLite tomava um bloqueio exclusivo
para gravar um registro grande de estágios (`database is locked`). O mesmo
encaminhamento da conexão é aplicado à leitura e restauração dos snapshots.

Depois de confirmado o cancelamento, os resultados podem ser removidos pela
interface, desde que não exista outra execução ativa no projeto. A remoção e a
nova execução mantêm suas verificações habituais de acesso e integridade.
Cliques repetidos em cancelar não criam eventos terminais duplicados.

`tests/test_restarted_cancellation.py` reproduz o estado de bloqueio exclusivo,
fecha e recria o serviço, cancela a execução persistida, verifica o histórico,
remove os resultados e submete uma nova execução com os motores simulados.
Os testes de interface também cobrem os botões de cancelamento na configuração,
na seleção de arquivos e durante execução.
