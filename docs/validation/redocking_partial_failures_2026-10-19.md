# Redocking com falhas parciais — 19/10/2026

## Comportamento

A preparação usada pelo redocking isola falhas por par. Receptores preparados são
reutilizados por seus ligantes; falhas do receptor são registradas para todos os
pares dependentes. Uma falha individual de ligante não invalida esse receptor.
Os arquivos parciais do par são limpos e os demais pares continuam. A preparação
estrita continua sendo o padrão para outros consumidores de `Docking`.

O wrapper executa Vina apenas para os pares preparados e validados. Ao terminar,
combina os RMSDs efetivamente calculados com as linhas que falharam em
`pdb_codes.csv`. As novas colunas são `STATUS`, `FAILED_STEP` e `ERROR`; uma falha
recebe `preparation_failed` e RMSD vazio. `redocking_failures.csv` reúne os erros.
Se todos falharem na preparação, Vina não é iniciado e a etapa termina com aviso.
Erros globais de instalação, configuração e armazenamento continuam sendo erros,
pois não são defeitos de um complexo específico. Falhas do cálculo Vina ou do
RMSD permanecem sob o tratamento anterior; esta mudança isola a preparação.

A interface mantém as falhas em linhas vermelhas, com estado textual e RMSD
**—**. **Ver erro** abre o motivo e os arquivos de diagnóstico. Pose 3D e
contatos ficam desativados nessas linhas. A contagem de falhas é agregada nos
lotes individuais, preservada no reaproveitamento e exibida na etapa e na
notificação ao encerrar. Pares com falha não são oferecidos como preparados
válidos nas próximas seleções.

## Exclusão explícita

A ação **Excluir falha e complexo de origem** exige confirmação, permissão de
editor/proprietário, execução encerrada e uma linha que ainda represente falha.
A origem é resolvida pelas entradas/procedência do lote; a busca não remove
arquivos de outras coleções apenas por terem o mesmo nome.

A exclusão remove o PDB de origem, a cópia da execução e cópias de conteúdo
verificado no cache de entradas. Atualiza metadados, artefatos, configurações de
pares, integridade e histórico. Arquivos importados também saem do cadastro de
assets e do bloco de importação. Uma falha durante a transação restaura os bytes
e os registros; uma origem modificada desde a execução é recusada.

Como o PDB pode ser compartilhado por vários ligantes, todos os pares desse PDB
saem da próxima execução do mesmo bloco. Resultados históricos bem-sucedidos e
seus arquivos preparados são preservados. A confirmação explica esse efeito.
A recuperação original recebe `excluded_pdb_ids`, editável no formulário, para
não baixar o complexo novamente em uma reexecução. As cópias do histórico para
restauração não são entradas ativas. Não há exclusão automática ao ocorrer erro.

## Validação

**106 testes passaram**, incluindo os cenários novos de falha de ligante,
falha de receptor, todos os pares com falha, reutilização de receptor, limpeza de
arquivos parciais, manutenção de RMSDs válidos, tabela destacada, diagnóstico,
confirmação de exclusão, preservação de resultados do mesmo PDB, isolamento de
coleções, assets importados, exclusão de recuperação, acesso negado, execução
ativa, origem alterada e restauração após falha transacional.

```bash
PATH=/home/michel/progs/miniconda3/envs/BioMolExplorer/bin:$PATH \
PYTHONPATH=src:tests MPLCONFIGDIR=/tmp/biomol-mpl \
/home/michel/progs/miniconda3/envs/BioMolExplorer/bin/python -m unittest \
  test_redocking_partial_failures test_redocking_results \
  test_redocking_configuration test_redocking_execution_review \
  test_preparation_failures test_flexible_retrieval test_retrieval_usability \
  test_localization test_pipeline_execution test_input_processing -q
```

Os testes novos usam armazenamento, preparação de scripts, metadados e controles
Flet reais, com respostas de Chimera/Open Babel/Vina substituídas para injetar
falhas controladas. Não avaliam qualidade científica, convergência ou RMSD com
motores reais. `git diff --check` passou. A documentação fonte foi atualizada e
os HTMLs correspondentes foram regenerados pelo renderizador do projeto.
