# Preparação molecular com Chimera

> Histórico da migração nativa. Desde 2026-10-26, docking e redocking voltaram ao preparo exclusivo pelo Chimera. Consulte a [configuração atual](molecular_preparation.md). As descrições abaixo documentam a implementação anterior.

## Caminho ativo

A preparação de receptores e ligantes usa UCSF Chimera. Os controles expõem
remoção de solvente e hidrogênios, adição de hidrogênios, cargas `gas` ou `am1`,
minimização e cofatores. Receptor e ligante possuem configurações independentes.
O pH é aplicado na conversão com Open Babel; não é um argumento dos templates
Chimera. Este guia descreve o código atual, sem instruções para backends alternativos.

## Pontos de entrada

| Operação | Implementação | Uso |
| --- | --- | --- |
| Preparação independente | `wrappers.redocking.prepare_structures` | Receptores e candidatos para Vina, DOCK6 ou ambos |
| Redocking | `wrappers.redocking.perform_redocking` | Pares explicitamente selecionados, preparação e RMSD com Vina |
| Docking Vina | `wrappers.docking.perform_docking_vina` | Candidatos e receptor, com reutilização dos preparados compatíveis |
| Docking DOCK6 | `wrappers.docking.perform_docking_dock6` | Candidatos, receptor, esferas, grades, docking e footprint |
| Consenso | `wrappers.docking.generate_consensus` | Combina scores de pares composto/receptor presentes nos dois motores |

`biomolexplorer.redocking_config` valida pares e aplica as opções aos templates.
`biomolexplorer.redocking_preparation` isola falhas por par e reutiliza receptores.
`biomolexplorer.docking_preparation` prepara candidatos e valida os formatos do receptor.
O backend de docking é implementado em `src/caad/docking.py`.

## Recursos e resultados

Os templates ativos de preparação estão em `src/biomolexplorer/resources/chimera/`:

| Template | Responsabilidade |
| --- | --- |
| `prepare_complex.template` | Selecionar cadeia e cofatores para o complexo de trabalho |
| `prepare_receptor.template` | Preparar o receptor e produzir MOL2 e PDB sem H |
| `prepare_ligand.template` | Separar e preparar o ligante de referência, produzindo PDB e MOL2 |
| `prepare_better_conform.template` | Preparar uma conformação para o protocolo DOCK6 |

Open Babel produz formatos como PDBQT. Centros e metadados identificam o sítio
utilizado. O receptor reutilizado deve conter os formatos exigidos pelo motor:
Vina requer PDBQT; DOCK6 também usa MOL2 e PDB sem H. A superfície SES é gerada
por `biomolexplorer.molecular_surface`, com NumPy e SciPy, e gravada no formato
consumido pelo `sphgen`. Não é necessário instalar o executável DMS.

## Executar e diagnosticar

```text
Conecte receptor e compostos ao bloco Preparar para docking.
Escolha o motor e revise as opções de receptor e ligante.
Execute e confira os arquivos preparados, centros e metadados.
Conecte as saídas ao motor correspondente e avalie seus resultados.
```

Chimera e os executáveis do protocolo precisam estar acessíveis ao worker.
Falhas de comando preservam informações de diagnóstico; consulte os logs do bloco
antes de alterar os parâmetros. O cancelamento encerra os processos científicos.
A existência de arquivos não constitui validação química: avalie estruturas,
RMSD e scores com referências adequadas ao estudo.

Consulte [instalação](installation.md), [redocking](redocking_configuration.md),
[manual](user_manual.md) e [superfícies moleculares](dms_migration.md).
