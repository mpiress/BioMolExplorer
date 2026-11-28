# Contrato das portas e remoção do reparo de cargas: C5 e C8 — 28/11/2026

Passo 3 e último dos três combinados após a
[triagem](../test_suite_triage_2026-11-27/README.md). O passo 1 está nas
[decisões aplicadas a testes e documentação](../test_suite_decisions_2026-11-27/README.md)
e o passo 2 nas [correções de código](../test_suite_code_fixes_2026-11-27/README.md).
Branch `fix/test-decisions-pr3`, a partir de `a87b39f`.

As duas pendências deste passo eram decisões de contrato, não defeitos. Cada
uma foi investigada antes de ser decidida, e a investigação mudou o quadro das
duas.

| Medida | Passo 2 | Agora |
| --- | --- | --- |
| `unittest discover -s tests` | OK, 1.674 testes em 734 s | **OK**, 1.670 testes em 741 s |
| Falhas esperadas rastreadas | 6 (C5 ×2, C8, I1 ×2, I2) | **3** (I1 ×2, I2) |
| Identidades de cache deslocadas | — | **0 de 26** (ver abaixo) |

Os quatro testes a menos são os de `charge_commands`, removidos com ele.

## C5 — O contrato científico das portas

### O que a investigação encontrou

As duas falhas esperadas de C5 divergiam em **uma célula cada**, não em uma
matriz inteira:

| Teste | Única divergência |
| --- | --- |
| `test_connection_matrix` | Uma importação sem tipo (`other`) alimentando a recuperação ZINC |
| `test_redocking_receptor_outputs` | Poses do redocking alimentando `dynamics_dock6` e `pose_validation` |

### Lista de downloads ZINC

**O código está certo; o teste e a documentação estavam errados.** O tipo
`other` não passa por nenhuma validação (`input_validation.validate_file`
retorna imediatamente), enquanto `zinc_urls` lê e interpreta a lista com
`zinc_retrieval.read_download_list`. Aceitar `other` faria a lista chegar à
recuperação sem leitura. O próprio padrão esperado de `other` já orienta:
*"Para listas de download ZINC, escolha o tipo Lista de downloads ZINC."*

A tabela esperada do teste foi corrigida, e com ela
`docs/pipeline_validation.md`, que listava `zinc_urls` **ou** `other` nos dois
idiomas.

### Saídas do redocking

**Decisão: o redocking entrega apenas os receptores que preparou e validou.**
É assim que docking, Amber Score e dinâmica os reaproveitam sem prepará-los de
novo. O ligante cristalográfico reposicionado é o controle da medida de RMSD
daquele protocolo, não um conjunto de dados para outro bloco.

O contrato tipado era o elemento divergente, e divergia de três lugares que já
diziam a mesma coisa:

| Onde | O que já dizia |
| --- | --- |
| `artifact_choices.redocking_receptors` | *"Public redocking handoffs expose receptors, never crystal ligands."* |
| `docs/pose_validation.md` §1 | As poses a validar vêm de Vina, DOCK6 ou consenso; do redocking vem o receptor |
| `docs/pose_validation.md` §6.4 | A exatidão contra o ligante cristalográfico está fora do bloco: *"o redocking calcula o RMSD"* |
| `docs/pipeline_validation.md` | Saídas do redocking: estruturas, metadados de avaliação e preparados |

Na prática a aresta já era um beco sem saída: `OUTPUT_PORTS['redocking']`
declarava `redocking_results: ('molecular_poses',)`, o canvas aceitava a
ligação, e o formulário de seleção não oferecia nenhum arquivo — o que o teste
vizinho `test_public_choices_contain_only_receptors_and_hide_crystal_ligands`,
que sempre passou, já verificava. Como `base_selected_mols` é obrigatório, o
bloco ficava aparentemente configurado e só falhava na execução.

**Correção.**

| Arquivo | Mudança |
| --- | --- |
| `src/biomolexplorer/artifact_contracts.py` | `INTERNAL_PORTS` declara as portas que publicam resultados sem serem entrada de outro bloco, hoje `redocking.redocking_results`. `output_types(operation,connectable=True)` as exclui. |
| `src/biomolexplorer/flow.py` | Os tipos conectáveis do canvas passam a vir desse cálculo. |

Optou-se por declarar a porta como interna em vez de retipá-la, porque
`OUTPUT_PORTS` tem dois papéis: tipar a conexão e mapear o papel de cada
artefato na proveniência. Retipar mudaria o `artifact_type` gravado nos
registros das poses (`molecular_poses`) e criaria divergência entre os
registros já persistidos e os novos. Com a porta interna, a proveniência,
`infer_output_port`, `ROLE_PORTS` e `operation_contract` seguem intactos: só os
tipos conectáveis mudam.

Os arquivos das poses continuam publicados, baixáveis e listados com seu RMSD.
`docs/pipeline_validation.md` ganhou a seção **Portas internas** e
`docs/pose_validation.md` a nota correspondente, nos dois idiomas.

**Verificação.** `test_redocking_receptor_outputs` passou a percorrer também
`pose_validation`, e a asserção final confirma que a entrega dos receptores
continua aberta para `dynamics_dock6`, `pose_validation`, `docking_vina`,
`docking_dock6` e `prepare_structures` — a razão de existir do bloco.

**Migração.** Um projeto que já tenha essa aresta inútil passa a recusá-la na
validação, com a mensagem de dados incompatíveis, em vez de aceitá-la e falhar
na materialização. A aresta não é construível pela interface desde esta
mudança, e não produzia resultado antes dela.

## C8 — Reparo automático de cargas do Chimera

**Decisão: remover o código morto.**

A investigação acrescentou o dado que faltava à triagem. O reparo nunca foi
religado após `9bdb316`, mas, mais que isso, **o redesenho em duas fases
detecta a condição que ele contornava**: `chimera_preparation.ligand_charging`
lê a carga formal do grafo já hidrogenado e recusa explicitamente um grafo de
camada aberta (*"O ligante hidrogenado contém um grafo molecular de camada
aberta"*) — exatamente o elétron desemparelhado que fazia o AM1 falhar. O
relatório de 13/11/2026 já registrava a expectativa de retry de carga como
falha preexistente.

Religar também não seria restaurar a chamada antiga: o `addh` passou a viver no
script da primeira fase, então o retry teria de reexecutar as duas fases com os
tipos corrigidos. Seria código novo no caminho científico de execução para
recuperar, pelo censo do passo 1, **6 dos 809 compostos excluídos** no projeto
Alzheimer (menos de 1%).

**Correção.** `charge_commands` foi removido de
`src/biomolexplorer/chimera_charge_repair.py`, com os quatro testes que o
exercitavam (três passavam, um era a falha esperada). `local_definition`
permanece — é o que fixa a identidade de cache da definição química local — e a
via CCD de reparo do MOL2 (`chimera_mol2.canonicalize_component`) continua
coberta por seus próprios testes. As variáveis `BIOMOL_CCD_PATH` e
`BIOMOL_PREPARATION_CACHE` continuam documentadas no módulo, porque
`local_definition` as usa.

## Identidades de cache: nenhuma deslocada

Remover código exigiu cuidado próprio. `chimera_charge_repair.py` entra em duas
identidades de cache: na lista explícita de
`ligand_cache.preparation_identity` e, por `caad/docking.py`, no fecho de
importação de `stage_cache.operation_digest`. As 26 identidades foram medidas
antes e depois:

| Identidade | Deslocada pela mudança |
| --- | --- |
| `operation_digest` de `redocking`, `prepare_structures`, `docking_vina`, `docking_dock6`, `consensus` | Sim, por C8 |
| Identidade do adaptador de preparação de ligantes | Sim, por C8 |
| `implementation_digest` e as 19 demais | Não |
| Todas as 26, por C5 | **Não** — a mudança de contrato é de tipagem do editor e do validador |

As seis identidades deslocadas por C8 reduzem-se a quatro pares de digests
distintos, registrados como equivalência auditada em
`resources/processing_compatibility.json`
(`audit_2026_11_28_charge_repair_removal`). O fundamento é direto: o gerador de
comandos IDATM não tinha chamador em `src/`, logo nenhum resultado de
preparação, atribuição de cargas, docking, redocking ou consenso pode diferir.

Medido depois do registro: **as 26 identidades voltaram ao valor anterior à
mudança**. Preparações, redockings e docagens já calculados continuam
reaproveitáveis.

## Pendências que seguem rastreadas

Encerram-se aqui as decisões C1–C8. Restam as duas investigações, ambas da
PHASE-07 (#14):

| Marcador | Teste | Resolve em |
| --- | --- | --- |
| I1 | `test_result_persistence`, dois testes de migração do cache v1 | PHASE-07 (#14) |
| I2 | `test_analysis_origins` | PHASE-07 (#14) |

Fora de #6, a investigação aberta no passo 1 continua: a topologia ou os tipos
atômicos MOL2 respondem por 694 compostos excluídos no projeto Alzheimer —
a causa dominante, cem vezes maior que a de C8.

## Reprodução

```sh
conda activate BioMolExplorer
export PYTHONPATH=src:tests PYTHONDONTWRITEBYTECODE=1
python -m unittest discover -s tests
```

O ambiente precisa estar ativado, e não apenas o seu interpretador invocado por
caminho: `test_molecular_core` e `test_ligand_preparation_performance` chamam
`obabel` por nome, então sem o `bin` do ambiente no `PATH` os dois acusam
`FileNotFoundError: 'obabel'` — uma falha de ambiente, não do código.

A medição das identidades de cache usou apenas leitura, por
`stage_cache.implementation_digest`, `stage_cache.operation_digest` e os
digests dos arquivos do adaptador de `ligand_cache.preparation_identity`.
