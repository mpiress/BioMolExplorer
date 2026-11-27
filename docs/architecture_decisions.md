# Registros de decisões de arquitetura

Os ADRs (Architecture Decision Records) documentam as decisões da arquitetura de
artefatos científicos descrita em #6 e nas sub-issues PHASE-00 a PHASE-08
(#7 a #15). Cada decisão registra contexto, alternativas, solução, justificativa,
compromissos, implicações científicas e de migração, requisitos de validação e
desvios da especificação de referência.

## Processo

- Um ADR é criado na fase indicada abaixo, com status **Proposta**, a partir do
  [modelo](adr/0000-template.md).
- A decisão passa a **Aceita** quando os critérios AC associados têm evidência
  verificável. Uma decisão substituída nunca é apagada: recebe o status
  **Substituída por ADR-XXX**.
- Os números são permanentes e não são reutilizados.
- Cada ADR tem uma versão em português (`docs/adr/`) e uma em inglês
  (`docs/en/adr/`).

## Índice

| ADR | Decisão | Fase | Status |
| --- | --- | --- | --- |
| ADR-001 | Modelo de identidade dos artefatos científicos | PHASE-01 | Não iniciado |
| ADR-002 | Tecnologia do catálogo persistente de metadados | PHASE-01 | Não iniciado |
| ADR-003 | Estratégia de armazenamento endereçado por conteúdo | PHASE-02 | Não iniciado |
| ADR-004 | Modelo de assinatura do cache computacional | PHASE-03 | Não iniciado |
| ADR-005 | Rastreamento de dependências e invalidação | PHASE-04 | Não iniciado |
| ADR-006 | Armazenamento de conjuntos de dados analíticos | PHASE-05 | Não iniciado |
| ADR-007 | Interfaces unificadas de acesso a dados | PHASE-05 | Não iniciado |
| ADR-008 | Publicação transacional e recuperação | PHASE-02 | Não iniciado |
| ADR-009 | Migração de projetos e compatibilidade retroativa | PHASE-07 | Não iniciado |
| ADR-010 | Ciclo de vida e retenção de artefatos | PHASE-07 | Não iniciado |
| ADR-011 | Recuperação de dados para visualização | PHASE-06 | Não iniciado |
| ADR-012 | Autoridade dos metadados científicos e procedência | PHASE-07 | Não iniciado |

## Evidência de partida

A avaliação da PHASE-00 está em
[docs/validation/artifact_architecture_assessment_2026-11-26](validation/artifact_architecture_assessment_2026-11-26/README.md).
Ela contém o inventário, as medições de referência, o registro de requisitos e o
plano de implementação que cada ADR deve considerar.
