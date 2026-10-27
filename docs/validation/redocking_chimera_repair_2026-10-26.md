# Redocking com Chimera: correção e validação em 26/10/2026

O projeto analisado contém 35 pares em 34 PDBs. Na execução original, dois pares concluíram e 33 falharam. As configurações salvas foram mantidas: receptor com cargas `gas`, ligantes com `am1`, hidrogênios e minimização habilitados para os dois papéis; Vina com pH 7,4, caixa 24 × 24 × 24 Å, exhaustiveness 20 e dez modos. As opções históricas `auto` são normalizadas para Chimera.

## Causas e correções

- `minimize ... prep false` preservava as cargas, mas também pulava a limpeza estrutural implícita do Chimera. Dock Prep agora resolve conformações alternativas antes de AddH/AddCharge e reconstitui cadeias laterais do receptor. O ligante recebe limpeza própria de alternates, sem alteração de cargas nessa etapa.
- A geometria cristalográfica podia produzir tipagem IDATM inadequada e hidrogênios insuficientes, causando falhas AM1 com número ímpar de elétrons. Definições CCD disponíveis orientam carbonos, oximas e nitronas; o Chimera continua adicionando hidrogênios, calculando cargas e minimizando. O método escolhido permanece igual.
- WriteMol2 pode inferir novamente ordens de ligação incompatíveis com a topologia já parametrizada. A exportação corrige ligações usando a definição química, mantendo coordenadas, cargas e o estado de fósforo preparado pelo Chimera. Referências incompletas preservam os átomos observados, com escopo e ausências registrados em `.export.json`.
- O redocking antes chamava apenas Vina. `redocking_engine` agora seleciona Vina ou DOCK6, com controles específicos, preflight e assinaturas de cache correspondentes. DOCK6 publica poses pelo contrato existente e calcula RMSD de átomos pesados com simetria, no referencial fixo do receptor.
- Os logs mantêm a saída completa do Chimera, inclusive quando retorna código zero com erro. Relatórios antigos de falhas são removidos após uma execução bem-sucedida.

## Execuções reais

**Vina: 35/35 pares concluídos**, com poses e RMSD finito. **MOL2: 35/35 válidos**, conectados e sem elétrons radicalares na leitura de validação; a carga formal é consistente com o inteiro da soma das cargas parciais. A normalização de exportação preservou as coordenadas e as cargas do preparo.

As estruturas e configurações foram lidas do projeto; os ensaios usaram cópias em `/tmp`. As definições químicas da validação completa vieram do cache local do projeto. As entradas e os resultados originais foram preservados.

| Par | RMSD Vina (Å) | Escopo |
|---|---:|---|
| 1P0I_BUA_606A | 1.5674 | Componente completo |
| 1XLU_ISP_1001A | 2.2625 | Fragmento observado |
| 1XLW_DEP_1001A | 8.8275 | Componente completo |
| 2WIG_TC3_601A | 2.4040 | Fragmento observado |
| 2WIJ_TN7_601A | 2.4843 | Fragmento observado |
| 2WIK_TC5_601A | 2.3033 | Fragmento observado |
| 2WSL_EFS_1544A | 2.4970 | Fragmento observado |
| 2XQF_VX_1530A | 9.7461 | Fragmento observado |
| 4BDS_THA_701A | 4.4840 | Componente completo |
| 5NN0_92H_627A | 5.9065 | Componente completo |
| 6QAA_HUN_613A | 4.6289 | Componente completo |
| 6ZWI_QRH_604A | 6.5644 | Fragmento observado |
| 7AWG_S6Q_605A | 6.1394 | Componente completo |
| 7BO3_IA4_610A | 3.5662 | Componente completo |
| 7QBR_AI6_604A | 1.4232 | Componente completo |
| 7QHD_C0I_606A | 2.6032 | Componente completo |
| 8AI7_M8X_612F | 8.4549 | Componente completo |
| 8AI7_MKF_198F | 10.4314 | Fragmento observado |
| 8QTX_WXN_605A | 6.0190 | Componente completo |
| 4M0E_1YL_604A | 0.1302 | Componente completo |
| 5HF5_DEP_604A | 8.2636 | Componente completo |
| 5HF9_HI6_607A | 7.0167 | Componente completo |
| 6NTO_L2Y_605A | 5.9864 | Fragmento observado |
| 6O5V_LND_601A | 6.5624 | Componente completo |
| 6WVO_HI6_604A | 4.4076 | Componente completo |
| 2QP8_SC7_1A | 1.3654 | Componente completo |
| 3VF3_0GS_401A | 0.6260 | Componente completo |
| 4B78_KGG_1385A | 0.5676 | Componente completo |
| 4DJX_0KQ_501A | 6.8221 | Componente completo |
| 4L7G_1W0_501A | 7.3150 | Componente completo |
| 5HU1_66F_501A | 7.4715 | Componente completo |
| 6EJ2_B7E_901A | 0.5913 | Componente completo |
| 6EQM_BUH_401A | 0.9850 | Componente completo |
| 6UWV_QK7_401A | 1.1046 | Componente completo |
| 7MYI_ZQS_401A | 1.2948 | Componente completo |

Execuções DOCK6 com busca flexível, usando receptores e ligantes preparados com as mesmas opções salvas:

| Par | RMSD de átomos pesados (Å) | Grid Score |
|---|---:|---:|
| 1P0I_BUA_606A | 1.6494 | -20.57712 |
| 2WIG_TC3_601A | 4.4552 | -25.46104 |
| 7QBR_AI6_604A | 7.1134 | -70.36139 |

## Regressão, interface e limites

Os 148 testes de regressão passaram. Cobrem configuração, preparação, curadoria e seleção, caches incrementais, isolamento de falhas, entrega de resultados, leitura MOL2, RMSD e localização. O fluxo `execute_operation` também concluiu 1XLU/ISP com preparo novo e reutilização desabilitada: um sucesso, nenhuma falha e nenhuma exclusão.

Os campos dos motores foram renderizados em larguras 560, 1024 e 1440 px, com parâmetros Vina e DOCK6. A documentação HTML foi regenerada e validada. O teste visual usa um fixture de configuração; a seleção real dos 35 pares foi validada no backend.

Esta é uma validação de execução e artefatos. RMSD alto continua representando resultado científico a avaliar; nenhum limiar de acurácia foi usado para converter uma execução em sucesso. Nove referências são fragmentos observados, com átomos ausentes da definição completa. Elas podem diferir do ligante livre integral. O DOCK6 foi executado em três pares, não em todos os 35.

Os valores, parâmetros e hashes estão em [redocking_chimera_repair_2026-10-26.json](redocking_chimera_repair_2026-10-26.json). Logs temporários: `/tmp/biomol-redocking-final-regression-148.log`, `/tmp/biomol-redocking-final-canonical-validation.log`, `/tmp/biomol-redocking-phosphoryl-summary.log`, `/tmp/biomol-redocking-backend-accepted.log` e `/tmp/biomol-redocking-dock6-ai6.log`.
