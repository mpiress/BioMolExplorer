# Recuperação do pipeline: redocking, seleções ADMET e finalização SQLite

Data: 2026-10-27. Nenhuma biblioteca foi instalada. Banco e arquivos originais foram consultados somente para leitura; ensaios científicos e exportações foram feitos em cópias em `/tmp`.

## Evidências do incidente

Na execução `22b178400b424daea5d27f563e7cd344`, o Vina terminou o cálculo de **5NN0/92H/627/A** com retorno zero. A etapa falhou depois, em `heavy_atom_rmsd` → `validate_mol2`, com carbono de valência 6. Os lotes anteriores tinham sido processados; o erro abortou o restante do redocking.

O MOL2 era exportado após salvar/reabrir PDB. Esse intermediário perdeu ordens de ligação e cargas, gerando `NO_CHARGES` e carbonos incorretamente classificados. Apenas eliminar esse intermediário não bastava: Chimera clássico reclassifica tipos IDATM pela geometria mesmo ao carregar MOL2. O ensaio nativo identificou CAM/CAT como C2 apesar da tipagem C.3 do arquivo. O reparo precisa ocorrer antes de adicionar H e calcular cargas.

A conversão Open Babel com `-p` também apagava nomes atômicos, mesmo com `-xn`. A pose ficava com nomes genéricos C/N/O, impossibilitando o mapeamento confiável do novo RMSD.

**PubChem** e **Generate fingerprints** falharam separadamente em `validate_reference`: suas seis referências persistidas usavam nomes como `B04__admet__bbb-pos.csv` e `B04__admet__rules-lipinski-veber.csv`. O ADMET regenerado passou a usar nomes contendo alvo e MOL/SIM. Além da mudança de identidade, a implementação não exportava os arquivos independentes de regras e alertas. O histórico também marcava as saídas antigas de ADMET como removidas; elas não estavam disponíveis para adoção como entradas íntegras.

Ao finalizar com erro, o log registrou **database is locked**. `_update` mantinha uma transação de escrita e `read_snapshot` abria outra conexão para consultar a pasta do projeto. Quando a escrita grande obtinha bloqueio exclusivo, a própria finalização aguardava a transação que ainda não podia concluir. A execução permanecia marcada como ativa.

## Correções

1. A exportação MOL2 conserva o modelo Chimera, sem reabrir o PDB. O preparo obtém a definição CCD do componente, reutilizando o cache local, conserva os átomos pesados observados e suas coordenadas e fixa tipos IDATM derivados do grafo validado antes de H/cargas/minimização. Átomos ausentes são registrados como fragmentos observados, sem inventar átomos pesados. O MOL2 preparado continua sendo validado estritamente.
2. Nomes atômicos são restaurados após a conversão para PDBQT por correspondência única de elemento/coordenadas, **antes** do docking. Coordenadas e campos de carga do arquivo convertido permanecem iguais. Não há tentativa de inferir esse mapeamento em poses já calculadas. O RMSD conserva átomos pesados, simetria e ausência de ajuste geométrico.
3. Uma falha recuperável Vina/RMSD passa a registrar o par e permitir os demais. Erros globais de ferramentas continuam interrompendo a etapa. A execução Vina recebe somente os scripts do par atual, evitando repetir cálculos de pares anteriores no mesmo diretório.
4. ADMET esquema 4 declara os papéis de suas saídas. BBB± e HIA± são partições independentes do conjunto avaliado; a seleção global não elimina o lado negativo. Regras e alertas voltam a ser exportados. A migração de nomes antigos exige identidade antiga válida, mesmo produtor, mesmo lote e um único papel correspondente. Identidades forjadas, conjuntos ambíguos, regras diferentes e arquivos ausentes continuam sendo recusados.
5. Versões internas de preparo e exportação impedem adotar resultados antigos incompletos como o protocolo corrigido. Exclusão isolada do redocking conserva o cache dos ramos independentes. Resultados removidos, modificados ou com entradas científicas diferentes não são declarados reaproveitados.
6. A leitura do snapshot na finalização usa a mesma conexão/transação SQLite. A regressão verifica a persistência do estado final e proíbe a abertura desse leitor aninhado.

## Validação

- **161 testes** passaram na suíte ampliada em 78,093 s: configuração, métodos, cache incremental e de etapas, referências de artefatos, ADMET, execução, portabilidade, falhas parciais, DOCK6 e idiomas. Duas asserções antigas de testes foram atualizadas para o manifesto de nomes ADMET e a posição da coluna Target, preservando a verificação dos dados e diagnósticos.
- **36 testes** adicionais de recuperação/configuração/reaproveitamento passaram após restringir os scripts Vina ao par atual. Outros **36 testes** passaram após o marcador de protocolo e cache (`/tmp/biomol-recovery-cache-tests.log`).
- As **seis referências reais** dos consumidores foram resolvidas. As **três seleções reais de fingerprints** foram processadas: MOL BBB+ (168 linhas), MOL regras Lipinski/Veber (210 linhas) e SIM BBB+ (203 linhas). PubChem foi validado até a resolução e leitura das três entradas, sem novas consultas de rede. [Evidências dos consumidores](pipeline_recovery_2026-10-27/consumers.json).
- As entradas CCD dos **35 pares configurados** foram construídas e validadas; nove são fragmentos observados. Essa verificação não equivale à execução de docking de todos os pares. [Escopo de cada par](pipeline_recovery_2026-10-27/ccd_inputs.json).
- **5NN0/92H**: novo preparo Chimera com AM1 para o ligante; MOL2 preparado com 71 átomos, zero elétrons radicais e carga formal total zero. Vina executou com caixa 24 Å, exhaustiveness 20 e dez modos; nenhum par falhou. RMSD cristalográfico de átomos pesados: **5,932488 Å**. Preparo do ligante e execução duraram 96,698 s. O receptor foi reutilizado. [Vina](pipeline_recovery_2026-10-27/vina.json).
- **Mesmo par com DOCK6 6.11**: preparação reutilizada, busca flexível, 20 orientações, uma pose e distância de seleção 10 Å; nenhum par falhou. RMSD **5,775106 Å**, Grid Score **−73,005066**, duração 188,010 s. [DOCK6](pipeline_recovery_2026-10-27/dock6.json).

Os RMSDs altos são resultados a avaliar, não aprovação científica do protocolo. A validação nativa não cobre a execução dos 35 pares nem toda a suíte do repositório. Não foi iniciada uma nova execução no projeto original nem alterado seu histórico. É necessário reiniciar a aplicação para carregar o código e executar novamente com reaproveitamento habilitado.
