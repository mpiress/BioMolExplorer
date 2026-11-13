# AMBER Score: saneamento estrutural e compatibilidade com DOCK6 6.11

## Escopo e arquitetura

O suporte é independente do identificador PDB. Não há ramificação por 2WSL no código de produção. Esse alvo é somente o caso de aceitação molecular utilizado nesta revisão.

O caminho de execução é: entrada imutável → cópia saneada → normalização de variantes explícitas → seleção do perfil completo da proteína → parâmetros de componentes retidos → preparação nativa do DOCK6 → validação de topologias → seleção do kernel compatível, quando necessária → refinamento e análise → publicação/cache do par.

Os módulos principais são `structure_sanitation.py`, `md_residue_support.py`, `md_cofactor_support.py` e `md_torsion_support.py`. A integração ocorre em `md_engine.run_pair` e `prepare`. O Redocking aplica a mesma política antes de executar o script de preparo do receptor. A geração registra a seleção; a verificação adiada é necessária porque o PDB do complexo pode ainda não existir nessa fase. As assinaturas públicas existentes permanecem compatíveis; `residue_parameter_policy` é uma opção adicional com valor padrão.

## Saneamento e classificação biológica

`HETATM` não determina se um resíduo pertence à proteína. O saneamento usa conectividade explícita `LINK`/`CONECT`, continuidade peptídica por geometria C–N, cofatores reconhecidos e escolhas explícitas do usuário. O relatório registra uma decisão por resíduo.

- Águas, aditivos reconhecidos e aminoácidos HETATM isolados são removidos somente quando não estão conectados nem foram explicitamente preservados.
- PTMs ligados à cadeia, cofatores reconhecidos e metais são preservados. Espécies não classificadas são conservadoramente retidas para diagnóstico de parametrização.
- No Redocking, o ligante de referência selecionado é excluído da cópia do receptor; ligação covalente desse ligante exige um protocolo explicitamente compatível.
- AltLoc é escolhido de forma consistente por resíduo, pela soma das ocupâncias, com desempate determinístico. Átomos compartilhados sem altLoc são mantidos. Duplicatas remanescentes geram erro.
- Coordenadas não finitas, seriais duplicados e múltiplos modelos sem seleção são rejeitados antes das ferramentas científicas.
- Aliases inequívocos de proteínas padrão são normalizados: OT1/O, OT2/OXT, HN/H e CD/CD1 em ILE. Cada mudança fica registrada.
- A normalização preexistente de HIS usa hidrogênios explícitos; CYS/CYX usa dissulfetos explícitos. Gaps C–N recebem TER na cópia de execução. Não são inventados loops nem estados de protonação a partir de um pH presumido.

O problema anterior do Redocking era a dependência da seleção ampla `protein` do Chimera e de uma limpeza focada em aminoácidos HETATM isolados. Essa abordagem não distinguia suficientemente o papel dos demais heterogêneos. A nova política é aplicada antes dessa seleção e acrescenta os componentes que devem permanecer. O DockPrep não recebe autorização para mutar silenciosamente resíduos modificados.

No 2WSL, o CSS é uma modificação covalente da cadeia proteica. A sua permanência era correta; o erro estava na ausência do modelo de parâmetros e na incompatibilidade das torções. Removê-lo como contaminante ou substituí-lo por CYS perderia um átomo de enxofre e alteraria a química.

## Perfis de parâmetros

`auto_curated` mantém ff94 nativo para proteínas padrão. Quando há um PTM suportado, toda a proteína usa ff03, compatível com os modelos Forcefield_PTM publicados. Não são misturadas cargas de backbone ff94 e ff03.

A distribuição inclui 32 modelos: 0AF, 2MR, ALY, BTK, CGU, CSO, CSP, CSS, DA1, DA2, DAH, DHP, FCY, GCY, HTR, HY3, HYP, LYZ, M3L, MLY, MLZ, ORN, PCA, PRK, PTN, PTR, PYH, SEN, SEP, SPC, TON e TPO. São modelos de estados químicos específicos: CSS neutro; SEP, TPO e PTR com carga −2. A seleção não constitui previsão de protonação fisiológica.

O manifesto contém átomos, tipos, cargas, ligações, carga líquida, hashes, origem e licença. A integridade dos arquivos é verificada antes do uso. Dois erros de nomenclatura na fonte publicada — hidrogênios duplicados em DA1 e SEN — receberam nomes únicos. Os índices e substituições estão em `manifest.json/adaptations`; cargas, coordenadas e parâmetros permanecem iguais.

Átomos pesados ausentes de cadeias laterais de modelos conhecidos podem ser completados pelo LEaP a partir do template. O backbone N/CA/C/O deve existir. O relatório lista os átomos reconstruídos. Após o preparo são verificados inventário completo, ligações internas, tipos AMBER, carga do modelo e preservação das coordenadas originais dentro de 0,002 Å. Modelos internos não são aplicados automaticamente a terminais ou crosslinks incompatíveis.

`native_ff94` permite manter explicitamente o perfil nativo; não fornece, por si só, parâmetros para um PTM ausente. `ff03_ffptm` força um perfil uniforme ff03.

Perfis adicionais validados podem ser fornecidos em `BIOMOL_AMBER_PTM_PROFILES`, apontando para uma pasta com `manifest.json`, `<RES>.prepin` e `<RES>.frcmod`. Usam o mesmo esquema da biblioteca distribuída, precisam declarar ff03, origem, licença, grafos, cargas e hashes, e não podem sobrescrever modelos distribuídos. O conteúdo desses perfis participa da identidade do cache. Esse mecanismo permite incorporar, por exemplo, um modelo de MSE validado sem fingir que selênio é enxofre.

## Parametrização dinâmica de componentes isolados

Componentes orgânicos isolados sem template nativo podem ser parametrizados automaticamente com o grafo do CCD, RDKit e utilitários `antechamber`/`parmchk2` do Chimera, usando GAFF e AM1-BCC. A resolução reutiliza o mecanismo central `run_amber_tool`; não exige ativação manual de Conda nem altera o ambiente DOCK6.

O fluxo verifica nomes/elementos, estereoquímica, estado de carga explícito, hidrogênios observados, ausência de radicais, inventário, ordem de átomos/coordenadas exportadas, soma das cargas e ausência de termos que precisem de revisão no frcmod. Instâncias do mesmo componente/estado compartilham um único modelo durante o preparo. Definições locais são procuradas antes de uma consulta ao CCD; sem definição química confiável, o processo falha explicitamente, sem inferir ordens de ligação pela distância.

O cache é endereçado por definição química, estado de carga, identidade molecular, geometria de parametrização e executáveis. Arquivos intactos já instalados no diretório de execução são recuperados antes de chamar ferramentas novamente. O manifesto de execução registra origem, licença declarada, hashes, método, ferramentas e chave de cache.

**Limite científico:** não existe uma parametrização universal segura para qualquer PTM, metal ou cofactor. PTMs covalentes desconhecidos, MSE sem perfil validado, centros metálicos não suportados, estados de spin e crosslinks especiais continuam exigindo modelos apropriados. GAFF de uma molécula destacada não substitui parâmetros de um resíduo covalentemente incorporado à proteína. O sistema diagnostica esses casos em vez de remover átomos ou inventar parâmetros. Não há previsão automática geral de pKa/tautômeros nem reconstrução de segmentos proteicos ausentes.

## Kernel AMBER local

O kernel original do DOCK6 6.11 trata apenas determinadas periodicidades e fases de torção. Alguns modelos publicados usam periodicidades maiores e fases arbitrárias. São examinados os parâmetros de torção das topologias do receptor, ligante e complexo antes da execução.

Quando necessário, uma cópia local dos fontes é compilada com suporte à expressão completa:

\[
E(\phi)=k\,[1+\cos(n\phi-\delta)]
\]

O gradiente é calculado por diferenciação automática em modo direto para as 12 coordenadas envolvidas. São preservados os termos múltiplos; não há truncamento de periodicidades, arredondamento de fases ou ajuste aproximado dos parâmetros. Termos suportados continuam no caminho nativo otimizado. O detector reconhece o arredondamento nativo de π do LEaP (3,141594), com tolerância de 1e-5 rad, evitando compilações desnecessárias para topologias padrão.

A adaptação requer os fontes da distribuição licenciada DOCK6 6.11 e um compilador C/C++, make, flex e zlib, declarados em `requirements.yml`. A construção ocorre em diretório temporário, com `make -j1`. O executável é armazenado no cache e copiado para `.amber-support/dock6` na execução. A instalação original não é editada e nenhum fonte DOCK6 é redistribuído pelo projeto. Layouts de fontes não reconhecidos produzem erro explícito.

O teste numérico compara 100 energias com o ângulo diedro calculado independentemente pelo RDKit e 1.200 derivadas com diferenças finitas. Essa verificação cobre a expressão e o gradiente; não substitui validação científica de todos os sistemas e parâmetros possíveis.

## Persistência e desempenho

Resultados finais usam o cache verificado por par; preparos de receptor, ligante e complexo têm checkpoints separados. A identidade de preparo é separada de configurações de análise/refinamento que não alteram a parametrização. Mudanças de temperatura podem reaproveitar topologias. O ligante não recebe dependência desnecessária do perfil PTM do receptor.

Arquivos derivados antigos são removidos da cópia de trabalho antes de uma tentativa de restauração/preparo. Assim, falhas não podem ser confundidas com topologias antigas de outro perfil. Logs `.out`/`.log` vazios, produzidos legitimamente quando o preparo é restaurado, não são publicados como objetos de cache. Antes, a exigência de conteúdo nesses diagnósticos podia transformar um refinamento concluído em falha de publicação. Artefatos científicos vazios continuam sendo rejeitados, e o índice do par só é publicado após a validação/gravação dos objetos. Resultados históricos e dados originais não são reescritos. Uma nova política científica pode exigir novo cálculo, sem invalidar a existência histórica dos artefatos anteriores.

## Rastreabilidade e licenças

`residue_support.json`, `protocol.json`, o relatório de saneamento e os manifests de parâmetros preservam a origem das decisões. Os dados Forcefield_PTM são distribuídos sob **CC BY-NC 4.0**, separados da licença do código. Essa condição também aparece na interface e em `LICENSE`/`NOTICE.md`.

Referências primárias:

- Khoury et al., Forcefield_PTM, DOI [10.1021/ct400556v](https://doi.org/10.1021/ct400556v).
- [Suplemento original e metadados de distribuição](https://acs.figshare.com/articles/journal_contribution/2344720).
- [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).
- [Chemical Component Dictionary](https://www.wwpdb.org/data/ccd) e [política de uso wwPDB](https://www.wwpdb.org/about/usage-policies).

## Validação desta revisão

O utilitário `scripts/validate_amber_support.py --dock6-root /caminho/dock6 --output /tmp/validacao` reproduz a validação da biblioteca; `--pairs-json` permite validar pares exportados e o reaproveitamento completo.

As evidências e os tempos finais estão em `validation/amber_support_2026-11-13.json`. Os artefatos completos de execução permanecem nos diretórios de validação em `/tmp`, sem modificar o projeto Alzheimer do usuário.

O caso molecular usa 2WSL e a pose CHEMBL92663 correspondente ao mesmo alvo, vinda do consenso. O protocolo mantém 100 ciclos antes do refinamento, 1.000 passos internos, 10.000 ciclos após, 300 K e GB5. O fluxo executa saneamento, parametrização, validação topológica, construção/recuperação do kernel, refinamento e análise. A repetição é testada com `run_pair` bloqueado para comprovar a recuperação do cache.

Também foram geradas topologias nativas de tripeptídeos para todos os 32 modelos, verificando inventário e cargas, e foi executada uma parametrização GAFF/AM1-BCC de componente orgânico sintético isolado com repetição sem ferramenta. Isso não equivale a simulações completas de todos os PTMs ou cofatores.

A validação de integridade compara 185 arquivos de entrada, resultados e caches já existentes, incluindo 4M0E e 2QMG, e o hash do executável DOCK6 original. Dois testes antigos de Redocking apresentam falhas reproduzidas com `caad/docking.py` do HEAD anterior: expectativa de retry de carga e texto do erro de método/motor. Essas falhas são preexistentes e não são contabilizadas como aprovação global. A interface foi coberta por testes de configuração/localização; não foi realizada aceitação interativa completa do desktop.

Um AMBER Score finito comprova a execução numérica, não afinidade favorável. Valores positivos elevados devem ser examinados quanto à pose, contatos e convergência. Esta revisão não certifica a adequação biológica de todas as poses nem a prontidão global da plataforma.


### Resultados finais no pacote instalado

| Verificação | Resultado |
|---|---|
| Suíte focada final | 108 testes aprovados em 103.76 s |
| Bibliotecas PTM | 32/32 topologias com inventário e cargas validados |
| Expressão de torção | 100 energias e 1.200 derivadas verificadas |
| 2WSL + CHEMBL92663 | 1 par concluído, 0 falhas; 23.89 s com preparos reaproveitados |
| Repetição do mesmo par | 1 resultado de cache; 0.135 s; cálculo bloqueado no teste |
| AMBER Score do par | 2027.359375; valor positivo requer avaliação científica da pose |
| Perfil do 2WSL | ff03 + modelo PTM; CSS preservado e validado |
| 1LGU padrão | ff94 nativo, sem necessidade do kernel adaptado |
| Componentes orgânicos isolados | GAFF/AM1-BCC executado; repetição sem chamar ferramentas |
| Histórico 4M0E/2QMG e demais artefatos auditados | 185 arquivos intactos |
| Instalação do aplicativo | 347 arquivos conferidos contra o wheel; sem diferenças |
| Dependências instaladas | pip check sem dependências quebradas |

O pacote foi aplicado no ambiente Python do BioMolExplorer. A aplicação que já estiver aberta precisa ser reaberta para carregar os módulos atualizados. A instalação original do DOCK6 conserva o SHA-256 `62ddcc919069b6d7af3872b6b06ad344a75ca1d6be2a34013df265579d92f4f2`. O executável local de compatibilidade deste teste tem SHA-256 `57a4b0ea566d4712297301889d6b032eb35938529065ceee51c2e70262bd0322`. Nenhuma edição manual do 2WSL ou da pose foi necessária.
