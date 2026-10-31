> Registro histórico. As instruções de ambiente/protocolos desta validação foram substituídas pela [revisão de Amber Score nativo de 2026-10-29](native_amber_score_2026-10-29/README.md). Não representam as opções atuais da aplicação.

# Migração da preparação molecular — 2026-10-23

## Resultado implementado

O novo caminho é funcional e preferencial no modo **Automático** para candidatos compatíveis sem minimização. Vina usa RDKit/Meeko; DOCK6 preserva a família de cargas Antechamber/Gasteiger ou AM1-BCC selecionada. O modo automático **não** remove a minimização solicitada: esse protocolo, receptores e referências cristalográficas continuam no fallback validado. A opção explícita `legacy_chimera` permanece disponível. Ambientes ainda sem Meeko/Gemmi conservam o fallback automático, com orientação para atualizar dependências.

Diferença entre RDKit-Gasteiger e Antechamber-Gasteiger não é critério de interrupção. A comparação metodológica foi efetivamente executada, os problemas identificados foram corrigidos e o encaminhamento automático foi alterado no código.

## A. Arquivos desta entrega

O workspace já continha alterações anteriores de docking/dinâmica. Esta relação descreve os arquivos criados ou editados nesta migração, não atribui mudanças anteriores à entrega atual.

- Núcleo novo: `src/biomolexplorer/molecular_core.py`, `molecular_adapters.py`, `preparation_tools.py`.
- Integração/validação/cache: `src/biomolexplorer/docking_preparation.py`, `docking_data.py`, `ligand_cache.py`, `dock6_validation.py`, `redocking_config.py`, `stage_cache.py`, `md_engine.py`, `md_inputs.py`.
- Interface/traduções: `src/biomolexplorer/ui/preparation_settings.py`, `src/biomolexplorer/resources/i18n/en.json`.
- Dependências: `requirements.yml`, `requirements-ambertools.yml`.
- Testes: `tests/test_molecular_core.py`, `test_docking_preparation.py`, `test_ligand_preparation_performance.py`, `test_preparation_failures.py`.
- Scripts: `scripts/benchmark_molecular_preparation.py`, `compare_antechamber_migration.py`, `validate_preparation_migration.py`, `build_docs.py`.
- Fontes de documentação: `README.md`, `docs/molecular_preparation.md`, `docs/en/molecular_preparation.md`, `docs/installation.md`, `docs/en/installation.md`, `docs/user_manual.md`, `docs/en/user_manual.md`, `index.html`.
- Relatórios/dados: este arquivo, `preparation_inventory_2026-10-23.md`, `preparation_benchmark_2026-10-23.json`, `preparation_dock6_benchmark_2026-10-23.json`, `preparation_scientific_2026-10-23.json`, `preparation_methodological_redocking_2026-10-23.json`, `antechamber_equivalence_2026-10-23.json`, todos em `docs/validation/`.
- HTML gerado, nos dois diretórios `docs/pt/` e `docs/en/`: `installation.html`, `user_manual.html`, `frontend.html`, `backend_usage.html`, `architecture.html`, `projects.html`, `pipeline_validation.html`, `chimera_migration.html`, `dms_migration.html`, `retrieval.html`, `redocking_configuration.html`, `logging.html`, `molecular_viewer.html`, `admet_research.html`, `molecular_preparation.html`. Os portais `docs/index.html` e `docs/pt/index.html` também foram regenerados para incluir o novo guia.

## B. Comportamento legado determinado no código instalado

Inventário detalhado: [preparation_inventory_2026-10-23.md](preparation_inventory_2026-10-23.md).

Foi inspecionado o código do **Chimera 1.17.3**, especialmente `share/AddCharge/__init__.py`, `share/AmberInfo/__init__.py` e o launcher `bin/chimera`:

- Proteínas padrão: cargas/templates do modelo AMBER selecionado (`14sb` nos templates de docking). `gas`/`am1` não significam calcular a proteína inteira como ligante.
- Ligantes e resíduos não padrão: `_nonStdCharge` constrói MOL2 com H presentes, preserva a ordem de entrada, estima a carga líquida e chama o Antechamber incluído em `bin/amber18`.
- `gas`: `-c gas`. `am1`: `-c bcc` e `-ek "qm_theory='AM1',"`. Ambos usam `-fi mol2 -fo mol2 -nc <carga> -j 5 -s 2 -dr n`, com tipagem GAFF para não padrão.
- Íons monoatômicos têm tratamento específico. Resíduos covalentemente ligados podem receber átomos adicionais/capping e compensação de carga; esse caso não é equivalente a carregar uma pequena molécula isolada.
- O launcher fornece bibliotecas próprias, incluindo `libgfortran.so.3`. O novo adaptador fornece essa biblioteca somente ao subprocesso quando usa o bundle legado, sem contaminar o ambiente Python.
- O preparo anterior de candidatos iniciava Open Babel, Chimera e conversões de exportação por composto não reutilizado. A minimização opcional fazia parte da metodologia e foi preservada, não substituída por um minimizador diferente.

## C. Arquitetura e correções

`MolecularPreparationCore → PreparedMolecule → AutoDockAdapter / Dock6Adapter / AmberAdapter` separa grafo/coordenadas, modelo de carga, tipagem e formato.

- O núcleo trabalha com RDKit Mol e registra formal charge, cargas por átomo, ordens de ligação, estereoquímica, mapeamento/nomes e origem disponível. Recusa coordenadas não finitas ou 2D como se fossem pose 3D.
- Conformeros iniciais de SMILES reutilizam a rotina existente de semente 42/UFF; essa etapa prévia não foi removida para ganhar velocidade e não é confundida com a minimização AMBER opcional.
- `GasteigerEngine`: RDKit dentro do processo. `AntechamberGasteigerEngine` e `AM1BCCEngine`: Antechamber direto, carga formal correta e argumentos metodológicos do legado. AM1-BCC inclui a geometria no cache.
- Vina/Vinardo em política rápida usa Gasteiger/RDKit, registrando método solicitado e efetivo. A opção **respeitar método selecionado** permite AM1-BCC explícito. Configurações como AD4 não recebem a troca rápida. DOCK6/both sempre preservam a seleção Gasteiger/AM1-BCC via Antechamber.
- Meeko lê cargas prontas, tipa AutoDock e pode fornecer PDBQT/string diretamente ao Vina. MOL2 para DOCK6 é gerado do grafo preparado, não de uma conversão PDBQT intermediária.
- Corrigido exportador que marcava ligações não cíclicas de carboxilato como `ar`: isso alterava a correção BCC. Ordens 1/2 e aromaticidade real agora vêm do grafo RDKit; ligações amida permanecem reconhecidas.
- Nomes de halogênios preservam o símbolo completo e o mapeamento não depende de posição de linha não verificada. IDs, conectividade e coordenadas são conferidos após Antechamber.
- O leitor Amber18 não aceita a extensão `UNITY_ATOM_ATTR`. Ela é removida apenas da entrada transitória do Antechamber; a carga formal continua explicitamente em `-nc` e no grafo interno.
- O leitor MOL2 RDKit falhava em certos SYBYL N.4/O.co2. Esses casos usam decodificação Open Babel → SDF em memória, mantendo validação estrita de tipos conhecidos, contagens, conectividade e números finitos.
- Cargas são validadas contra a carga formal com tolerância de 0,001 e epsilon numérico de 10⁻¹⁰; não são renormalizadas arbitrariamente para esconder erros.
- O cache relacional existente distingue as duas implementações de Gasteiger. Cache de exportações e fingerprints de bloco incluem os parâmetros/versões relevantes. A tabela nativa deixa de ser inteiramente reescrita por molécula; há flush final/em interrupção e checkpoints individuais persistentes.
- `AmberAdapter` suporta GAFF/GAFF2 explícitos, `-c rc`, parmchk2 e tleap, validando preservação de cargas e topologia. Parametrização usa diretório isolado para não aceitar topologia antiga quando tleap falha sem código de saída útil.
- A dinâmica só reutiliza AM1-BCC nativo após verificar hash do MOL2 e vetor de cargas no índice compartilhado; não confia apenas em metadados declarados pelo CSV. O protocolo ff94/GAFF, minimização, GB, passos e trajetória existentes não foram alterados.

## D. Dependências

Mantidos RDKit, Open Babel, Vina, NumPy, SciPy, Biopython, DOCK6 externo e os demais componentes já utilizados. Adicionados **Meeko >=0.7.1,<0.8** e **Gemmi >=0.7,<0.8** ao `requirements.yml`; Gemmi é importado pela API Meeko utilizada. Os dois foram instalados e testados somente em `/tmp/biomol-preparation-deps`, sem modificar o ambiente principal.

AmberTools permanece isolado em `requirements-ambertools.yml`, também utilizável para cargas/parametrização, não só trajetórias. A separação evita o conflito de bibliotecas nativas já identificado com Open Babel/PyMOL. Não foram adicionados OpenFF, OpenMM ou PDBFixer sem utilização real. Nenhum pacote existente foi removido por suposição.

O dry-run Conda do requirements atualizado resolveu com sucesso; pacotes pip foram instalados/importados separadamente para os testes. Um dry-run Conda não executa, por si só, toda a instalação pip. A atualização efetiva do ambiente principal continua necessária para habilitar Meeko/Gemmi fora dos testes.

## E. Compatibilidade e ativação

Nomes de outputs, colunas prepared_pdbqt/prepared_mol2, seleção por papéis, opções de pH/H/cofatores e aliases antigos continuam. Não houve alteração de exhaustiveness, sampling, passos ou cutoffs.

O padrão **auto** já está implementado, não é somente uma proposta: sem minimização e com ferramentas disponíveis, usa o novo caminho. Minimização AMBER, preparo de receptor e referências de redocking mantêm o legado necessário. A opção Python explícita não ignora um pedido de minimização incompatível. Entradas somente PDBQT sem grafo de origem recebem orientação para fornecer SDF/MOL2 ou selecionar o legado, em vez de inventar conectividade silenciosamente.

## F. Testes e integração nativa

Regressão ampla: **170 aprovados, 0 falhas, 0 ignorados, 91,172 s**. Abrange núcleo, dinâmica, métricas/painel, localização, preparo/falhas, handoffs, cache incremental dos motores/redocking, pipeline e canvas. Rodadas direcionadas complementares após os últimos ajustes estão registradas no histórico de execução; não devem ser somadas como testes independentes.

Rodada direcionada final de núcleo, opções/traduções, métricas, dinâmica e desempenho/reuso: **64 aprovados, 0 falhas, 0 ignorados, 22,328 s**. Inclui bibliotecas MOL2 indexadas, alteração do modelo de carga mesmo em entrada já preparada, parametrização isolada e reuso entre códigos distintos. Há sobreposição com a regressão ampla.

As ferramentas estavam disponíveis localmente, por isso os testes nativos não foram ignorados. Em máquinas sem Meeko ou AmberTools de validação, os testes condicionais são ignorados com motivo explícito. Não foi reexecutada a suíte completa de todos os subsistemas do repositório.

Incluídos: neutro, cátion, ânion, aromático, heterociclo, halogenado, estereocentro, rejeição de radical/metal/multifragmento/NaN, política de scoring e caches distintos por engine. Parametrização Antechamber/parmchk2/tleap, reuso de GAFF/cargas, AMBER Score e trajetória curta sander/cpptraj foram executados nativamente. Casos de falha parcial e adição de compostos preservam sucessos e reprocessam somente o necessário.

## G. Benchmark medido

Preparação sem docking, minimização explicitamente desativada nos dois caminhos, um worker. Sete grafos representativos são repetidos sob IDs distintos para expor o reuso de cargas; **não são 1000 moléculas quimicamente únicas**. As duas séries foram executadas no mesmo host e podem sofrer carga concorrente; não constituem estimativa estatística de desempenho universal.

| Destino / lote frio | Legado | Nativo | Ganho observado | Chamadas diretas legado → nativo |
| --- | --- | --- | --- | --- |
| Vina / 1 | 1,3510 s | 0,5161 s | 2,62× | 3 → 1 (lscpu de inicialização NumPy; zero Chimera) |
| Vina / 10 | 7,0261 s | 0,4877 s | 14,41× | 30 → 0 |
| DOCK6 / 1 | 1,2849 s | 0,2663 s | 4,82× | 2 → 1 Antechamber |
| DOCK6 / 10 | 6,6537 s | 0,7630 s | 8,72× | 20 → 7 Antechamber |

Vina nativo: 100 em 0,5702 s; 1000 em 1,4921 s; reexecução de 1000 em 0,9750 s. Além de cargas, exportações nativas idênticas são compartilhadas entre códigos diferentes: o lote de 1000 usa sete estruturas e 993 hits de exportação dentro da própria primeira execução. Os 100/1000 legados não foram medidos; não há extrapolação apresentada como resultado. DOCK6 nativo de 10 com cache: 0,0309 s e nenhuma chamada externa direta.

Contagens são lançamentos **diretos pelo Python**: os processos internos do Antechamber/Chimera não são contados como se tivessem sido observados. `observed_tool_work_files` é um limite inferior, não todos os arquivos do SO. Os JSONs incluem timings de leitura, Gasteiger e exportadores; timings aninhados não são somáveis e protonação/cache/embedding não foram isolados artificialmente. Os comparativos pareados incluem tempos por método, com/sem overhead Chimera.

O receptor conservou o mesmo backend; não se inventa speedup de receptor. No caso 1lgu, o preparo legado medido levou 0,8775 s, preservou 162 resíduos e cadeia A e passou de 1290 para 2600 átomos pela adição de H. Os checkpoints existentes continuam evitando repeti-lo por ligante.

Dados: [Vina](preparation_benchmark_2026-10-23.json), [DOCK6](preparation_dock6_benchmark_2026-10-23.json).

## H. Equivalência metodológica e diagnóstico de scores

Foram comparadas **oito moléculas × dois métodos = 16 pares**, com a mesma versão Antechamber do Chimera, geometria/H/protonação, carga formal, nomes, receptor e parâmetros de score. Casos: etanol, metilamônio, acetato, piridina, benzeno, cloroetano, etilenoglicol e um estereocentro fluorado. Todos os 16 passaram nas tolerâncias, com tipos SYBYL iguais.

| Métrica | Limite estabelecido | Maior diferença observada |
| --- | --- | --- |
| Carga por átomo, máxima | 0,001 | 0,000500 |
| Carga por átomo, média | 0,0002 | 0,00009091 |
| Soma de cargas, diferença | 0,001 | 0,001000 (limite inclusivo com epsilon numérico) |
| Continuous Score | 0,05 kcal/mol | 0,000194 kcal/mol |
| Componente eletrostático | 0,02 kcal/mol | 0,000195 kcal/mol |
| Componente vdW | 0,02 kcal/mol | 0 |
| RMSD das coordenadas fornecidas | 0,002 Å | 0 |

São comparações de **score contínuo da pose fornecida**, não buscas DOCK6 completas. Diferenças pequenas são compatíveis com a precisão de escrita do MOL2 legado. O controle inicial com ureia falhou no AM1-BCC legado do Chimera e não foi falsamente contado como par validado; o caso mult heteroátomo da série final é etilenoglicol.

No redocking Vina de **1lgu/AAH**, seed 2026, CPU 1, exhaustiveness 20, 10 poses e caixa 24 Å: score **−5,579** e RMSD pesado **2,78612 Å** nos dois pipelines. RMSD é da primeira pose por score, sem alinhamento nem correção de simetria, contra a geometria cristalográfica preservada na preparação. Isso verifica consistência do caso, não qualidade universal de redocking.

No mesmo caso com H regenerados por preparadores diferentes, DOCK6 mudou de 7,930379 para 7,259476 kcal/mol: vdW 8,880860 → 8,164497; eletrostática −0,950482 → −0,905021. Cargas da mesma família ficaram próximas; os H apresentaram RMSD de aproximadamente 0,04048 Å por melhor correspondência geométrica, sem alinhamento. Portanto essa diferença **não** foi atribuída automaticamente à carga nem considerada motivo para interromper a migração. O controle de geometria idêntica acima elimina essa diferença estrutural e confirma a equivalência da atribuição direta.

Dados: [comparativos pareados](antechamber_equivalence_2026-10-23.json), [redocking metodológico](preparation_methodological_redocking_2026-10-23.json). O arquivo [comparação experimental anterior](preparation_scientific_2026-10-23.json) preserva o resultado RDKit-Gasteiger versus Antechamber-Gasteiger; não representa o caminho DOCK6 padrão atual.

## I. Limitações explícitas

- Receptores, minimização AMBER validada e referências cristalográficas não foram reescritos como um protocolo Python diferente. Eles são fallbacks conhecidos, não funções removidas ou prometidas como migradas.
- Complexos metálicos, covalência com proteína, radicais, cofatores sem parâmetros e sistemas altamente carregados continuam exigindo avaliação/ação explícita. Não houve substituição nativa validada de receptores com metais/cofatores.
- Protonação Open Babel não é AddH Chimera nem enumeração exaustiva de microestados. Comparações devem controlar H e geometria antes de concluir sobre diferenças de carga/score.
- Não foi executada uma campanha sistemática de redocking DOCK6 nem uma longa MD. Os testes curtos validam integração e números finitos, não estabilidade termodinâmica ou afinidade.
- OpenFF/OpenMM/PDBFixer não foram acrescentados só para preencher uma arquitetura conceitual. Meeko/Gemmi devem ser instalados pelo requirements no ambiente real antes de contar com o caminho Vina nativo.
