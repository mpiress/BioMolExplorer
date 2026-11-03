# Preparo concorrente de candidatos com Chimera — 02/11/2026

## Diagnóstico e escopo

Os caminhos ativos Vina e DOCK6 chamam `docking_data.prepare_ligands`, que usa `docking_preparation.prepare_candidates`. Este último já possuía um ThreadPoolExecutor, fila limitada e diretórios temporários independentes. Portanto, não se atribui todo o ganho serial/paralelo a uma capacidade que antes não existia.

A alteração acrescenta admissão conservadora de recursos e coordenação entre workers, limita as threads numéricas do Chimera e reduz regravações do CSV crescente. O objetivo é tornar a concorrência existente previsível e segura, principalmente com vários blocos/projetos e grandes bibliotecas, preservando ciência e resultados. Preparação de receptores/ligantes cristalográficos por outros caminhos e sessões Chimera manuais não receberam uma nova implementação de paralelismo nesta tarefa.

## Implementação

- `preparation_resources.py`: considera afinidade e quotas de CPU, memória disponível e limites cgroup-v2 herdados. Usa metade das CPUs disponíveis, estimativa de 1 GiB por processo e reserva de pelo menos 1 GiB ou 25% da memória disponível; limita o teto a quatro e respeita o orçamento configurado do worker. Memória pequena/desconhecida implica execução serial.
- Leases de arquivo Linux limitam Chimera de candidatos entre workers do mesmo usuário. O lock permanece durante o comando nativo e é liberado ao concluir, falhar ou encerrar seu grupo de processos. Processos de outros usuários/manuais ficam fora desse controle. Trata-se de admissão por estimativa, não de limite rígido de RAM ou garantia absoluta contra OOM.
- `BIOMOL_PREPARATION_WORKERS` pode reduzir o paralelismo; não ultrapassa o orçamento do worker. A configuração padrão do workspace continua sendo dois workers. Cargas, pH, hidrogênios, minimização e validação molecular não foram alterados.
- O índice `compounds.csv` tem gravações periódicas coalescidas com intervalo mínimo de cinco segundos e gravação forçada ao finalizar ou tratar erro fatal. Arquivos nativos/manifests de cache continuam duráveis por molécula; o CSV não é um checkpoint por item. Isso reduz o custo quadrático de regravar o índice em lotes com muitos hits de cache.
- A mesma conformação preparada continua fornecendo MOL2/PDBQT quando ambos são solicitados. A mudança não cria outro método de cargas nem substitui falhas por aproximações.
- Equivalências exatas em `processing_compatibility.json` preservam identidades anteriores somente para os hashes auditados desta alteração de escalonamento/índice. Alterações futuras de código científico não recebem equivalência automática.

## Medições e limites

Resultados nativos e validações finais estão nos arquivos de evidência deste diretório. Toda a execução experimental usou cópias em `/tmp`. Nenhum pipeline de pesquisa foi disparado, cancelado ou reiniciado; arquivos originais foram apenas lidos.

Os testes com moléculas pequenas e amostras BBB+ reais não constituem validação de toda a biblioteca ou de todas as químicas. Paralelismo não corrige falhas de tipagem de ligações MOL2 ou falhas ANTECHAMBER/SQM do Chimera. As falhas observadas são registradas e não devem ser interpretadas como falhas causadas pelo paralelismo sem comparação do protocolo serial.

### Resultados quantitativos

| Caso | Série | Dois workers | Reúso/cache | Escopo |
| --- | --- | --- | --- | --- |
| 4 moléculas pequenas, Chimera AM1-BCC + minimização, saída MOL2/PDBQT | 4,088 s | 2,312 s | 0,492 s | Quatro aceitas; cargas e coordenadas idênticas entre série/paralelo |
| Primeira amostra BBB+ real | 92,508 s | 58,435 s | 38,299 s | Uma aceita em todos os modos; três falhas químicas reavaliadas |
| Segunda amostra BBB+ real | 74,749 s | 42,122 s | 23,209 s | Duas aceitas em série/paralelo; três na tentativa com cache, que reprocessou as falhas |
| 1.000 entradas de cache, microbenchmark sintético de I/O | Versão anterior: 8,517 s e 1.002 índices gravados | Versão atual: 0,857 s e dois índices gravados | Sem comandos nativos nessa comparação | Isola o custo de índices/cache, não mede o algoritmo químico |

O ganho série/paralelo foi 1,77× nas moléculas pequenas; os testes com dados reais mediram aproximadamente 1,58× e 1,77× no tempo de avaliação dos lotes. Isso não significa aceleração equivalente de toda a aplicação, nem que a versão anterior fosse sempre serial: ela já utilizava concorrência. A melhoria comparativa diretamente atribuível à alteração do índice aparece no microbenchmark sintético de cache.

Nas amostras reais, os scripts experimentais exigiam inicialmente quatro preparos bem-sucedidos e terminaram com a asserção de aceitação não satisfeita. Esses experimentos **não são declarados aprovados integralmente**. Suas medições e conjuntos aceitos foram recuperados e conferidos separadamente em `native-real-samples.json`. Duas falhas foram reproduzidas com `docking_preparation.py` intacto do HEAD anterior: tipagem/topologia MOL2 de CHEMBL4296710 e ANTECHAMBER/SQM de CHEMBL176366. Esses problemas não foram corrigidos alterando cargas ou algoritmos. O novo sucesso de uma molécula reavaliada na etapa warm reforça que essa etapa não é uma medição de cache puro: falhas não são cacheadas como êxitos.

Cargas e coordenadas dos compostos aceitos coincidiram também entre série/paralelo nas duas amostras reais. Os testes não verificam validade farmacológica, acurácia de afinidade, todas as protonações ou toda a biblioteca. Não foram executados novos cálculos de docking/scoring Vina/DOCK6 nativos; foram verificados seus formatos de preparo e as integrações por regressão. A métrica `max_child_rss_kib` dos JSONs é o maior RSS individual de subprocesso observado por getrusage, não o pico de RAM agregado do conjunto de processos.

### Testes e integridade

- 66 testes de regressão aprovados: recursos/cgroups, concorrência entre processos, liberação após SIGTERM, cache/tamper, falhas parciais, formatos de preparo, contexto/ordem das moléculas e integração incremental Vina/DOCK6. A fixture de comandos nativos simulados foi isolada das leases reais para evitar que um benchmark concorrente alterasse sua expectativa de paralelismo.
- Mais 25 testes finais de recursos/cache/resultados publicados foram executados com sucesso após o ajuste final de compatibilidade (`final-cache-and-published-25-tests.log`); há sobreposição de casos com a bateria de 66.
- Lote de 300 entradas em cache preservou todas as linhas/ordem sem novo comando Chimera, com gravações de índice coalescidas.
- Identidade real do preparo comparada com o código anterior: preservada; versões efetivas das 16 operações e da normalização de compostos: preservadas. O cache científico não foi invalidado só pela mudança de escalonamento.
- Arquivos BBB+ originais comparados por SHA-256: íntegros, sem escrita no projeto. Nenhum resultado original removido ou pipeline do usuário reiniciado.
- As demonstrações pequenas e reais usaram o Chimera 1.17.3 instalado e Open Babel no ambiente científico existente. Nenhuma biblioteca foi instalada. A regressão não representa aprovação global de release.

## Arquivos alterados

- `src/biomolexplorer/preparation_resources.py` (novo): política de admissão, leases e ambiente do processo.
- `src/biomolexplorer/docking_preparation.py`: integração dos limites e coalescência do índice.
- `src/biomolexplorer/ligand_cache.py`, `stage_cache.py` e `resources/processing_compatibility.json`: equivalências exatas de versões e separação da política de recursos do contrato científico.
- `tests/test_preparation_resources.py` (novo) e `tests/test_ligand_preparation_performance.py`: limites, concorrência, recuperação e I/O.
- Guias de preparo PT/EN e HTML correspondente, mais este relatório.

Para os dois motores, preferir uma preparação com saída **Vina e DOCK6** já permite compartilhar a conformação e o cálculo de cargas quando os requisitos são equivalentes. Para restringir a carga do equipamento, `BIOMOL_PREPARATION_WORKERS=1` força o modo serial sem mudar o protocolo científico.

### Como repetir a regressão focada

```bash
PYTHONPATH=src:tests python -m unittest test_preparation_resources test_ligand_preparation_performance test_docking_preparation test_preparation_failures test_chimera_preparation test_docking_incremental test_vina_reuse test_dock6_incremental
PYTHONPATH=src:tests python -m unittest test_preparation_resources test_ligand_preparation_performance test_published_results
python scripts/validate_docs.py
```

Os testes nativos requerem Chimera/Open Babel instalados e não são disparados por esses comandos de regressão simulada. Os tempos do relatório são medições locais, não estimativas garantidas para outras bibliotecas.
