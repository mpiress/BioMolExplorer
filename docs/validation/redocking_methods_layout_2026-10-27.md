# Redocking: configuração comum, motores e organização visual

Data: 2026-10-27. Nenhuma biblioteca foi instalada. Ensaios nativos foram executados em cópias temporárias; os dados originais dos projetos foram preservados.

## Contratos implementados

- Uma configuração de preparo de receptores e uma de ligantes, ambas comuns a todos os pares; `preparation_options` é a configuração canônica. `preparation_pairs` legado continua suportado e é preenchido para consumidores existentes. Ao abrir opções antigas heterogêneas, a interface informa que salvar aplicará a configuração comum.
- Tabela de curadoria somente de leitura: PDB, Ligand, Residue number, Chain e Actions. Editar abre um rascunho em popup; somente Salvar altera a seleção. Cancelar descarta alterações. Remover retira a linha da seleção. Identidades duplicadas e pares incompatíveis com os metadados são recusados; a resolução pertence à identidade validada.
- Receptores/ligantes com a mesma organização vertical. Cofatores são uma seção expansível; o campo só aparece após habilitar sua inclusão. Painéis responsivos, espaçamento entre campos e textos de ajuda com quebra de linha. Opções de Vina e DOCK6 aparecem conforme a escolha, incluindo as subseções avançadas.
- Textos dos novos painéis, dicas e ações nos dicionários inglês/português da sessão. Corrigida também a tradução parcial de “Remover solvente”. Rótulos acessíveis identificam as ações por ícone.
- Vina usa PDBQT; DOCK6 6.11 usa MOL2 e receptor sem H, sem exigir conversões PDBQT desnecessárias. Seletores aceitam receptores MOL2 quando não existe PDBQT e preferem o PDBQT se ambos existirem. Um consumidor Vina pode exportar o formato rígido em uma cópia de trabalho, sem repetir o preparo Chimera nem reaplicar cargas.
- DOCK6 reutiliza superfície/esferas do receptor, mas seleciona o sítio e calcula grades por ligante. Cada execução copia somente os arquivos daquele receptor. Poses adicionais só geram arquivo separado quando existem múltiplas poses; scripts/logs nativos são preservados mesmo em falhas.
- Cache separa motores, parâmetros, ferramentas, templates, configurações e integridade de arquivos. Mudanças de semente ou binário DOCK6 invalidam as poses correspondentes.
- RMSD novo: átomos pesados, correção de simetria sem ajuste geométrico, referência nas coordenadas cristalográficas quando o PDB está disponível. Fallback de referência preparada e cálculo Vina legado são identificados nos metadados e no CSV; resultados antigos não são recalculados ao serem consultados.

## Verificações automatizadas

122 testes passaram em 35,034 s, incluindo os novos testes em `tests/test_redocking_methods.py` e as suítes de configuração, curadoria, reaproveitamento, tradução, falhas DOCK6, integração de docking, resultados e orientação do pipeline.

```bash
PATH=/home/michel/progs/miniconda3/envs/BioMolExplorer/bin:$PATH \
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/biomol-mpl PYTHONPATH=src:tests \
/home/michel/progs/miniconda3/envs/BioMolExplorer/bin/python -m unittest \
  test_redocking_methods test_redocking_configuration test_redocking_incremental test_localization \
  test_docking_configuration test_docking_handoffs test_docking_preparation test_docking_interoperability \
  test_dock6_failures test_redocking_execution_review test_curated_redocking test_redocking_results test_pipeline_guidance
```

O escopo acima é focado. Não foi executada toda a suíte do repositório e não constitui aprovação de release.

## Interface renderizada

Flet real em Chrome local com recursos locais, inglês em 390/900/1440 px e português em 1440 px, altura de 1000 px. Foram inspecionadas tabela, preparação comum, método e opções avançadas dos dois motores. O roteiro abriu **Edit pair** e fechou com **Cancel**; não registrou exceções JavaScript. Os testes de unidade verificaram separadamente salvar, cancelar, duplicatas, remoção, cofatores, idioma e preservação das configurações.

- [Tabela em inglês](redocking_layout_2026-10-27/pairs-en.png)
- [Preparação alinhada em inglês](redocking_layout_2026-10-27/preparation-en.png)
- [Popup de edição](redocking_layout_2026-10-27/edit-en.png)
- [Opções DOCK6 em tela estreita](redocking_layout_2026-10-27/dock6-mobile-en.png)
- [Resumo do navegador](redocking_layout_2026-10-27/browser.json)

## Execuções nativas

| Entrada | Preparação e busca | Resultado observado |
|---|---|---|
| Tutorial DOCK6 1ABE | MOL2 preparado, rigid, 50 orientações | DOCK6 6.11 concluiu; RMSD 0,196376 Å; Grid Score -37,541985; referência preparada |
| Par previamente configurado 2WSL/EFS/1544/A | PDB bruto → Chimera, gas para ambos, sem minimização; DOCK6 rigid, 50 orientações | Concluiu sem falhas; RMSD 6,797907 Å; Grid Score -19,606789; referência cristalográfica |
| Mesmo par e geometria Chimera | Conversão de formato para Vina; exhaustiveness=1, num_modes=2, caixa 18 Å | Vina concluiu sem falhas; RMSD 1,814364 Å; referência cristalográfica |
| Mesmo par preparado | DOCK6 flex, 30 orientações, 2 poses | Concluiu sem falhas; RMSD 3,357124 Å; Grid Score -20,492163; referência cristalográfica |
| Mesmo par, revisão final | DOCK6 flex, 10 orientações, 1 pose | Concluiu sem falhas em 161,087 s; RMSD 3,691323 Å; Grid Score -20,075611; referência cristalográfica |

Esses ensaios comprovam execução, formatos, referência e persistência nos casos observados. Usaram busca reduzida; RMSDs altos não foram classificados como aprovação. Não validam o protocolo de todos os alvos/ligantes, convergência, qualidade das cargas ou superioridade de um motor. O custo observado do caso final é específico deste receptor, superfície e equipamento.
