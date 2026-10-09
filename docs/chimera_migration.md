# Avaliação da substituição do UCSF Chimera

[Documentação](README.md) · Português · [English](en/chimera_migration.md)

## Decisão e escopo

A revisão de 5 de outubro de 2026 não autoriza retirar o Chimera do backend atual: não foi demonstrada uma substituição completa que preserve os protocolos científicos, as entradas aceitas e os templates editáveis. O backend, os templates e a instalação do Chimera foram mantidos. Nenhuma nova funcionalidade científica foi ativada. Isso não significa que uma implementação Python seja impossível; significa que sua equivalência ainda precisa ser implementada e validada antes de substituir o fluxo existente.

A inspeção abrangeu as referências ao Chimera no código Python, workflows, recursos, instalador, interface, testes e documentação, incluindo a busca por notebooks. Foram consultadas fontes oficiais das alternativas. Não houve comparação experimental entre motores: `chimera` não foi encontrado no PATH desta sessão, as bibliotecas candidatas OpenMM/PDBFixer/OpenFF/ParmEd não estão instaladas no ambiente científico inspecionado e não foram encontrados arquivos de referência PDB/MOL2/PDBQT versionados pela busca nos arquivos do projeto. Os testes com motores simulados não resolvem essa lacuna.

No redocking, as seleções e opções efetivas são geradas para cada par validado. O template base ainda remove solvente e hidrogênios quando usado diretamente pela preparação independente; no redocking, essas remoções são transferidas para as etapas de receptor e ligante para respeitar as opções de cada um.

## Inventário das funções

Os recursos reais estão em `src/biomolexplorer/resources/chimera/`. As referências legadas a `src/scripts/chimera/` são resolvidas por `biomolexplorer.paths.resolve_path`, inclusive para a cópia de recursos de cada worker.

| Recurso ou ponto de entrada | Trabalho realizado | Contrato relevante |
| --- | --- | --- |
| `prepare_complex.template` | Abre PDB, conserva a cadeia selecionada e os cofatores explicitamente declarados, inverte a seleção e remove o restante | `{PDB}_{CHAIN}.complex.pdb`; nenhum cofator é incluído por padrão; solvente e hidrogênios são tratados nas etapas de receptor/ligante |
| `prepare_receptor.template` | Remove ligantes, conserva a seleção `protein`, adiciona hidrogênios, atribui cargas com `chargeModel 14sb method gas`, minimiza e grava MOL2; reabre esse MOL2 e remove H | `{PDB}_{CHAIN}.dockprep.mol2` e `{PDB}_{CHAIN}.noH.pdb` |
| `prepare_ligand.template` | Isola o resíduo por número/cadeia, adiciona H, atribui cargas pelo método escolhido, minimiza, grava PDB e o reabre para exportar MOL2 | `{PDB}_{LIGAND}_{RESNUM}{CHAIN}.lig.pdb` e `.lig.mol2`; o ciclo de reabertura deve ser avaliado ao comparar a preservação de atributos |
| `prepare_better_conform.template` | Abre a primeira pose Vina extraída em PDB, remove H/solvente, adiciona H, atribui cargas, minimiza e exporta MOL2 | `.lig.mol2` usado no refinamento DOCK6; não é apenas conversão de formato |
| `prepare_md.template` | Remove solvente/H, adiciona H e escreve PDB | Recurso disponível no catálogo/editor, mas sem chamada de execução identificada nos fluxos atuais; não implementa uma simulação MD |
| `Docking.prepare_on_chimera` | Executa `chimera --nogui --silent` com o arquivo `.com`, propaga falhas e remove o script apenas em caso de sucesso | Existe em `caad/docking.py` e no módulo duplicado `caad/redocking.py` |
| `wrappers.docking.perform_consensus` | Quando falta `pdb_code`, orienta refinamento manual de loops no Chimera e interrompe o fluxo | Ação humana indicada, sem implementação automática de refinamento |

O wrapper ativo de redocking importa `Docking` e `DockVina` de `caad.docking`; `caad.redocking` mantém outra implementação que também precisaria ser migrada para consumidores desse módulo. `prepare_for_docking` executa complexos, receptores e ligantes nessa ordem, paralelizando cada grupo. No redocking, receptor e ligante usam a mesma cadeia escolhida no formulário de pares. Opções do ligante são compartilhadas com sua conformação, e cofatores são preservados conforme a seleção explícita. Preservar o fluxo inclui essas escolhas, nomes de arquivos e metadados.

`catalog.py` publica os cinco templates para preparação, redocking, Vina e DOCK6. `templates.py` aceita comandos científicos `open`, `delete`, `select`, `write`, `close`, `addh`, `addcharge` e `minimize`; não limita todos os argumentos à configuração original. `ui/guided.py` permite desligar adição de H, minimização e remoções, além de editar métodos de cargas. Uma migração que implemente apenas os cinco templates originais perderia comportamentos configuráveis e precisaria de migração das configurações salvas.

## Dependências que o Chimera não executa

Open Babel gera PDBQT e aplica o parâmetro `pH` posteriormente; os templates Chimera não recebem esse parâmetro. PyMOL calcula centros de massa e participa das análises estruturais. AutoDock Vina realiza docking/redocking. O porte Python do algoritmo DMS calcula a superfície a partir de `.noH.pdb`. DOCK6 e auxiliares realizam esferas, caixas, grades, minimização, docking e footprint. Consenso e tabelas são processados pelo código Python.

Assim, retirar o Chimera não eliminaria todas as soluções de terceiros: bibliotecas Python também são dependências externas, e DOCK6/Vina/Open Babel continuariam no fluxo atual.

## Alternativas pesquisadas

As opções abaixo são candidatas para uma implementação futura, não dependências ou recursos já ativados. Versões publicadas, APIs de documentação e versões resolvidas pelo Conda devem ser distinguidas; não foi escolhido um conjunto de versões sem instalar e validar o ambiente completo.

| Biblioteca | Aplicação indicada | Limite para esta substituição |
| --- | --- | --- |
| Biopython, já declarado | Ler PDB/mmCIF e selecionar cadeias/resíduos | Requer política explícita para modelos, altlocs, inserções, resíduos modificados e classificação `protein`/`ligand`/`solvent`; não faz cargas ou minimização |
| RDKit, já declarado | Identidade química, ligação com SMILES/SDF, sanitização, estereoquímica, Gasteiger e conformeros | MMFF/UFF não reproduzem o protocolo Amber do Chimera; inferência de química a partir de coordenadas não garante equivalência |
| Open Babel, já declarado | API Python para conversões, percepção química e Gasteiger | Não demonstrado como substituto dos métodos e parâmetros atuais; exportar MOL2 não prova equivalência química |
| PDBFixer | Reparar átomos/resíduos ausentes e preparar proteínas | Reparo ou conversão de resíduos altera a estrutura; deve ser opção explícita com proveniência e validação |
| OpenMM | Hidrogênios, parametrização Amber e minimização de proteínas; futura dinâmica | Precisa de templates e parametrização para casos não padrão; minimizador padrão é diferente do Chimera |
| OpenFF Toolkit + openmmforcefields | Parametrização de ligantes e integração com OpenMM | Precisa de ordens de ligação/cargas formais conhecidas; Sage/SMIRNOFF introduz outro protocolo, enquanto GAFF exige backend adequado |
| AmberTools | AM1-BCC e GAFF pela rota Antechamber | Mantém executáveis compilados; uma API Python que os chama não atende à eliminação dessa classe de dependência |
| ParmEd | Transporte de parâmetros/cargas e leitura/escrita MOL2 | Escrita de MOL2 não garante que os tipos atômicos resultantes sejam os esperados pelo DOCK6 |
| Meeko | Preparar entradas AutoDock e recuperar poses com informação química | Melhoria possível na ligação Vina → preparação; não cobre sozinho receptor Amber, minimização e MOL2 para DOCK6 |

[OpenMM Modeller](https://docs.openmm.org/latest/userguide/application/03_model_building_editing.html) oferece hidrogênios com estados dependentes de pH e definições adicionais. Isso permite controles futuros, mas não comprova reproduzir as escolhas locais de [Chimera AddH](https://www.rbvi.ucsf.edu/chimera/docs/UsersGuide/midas/addh.html).

[PDBFixer](https://github.com/openmm/pdbfixer) oferece reparo de átomos, loops e resíduos. Seu uso para preencher estruturas pode ser útil, mas modelar um loop ausente não comprova equivalência ao refinamento manual indicado no wrapper.

[RDKit](https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html) oferece MMFF/UFF e [Open Babel](https://open-babel.readthedocs.io/en/latest/Charges/charges.html) oferece modelos de cargas. São candidatos para tarefas delimitadas; usar outro campo de força não deve ser apresentado como preservação do protocolo existente.

[OpenFF](https://docs.openforcefield.org/en/latest/faq.html) explica a necessidade de identidade química completa e a ambiguidade de PDB. [openmmforcefields](https://github.com/openmm/openmmforcefields) integra parametrização de moléculas pequenas a OpenMM. Minha avaliação é que essas ferramentas compõem uma boa arquitetura futura, mas não resolvem automaticamente a identidade química das entradas atuais.

[ParmEd MOL2](https://parmed.github.io/ParmEd/html/api/parmed/parmed.formats.mol2.html) pode transportar estruturas e [Meeko](https://github.com/forlilab/Meeko) oferece preparação e exportação AutoDock. Uma futura migração pode conservar uma molécula de referência e um mapeamento de átomos durante a recuperação das poses.

## Por que a substituição completa não foi executada

1. **Identidade química insuficiente.** `MolConverter.extract_pdb_to_pdbqt` extrai a região `MODEL 1`/`MODEL 2`, mantém registros ATOM/HETATM/CONECT e trunca linhas a 66 caracteres. O caminho para `recover_better_conforms_of_vina` não passa uma molécula de referência com ordens de ligação e cargas formais. Exigir SDF/SMILES/CCD ou metadados extras sem alternativa validada restringiria entradas atualmente aceitas.
2. **Dois métodos de cargas.** A API expõe `gas` e `am1`; não seria completo oferecer apenas Gasteiger. `am1` é a abreviação Chimera para a rota AM1-BCC, não uma autorização para trocar por cargas MMFF ou predições neurais. Mesmo Gasteiger exige comparação de percepção de ligações, protonação e cargas por átomo entre implementações. O receptor usa ff14SB para resíduos padrão; `method gas` não significa substituir todas as cargas proteicas por Gasteiger.
3. **Minimização diferente.** O [Chimera minimize](https://www.rbvi.ucsf.edu/chimera/docs/UsersGuide/midas/minimize.html) usa MMTK, parâmetros Amber/Antechamber e uma sequência de descida íngreme e gradiente conjugado; também chama preparação estrutural. O [OpenMM LocalEnergyMinimizer](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.LocalEnergyMinimizer.html) usa L-BFGS. As diferenças não provam dano, mas exigem comparar geometrias e resultados antes de afirmar preservação.
4. **Contratos downstream.** DOCK6 consome cargas e tipos de MOL2 para grades/scores/footprints; O gerador nativo de superfícies consome o PDB sem H; centros calculados após minimização definem a caixa Vina. Preservar apenas extensões e nomes não preserva esses dados.
5. **Configurações e ação manual.** Seletores/argumentos dos templates editáveis e a orientação de refinamento de loops precisam de compatibilidade definida. Não foi implementado um interpretador equivalente da linguagem científica Chimera.
6. **Ausência de comparação de referência.** Não há nesta revisão resultados pareados que demonstrem cobertura dos dois métodos e das entradas/configurações aceitas. Aprovação de testes de interface ou arquivos não substitui essa comparação.

A identificação dos métodos é respaldada pela documentação de minimização e pela [orientação da equipe Chimera sobre `method gas`](https://mail.cgl.ucsf.edu/mailman/archives/list/chimera-users%40cgl.ucsf.edu/thread/MRTULUWKFOL6TJWDAQVU2WFPWI3UQSH2/). A página específica de `addcharge` não pôde ser recuperada nesta sessão por restrição de acesso do site.

## Ambiente e validação

`requirements.yml` declara RDKit, Biopython, Open Babel, PyMOL e Vina. Como não há migração ativada, essa lista foi preservada; não foram adicionados OpenMM, PDBFixer, OpenFF, ParmEd, Meeko ou AmberTools sem uso pelo código. A revisão foi anotada nos arquivos. Não foi feita atualização geral de versões nem resolução/instalação Conda, que seriam mudanças independentes e exigiriam sua própria verificação.

Os testes existentes `test_docking_handoffs`, `test_workspace` e `test_stage_default_isolation` verificam contratos de preparação/docking, templates, isolamento e fluxo do workspace. Os motores científicos são simulados nos casos pertinentes. Na revisão, os três módulos passaram: 31 testes em 18,151 segundos, executados com o Python do ambiente BioMolExplorer e `PYTHONPATH=src:tests`. Esse resultado é verificação de software, sem alegação de equivalência química. As páginas HTML foram regeneradas e `git diff --check` não apontou problemas de whitespace.

## Critérios para uma migração futura

1. Criar referências pareadas Chimera/Python para proteínas, terminais, histidinas, pontes dissulfeto, cadeias múltiplas, FAD, ligantes aromáticos/carregados, resíduos modificados e poses importadas, além dos casos de uso reais do projeto.
2. Preservar uma referência química e um mapeamento de átomos das poses, com funcionamento offline quando necessário. Não inferir silenciosamente uma química diferente.
3. Implementar os dois métodos de cargas, parametrização de receptor/ligante, seleções e todas as opções editáveis; registrar versões, campo de força, protonação e alterações estruturais.
4. Definir tolerâncias científicas antes dos testes: identidade e contagem de átomos, ordens de ligação, carga total/por átomo, tipos MOL2, geometria, centros, RMSD, scores Vina/DOCK6, footprint e classificação final. Igualdade byte a byte não é requisito apropriado, mas resultados divergentes precisam de avaliação explícita.
5. Executar os motores reais e os contratos do workspace nos dois modos, com as mesmas entradas. Validar instalação limpa compatível com Python 3.12, atualizar dependências efetivamente utilizadas e migrar templates/configurações salvas.
6. Somente após aprovação desses critérios remover executável, instalador, templates antigos e instruções Chimera, inclusive o módulo duplicado.

Funcionalidades possíveis nesse processo incluem relatório de preparação por átomo/resíduo, reparo opcional de estruturas, estados de protonação explícitos e recuperação de poses com química preservada. Elas são propostas, sem implementação nesta revisão, porque também mudam o protocolo e precisam de validação.
