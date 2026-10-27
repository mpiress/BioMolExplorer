# Migração do caminho principal sem Chimera

> Histórico da migração nativa. Desde 2026-10-26, docking e redocking voltaram ao preparo exclusivo pelo Chimera. Consulte a [configuração atual](molecular_preparation.md). As descrições abaixo documentam a implementação anterior.

Atualização: 23 de outubro de 2026.

## Estado atual

Seleção/separação estrutural, extração cristalográfica com topologia confiável, preparação padrão de receptores ff14SB, cargas de pequenas moléculas e minimização sander agora seguem Python/Gemmi/AmberTools. Vina usa Meeko/AutoDockAdapter; DOCK6 recebe MOL2 com cargas preservadas. Chimera não é dependência desse caminho principal.

O backend legacy_chimera permanece explícito para protocolos antigos, templates científicos personalizados e sistemas especiais que tenham preparação legada validada. Não é escolhido automaticamente quando faltam AmberTools, Meeko/Gemmi ou uma definição química confiável. Esses casos produzem diagnóstico claro.

## Operação, metodologia e mudança científica

| Categoria | Substituição | Interpretação |
| --- | --- | --- |
| STRUCTURE_SELECTION | Gemmi/Python, modelo e altloc explícitos | Equivalência de operação: átomos/identificadores/coordenadas preservados |
| CRYSTAL_LIGAND_PREPARATION | CCD ou template SDF/MOL2 + mapeamento por nomes/elementos | A extração não desloca átomos; não inventa bond orders |
| RECEPTOR_PREPARATION | LEaP, ff14SB e estados de protonação registrados | Campo de força preservado; H, átomos reparados e decisões automáticas são mudanças documentadas |
| CHARGE_ASSIGNMENT | Antechamber Gasteiger/AM1-BCC; RDKit onde permitido no Vina sem minimização | Métodos/implementações distintos não compartilham cache indevidamente |
| MINIMIZATION | sander, 100 SD + 10 CG, ff14SB/GAFF, vácuo | Mudança MMTK → sander controlada, sem exigir identidade numérica |
| FORMAT_CONVERSION | Adaptadores Meeko/Open Babel em processo | Nenhum processo externo apenas para seleção/exclusão de átomos |
| LEGACY_FALLBACK | Templates Chimera somente por escolha explícita | Compatibilidade excepcional, sem reintroduzir Chimera no padrão |

O módulo histórico caad.redocking delega os adaptadores de preparação ao módulo compartilhado caad.docking. A recuperação histórica de poses segue o núcleo nativo quando há grafo Meeko explícito. PDBQT sem química confiável exige um arquivo químico correspondente.

## Integridade e casos especiais

Metais/coordenados, ligantes covalentes, radicais, resíduos sem parâmetros e estados de carga não suportados não são forçados pelo caminho geral. Cofatores orgânicos são tratados separadamente. Histidinas/variantes protonadas explícitas são conservadas; escolhas automáticas de pH e tautômero são registradas e não equivalem a uma previsão de pKa local.

Também foi corrigida uma perda de cargas no template legado: escrever PDB e reabri-lo antes do MOL2 descartava as cargas atribuídas. O MOL2 agora é escrito diretamente da molécula carregada. Projetos que dependiam acidentalmente de cargas zero precisam ter seus resultados reavaliados, não adotados como referência química correta.

## Validação

O sistema DOCK6 1LGU/AAH é comparado com Chimera para seleção, receptor e redocking Vina; há scores contínuos pareados e execução real de superfície/esferas/caixa/grades/docking flexível DOCK6. A minimização mede energia, geometria, RMSD e clashes, sem alegar convergência ou equilibração. Cache frio/quente e ausência de chamadas Chimera no padrão são verificados.

Resultados e limites: [dados reais](validation/native_structures_2026-10-23.json). Guia de instalação/configuração: [preparação molecular](molecular_preparation.md). Os [relatórios anteriores](validation/chimera_remaining_2026-10-23.md) descrevem o estado anterior à migração principal e não a dependência atual.

Uma diminuição de energia intramolecular não garante melhoria do score de uma pose fixa. A minimização separada de receptor e ligante a vácuo altera a interface; diferenças eletrostáticas, vdW e poses precisam de análise independente. Os testes atuais não demonstram ausência de regressão sistemática para todos os alvos.

Fontes primárias: [CCD wwPDB](https://www.wwpdb.org/data/ccd), [sander Amber](https://ambermd.org/tutorials/basic/tutorial5/index.php). O bloco de dinâmica continua independente de Chimera e mantém seu próprio protocolo DOCK6/AmberTools.
