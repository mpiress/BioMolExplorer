# Preparação molecular para docking e redocking

## Ferramenta utilizada

UCSF Chimera é a única ferramenta de preparação nos blocos Preparar para docking,
Redocking, Docking Vina e Docking DOCK6. As opções históricas `auto` e `python`
são normalizadas para `legacy_chimera`, inclusive em projetos já salvos e na
retomada de execuções. Não existe seleção de AmberTools ou Meeko nesses blocos.
ChimeraX não substitui o executável `chimera` nesses scripts.

## Configuração

Uma seção define o preparo de todos os receptores e outra define o de todos os
ligantes: remoção de solvente/hidrogênios, adição de hidrogênios, minimização e
método de cargas. Identidade do par, cadeia e cofatores continuam vinculados às
estruturas selecionadas. As opções não são apagadas durante a migração.

- `gas`: cargas Gasteiger pelo Chimera.
- `am1`: cargas AM1-BCC pelo Chimera.
- Minimização: comando de minimização do Chimera.
- Proteínas padrão usam os templates de cargas do Chimera; a escolha gas/am1
  controla os componentes não padrão, incluindo ligantes e cofatores.

Chimera pode usar utilitários fornecidos em sua própria distribuição para
AM1-BCC. O aplicativo não encaminha docking/redocking ao ambiente separado
BioMolExplorer-AmberTools. Esse ambiente continua disponível para os protocolos
de dinâmica que o exigem.

## Entrada e exportação

Moléculas de ChEMBL, PubChem, ZINC e arquivos próprios usam o mesmo preparo.
RDKit pode gerar uma conformação inicial de SMILES; Open Babel faz conversões
de formato. A atribuição final de cargas e a minimização ficam no Chimera.
MOL2 e PDBQT são exportados da conformação preparada, preservando o método de
cargas escolhido. Vina e DOCK6 continuam como motores de docking.

Redocking valida PDB, ligante, resíduo, cadeia e cofatores antes dos cálculos.
Chimera/Open Babel precisam estar disponíveis no ambiente do worker. No redocking,
escolha Vina ou DOCK6; apenas o motor escolhido é exigido, além das ferramentas
de preparação. A instalação DOCK6 é definida pelo administrador.
O diretório do interpretador científico é incluído no PATH do worker.

## Reutilização e diagnóstico

Receptores explicitamente fornecidos como preparados podem ser reutilizados.
Cargas, pH, minimização, estrutura, ferramentas e parâmetros participam das
assinaturas de cache. Preparos nativos antigos não são reclassificados como
Chimera para evitar reutilização de resultados incompatíveis. Falhas de pares
continuam isoladas e registradas; falta de ferramenta gera diagnóstico explícito.

Antes de AddH/AddCharge, o Chimera executa a limpeza estrutural do Dock Prep:
remove conformações alternativas de receptores e ligantes e reconstitui cadeias
laterais incompletas do receptor. A minimização usa `prep false` após essas
etapas, preservando as cargas escolhidas.

Definições CCD locais disponíveis orientam a hibridização de carbonos e a
tipagem de oximas/nitronas antes da protonação. Uma falha de cargas pode
consultar a definição do componente e repetir o preparo no Chimera, com o mesmo
método de cargas. O aplicativo não troca AM1 por gas para contornar falhas.

A exportação MOL2 pode corrigir ordens de ligação usando a definição do
componente, preservando coordenadas e cargas calculadas. Conformações
alternativas são resolvidas antes do preparo; átomos ausentes não são inventados.
Quando a referência contém apenas o fragmento observado, o arquivo interno
`.export.json` registra o escopo e os nomes dos átomos ausentes.
