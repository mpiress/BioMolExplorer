# Controles e interações do visualizador — 10/10/2026

O visualizador de docking agora utiliza os campos de estilo do ligante e
hidrogênios já presentes no HTML. O campo de representação principal controla
apenas o receptor. Os mesmos controles funcionam na abertura de estruturas
isoladas. Hidrogênios ausentes no arquivo não são gerados.

As interações da pose selecionada são desenhadas como linhas tracejadas nas
coordenadas originais, com etiquetas e uma lista clicável de resíduos, tipos e
distâncias. Há filtros por tipo e um controle geral de visibilidade. Ocultar
receptor, pose ou ligantes também oculta as interações. A representação por fitas
mostra em bastões os resíduos envolvidos. Referências cristalográficas/preparadas
continuam sendo camadas de comparação; a classificação descreve apenas a pose.

A classificação implementada é geométrica e limitada: π–π paralelo/em T,
ligação de hidrogênio com H explícito, contato hidrofóbico e contato van der Waals.
Os limites de π–π, hidrofobicidade e ligação de hidrogênio seguem os
[parâmetros publicados do PLIP](https://github.com/pharmai/plip/blob/master/plip/basic/config.py).
Esta implementação não executa o PLIP completo. O contato van der Waals usa os
raios do RDKit, com separação entre superfícies de −0,4 a +0,5 Å; não mede energia
nem demonstra afinidade. Mantém-se o par mais próximo por resíduo para contatos
hidrofóbicos e van der Waals. Uma classificação ausente não prova ausência da
interação: ligações de hidrogênio exigem H explícitos e anéis exigem topologia
química válida. PDBQT depende do SMILES associado para recuperar ordens de ligação.
Quando esses dados são insuficientes, a tela informa a limitação e mantém os
contatos por distância disponíveis.

Validação: testes de geometria com anéis paralelos/perpendiculares, ligação de
hidrogênio com orientação correta/incorreta, exclusão de sobreposição e distância
excessiva, topologia Vina com H explícito e topologia incompatível. Também foram
executados os testes de cenas/autorização, visualizador, resultados de redocking
e visualizações de resultados.

Validação WebGL com Chrome headless e duas cenas sintéticas equivalentes
(benzeno/PHE, MOL2 DOCK6 e PDBQT Vina): estilos independentes, ocultação/exibição de
H, três interações 3D, filtros por tipo, visibilidade, rotação, zoom e exportação
PNG. Nenhuma requisição externa nem erro JavaScript observado. Essas cenas
validam o fluxo e a geometria controlada, não a qualidade científica de uma
simulação real.

Reiniciar o aplicativo e abrir uma nova janela do visualizador carrega os novos
scripts; as URLs dos assets já são versionadas pelo conteúdo.

## Traçados e identificação ao passar o mouse

Cada relação usa um traçado 3D de segmentos cilíndricos finos, com intervalos,
para manter espessura visível e permitir detecção do mouse na geometria. Relações
com os mesmos extremos recebem um pequeno desvio no ponto intermediário, sem
mover os extremos moleculares. A distância informada continua sendo a distância
química original, não o comprimento do traçado de apresentação.

Paleta: ligação de hidrogênio verde; π–π paralelo rosa; π–π em T coral;
contato hidrofóbico lilás; contato van der Waals verde-lima. A legenda mostra uma
amostra tracejada de cada cor, inclusive para tipos ausentes na pose. O tooltip
usa a mesma cor na borda e mostra tipo, resíduo, cadeia e distância. Ele desaparece
ao sair do traçado ou ocultar/filtrar as interações.

O validador WebGL projeta os centros dos segmentos na tela e envia movimentos
reais do mouse para verificar a identificação de cada tipo presente e o
fechamento do tooltip ao sair. A implementação usa os callbacks de hover das
[formas do 3Dmol](https://3dmol.org/doc/ShapeSpec.html), conferidos também na versão
embarcada no projeto.
