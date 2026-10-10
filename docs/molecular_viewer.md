# Visualizador molecular: poses e interações em 3D

[Documentação](README.md) · Português · [English](en/molecular_viewer.md)

O explorador WebGL local permite inspecionar estruturas, comparar a pose calculada
com seu receptor e identificar interações geométricas. Os arquivos moleculares
permanecem na aplicação; a visualização não envia estruturas a serviços externos.

## Abrir a estrutura ou o resultado

| Origem | Ação e conteúdo |
| --- | --- |
| PDB recuperado ou arquivo molecular | **Visualizar estrutura 3D** abre o arquivo selecionado |
| Tabela de compostos | **3D** gera um conformero local do SMILES; não é uma pose de docking |
| Redocking concluído | **3D** sobrepõe a melhor pose Vina, receptor e referência disponível |
| Docking Vina ou DOCK6 | **3D** sobrepõe a melhor pose e seu receptor; uma referência preparada pode acompanhar o resultado |
| Consenso | **3D Vina** e **3D DOCK6** abrem as poses de cada motor quando seus arquivos e receptores estão disponíveis |

PDB, PDBQT, MOL2 e SDF são aceitos pelo componente de estruturas. Na abertura de
um arquivo com vários modelos, **Modelo** permite escolher o modelo. Na cena de
resultados, a pose é selecionada automaticamente pelo menor score; sem score,
usa-se a primeira pose, identificada como tal. Essa seleção não oferece o campo
**Modelo**. Para consultar outras poses, abra o arquivo correspondente da lista
**Ver simulação** ou baixe-o.

O consenso pode calcular scores sem anexos estruturais. Para visualizar uma pose
com seu receptor, ambos precisam continuar disponíveis e associados ao resultado.
Arquivos ausentes ou ambíguos são informados; não se escolhe outro receptor.

## Navegar e escolher os estilos

Arraste para girar, use a roda do mouse para zoom e o botão direito para deslocar.
**Recentrar** recupera a visão do conjunto; **Tela inteira** amplia a área e
**Salvar imagem** exporta a cena como PNG.

1. Em **Estilo do receptor**, escolha fitas, bastões, esferas ou linhas. Para receptores MOL2, a cena de docking inicia em bastões.
2. Em **Estilo do ligante**, escolha bastões, esferas ou linhas, independentemente do receptor. A opção também controla os ligantes de referência visíveis.
3. Use **Colorir por** e **Cadeia** para explorar o receptor. A pose aparece em ciano e a referência em dourado.
4. Marque **Hidrogênios do ligante** para mostrar os H presentes no arquivo; desmarque para ocultá-los. O controle não gera átomos nem muda a preparação ou os scores.
5. Use **Ligantes**, **Água** e as caixas de camadas para escolher os elementos exibidos. Clique em um átomo para consultar seus detalhes.

Vina e DOCK6 podem preservar conjuntos diferentes de hidrogênios. H ausentes na
saída Vina não aparecem ao marcar o campo. Ocultar H é uma escolha de apresentação:
a classificação continua usando os átomos disponíveis no arquivo original.

## Comparar pose, receptor e referência

As camadas identificam **Receptor utilizado**, **Melhor pose** ou **Primeira pose
sem score**, e a referência disponível. No redocking, o ligante cristalográfico
aparece em dourado; o restante do complexo cristalográfico pode ser ativado como
camada de comparação. Sem o PDB original, uma referência preparada disponível
aparece com esse nome.

As coordenadas originais são preservadas, sem realinhamento automático. A tabela
**Resíduos** compara as menores distâncias entre átomos pesados da pose e da
referência; pode ser baixada em CSV. No visualizador, **Resíduos próximos** permite
selecionar e destacar um resíduo. O limite inicial é 4 Å, ajustável entre 2 e 8 Å.
Essa distância filtra a lista de proximidade, não os critérios da classificação
química das interações.

## Identificar as interações

Quando a topologia química e os resíduos permitem a classificação, o painel
**Mostrar interações em 3D** oferece filtros por tipo. Cada relação recebe um
traçado tracejado entre os átomos envolvidos, ou entre os centros dos anéis para
π–π. Traçados com os mesmos extremos recebem um pequeno desvio visual para
permitir sua identificação; os extremos e a distância química permanecem iguais.

| Cor da legenda e do traçado | Tipo identificado |
| --- | --- |
| Verde | Ligação de hidrogênio |
| Rosa | π–π paralelo |
| Coral | π–π em T |
| Lilás | Contato hidrofóbico |
| Verde-lima | Contato van der Waals geométrico |

Passe o mouse sobre um segmento tracejado: o tooltip mostra **tipo, resíduo,
cadeia e distância em Å**, com borda na mesma cor da relação. A legenda mantém as
cores mesmo quando um tipo não ocorre na pose. A lista de interações mostra os
mesmos dados; clique em uma linha para destacar e aproximar o resíduo.

Os filtros permitem ocultar tipos individuais. Ocultar a pose, o receptor, os
ligantes ou todas as interações também remove os traçados e o tooltip. No modo
por fitas, os resíduos envolvidos aparecem adicionalmente em bastões. As
interações classificadas descrevem a pose selecionada, não a referência de
comparação.

## Interpretar os resultados

A classificação é geométrica e limitada aos cinco tipos da legenda. Usa química
do ligante e posições originais; não calcula energias, afinidade nem comprova
ligação experimental. Ligações de hidrogênio exigem H explícitos e orientação
adequada; anéis aromáticos exigem topologia válida. Poses PDBQT dependem do SMILES
associado para recuperar ordens de ligação. Sem essas informações, a tela explica
a indisponibilidade e mantém os contatos por distância. Ausência de uma relação
na cena não prova que ela não exista.

Os limites de hidrofobicidade, ligação de hidrogênio e π–π seguem os
[parâmetros do PLIP](https://github.com/pharmai/plip/blob/master/plip/basic/config.py),
mas a ferramenta não executa o classificador PLIP completo. O contato van der
Waals usa proximidade em relação aos raios atômicos do RDKit. Para contatos
hidrofóbicos e van der Waals, exibe-se o par mais próximo por resíduo.

**Footprint** é outra análise: em resultados DOCK6 e no consenso com um PDF
associado, apresenta energias de van der Waals e eletrostáticas por resíduo,
comparando a referência com a pose final. O contato van der Waals desenhado em
3D não é essa energia. Consulte o [manual de docking e consenso](user_manual.md#vina-dock6-e-consenso).

## Acesso, requisitos e atualização

Leitores podem consultar e baixar resultados. A autorização é conferida ao abrir
as estruturas e buscar a cena. O visualizador exige navegador com WebGL e
aceleração gráfica; o limite de abertura é 32 MB. No desktop, abre uma página
local; na versão web, usa a mesma origem da aplicação.

Após atualizar o código, reinicie o aplicativo e abra uma nova janela do
visualizador. Janelas já abertas mantêm os scripts anteriores. As URLs dos
arquivos de interface são versionadas pelo conteúdo. Se houver falha, volte ao
projeto e abra o visualizador novamente.

Os testes de geometria e os movimentos reais do mouse no Chrome estão descritos
na [validação do visualizador](validation/viewer_controls_2026-10-10.md).
