# Configuração do redocking

[Documentação](README.md) · Português · [English](en/redocking_configuration.md)

## Selecionar e preparar os pares

Após recuperar os PDBs, selecione os arquivos de entrada e abra a configuração do redocking. Em **Entradas**, escolha explicitamente cada par PDB / ligante / resíduo. Nenhum par é selecionado automaticamente. É obrigatório selecionar pelo menos um antes de aplicar a configuração.

Cada seleção mostra a cadeia comum ao receptor e ao ligante, os cofatores e as opções de preparação. A resolução acompanha os metadados do par e não precisa ser digitada. Para PDBs próprios sem metadados, use **Informar par manualmente**.

Ative **Considerar cofatores como parte do receptor** e informe os códigos, separados por vírgula (por exemplo, `FAD, MG`). Sem cofatores, a seleção do complexo usa `select #0:.{chain}`. Os cofatores declarados também são preservados na preparação do receptor, inclusive fora da cadeia selecionada. O ligante de redocking não pode ser declarado como cofator.

Configure remoção de solvente, remoção dos hidrogênios existentes, adição de hidrogênios, minimização e método de cargas para receptor e ligante. A configuração do ligante é compartilhada entre preparação e conformação. O pH continua sendo definido uma vez para o bloco.

Pares que usam o mesmo PDB e a mesma cadeia compartilham o arquivo de receptor preparado; devem usar os mesmos cofatores e as mesmas opções do receptor. As opções de seus ligantes podem ser diferentes. No processamento individual, cada execução recebe somente os pares e as configurações da estrutura correspondente.

Quando os pares já foram configurados, o pipeline usa essa seleção sem solicitar novamente. Antes dos cálculos, valida a presença do receptor, ligante, resíduo, cadeia e cofatores na estrutura. Complexos preparados devem conter os arquivos PDBQT e centros válidos. Chimera, Open Babel e Vina precisam estar instalados no ambiente de execução; com entradas já preparadas, somente Vina é exigido. Verbosidade permanece desativada.

## Configuração pelo backend

Consulte o executável e a mensagem no log de execução. Scripts Chimera que falharam ficam preservados. Centro do ligante inválido ou arquivos preparados incompletos impedem o docking; corrija a configuração ou instalação e tente novamente.

Informe, por exemplo, `pdb_codes=[["4M0E", "1YL", 604, "A", 2.0]]`. `preparation_pairs` usa a chave `4M0E|1YL|604|A`, contendo `cofactors`, `receptor` e `ligand`. As opções permitidas são `remove_solvent`, `remove_hydrogens`, `add_hydrogens`, `minimize` e `charge_type` (`gas` ou `am1`). Opções omitidas usam os padrões. O quinto elemento do registro é a resolução. Projetos legados com `ligand_chain` passam a usar a cadeia comum `CHAIN`.

Veja [validação do pipeline](pipeline_validation.md) para cobertura e limites da execução nativa.

A seleção de um subconjunto não remove PDBs nem arquivos preparados de outras estruturas da coleção. Quando o par não define seu método de cargas, a preparação do ligante conserva o método global do bloco.

## Falhas de preparação com Chimera clássico

Um processo Chimera pode terminar com código zero mesmo quando um comando `.com` falhou. O pipeline verifica essas mensagens e os arquivos esperados antes de avançar. Cada saída deve existir, estar preenchida e conter átomos; o mesmo vale para a conversão pelo Open Babel. Saídas inválidas interrompem a etapa e preservam o script Chimera para diagnóstico.

Os caminhos de `open` e `write` são passados sem aspas de shell: no Chimera clássico, essas aspas fazem parte do nome do arquivo. Caminhos com espaços foram verificados com o executável real. A saída padrão e os erros das ferramentas são registrados nos logs científicos.

## Consultar resultados de redocking

Ao expandir uma etapa concluída em **Execuções**, a tabela mostra PDB, ligante, resíduo, cadeia e **RMSD (Å)**. **Ver simulação** abre um popup com os arquivos correspondentes àquela seleção. A lista distingue o receptor preparado, o ligante de referência, as poses e os metadados. Arquivos compartilhados pelo receptor ou pela coleção acompanham as simulações correspondentes.

Use o botão de download de cada linha para baixar um arquivo, ou **Baixar todos (ZIP)** para obter o conjunto da simulação. O ZIP conserva as pastas, inclusive quando o ligante de referência e a saída do Vina têm o mesmo nome. Arquivos de outras simulações ficam fora desse conjunto.

O botão **Visualizar estrutura 3D** abre o visualizador no navegador, como na consulta dos PDBs. Ele aceita PDB, PDBQT e MOL2 para examinar receptor, ligante e poses antes do download. Ao abrir um arquivo com várias poses, use **Modelo** no visualizador. A cena aberta pelo botão **3D** da tabela seleciona automaticamente a pose de menor score, ou a primeira sem score; essa cena não oferece troca de modelo. As legendas e ações acompanham o idioma selecionado. Leitores também podem consultar e baixar os resultados.

A tabela resume os valores finitos e não negativos da coluna `RMSD` de `pdb_codes.csv`, com três casas decimais na apresentação; o arquivo conserva a precisão original. O RMSD compara o ligante de referência com o resultado do redocking. A interface não aplica um limiar automático de aprovação científica. A identificação PDB / ligante / resíduo / cadeia permite distinguir cada simulação. Etapas com falha ou ainda em execução mantêm a lista geral de arquivos para diagnóstico.

[Uso do backend](backend_usage.md) inclui um exemplo completo do dicionário de parâmetros e a organização dos arquivos de saída.

## Reutilizar o receptor no docking de candidatos

Conecte o redocking a **Receptor PDB (retrieval ou redocking)** do bloco **Preparar para docking**. Escolha o receptor `.dockprep.pdbqt` no campo **Receptor que deseja utilizar**; arquivos exclusivos do ligante não aparecem nessa seleção. O receptor, os formatos complementares e os centros são reutilizados, com as opções de preparo do receptor desabilitadas. Os candidatos vêm de uma ou mais fontes ChEMBL, PubChem, ZINC ou arquivos próprios e continuam sendo preparados. Escolha Vina, DOCK6 ou ambos como saída e conecte receptor e compostos a cada motor. Consulte o [manual do usuário](user_manual.md#preparacao-e-redocking) para processamento individual e merge.

## Estilos e interações em 3D

A cena **3D** compara pose, receptor e referência nas coordenadas originais.
Escolha estilos independentes para ligante e receptor e use **Hidrogênios do
ligante** para mostrar apenas os H existentes no arquivo. Interações classificadas
recebem linhas tracejadas coloridas, filtros e identificação pelo mouse. PDBQT
exige SMILES associado para recuperar topologia; resultados de redocking sem essa
informação continuam oferecendo contatos por distância, com a indisponibilidade
química indicada. Consulte o [guia do visualizador](molecular_viewer.md) para a
legenda, requisitos e limites da classificação.
