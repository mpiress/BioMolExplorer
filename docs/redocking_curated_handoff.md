# Curadoria e configuração do redocking

A escolha das estruturas transporta os pares receptor/ligante para a configuração
principal. Os resíduos e cadeias salvos pela curadoria têm prioridade sobre a
cópia da configuração da execução que estava aguardando entradas. Na ausência
de uma revisão, são preservados os pares explicitamente configurados ou usados
os metadados das estruturas escolhidas. Estruturas desmarcadas não deixam pares
nem opções de preparo órfãs. Uma revisão sem ligantes permanece vazia.

A configuração principal apresenta os pares já selecionados e não exige sua
seleção novamente. A seleção manual continua disponível quando ainda não há
pares associados às entradas. A validação de cadeia, resíduo, cofatores e
arquivos continua obrigatória antes de executar.

O preparo possui duas seções gerais: **Preparação de todos os receptores** e
**Preparação de todos os ligantes**. Hidrogênios, solvente, cargas e minimização
são aplicados a todos os pares. Identidades estruturais, cadeias especiais e
cofatores permanecem vinculados aos pares. Para configurações antigas com
opções diferentes, a interface informa que parte do primeiro par; a confirmação
aplica os valores gerais exibidos a todo o conjunto.

As opções continuam serializadas em `preparation_pairs` para compatibilidade com
os motores científicos e suas assinaturas de cache. Resultados existentes só
são reutilizados quando suas entradas e opções de preparo continuam válidas.

## Confirmação do formulário

Ao aplicar o redocking, converter uma conexão `auto` em referências explícitas
às estruturas não representa uma desconexão. O manipulador do botão preserva
os pares e o preparo devolvidos pelo formulário; não restaura os parâmetros
padrão antes de salvar ou chamar `PipelineService.resume`. A reconciliação das
estruturas removidas permanece a cargo da seleção curada e da validação.

Os testes em `test_redocking_apply.py` exercitam o botão de aplicação até a
fronteira do serviço e a gravação no projeto, incluindo a troca de referências
automáticas por seletores e identidades de artefatos.

Ao concluir a configuração, toda a sequência de janelas de seleção, curadoria e
configuração daquele bloco é encerrada. Seus botões antigos são invalidados e
o rascunho da seleção é liberado. Falhas de validação mantêm a configuração
aberta para correção; concluir o fluxo não fecha janelas de outros blocos.
