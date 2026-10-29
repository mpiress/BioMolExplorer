# Clear All e recuperação dos pares de redocking — 2026-10-27

## Escopo e dados

A implementação adiciona **Clear All / Limpar tudo** ao lado de Refresh Results.
Por escolha do usuário, preserva as conexões do pipeline e limpa as seleções de
arquivos e pares. Foram examinados, em modo somente leitura, os registros e os
artefatos da execução `fcd0aef60bd84a39a7ac7208c7d73351` do projeto Alzheimer.
O redocking dessa execução registrava 24 falhas de preparo entre 35 pares. O banco
e os arquivos originais não foram alterados pelos testes. Nenhuma biblioteca foi
instalada.

## Clear All

- Limpeza auditada em transação: remove resultados calculados publicados e parciais,
  entradas materializadas e seleções derivadas de todas as execuções do projeto.
- Mantém arquivos importados/providos, conexões, parâmetros científicos comuns,
  consultas de recuperação, cache de definições CCD, logs e histórico.
- Invalida manifestos, lotes e chaves de reutilização dos resultados removidos;
  a interface mostra **Results cleared**, sem apresentar os blocos limpos como
  concluídos ou reutilizados. As datas e a duração histórica da execução permanecem.
- Bloqueia limpeza durante execução ativa ou por colaboradores sem edição.
- Confirmação com progresso e bloqueio de cliques repetidos; restaura arquivos,
  seleções e estado do banco se a transação falhar.

Chrome/Flet reais: inglês em 1200 px e português em 600 px. Cancelar preservou os
resultados; confirmar removeu os calculados, preservou os importados e as conexões,
resetou as seleções e atualizou a tela. Nenhum erro JavaScript foi observado.
Evidências: [relatório](clear_all_redocking_2026-10-27/ui-report.json),
[confirmação em tela estreita](clear_all_redocking_2026-10-27/pt-confirmation.png),
[progresso](clear_all_redocking_2026-10-27/en-progress.png),
[resultados limpos](clear_all_redocking_2026-10-27/en-cleared.png).
Uma verificação adicional em 360 px também passou:
[relatório](clear_all_redocking_2026-10-27/ui-360-report.json) e
[confirmação](clear_all_redocking_2026-10-27/pt-360-confirmation.png).

## Causas corrigidas no preparo científico

1. Oxigênios de hidroxilas conjugadas podiam ser tipados como carbonilas por usar
   hibridização em vez da ordem real de ligação. O Chimera deixava de adicionar
   hidrogênios necessários, causando falhas de AM1-BCC. BUA reproduziu o erro.
2. Nitrogênios aromáticos de pirrol substituídos precisavam do tipo nativo `Npl`,
   em vez de `N2`. S6Q reproduziu o erro.
3. A heurística de carga líquida do Chimera não correspondia ao grafo observado
   de alguns ligantes. ISP/TC3/DEP/TN7/TC5 incluem fragmentos fosforilados; 0GS
   reproduziu uma atribuição incorreta de +2 em um grafo neutro. A carga agora é
   derivada do grafo validado após a adição nativa de hidrogênios e informada
   explicitamente à mesma rotina nativa AM1-BCC/Gasteiger escolhida pelo usuário.
4. `minimize` executava DockPrep novamente, podendo adicionar hidrogênios e
   substituir cargas. Agora minimiza o grafo já preparado, com `prep false`, mantendo
   o campo de força, o método escolhido e a minimização. Há validação de identidade
   química e soma de cargas após o preparo, com relatório `.lig.preparation.json`.
5. O receptor 2XQF incluía uma glicina livre de cristalização (HETATM GLY1531)
   como proteína, provocando erro MMTK em H1. Aminoácidos HETATM isolados, sem
   ligação peptídica e não selecionados explicitamente como cofatores, são removidos.
   Resíduos peptídicos e cofatores declarados permanecem.
6. O erro nativo FindHBond `accPhiPsi() takes exactly 8 arguments (9 given)`
   impedia QK7 e ZQS. Apenas nesse erro há uma repetição com AddH geométrico nativo;
   cargas e minimização permanecem. A recuperação fica registrada.
7. RDKit `RemoveHs` conserva alguns hidrogênios de valência não convencional.
   O RMSD de átomos pesados passa a usar `RemoveAllHs`, evitando falhas falsas
   de mapeamento em fragmentos fosforilados.

Não há troca automática de AM1-BCC por Gasteiger, remoção silenciosa da minimização
ou aceitação de topologia final sem sanitização. RMSD alto permanece um resultado
válido; não é convertido em falha. Erros reais conservam diagnóstico e isolamento
por par. A versão do preparo e os hashes do código invalidam caches incompatíveis.

## Validação automatizada

- 137 testes passaram na suíte de preparo, configuração, falhas parciais,
  handoffs, regressão do pipeline, interface, progresso, splash e traduções.
- Os 6 testes finais de Clear All passaram, incluindo um caso adicional de pares
  sem conexão, datas históricas, permissões, execução ativa e rollback.
- Outros 25 testes de cache, revisão de execução e falhas de preparo passaram.
- Os testes de carga incluem grafos reais ISP (+1) e 0GS (0), métodos am1/gas,
  manutenção de minimização, rejeição de soma de cargas incompatível e RMSD sem H.
- A compilação gerou 28 páginas de documentação e dois portais de idioma.

O caminho nativo Gasteiger também foi executado em 1P0I/BUA, mantendo minimização:
Vina e RMSD concluíram (1,779 Å), com grafo preservado, carga líquida 0 e soma de
cargas parciais 0,0001 (arredondamento do MOL2).
[Registro](clear_all_redocking_2026-10-27/vina-gasteiger.json).

## Limites de interpretação

Os testes de Vina usam os 35 pares já selecionados, os arquivos cristalográficos e
as opções salvas (receptor Gasteiger, ligante AM1-BCC, hidrogênios e minimização
habilitados, caixa 24 Å, exhaustiveness 20, 10 poses, seed 73). Receptores previamente
preparados válidos são reutilizados; o 2XQF requer novo preparo completo. As cópias
são independentes do projeto original. Fragmentos observados continuam identificados
como tais na proveniência CCD: executar docking não valida um mecanismo covalente,
um estado de protonação fisiológico ou afinidade experimental. Esta é validação
computacional dos casos configurados, não uma certificação de todos os PDBs.

## Resultado final dos testes nativos

**35/35 pares selecionados executaram com AutoDock Vina**, sem falhas de preparo
ou avaliação na rodada final. Todos os grafos hidrogenados e preparados foram
comparados e as somas de cargas validadas. Os RMSDs vão de 0,210 a 9,315 Å;
valores altos são mantidos como resultado de avaliação, sem falsa exclusão do par.
O receptor 2XQF passou por preparo e minimização completos novos (aproximadamente
13 minutos na máquina durante os testes concorrentes); seu ligante foi então
reprocessado pelo fluxo final, com RMSD 9,170 Å. Os outros receptores válidos foram
reutilizados. A preparação e o Vina de cada ligante foram realmente executados.

- [35 pares, opções, validações e hashes](clear_all_redocking_2026-10-27/vina-35-pairs.json)
- [Tabela CSV dos 35 pares](clear_all_redocking_2026-10-27/vina-35-pairs.csv)

**DOCK6 6.11: 2/2 pares executaram**, com busca flexível, 1000 orientações e
preparação completa de superfície, esferas e grid: 1P0I/BUA (RMSD 10,871 Å,
score −22,416943) e 6NTO/L2Y (RMSD 4,506 Å, score −41,664238). Os RMSDs altos
não são apresentados como execução impossível. A compatibilidade com DOCK6 foi
verificada nesses dois casos, não em todos os 35 pares.
[Registro DOCK6](clear_all_redocking_2026-10-27/dock6-two-pairs.json).

Após a revisão final de proveniência, **39 testes de cache, reutilização,
recuperação e falhas passaram**. O cache inclui agora os relatórios JSON do
preparo, com hashes, preservando o método e as validações em pares reutilizados.
[Log](clear_all_redocking_2026-10-27/cache-provenance-tests.log).

Versões medidas: Chimera clássico 1.17.3, AutoDock Vina `36dd023-mod`,
Open Babel 3.1.0 e DOCK6 6.11. Resultados foram obtidos em cópias temporárias,
sem limpar ou reexecutar o projeto original. Reinicie o processo da aplicação
para carregar o código atualizado.
