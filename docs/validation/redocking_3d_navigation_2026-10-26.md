# Redocking: abertura da visualização 3D — 2026-10-26

## Causa e correção

A abertura da aba ocorria após a preparação assíncrona da cena no servidor.
Nesse momento, o navegador podia bloquear a abertura por não haver mais um gesto
ativo do usuário. O lançador de URLs do Flet instalado não devolvia esse bloqueio
ao Python, deixando apenas o splash como feedback.

Os botões web agora possuem uma URL direta e abrem a aba no próprio clique.
A página do visualizador carrega primeiro; os dados moleculares são preparados
na requisição da cena, com nova validação de sessão, projeto, resultado e arquivos.
O mesmo ajuste foi aplicado aos botões de docking e comparação de resíduos.
A aplicação desktop conserva seu caminho de abertura pelo navegador externo.

## Validação

- Uma cópia temporária dos resultados existentes continha 35 pares. As 35
  requisições de cena retornaram HTTP 200, sem executar novamente o redocking.
  A cópia omitiu referências de proveniência a execuções históricas não copiadas;
  arquivos moleculares e poses foram preservados. O projeto original não foi alterado.
- No Chrome, com bloqueio de popups habilitado, um clique real no botão corrigido
  abriu a aba em aproximadamente 0,07 s. A cena renderizou receptor, pose e
  referência, respondeu à rotação e não apresentou erros JavaScript.
- O fluxo anterior foi reproduzido com atraso controlado de seis segundos na
  preparação: a aba foi bloqueada. Esse atraso foi uma condição de teste,
  não uma medição do tempo de preparação do usuário.
- Testes automatizados verificam navegação direta, carregamento adiado,
  autorização, sessão encerrada, arquivo removido, projeto trocado e manutenção
  do callback desktop. A renderização WebGL foi verificada em um dos 35 pares;
  os demais tiveram seus dados de cena validados.
- A suíte focada passou com 41 testes. Os 34 documentos HTML passaram na
  validação de arquivos locais, links, âncoras e identificadores.

Não houve mudanças nos cálculos científicos de docking ou redocking.
