# Revisão do editor de pipeline — 18/10/2026

A revisão encontrou problemas no editor compartilhado pelas operações, que também
atingem a inclusão e configuração de `retrieve_structures` (PDB). Adicionar esse
bloco não inicia a recuperação de estruturas: `FlowCanvas.add()` apenas cria a
configuração local e agenda sua persistência. Não foi identificada consulta ao
PDB ou renderização molecular nesse caminho.

## Causas e correções

1. **Cancelamento de uma gravação já iniciada.** `schedule_save()` cancelava a
   tarefa anterior mesmo depois de ela chamar `save_pipeline()` por
   `asyncio.to_thread()`. O cancelamento interrompe a espera, mas a thread pode
   concluir a gravação. A interface perde a resposta com a nova revisão. Uma
   edição seguinte do bloco recém-criado usa a base antiga e pode gerar o erro
   de outro colaborador, mesmo na mesma sessão. Agora apenas o temporizador de
   600 ms pode ser cancelado; as gravações iniciadas concluem sua confirmação,
   serializadas pelo lock existente. O temporizador também pertence à sessão
   e ao projeto que o criou.
2. **Diálogo invalidado pelo próprio salvamento.** A confirmação de salvamento
   substituía `self.current`. O diálogo mantinha a referência anterior e sua
   proteção contra troca de projeto impedia aplicar a configuração. Agora a
   identidade do projeto e dos blocos inalterados é preservada.
3. **Trabalho excessivo em cada movimento.** O arraste recriava a camada de
   conexões e as opções do seletor; chamava `page.update()` para a página toda,
   percorrendo e traduzindo novamente seus controles. Agora atualiza apenas o
   bloco movido e a camada persistente de conexões. A grade só é substituída
   quando o tamanho do canvas muda, as opções persistem durante o arraste e
   a seleção altera as bordas sem reconstruir os blocos. Eventos de movimento
   têm intervalo de 32 ms. O cálculo de conexões usa posições indexadas por ID,
   eliminando buscas lineares por origem.
4. **Tradução repetida.** Textos sem correspondência direta voltavam a percorrer
   os padrões regulares em toda atualização. Um cache por tradutor/sessão,
   limitado a 4096 entradas, reutiliza os resultados e é limpo ao trocar idioma.
5. **Reconstrução desnecessária após salvar.** Confirmações que não modificam
   o pipeline deixam o canvas existente montado. Alterações recebidas do
   armazenamento continuam a provocar atualização.
6. **Edições durante a gravação e respostas atrasadas.** A confirmação agora
   combina alterações locais posteriores ao envio com alterações remotas já
   incorporadas no salvamento. A consulta periódica e o redesenho não aplicam
   respostas anteriores à revisão/data já confirmada na interface. Conflitos
   reais no mesmo campo continuam sujeitos à validação de colaboração.

A inspeção também verificou o formulário PDB, a descoberta de artefatos, a
sincronização periódica e os contratos de portas. `catalog.operation_fields()`
já tinha cache: a leitura de assinaturas Python não foi tratada como causa da
lentidão recorrente. A tela ainda consulta até 100 execuções para descobrir
artefatos ao ser aberta; isso pode custar mais em projetos com histórico grande,
mas não é executado em cada movimento. Não foi introduzido cache de permissões
nem de resultados científicos.

## Medição reproduzível

```bash
PYTHONPATH=src python scripts/benchmark_pipeline_editor.py \
  --baseline 10091e3b28a1943e3bf24c16cb6bcdd6882376ff
```

O benchmark compara o canvas e a localização da revisão indicada com os arquivos
atuais, usando os demais módulos atuais em ambos os casos. São 50 eventos de
movimento, após aquecimento, com blocos PDB sem conexões e página simulada.
Mede a mediana do processamento Python por evento. Não mede transporte,
serialização real do Flet, pintura Flutter/GPU, memória ou FPS do navegador.
Os números dependem da máquina. A revisão foi feita com Flet 1.0.3.

| Blocos PDB | Antes (ms/evento) | Depois (ms/evento) |
| --- | ---: | ---: |
| 1 | 1,357 | 0,143 |
| 20 | 13,372 | 0,293 |
| 80 | 53,690 | 0,777 |

## Regressões cobertas

`tests/test_pipeline_editor_regressions.py` verifica a gravação interrompida por
nova edição, aplicação do diálogo PDB aberto durante o salvamento, preservação
de alterações remotas e locais, rejeição de uma resposta periódica antiga,
atualização restrita durante arraste e invalidação do cache de idioma.
As suítes existentes exercitam permissões, colaboração, conexão por arraste,
formulários, navegação, zoom e preservação das entradas.

Validação: **90 testes passaram** com o comando abaixo; `git diff --check`
também passou.

```bash
PYTHONPATH=src:tests python -m unittest \
  test_pipeline_editor_regressions test_canvas_drag_regression \
  test_ui_permissions test_localization test_pipeline_guidance test_guided_ui \
  test_ui_layout test_frontend_navigation test_zoom_visualizations \
  test_project_portability test_guided_input_preservation -q
```

A validação realizada usa controles Flet e armazenamento reais, com página e
transporte simulados. Falta confirmar visualmente a fluidez e o consumo de
memória no aplicativo em execução com o projeto afetado.
