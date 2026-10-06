# Logs e diagnóstico

[Documentação](README.md) · Português · [English](en/logging.md)

O sistema combina logs de texto para leitura direta, eventos JSONL para filtros e um resumo por job. Interface, pipeline, supervisor, worker e módulos científicos usam o mesmo formato. Mensagens técnicas e saídas dos executáveis conservam o idioma de origem; os controles e mensagens apresentados na interface continuam acompanhando o idioma da sessão.

## Localizar a causa de uma falha

1. Em **Execuções**, identifique a etapa que falhou e abra **Ver log da etapa**. O `execution.log` reúne a saída do worker, inclusive os eventos científicos e seus tracebacks.
2. Use o `job_id` da etapa para abrir `logs/jobs/<job_id>/diagnostic.json`. Consulte `status`, `error_code`, `error` e `action`.
3. No mesmo diretório, procure `tool.failed` em `events.jsonl`. O `command_id` liga início, saída e resultado de um comando; `configuration` identifica o script Chimera ou a configuração Vina e `cwd` informa onde ele foi executado.
4. Consulte o evento anterior `tool.output` para a mensagem original da ferramenta. Para falhas Python, `source` e `exception.traceback` mostram o local e a cadeia de exceções. Os eventos posteriores `worker.failed`, `job.failed` e `stage.failed` descrevem a propagação da mesma falha entre camadas.

Não existe um `tool.failed` para toda falha: erros de validação podem ocorrer antes de iniciar ferramentas. Nesse caso, consulte o evento de erro da etapa ou do worker e o caminho/registro indicado. Um pipeline que falha antes de criar um job registra a falha em `logs/backend.log` e `logs/events.jsonl`.

## Arquivos e retenção

O padrão é `logs/` na raiz do checkout; numa instalação sem checkout, é `logs/` no diretório de execução. Configure `BIOMOL_LOG_DIR` ou `biomolexplorer-ui --log-dir /caminho/logs` para mudar a raiz.

| Arquivo | Finalidade |
| --- | --- |
| `frontend.log` | Ações e falhas da interface |
| `backend.log` | Pipeline, supervisão e ciclo dos jobs |
| `errors.log` | Eventos ERROR/CRITICAL, com causa e traceback quando disponíveis |
| `events.jsonl` | Um objeto JSON por evento; inclui contexto e campos de diagnóstico |
| `jobs/<job_id>/diagnostic.json` | Resumo mais recente do job, gravado atomicamente |
| `jobs/<job_id>/` | Eventos e arquivos científicos exclusivos daquele worker |
| `jobs/<job_id>/execution.log` | Cópia do log de execução; o original continua no diretório privado do job |
| `docking.log`, `complex.log`, `bioactivities.log`, etc. | Arquivos científicos existentes, agora com o formato comum |

Os arquivos gerenciados de texto e JSONL rotacionam em aproximadamente 5 MiB, com três cópias anteriores (`.1` a `.3`). Um evento individual muito grande pode superar esse tamanho. A escrita e a rotação são sincronizadas entre threads e processos no ambiente Linux. Arquivos `.lock` são auxiliares dessa sincronização. O resumo é substituído a cada atualização; `execution.log` e os diretórios antigos de jobs não têm limpeza automática por idade.

Os diagnósticos incluem caminhos e identificadores do estudo. Credenciais não devem ser registradas: parâmetros completos e comandos completos não são serializados nos novos eventos de execução. O formatador também mascara padrões comuns como `token=...` e `password=...`; essa proteção não identifica todo conteúdo sensível possível. Revise os arquivos antes de compartilhá-los.

## Filtrar rapidamente

Execute no ambiente Python do projeto:

```bash
# Erros e avisos do processo principal
python -m biomolexplorer.log_report

# Erros de um worker; substitua JOB_ID pelo identificador real
python -m biomolexplorer.log_report --job JOB_ID --level ERROR

# Histórico de comandos, incluindo saídas e duração
python -m biomolexplorer.log_report --job JOB_ID --level INFO --limit 100

# Eventos do pipeline para uma execução/etapa
python -m biomolexplorer.log_report --run RUN_ID --stage STAGE_ID

# Saída estruturada para análise externa
python -m biomolexplorer.log_report --job JOB_ID --json
```

`--directory` seleciona uma raiz diferente de logs; com `--job`, o leitor procura a subpasta `jobs/<job_id>` nessa raiz. `--operation redocking` filtra pela operação. O leitor percorre os arquivos JSONL rotacionados e conserva apenas os últimos eventos correspondentes, sem carregar o histórico inteiro na memória. Logs antigos, anteriores ao JSONL, continuam disponíveis como texto.

## Códigos de diagnóstico

| Código | Significado e verificação |
| --- | --- |
| `TOOL_NOT_FOUND` | Executável não localizado; confira instalação e PATH do worker |
| `TOOL_START_FAILED` | Não foi possível iniciar; confira permissões e diretório |
| `TOOL_TIMEOUT` | Comando excedeu seu limite; confira a saída, entrada e `BIOMOL_COMMAND_TIMEOUT` |
| `TOOL_EXIT_FAILED` | Ferramenta terminou com código diferente de zero; consulte stdout/stderr |
| `TOOL_REPORTED_ERROR` | Chimera relatou erro mesmo retornando zero; confira o script preservado |
| `PREPARED_OUTPUT_INVALID` | Preparação produziu arquivo ausente, vazio ou sem átomos; confira a ferramenta e o script |
| `EXECUTION_TIMEOUT` | Job excedeu o limite do supervisor; confira o último comando |
| `INPUT_NOT_FOUND` | Arquivo de entrada não localizado |
| `ACCESS_DENIED` | Acesso negado pelo sistema de arquivos |
| `VALIDATION_FAILED` | Entrada ou configuração inválida |
| `WORKER_EXIT_FAILED` | Worker encerrou sem diagnóstico válido; consulte `execution.log` |
| `UNEXPECTED_ERROR` | Exceção sem categoria específica; consulte a causa e o código de origem |

Os códigos são estáveis para filtros. A orientação `action` indica a próxima verificação, sem decidir automaticamente a validade científica da entrada. Cancelamento pelo usuário é um estado próprio, não um erro de ferramenta.

## Contrato para desenvolvimento

Use `get_logger` ou o adaptador existente `LoggerManager`; não adicione handlers por conta própria. Imports não criam arquivos. `log_context` mantém o contexto por thread/tarefa; o supervisor transmite explicitamente projeto, execução, etapa, job e operação ao worker. O contexto molecular é acrescentado na execução dos scripts Chimera e Vina.

```python
from biomolexplorer.diagnostics import get_logger, log_context, event

logger = get_logger("analysis")
with log_context(job_id="example-job", operation="redocking", pair="4M0E|1YL|604|A"):
    event(logger, "analysis.started", "Starting reference analysis")
```

O JSONL usa `schema_version=1`, horário UTC com milissegundos, nível, componente, evento, processo, thread e código de origem. IDs, duração, ferramenta, configuração, código de saída, causa e ação são opcionais conforme o evento. Falhas na escrita de diagnóstico não devem invalidar resultados científicos.

`tests/test_structured_diagnostics.py` verifica isolamento do contexto, máscara de credenciais, ausência de duplicação, causas encadeadas, códigos de comandos, filtros e rotação compartilhada por processos. Os testes de workers e pipelines verificam o fluxo real entre processos. Veja também [validação do pipeline](pipeline_validation.md).
