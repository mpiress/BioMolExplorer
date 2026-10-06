# Logs and diagnostics

[Documentation](../README.md) · English · [Português](../logging.md)

The system combines readable text logs, filterable JSONL events and a per-job summary. The interface, pipeline, supervisor, worker and scientific modules share the same format. Technical messages and executable output retain their original language; interface controls and displayed messages continue following the session language.

## Find the cause of a failure

1. In **Runs**, identify the failed stage and open **View stage log**. Its `execution.log` contains worker output, including scientific events and tracebacks.
2. Use the stage's `job_id` to open `logs/jobs/<job_id>/diagnostic.json`. Review `status`, `error_code`, `error` and `action`.
3. In that directory, find `tool.failed` in `events.jsonl`. The `command_id` connects a command's start, output and result; `configuration` identifies the Chimera script or Vina configuration and `cwd` identifies its working directory.
4. Review the preceding `tool.output` for the tool's original message. For Python failures, `source` and `exception.traceback` identify the location and exception chain. Subsequent `worker.failed`, `job.failed` and `stage.failed` events describe propagation through those layers.

Not every failure has a `tool.failed` event: validation can fail before tools start. In that case, inspect the stage or worker error and the reported path/record. A pipeline failing before job creation records the error in `logs/backend.log` and `logs/events.jsonl`.

## Files and retention

The default is `logs/` at the checkout root; installations without a checkout use `logs/` in the working directory. Set `BIOMOL_LOG_DIR` or use `biomolexplorer-ui --log-dir /path/to/logs` to change the root.

| File | Purpose |
| --- | --- |
| `frontend.log` | Interface actions and failures |
| `backend.log` | Pipeline, supervision and job lifecycle |
| `errors.log` | ERROR/CRITICAL events, with causes and tracebacks when available |
| `events.jsonl` | One JSON object per event, including context and diagnostic fields |
| `jobs/<job_id>/diagnostic.json` | Most recent job summary, written atomically |
| `jobs/<job_id>/` | Events and scientific files isolated to that worker |
| `jobs/<job_id>/execution.log` | Copy of the execution log; the original remains in the private job directory |
| `docking.log`, `complex.log`, `bioactivities.log`, etc. | Existing scientific files using the common format |

Managed text and JSONL files rotate at approximately 5 MiB, retaining three earlier copies (`.1` through `.3`). An individual large event can exceed this size. Writes and rotation are synchronized across threads and processes in the Linux environment. `.lock` files support that synchronization. The summary is replaced on updates; `execution.log` and old job directories have no automatic age-based cleanup.

Diagnostics contain study paths and identifiers. Credentials must not be logged: new execution events do not serialize complete parameters or commands. Formatters also mask common patterns such as `token=...` and `password=...`; this protection cannot identify every kind of sensitive content. Review files before sharing them.

## Filter quickly

Run inside the project's Python environment:

```bash
# Main-process errors and warnings
python -m biomolexplorer.log_report

# Worker errors; replace JOB_ID with the actual identifier
python -m biomolexplorer.log_report --job JOB_ID --level ERROR

# Command history, including output and duration
python -m biomolexplorer.log_report --job JOB_ID --level INFO --limit 100

# Pipeline events for a run/stage
python -m biomolexplorer.log_report --run RUN_ID --stage STAGE_ID

# Structured output for external analysis
python -m biomolexplorer.log_report --job JOB_ID --json
```

`--directory` selects another log root; with `--job`, the reader looks for `jobs/<job_id>` under that root. `--operation redocking` filters by operation. The reader scans rotated JSONL files and retains only the latest matching events without loading the entire history into memory. Older logs created before JSONL support remain available as text.

## Diagnostic codes

| Code | Meaning and verification |
| --- | --- |
| `TOOL_NOT_FOUND` | Executable unavailable; check installation and worker PATH |
| `TOOL_START_FAILED` | Could not start; check permissions and working directory |
| `TOOL_TIMEOUT` | Command exceeded its limit; review output, input and `BIOMOL_COMMAND_TIMEOUT` |
| `TOOL_EXIT_FAILED` | Tool returned a nonzero exit code; inspect stdout/stderr |
| `TOOL_REPORTED_ERROR` | Chimera reported an error despite returning zero; inspect the retained script |
| `PREPARED_OUTPUT_INVALID` | Preparation output missing, empty or without atoms; inspect the tool output and script |
| `EXECUTION_TIMEOUT` | Job exceeded its supervisor limit; inspect the last command |
| `INPUT_NOT_FOUND` | Input file unavailable |
| `ACCESS_DENIED` | Filesystem access denied |
| `VALIDATION_FAILED` | Invalid input or configuration |
| `WORKER_EXIT_FAILED` | Worker exited without a valid diagnostic; inspect `execution.log` |
| `UNEXPECTED_ERROR` | Unclassified exception; inspect its cause and source location |

Codes are stable for filtering. The `action` guidance suggests the next check without automatically judging scientific validity. User cancellation is a separate status, not a tool error.

## Development contract

Use `get_logger` or the existing `LoggerManager` adapter rather than adding handlers yourself. Imports do not create files. `log_context` scopes context to the thread/task; the supervisor explicitly passes project, run, stage, job and operation to the worker. Molecular context is added during Chimera and Vina script execution.

```python
from biomolexplorer.diagnostics import get_logger, log_context, event

logger = get_logger("analysis")
with log_context(job_id="example-job", operation="redocking", pair="4M0E|1YL|604|A"):
    event(logger, "analysis.started", "Starting reference analysis")
```

JSONL uses `schema_version=1`, UTC timestamps with milliseconds, severity, component, event, process, thread and source location. IDs, duration, tool, configuration, exit code, cause and action are optional depending on the event. Diagnostic write failures must not invalidate scientific results.

`tests/test_structured_diagnostics.py` covers context isolation, credential masking, duplicate prevention, chained causes, command codes, filters and shared process rotation. Worker and pipeline tests exercise the actual process boundary. See [pipeline validation](pipeline_validation.md).
