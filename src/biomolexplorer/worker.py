"""Child process entry point; JSON is the boundary to scientific adapters."""
import json
import sys
import time
from pathlib import Path
from .operations import execute_operation
from .diagnostics import configure_logging, diagnose_exception, event, write_summary
from .progress import report_progress


def main():
    logger=configure_logging('backend', console=True)
    started = time.monotonic()
    request_path, result_path = map(Path, sys.argv[1:])
    try:
        request = json.loads(request_path.read_text())
        event(logger, 'worker.started', f'Worker iniciado; python={sys.executable}', operation=request['operation'])
        write_summary('running', operation=request['operation'])
        report_progress('Preparando a execução da etapa no Python…')
        result = execute_operation(**request).to_dict()
        report_progress('Etapa concluída. Organizando os resultados…')
        exit_code = 0
        duration = round((time.monotonic()-started)*1000)
        event(logger, 'worker.succeeded', 'Scientific worker completed', duration_ms=duration, artifacts=len(result.get('artifacts', [])))
        write_summary('succeeded', duration_ms=duration, artifacts=len(result.get('artifacts', [])))
    except Exception as exc:
        diagnosis = diagnose_exception(exc)
        duration = round((time.monotonic()-started)*1000)
        logger.exception('Scientific worker failed: %s', exc, extra={'event': 'worker.failed', 'duration_ms': duration, **diagnosis})
        write_summary('failed', error=exc, duration_ms=duration)
        result = {'error': f'{type(exc).__name__}: {exc}', 'diagnostic': diagnosis}
        exit_code = 1
    temporary = result_path.with_suffix('.tmp')
    temporary.write_text(json.dumps(result, allow_nan=False))
    temporary.replace(result_path)
    return exit_code


if __name__ == '__main__':
    sys.exit(main())
