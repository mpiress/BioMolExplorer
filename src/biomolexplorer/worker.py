"""Child process entry point; JSON is the boundary to scientific adapters."""
import json
import sys
import traceback
from pathlib import Path
from .operations import execute_operation
from .diagnostics import configure_logging
from .progress import report_progress


def main():
    logger=configure_logging('backend')
    request_path, result_path = map(Path, sys.argv[1:])
    try:
        request = json.loads(request_path.read_text())
        logger.info('Worker iniciado; python=%s operation=%s', sys.executable, request['operation'])
        report_progress('Preparando a execução da etapa no Python…')
        result = execute_operation(**request).to_dict()
        report_progress('Etapa concluída. Organizando os resultados…')
        exit_code = 0
    except Exception as exc:
        logger.exception('Falha no worker; request=%s',request_path)
        traceback.print_exc()
        result = {'error': f'{type(exc).__name__}: {exc}'}
        exit_code = 1
    temporary = result_path.with_suffix('.tmp')
    temporary.write_text(json.dumps(result, allow_nan=False))
    temporary.replace(result_path)
    return exit_code


if __name__ == '__main__':
    sys.exit(main())
