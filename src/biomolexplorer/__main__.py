"""Command-line access to the same services used by future Flet interfaces."""
import argparse
import json
import sys
from contextlib import redirect_stdout
from pathlib import Path
from .operations import OPERATIONS, execute_operation
from .diagnostics import configure_logging


def main():
    logger=configure_logging()
    parser = argparse.ArgumentParser(description='BioMolExplorer application services')
    parser.add_argument('operation', choices=sorted(OPERATIONS))
    parser.add_argument('--parameters', type=Path, required=True, help='JSON file with operation parameters')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        with redirect_stdout(sys.stderr):
            result = execute_operation(args.operation, json.loads(args.parameters.read_text()), args.output)
    except Exception as exc:
        logger.exception('Falha na operação CLI %s',args.operation)
        parser.exit(1, f'{type(exc).__name__}: {exc}\n')
    print(json.dumps(result.to_dict(), indent=2))


if __name__ == '__main__':
    main()
