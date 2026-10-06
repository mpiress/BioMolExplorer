"""Command-line access to the same services used by future Flet interfaces."""
import argparse
import json
import sys
from contextlib import redirect_stdout
from pathlib import Path
from .operations import OPERATIONS, execute_operation
from .diagnostics import configure_logging, event, diagnose_exception


def main():
    logger=configure_logging()
    parser = argparse.ArgumentParser(description='BioMolExplorer application services')
    parser.add_argument('operation', choices=sorted(OPERATIONS))
    parser.add_argument('--parameters', type=Path, required=True, help='JSON file with operation parameters')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    event(logger, 'cli.started', 'CLI operation started', operation=args.operation)
    try:
        with redirect_stdout(sys.stderr):
            result = execute_operation(args.operation, json.loads(args.parameters.read_text()), args.output)
    except Exception as exc:
        logger.exception('CLI operation failed: %s', exc, extra={'event': 'cli.failed', 'operation': args.operation, **diagnose_exception(exc)})
        parser.exit(1, f'{type(exc).__name__}: {exc}\n')
    event(logger, 'cli.succeeded', 'CLI operation completed', operation=args.operation, artifacts=len(result.artifacts))
    print(json.dumps(result.to_dict(), indent=2))


if __name__ == '__main__':
    main()
