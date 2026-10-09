"""Bounded scientific commands with correlated lifecycle and failure diagnostics."""
import logging
import os
import re
import shlex
import subprocess
import time
from contextlib import ExitStack
from uuid import uuid4
from .paths import resolve_path
from .diagnostics import diagnose_exception, event, get_logger, log_context, redact


class ScientificToolError(RuntimeError):
    def __init__(self, message, error_code, action):
        super().__init__(message)
        self.error_code = error_code
        self.action = action

    def __reduce__(self):
        # ProcessPoolExecutor must preserve the original scientific failure.
        return type(self), (str(self), self.error_code, self.action)


def run_command(command, cwd=None, check=True, timeout=None):
    arguments = shlex.split(command) if isinstance(command, str) else list(command)
    if not arguments:
        raise ValueError('Command cannot be empty')
    logger = get_logger('scientific_tools')
    tool = os.path.basename(str(arguments[0]))
    working_dir = str(resolve_path(cwd or '.'))
    configuration = next((str(a) for a in arguments[1:] if str(a).endswith(('.com', '.vina'))), None)
    with log_context(command_id=uuid4().hex[:12]), ExitStack() as stack:
        stdin = None
        if '<' in arguments:
            index = arguments.index('<')
            if index != len(arguments)-2:
                raise ValueError('Unsupported command redirection')
            input_path = resolve_path(cwd or '.')/arguments[-1]
            stdin = stack.enter_context(input_path.open('rb'))
            arguments = arguments[:index]
        started = time.monotonic()
        # Do not serialize the complete command: arguments can contain credentials/SMILES.
        event(logger, 'tool.started', 'Scientific command started', tool=tool, cwd=working_dir, configuration=configuration)
        try:
            result = subprocess.run(arguments, cwd=resolve_path(cwd) if cwd else None,
                stdin=stdin, capture_output=True, text=True, check=False,
                timeout=timeout if timeout is not None else float(os.environ.get('BIOMOL_COMMAND_TIMEOUT', '3600')))
        except (FileNotFoundError, subprocess.TimeoutExpired, OSError) as exc:
            elapsed = round((time.monotonic()-started)*1000)
            if isinstance(exc, FileNotFoundError):
                error = ScientificToolError(f'Ferramenta científica não encontrada: {arguments[0]}. Configure a instalação antes de executar.',
                    'TOOL_NOT_FOUND', 'Install the executable and verify PATH in the scientific worker environment.')
            elif isinstance(exc, subprocess.TimeoutExpired):
                error = exc
                error.error_code = 'TOOL_TIMEOUT'
                error.action = 'Review the last tool output, input size and command timeout before retrying.'
            else:
                error = ScientificToolError(f'Não foi possível iniciar a ferramenta {arguments[0]}: {exc}',
                    'TOOL_START_FAILED', 'Verify executable permissions and the working directory.')
            for channel in ('stdout', 'stderr'):
                output = getattr(exc, channel, None)
                if output:
                    if isinstance(output, bytes): output = output.decode('utf-8', errors='replace')
                    event(logger, 'tool.output', f'{tool} {channel}:\n{redact(output[-8000:])}', tool=tool, channel=channel)
            logger.error(str(error), extra={'event': 'tool.failed', 'tool': tool, 'cwd': working_dir, 'configuration': configuration,
                         'duration_ms': elapsed, **diagnose_exception(error)})
            if error is exc:
                raise
            raise error from exc
        elapsed = round((time.monotonic()-started)*1000)
        for channel, output in (('stdout', result.stdout or ''), ('stderr', result.stderr or '')):
            if output.strip():
                event(logger, 'tool.output', f'{tool} {channel}:\n{redact(output.rstrip())}', tool=tool, channel=channel)
        details = '\n'.join((result.stdout or '', result.stderr or '')).strip()
        hidden_failure = tool == 'chimera' and re.search(
            r'Error while sourcing|Traceback \(most recent call last\)|^Error(?:\s|:)', details, re.M)
        if result.returncode != 0 or hidden_failure:
            error = ScientificToolError(f'Falha na ferramenta {arguments[0]} (código {result.returncode}): {details[-2000:]}',
                'TOOL_REPORTED_ERROR' if hidden_failure else 'TOOL_EXIT_FAILED',
                'Inspect tool.output for the original cause, then verify input files and the generated configuration.')
            event(logger, 'tool.failed', str(error), level=logging.ERROR, tool=tool, cwd=working_dir,
                  configuration=configuration, returncode=result.returncode, duration_ms=elapsed, **diagnose_exception(error))
            if check or hidden_failure:
                raise error
            return False
        event(logger, 'tool.succeeded', 'Scientific command completed', tool=tool, cwd=working_dir,
              configuration=configuration, returncode=result.returncode, duration_ms=elapsed)
        return True
