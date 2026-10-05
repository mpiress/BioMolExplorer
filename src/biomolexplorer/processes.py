"""Bounded execution of scientific tools without a shell interpreter."""
import os
import shlex
import subprocess
from contextlib import ExitStack
from .paths import resolve_path


def run_command(command, cwd=None, check=True, timeout=None):
    arguments = shlex.split(command) if isinstance(command, str) else list(command)
    if not arguments:
        raise ValueError('Command cannot be empty')
    with ExitStack() as stack:
        stdin = None
        # Existing DOCK6 templates use showbox < file; handle this explicitly.
        if '<' in arguments:
            index = arguments.index('<')
            if index != len(arguments) - 2:
                raise ValueError('Unsupported command redirection')
            input_path = (resolve_path(cwd) if cwd else resolve_path('.')) / arguments[-1]
            stdin = stack.enter_context(input_path.open('rb'))
            arguments = arguments[:index]
        result = subprocess.run(arguments, cwd=resolve_path(cwd) if cwd else None,
            stdin=stdin, capture_output=True, text=True, check=check,
            timeout=timeout or float(os.environ.get('BIOMOL_COMMAND_TIMEOUT', '3600')))
        return result.returncode == 0
