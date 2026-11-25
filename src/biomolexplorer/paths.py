"""Path compatibility for scripts and explicit application workspaces."""
import os
from pathlib import Path

SOURCE_ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value, root=None):
    if value is None:
        raise ValueError('A path is required')
    path = Path(value).expanduser()
    text = path.as_posix()
    # Existing examples use /datasets as a project-relative path.
    if text == '/datasets' or text.startswith('/datasets/'):
        path = Path(text.lstrip('/'))
    if text.startswith('/src/scripts/') or text.startswith('src/scripts/'):
        resource_root = Path(os.environ.get('BIOMOL_RESOURCE_DIR', Path(__file__).parent / 'resources'))
        return (resource_root / text.removeprefix('/').removeprefix('src/scripts/')).resolve()
    if text.startswith('/src/') or text.startswith('src/'):
        return (SOURCE_ROOT / text.removeprefix('/').removeprefix('src/')).resolve()
    return (Path(root or os.environ.get('BIOMOL_WORKSPACE', Path.cwd())) / path).resolve()


def within(path, root):
    """Whether resolved ``path`` is ``root`` or inside it.

    Same result as ``Path.is_relative_to`` for resolved absolute paths, which
    is purely lexical, without building a Path for every parent: result
    verification calls this for tens of thousands of artifacts.
    """
    path, root = os.fspath(path), os.fspath(root)
    return path == root or path.startswith(root.rstrip(os.sep) + os.sep)


def directory(value):
    """Resolved directory string for legacy filename concatenation."""
    return str(resolve_path(value)) + os.sep


def worker_count():
    return max(1, int(os.environ.get('BIOMOL_CPU_WORKERS', min(4, os.cpu_count() or 1))))
