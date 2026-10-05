"""Validated configuration independent of scientific dependencies and UI."""
from dataclasses import dataclass
from pathlib import Path
import math
import os


@dataclass(frozen=True)
class AppConfig:
    workspace: Path
    max_jobs: int = 1
    cpu_workers: int = 2
    max_queued_jobs: int = 100
    job_timeout: float = 86400
    worker_python: Path | None = None
    resource_dir: Path | None = None

    def __post_init__(self):
        object.__setattr__(self, 'workspace', Path(self.workspace).expanduser().resolve())
        if self.resource_dir is not None:
            resource_dir = Path(self.resource_dir).resolve()
            if not resource_dir.is_dir():
                raise ValueError('resource_dir must point to a template directory')
            object.__setattr__(self, 'resource_dir', resource_dir)
        if self.worker_python is not None:
            worker = Path(self.worker_python).expanduser().resolve()
            if not worker.is_file() or not os.access(worker, os.X_OK):
                raise ValueError('worker_python must point to an installed Python executable')
            object.__setattr__(self, 'worker_python', worker)
        for name in ('max_jobs', 'cpu_workers', 'max_queued_jobs'):
            value = getattr(self, name)
            if type(value) is not int or value < 1:
                raise ValueError(f'{name} must be a positive integer')
        if type(self.job_timeout) not in (int, float) or not math.isfinite(self.job_timeout) or self.job_timeout <= 0:
            raise ValueError('job_timeout must be positive and finite')

    @property
    def state_dir(self):
        return self.workspace / '.biomolexplorer'
