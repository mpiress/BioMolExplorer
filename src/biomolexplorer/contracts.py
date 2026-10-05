"""Serializable results and operation metadata for CLI and future Flet UI."""
from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(frozen=True)
class OperationResult:
    operation: str
    artifacts: list[str] = field(default_factory=list)
    details: dict[str, Any] = field(default_factory=dict)

    def to_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class OperationSpec:
    module: str
    function: str
    required: tuple[str, ...]
    optional: tuple[str, ...] = ()
    enums: dict[str, str] = field(default_factory=dict)
