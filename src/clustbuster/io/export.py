"""Export protocol and errors implemented by the MVP export service."""

from pathlib import Path
from typing import Protocol

from clustbuster.models import ExportResult, Workspace


class ExportError(ValueError):
    """Raised when an export cannot be safely generated."""


class ExportCollisionError(ExportError):
    """Raised rather than overwriting an existing ClustBuster export namespace."""


class WorkspaceExporter(Protocol):
    def export(self, workspace: Workspace, destination: Path) -> ExportResult: ...

