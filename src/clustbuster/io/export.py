"""Export protocol implemented by the CSV and H5AD adapters in the MVP slice."""

from pathlib import Path
from typing import Protocol

from clustbuster.models import ExportResult, Workspace


class WorkspaceExporter(Protocol):
    def export(self, workspace: Workspace, destination: Path) -> ExportResult: ...

