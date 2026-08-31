"""Importer protocol shared by format-specific adapters."""

from pathlib import Path
from typing import Protocol

from clustbuster.models import ImportOptions, ImportResult, ProbeResult


class ObjectImporter(Protocol):
    def probe(self, path: Path) -> ProbeResult: ...

    def load(self, path: Path, options: ImportOptions) -> ImportResult: ...

