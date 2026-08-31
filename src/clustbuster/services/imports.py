"""Session-safe upload handling and H5AD import orchestration."""

from __future__ import annotations

import shutil
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from uuid import uuid4

from clustbuster.io.h5ad import H5adImporter
from clustbuster.models import ImportResult


class UploadError(ValueError):
    """Raised when an uploaded file cannot safely enter the session workspace."""


@dataclass(frozen=True, slots=True)
class SessionFiles:
    root: Path
    uploads: Path
    cache: Path
    state: Path
    exports: Path

    @classmethod
    def create(cls, workspace_root: Path) -> SessionFiles:
        root = workspace_root.resolve() / uuid4().hex
        paths = [root / name for name in ("uploads", "cache", "state", "exports")]
        for path in paths:
            path.mkdir(parents=True, mode=0o700, exist_ok=False)
        return cls(root=root, uploads=paths[0], cache=paths[1], state=paths[2], exports=paths[3])

    def cleanup(self) -> None:
        if self.root.exists():
            shutil.rmtree(self.root)


class ImportService:
    def __init__(self, max_upload_mb: int, importer: H5adImporter | None = None) -> None:
        self.max_upload_bytes = max_upload_mb * 1024 * 1024
        self.importer = importer or H5adImporter()

    def import_upload(
        self, upload: Mapping[str, Any], session_files: SessionFiles
    ) -> ImportResult:
        original_name = str(upload.get("name", "uploaded.h5ad"))
        source_path = Path(str(upload.get("datapath", "")))
        if not source_path.is_file():
            raise UploadError("The uploaded temporary file is unavailable; please upload it again")
        size = source_path.stat().st_size
        if size == 0:
            raise UploadError("The uploaded file is empty")
        if size > self.max_upload_bytes:
            limit_mb = self.max_upload_bytes // (1024 * 1024)
            raise UploadError(f"The upload exceeds the configured {limit_mb} MB limit")
        free_bytes = shutil.disk_usage(session_files.root).free
        if free_bytes < size * 2:
            raise UploadError("There is not enough session disk space to safely import this file")

        destination = session_files.uploads / f"{uuid4().hex}.h5ad"
        shutil.copyfile(source_path, destination)
        result = self.importer.load(destination)
        result.report.source_filename = original_name
        return result
