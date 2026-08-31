"""Session-safe upload handling and format-aware import orchestration."""

from __future__ import annotations

import shutil
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from uuid import uuid4

from clustbuster.io.base import ObjectImporter
from clustbuster.io.h5ad import H5adImporter
from clustbuster.io.seurat import SeuratImporter
from clustbuster.models import ImportOptions, ImportResult, ObjectFormat


class UploadError(ValueError):
    """Raised when an uploaded file cannot safely enter the session workspace."""


@dataclass(frozen=True, slots=True)
class SessionFiles:
    workspace_root: Path
    session_id: str
    root: Path
    uploads: Path
    cache: Path
    state: Path
    exports: Path

    @classmethod
    def create(cls, workspace_root: Path) -> SessionFiles:
        resolved_workspace_root = workspace_root.resolve()
        session_id = uuid4().hex
        root = resolved_workspace_root / session_id
        paths = [root / name for name in ("uploads", "cache", "state", "exports")]
        for path in paths:
            path.mkdir(parents=True, mode=0o700, exist_ok=False)
        return cls(
            workspace_root=resolved_workspace_root,
            session_id=session_id,
            root=root,
            uploads=paths[0],
            cache=paths[1],
            state=paths[2],
            exports=paths[3],
        )

    def cleanup(self) -> None:
        expected_root = self.workspace_root / self.session_id
        if self.root != expected_root or self.root.parent != self.workspace_root:
            raise RuntimeError("Refusing to clean a workspace outside its session boundary")
        if self.root.exists():
            shutil.rmtree(self.root)


class ImportService:
    def __init__(self, max_upload_mb: int, *, enable_seurat_import: bool = False) -> None:
        self.max_upload_bytes = max_upload_mb * 1024 * 1024
        self.importers: list[ObjectImporter] = [H5adImporter()]
        if enable_seurat_import:
            self.importers.append(SeuratImporter())

    def import_upload(
        self, upload: Mapping[str, Any], session_files: SessionFiles
    ) -> ImportResult:
        original_name = str(upload.get("name", "uploaded"))
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

        suffix = Path(original_name).suffix.lower()
        destination = session_files.uploads / f"{uuid4().hex}{suffix}"
        shutil.copyfile(source_path, destination)
        probes = sorted(
            ((importer.probe(destination), importer) for importer in self.importers),
            key=lambda item: item[0].confidence,
            reverse=True,
        )
        probe, importer = probes[0]
        if probe.format is ObjectFormat.UNKNOWN or probe.confidence <= 0:
            supported = "H5AD, Seurat RDS, and H5Seurat" if len(self.importers) > 1 else "H5AD"
            raise UploadError(f"The upload is not a recognized {supported} file")
        result = importer.load(destination, ImportOptions())
        result.report.source_filename = original_name
        return result
