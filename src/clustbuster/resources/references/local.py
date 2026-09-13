"""Validated, read-only local reference-matrix provider."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from clustbuster.models import (
    LoadedReference,
    ProviderStatus,
    ReferenceFilters,
    ReferenceSummary,
)
from clustbuster.resources.errors import (
    ProviderSchemaError,
    ProviderUnavailableError,
    ResourceNotFoundError,
)

REQUIRED_METADATA = {
    "schema_version",
    "reference_id",
    "name",
    "species",
    "matrix",
    "resource_version",
}


def _checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class LocalReferenceProvider:
    def __init__(self, root: Path) -> None:
        self.root = root.resolve()
        self._metadata: dict[str, dict[str, Any]] | None = None

    def _load_metadata(self) -> dict[str, dict[str, Any]]:
        if self._metadata is not None:
            return self._metadata
        if not self.root.is_dir():
            raise ProviderUnavailableError(f"Reference directory is unavailable: {self.root}")
        records: dict[str, dict[str, Any]] = {}
        for metadata_path in sorted(self.root.glob("*.json")):
            try:
                metadata = json.loads(metadata_path.read_text())
            except (OSError, json.JSONDecodeError) as exc:
                raise ProviderSchemaError(
                    f"Reference metadata is unreadable: {metadata_path.name}"
                ) from exc
            if not isinstance(metadata, dict):
                raise ProviderSchemaError(
                    f"Reference metadata must be an object: {metadata_path.name}"
                )
            missing = sorted(REQUIRED_METADATA - set(metadata))
            if missing:
                raise ProviderSchemaError(
                    f"Reference metadata {metadata_path.name} is missing: "
                    + ", ".join(missing)
                )
            reference_id = str(metadata["reference_id"])
            if reference_id in records:
                raise ProviderSchemaError(f"Duplicate reference identifier: {reference_id}")
            matrix_path = (self.root / str(metadata["matrix"])).resolve()
            if not matrix_path.is_relative_to(self.root):
                raise ProviderSchemaError(
                    f"Reference matrix escapes the configured directory: {reference_id}"
                )
            # Header-only discovery supports older sidecars without reading expression values.
            if (not metadata.get("cell_types") and matrix_path.is_file()
                    and matrix_path.suffix.casefold() in {".csv", ".tsv", ".txt"}):
                try:
                    header = pd.read_csv(
                        matrix_path, sep="," if matrix_path.suffix.casefold() == ".csv" else "\t",
                        nrows=0,
                    )
                    metadata["cell_types"] = [str(c) for c in header.columns if c != "gene"]
                except (OSError, pd.errors.ParserError, ValueError):
                    pass  # Validation reports the error when the reference is selected.
            metadata["_matrix_path"] = matrix_path
            records[reference_id] = metadata
        if not records:
            raise ProviderUnavailableError("No local reference metadata files were found")
        self._metadata = records
        return records

    @staticmethod
    def _summary(metadata: dict[str, Any]) -> ReferenceSummary:
        return ReferenceSummary(
            reference_id=str(metadata["reference_id"]),
            name=str(metadata["name"]),
            species=str(metadata["species"]),
            tissue=str(metadata.get("tissue") or "") or None,
            disease=str(metadata.get("disease") or "") or None,
            assay=str(metadata.get("assay") or "") or None,
            platform=str(metadata.get("platform") or "") or None,
            resource_version=str(metadata["resource_version"]),
            cell_types=tuple(str(name) for name in metadata.get("cell_types", [])
                             if str(name) != "gene"),
            metadata={key: value for key, value in metadata.items() if not key.startswith("_")},
        )

    def status(self) -> ProviderStatus:
        try:
            metadata = self._load_metadata()
        except (ProviderUnavailableError, ProviderSchemaError) as exc:
            return ProviderStatus(False, str(exc))
        versions = sorted({str(item["resource_version"]) for item in metadata.values()})
        return ProviderStatus(
            True,
            f"{len(metadata):,} local reference matrix available",
            ", ".join(versions),
        )

    def list_references(self, filters: ReferenceFilters) -> list[ReferenceSummary]:
        summaries = [self._summary(item) for item in self._load_metadata().values()]
        if filters.species:
            summaries = [
                item for item in summaries if item.species.casefold() == filters.species.casefold()
            ]
        if filters.tissue:
            summaries = [
                item
                for item in summaries
                if item.tissue and item.tissue.casefold() == filters.tissue.casefold()
            ]
        if filters.disease:
            summaries = [
                item
                for item in summaries
                if item.disease and item.disease.casefold() == filters.disease.casefold()
            ]
        return sorted(summaries, key=lambda item: item.name.casefold())

    def load_reference(self, reference_id: str) -> LoadedReference:
        metadata = self._load_metadata().get(reference_id)
        if metadata is None:
            raise ResourceNotFoundError(f"Unknown reference matrix: {reference_id}")
        matrix_path = Path(metadata["_matrix_path"])
        if not matrix_path.is_file():
            raise ProviderUnavailableError(
                f"Reference matrix file is unavailable: {matrix_path.name}"
            )
        suffix = matrix_path.suffix.casefold()
        try:
            if suffix == ".csv":
                table = pd.read_csv(matrix_path)
            elif suffix in {".tsv", ".txt"}:
                table = pd.read_csv(matrix_path, sep="\t")
            elif suffix == ".parquet":
                table = pd.read_parquet(matrix_path)
            else:
                raise ProviderSchemaError(
                    f"Unsupported reference matrix format: {matrix_path.suffix}"
                )
        except (OSError, ImportError, pd.errors.ParserError) as exc:
            raise ProviderUnavailableError("Reference matrix could not be read") from exc
        if "gene" not in table or table["gene"].isna().any() or not table["gene"].is_unique:
            raise ProviderSchemaError(
                "Reference matrix requires a unique, non-missing gene column"
            )
        reference_columns = [column for column in table.columns if column != "gene"]
        if not reference_columns:
            raise ProviderSchemaError("Reference matrix requires at least one cell-type column")
        numeric = table[reference_columns].apply(pd.to_numeric, errors="coerce")
        if numeric.isna().any().any() or not np.isfinite(numeric.to_numpy()).all():
            raise ProviderSchemaError("Reference matrix values must be finite numbers")
        checksum = _checksum(matrix_path)
        expected_checksum = str(metadata.get("checksum") or "")
        if expected_checksum and checksum != expected_checksum:
            raise ProviderSchemaError("Reference matrix checksum does not match its metadata")
        public_metadata = {
            key: value for key, value in metadata.items() if not key.startswith("_")
        }
        public_metadata["gene_count"] = len(table)
        public_metadata["cell_types"] = reference_columns
        return LoadedReference(
            summary=self._summary(metadata),
            matrix_path=matrix_path,
            checksum=checksum,
            metadata=public_metadata,
        )
