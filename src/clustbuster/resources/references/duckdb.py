"""Read-only reference-matrix provider for the sovereign-atlas catalog."""

from __future__ import annotations

import hashlib
from contextlib import closing
from pathlib import Path
from typing import Any

import duckdb
import numpy as np
import pandas as pd

from clustbuster.models import LoadedReference, ProviderStatus, ReferenceFilters, ReferenceSummary
from clustbuster.resources.errors import (
    ProviderSchemaError,
    ProviderUnavailableError,
    ResourceNotFoundError,
)

TABLE_NAME = "reference_matrices"
COLUMN_ALIASES = {
    "id": ("matrix_id", "reference_id"),
    "name": ("title", "name"),
    "species": ("species_scientific_name", "species"),
    "tissue": ("tissue",),
    "disease": ("condition", "disease"),
    "assay": ("assay",),
    "platform": ("platform",),
    "filename": ("stored_filename", "filename"),
    "sha256": ("sha256", "checksum"),
    "row_count": ("row_count",),
    "column_count": ("column_count",),
    "column_names": ("column_names",),
    "normalization": ("normalization",),
    "value_type": ("value_type",),
    "source_title": ("source_title",),
    "citation": ("citation",),
    "doi": ("doi",),
    "pmid": ("pmid",),
    "source_url": ("source_url",),
    "data_license": ("data_license",),
}
REQUIRED_FIELDS = {"id", "name", "species", "filename", "sha256"}


def _quote(identifier: str) -> str:
    return '"' + identifier.replace('"', '""') + '"'


def _checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class DuckDbReferenceProvider:
    def __init__(
        self,
        catalog_path: Path,
        files_root: Path,
        *,
        resource_version: str | None = None,
    ) -> None:
        self.catalog_path = catalog_path
        self.files_root = files_root.resolve()
        self._configured_version = resource_version
        self._columns: dict[str, str] | None = None
        self._validated: set[tuple[str, int, int]] = set()

    def _connect(self) -> duckdb.DuckDBPyConnection:
        if not self.catalog_path.is_file():
            raise ProviderUnavailableError(f"Reference catalog is unavailable: {self.catalog_path}")
        try:
            return duckdb.connect(str(self.catalog_path), read_only=True)
        except duckdb.Error as exc:
            raise ProviderUnavailableError("Reference catalog could not be opened") from exc

    def _schema(self) -> dict[str, str]:
        if self._columns is not None:
            return self._columns
        try:
            with closing(self._connect()) as connection:
                rows = connection.execute(f"DESCRIBE {_quote(TABLE_NAME)}").fetchall()
        except duckdb.Error as exc:
            raise ProviderSchemaError(
                f"Reference catalog does not expose the {TABLE_NAME} table"
            ) from exc
        available = {str(row[0]) for row in rows}
        resolved = {
            field: next((name for name in aliases if name in available), "")
            for field, aliases in COLUMN_ALIASES.items()
        }
        missing = sorted(field for field in REQUIRED_FIELDS if not resolved[field])
        if missing:
            raise ProviderSchemaError(
                f"Reference catalog is incompatible; missing fields: {', '.join(missing)}"
            )
        self._columns = {field: name for field, name in resolved.items() if name}
        return self._columns

    def status(self) -> ProviderStatus:
        try:
            self._schema()
            with closing(self._connect()) as connection:
                row = connection.execute(f"SELECT count(*) FROM {_quote(TABLE_NAME)}").fetchone()
                assert row is not None
                count = int(row[0])
        except (ProviderUnavailableError, ProviderSchemaError) as exc:
            return ProviderStatus(False, str(exc))
        except duckdb.Error as exc:
            return ProviderStatus(False, f"Reference catalog query failed: {exc}")
        return ProviderStatus(
            True,
            f"{count:,} reference matrix{' is' if count == 1 else 'es are'} available",
            self._configured_version,
        )

    def _select_expression(self, field: str) -> str:
        name = self._schema().get(field)
        return _quote(name) if name else "NULL"

    def _rows(self, filters: ReferenceFilters | None = None) -> list[tuple[Any, ...]]:
        columns = self._schema()
        fields = (
            "id",
            "name",
            "species",
            "tissue",
            "disease",
            "assay",
            "platform",
            "filename",
            "sha256",
            "row_count",
            "column_count",
            "column_names",
            "normalization",
            "value_type",
            "source_title",
            "citation",
            "doi",
            "pmid",
            "source_url",
            "data_license",
        )
        clauses: list[str] = []
        parameters: list[str] = []
        if filters:
            for field, value in (
                ("species", filters.species),
                ("tissue", filters.tissue),
                ("disease", filters.disease),
            ):
                if value:
                    name = columns.get(field)
                    if not name:
                        return []
                    clauses.append(f"lower(cast({_quote(name)} AS VARCHAR)) = lower(?)")
                    parameters.append(value)
        where = " WHERE " + " AND ".join(clauses) if clauses else ""
        select_fields = ", ".join(self._select_expression(field) for field in fields)
        name_column = _quote(columns["name"])
        sql = (
            f"SELECT {select_fields} FROM {_quote(TABLE_NAME)}{where} "
            f"ORDER BY lower(cast({name_column} AS VARCHAR))"
        )
        try:
            with closing(self._connect()) as connection:
                return connection.execute(sql, parameters).fetchall()
        except duckdb.Error as exc:
            raise ProviderUnavailableError("Reference catalog query failed") from exc

    def _summary(self, row: tuple[Any, ...]) -> ReferenceSummary:
        return ReferenceSummary(
            reference_id=str(row[0]),
            name=str(row[1]),
            species=str(row[2]),
            tissue=str(row[3]) if row[3] else None,
            disease=str(row[4]) if row[4] else None,
            assay=str(row[5]) if row[5] else None,
            platform=str(row[6]) if row[6] else None,
            resource_version=self._configured_version or f"sha256:{str(row[8])[:12]}",
        )

    def list_references(self, filters: ReferenceFilters) -> list[ReferenceSummary]:
        return [self._summary(row) for row in self._rows(filters)]

    def load_reference(self, reference_id: str) -> LoadedReference:
        row = next((item for item in self._rows() if str(item[0]) == reference_id), None)
        if row is None:
            raise ResourceNotFoundError(f"Unknown reference matrix: {reference_id}")
        filename = Path(str(row[7])).name
        matrix_path = (self.files_root / filename).resolve()
        if not matrix_path.is_relative_to(self.files_root):
            raise ProviderSchemaError("Reference matrix path escapes the configured directory")
        if not matrix_path.is_file():
            raise ProviderUnavailableError(f"Reference matrix file is unavailable: {filename}")

        expected_checksum = str(row[8])
        stat = matrix_path.stat()
        validation_key = (expected_checksum, stat.st_size, stat.st_mtime_ns)
        if validation_key not in self._validated:
            try:
                actual_checksum = _checksum(matrix_path)
            except OSError as exc:
                raise ProviderUnavailableError("Reference matrix could not be read") from exc
            if actual_checksum != expected_checksum:
                raise ProviderSchemaError("Reference matrix checksum does not match its catalog")
            self._validated.add(validation_key)

        try:
            if matrix_path.suffix.casefold() == ".csv":
                table = pd.read_csv(matrix_path)
            elif matrix_path.suffix.casefold() in {".tsv", ".txt"}:
                table = pd.read_csv(matrix_path, sep="\t")
            elif matrix_path.suffix.casefold() == ".parquet":
                table = pd.read_parquet(matrix_path)
            else:
                raise ProviderSchemaError(
                    f"Unsupported reference matrix format: {matrix_path.suffix}"
                )
        except (OSError, ImportError, pd.errors.ParserError) as exc:
            raise ProviderUnavailableError("Reference matrix could not be read") from exc
        if "gene" not in table or table["gene"].isna().any() or not table["gene"].is_unique:
            raise ProviderSchemaError("Reference matrix requires a unique, non-missing gene column")
        reference_columns = [column for column in table.columns if column != "gene"]
        if not reference_columns:
            raise ProviderSchemaError("Reference matrix requires at least one cell-type column")
        numeric = table[reference_columns].apply(pd.to_numeric, errors="coerce")
        if numeric.isna().any().any() or not np.isfinite(numeric.to_numpy()).all():
            raise ProviderSchemaError("Reference matrix values must be finite numbers")

        metadata_names = (
            "matrix_id",
            "title",
            "species",
            "tissue",
            "condition",
            "assay",
            "platform",
            "stored_filename",
            "sha256",
            "row_count",
            "column_count",
            "column_names",
            "normalization",
            "value_type",
            "source_title",
            "citation",
            "doi",
            "pmid",
            "source_url",
            "data_license",
        )
        metadata = dict(zip(metadata_names, row, strict=True))
        metadata["gene_count"] = len(table)
        metadata["cell_types"] = reference_columns
        return LoadedReference(
            summary=self._summary(row),
            matrix_path=matrix_path,
            checksum=expected_checksum,
            metadata=metadata,
        )
