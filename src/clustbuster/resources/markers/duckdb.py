"""Read-only MarkerCodex provider backed by its stable DuckDB consumer view."""

from __future__ import annotations

import hashlib
from contextlib import closing
from pathlib import Path
from typing import Any

import duckdb

from clustbuster.models import (
    CellTypeSummary,
    MarkerFacets,
    MarkerRecord,
    MarkerSet,
    ProviderStatus,
)
from clustbuster.resources.errors import (
    ProviderSchemaError,
    ProviderUnavailableError,
    ResourceNotFoundError,
)

VIEW_NAME = "marker_atlas"
COLUMN_ALIASES = {
    "cell_type": ("cell_type",),
    "gene": ("gene_symbol", "gene"),
    "species": ("species", "species_scientific_name"),
    "tissue": ("tissue",),
    "direction": ("marker_direction", "direction"),
    "confidence": ("confidence",),
    "verified": ("bbsr_verified", "verified"),
    "evidence": ("source_titles", "evidence"),
    "citation": ("source_identifiers", "citation"),
}
REQUIRED_FIELDS = {"cell_type", "gene", "species", "direction"}


def _quote(identifier: str) -> str:
    return '"' + identifier.replace('"', '""') + '"'


def _file_version(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return f"sha256:{digest.hexdigest()[:12]}"


class DuckDbMarkerProvider:
    """Query MarkerCodex without exposing its normalized internal schema."""

    def __init__(self, path: Path, *, resource_version: str | None = None) -> None:
        self.path = path
        self._configured_version = resource_version
        self._columns: dict[str, str] | None = None
        self._resource_version: str | None = None
        self._search_columns: tuple[str, ...] = ()

    def _connect(self) -> duckdb.DuckDBPyConnection:
        if not self.path.is_file():
            raise ProviderUnavailableError(f"Marker database is unavailable: {self.path}")
        try:
            return duckdb.connect(str(self.path), read_only=True)
        except duckdb.Error as exc:
            raise ProviderUnavailableError("Marker database could not be opened") from exc

    def _schema(self) -> dict[str, str]:
        if self._columns is not None:
            return self._columns
        try:
            with closing(self._connect()) as connection:
                rows = connection.execute(f"DESCRIBE {_quote(VIEW_NAME)}").fetchall()
        except duckdb.Error as exc:
            raise ProviderSchemaError(
                f"Marker database does not expose the {VIEW_NAME} consumer view"
            ) from exc
        available = {str(row[0]) for row in rows}
        self._search_columns = tuple(sorted(available))
        resolved = {
            field: next((name for name in aliases if name in available), "")
            for field, aliases in COLUMN_ALIASES.items()
        }
        missing = sorted(field for field in REQUIRED_FIELDS if not resolved[field])
        if missing:
            raise ProviderSchemaError(
                f"Marker view is incompatible; missing consumer fields: {', '.join(missing)}"
            )
        self._columns = {field: name for field, name in resolved.items() if name}
        return self._columns

    def _version(self) -> str:
        if self._configured_version:
            return self._configured_version
        if self._resource_version is None:
            try:
                self._resource_version = _file_version(self.path)
            except OSError as exc:
                raise ProviderUnavailableError("Marker database could not be checksummed") from exc
        return self._resource_version

    def status(self) -> ProviderStatus:
        try:
            self._schema()
            with closing(self._connect()) as connection:
                row = connection.execute(f"SELECT count(*) FROM {_quote(VIEW_NAME)}").fetchone()
                assert row is not None
                count = int(row[0])
            version = self._version()
        except (ProviderUnavailableError, ProviderSchemaError) as exc:
            return ProviderStatus(False, str(exc))
        except duckdb.Error as exc:
            return ProviderStatus(False, f"Marker database query failed: {exc}")
        return ProviderStatus(True, f"{count:,} MarkerCodex assertions available", version)

    def list_facets(self) -> MarkerFacets:
        columns = self._schema()
        species = _quote(columns["species"])
        tissue_name = columns.get("tissue")
        tissue_query = (
            f"SELECT DISTINCT trim(cast({_quote(tissue_name)} AS VARCHAR)) AS value "
            f"FROM {_quote(VIEW_NAME)} WHERE value <> '' ORDER BY value"
            if tissue_name
            else None
        )
        try:
            with closing(self._connect()) as connection:
                species_values = tuple(
                    str(row[0])
                    for row in connection.execute(
                        f"SELECT DISTINCT trim(cast({species} AS VARCHAR)) AS value "
                        f"FROM {_quote(VIEW_NAME)} WHERE value <> '' ORDER BY value"
                    ).fetchall()
                )
                tissue_values = (
                    tuple(str(row[0]) for row in connection.execute(tissue_query).fetchall())
                    if tissue_query
                    else ()
                )
        except duckdb.Error as exc:
            raise ProviderUnavailableError("Marker facets could not be read") from exc
        return MarkerFacets(species_values, tissue_values)

    def _where(
        self,
        *,
        species: str | None,
        tissue: str | None,
    ) -> tuple[list[str], list[Any]]:
        columns = self._schema()
        clauses: list[str] = []
        parameters: list[Any] = []
        if species:
            clauses.append(f"lower(cast({_quote(columns['species'])} AS VARCHAR)) = lower(?)")
            parameters.append(species)
        if tissue:
            tissue_column = columns.get("tissue")
            if not tissue_column:
                return ["FALSE"], []
            clauses.append(f"lower(cast({_quote(tissue_column)} AS VARCHAR)) = lower(?)")
            parameters.append(tissue)
        return clauses, parameters

    def search_cell_types(
        self,
        query: str,
        *,
        species: str | None = None,
        tissue: str | None = None,
        limit: int = 50,
    ) -> list[CellTypeSummary]:
        if limit < 1:
            raise ValueError("Marker search limit must be positive")
        columns = self._schema()
        cell_type = _quote(columns["cell_type"])
        gene = _quote(columns["gene"])
        species_column = _quote(columns["species"])
        tissue_column = columns.get("tissue")
        tissue_expression = (
            f"nullif(trim(cast({_quote(tissue_column)} AS VARCHAR)), '')"
            if tissue_column
            else "NULL"
        )
        clauses, parameters = self._where(species=species, tissue=tissue)
        if query.strip():
            clauses.append("(" + " OR ".join(
                f"contains(lower(cast({_quote(name)} AS VARCHAR)), lower(?))"
                for name in self._search_columns
            ) + ")")
            parameters.extend([query.strip()] * len(self._search_columns))
        where = " WHERE " + " AND ".join(clauses) if clauses else ""
        sql = f"""
            SELECT cast({cell_type} AS VARCHAR) AS cell_type_value,
                   cast({species_column} AS VARCHAR) AS species_value,
                   {tissue_expression} AS tissue_value, count(DISTINCT {gene})
            FROM {_quote(VIEW_NAME)}{where}
            GROUP BY 1, 2, 3
            ORDER BY 1, 2, 3
            LIMIT ?
        """
        parameters.append(limit)
        try:
            with closing(self._connect()) as connection:
                rows = connection.execute(sql, parameters).fetchall()
        except duckdb.Error as exc:
            raise ProviderUnavailableError("Marker search failed") from exc
        return [
            CellTypeSummary(str(row[0]), str(row[1]), str(row[2]) if row[2] else None, int(row[3]))
            for row in rows
        ]

    def get_markers(
        self,
        cell_type: str,
        *,
        species: str | None = None,
        tissue: str | None = None,
    ) -> MarkerSet:
        columns = self._schema()
        clauses, parameters = self._where(species=species, tissue=tissue)
        clauses.append(f"lower(cast({_quote(columns['cell_type'])} AS VARCHAR)) = lower(?)")
        parameters.append(cell_type)

        def expression(field: str, fallback: str = "NULL") -> str:
            name = columns.get(field)
            return _quote(name) if name else fallback

        sql = f"""
            SELECT cast({_quote(columns["cell_type"])} AS VARCHAR),
                   cast({_quote(columns["gene"])} AS VARCHAR),
                   cast({_quote(columns["species"])} AS VARCHAR),
                   {expression("tissue")}, cast({_quote(columns["direction"])} AS VARCHAR),
                   {expression("evidence")}, {expression("citation")},
                   {expression("confidence")}, {expression("verified")}
            FROM {_quote(VIEW_NAME)}
            WHERE {" AND ".join(clauses)}
            ORDER BY lower(cast({_quote(columns["gene"])} AS VARCHAR))
        """
        try:
            with closing(self._connect()) as connection:
                rows = connection.execute(sql, parameters).fetchall()
        except duckdb.Error as exc:
            raise ProviderUnavailableError("Marker set could not be read") from exc
        if not rows:
            raise ResourceNotFoundError(f"No marker set found for cell type: {cell_type}")

        def confidence(value: Any) -> tuple[float | None, str | None]:
            if value is None:
                return None, None
            try:
                return float(value), None
            except (TypeError, ValueError):
                return None, str(value)

        version = self._version()
        records: list[MarkerRecord] = []
        for row in rows:
            score, label = confidence(row[7])
            records.append(
                MarkerRecord(
                    cell_type=str(row[0]),
                    gene=str(row[1]),
                    species=str(row[2]) or None,
                    tissue=str(row[3]) or None,
                    direction=str(row[4]),
                    evidence=str(row[5]) if row[5] else None,
                    citation=str(row[6]) if row[6] else None,
                    confidence=score,
                    confidence_label=label,
                    verified=bool(row[8]) if row[8] is not None else None,
                    resource_version=version,
                )
            )
        return MarkerSet(records[0].cell_type, tuple(records))
