"""Validated, read-only CSV marker provider."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

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

REQUIRED_COLUMNS = {
    "cell_type",
    "gene",
    "species",
    "tissue",
    "direction",
    "evidence",
    "citation",
    "confidence",
    "resource_version",
}


class CsvMarkerProvider:
    def __init__(self, path: Path) -> None:
        self.path = path
        self._table: pd.DataFrame | None = None

    def _load(self) -> pd.DataFrame:
        if self._table is not None:
            return self._table
        if not self.path.is_file():
            raise ProviderUnavailableError(f"Marker catalog is unavailable: {self.path}")
        try:
            table = pd.read_csv(self.path, keep_default_na=False)
        except (OSError, pd.errors.ParserError) as exc:
            raise ProviderUnavailableError("Marker catalog could not be read") from exc
        missing = sorted(REQUIRED_COLUMNS - set(table.columns))
        if missing:
            raise ProviderSchemaError(
                "Marker catalog is missing required columns: " + ", ".join(missing)
            )
        if table.empty:
            raise ProviderSchemaError("Marker catalog contains no records")
        for column in REQUIRED_COLUMNS - {"confidence"}:
            table[column] = table[column].astype(str)
        if (table["cell_type"].str.strip() == "").any() or (
            table["gene"].str.strip() == ""
        ).any():
            raise ProviderSchemaError("Marker cell types and genes must not be empty")
        invalid_directions = sorted(set(table["direction"]) - {"positive", "negative"})
        if invalid_directions:
            raise ProviderSchemaError(
                "Marker direction must be positive or negative: "
                + ", ".join(invalid_directions)
            )
        confidence = pd.to_numeric(table["confidence"], errors="coerce")
        if confidence.isna().any() or ((confidence < 0) | (confidence > 1)).any():
            raise ProviderSchemaError("Marker confidence values must be between 0 and 1")
        table = table.copy()
        table["confidence"] = confidence
        self._table = table
        return table

    def status(self) -> ProviderStatus:
        try:
            table = self._load()
        except (ProviderUnavailableError, ProviderSchemaError) as exc:
            return ProviderStatus(False, str(exc))
        versions = sorted(set(table["resource_version"].astype(str)))
        return ProviderStatus(
            True,
            f"{len(table):,} marker records available",
            ", ".join(versions),
        )

    def list_facets(self) -> MarkerFacets:
        table = self._load()
        return MarkerFacets(
            species=tuple(sorted(filter(None, set(table["species"])), key=str.casefold)),
            tissues=tuple(sorted(filter(None, set(table["tissue"])), key=str.casefold)),
        )

    @staticmethod
    def _filter(
        table: pd.DataFrame,
        *,
        species: str | None,
        tissue: str | None,
    ) -> pd.DataFrame:
        selected = table
        if species:
            selected = selected.loc[
                selected["species"].str.casefold() == species.casefold()
            ]
        if tissue:
            selected = selected.loc[
                selected["tissue"].str.casefold() == tissue.casefold()
            ]
        return selected

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
        table = self._filter(self._load(), species=species, tissue=tissue)
        if query.strip():
            table = table.loc[
                table.astype(str).apply(
                    lambda column: column.str.contains(query.strip(), case=False, regex=False)
                ).any(axis=1)
            ]
        summaries: list[CellTypeSummary] = []
        grouped = table.groupby(["cell_type", "species", "tissue"], sort=True)
        for (cell_type, record_species, record_tissue), group in grouped:
            summaries.append(
                CellTypeSummary(
                    str(cell_type),
                    str(record_species) or None,
                    str(record_tissue) or None,
                    int(group["gene"].nunique()),
                )
            )
        return summaries[:limit]

    def get_markers(
        self,
        cell_type: str,
        *,
        species: str | None = None,
        tissue: str | None = None,
    ) -> MarkerSet:
        table = self._filter(self._load(), species=species, tissue=tissue)
        table = table.loc[table["cell_type"].str.casefold() == cell_type.casefold()]
        if table.empty:
            raise ResourceNotFoundError(f"No marker set found for cell type: {cell_type}")
        records = tuple(
            MarkerRecord(
                cell_type=str(row.cell_type),
                gene=str(row.gene),
                species=str(row.species) or None,
                tissue=str(row.tissue) or None,
                direction=str(row.direction),
                evidence=str(row.evidence) or None,
                citation=str(row.citation) or None,
                confidence=float(row.confidence),
                resource_version=str(row.resource_version) or None,
            )
            for row in table.itertuples(index=False)
        )
        return MarkerSet(records[0].cell_type, records)
