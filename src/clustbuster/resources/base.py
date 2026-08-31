"""Stable protocols for marker and reference providers."""

from typing import Protocol

from clustbuster.models import (
    CellTypeSummary,
    LoadedReference,
    MarkerSet,
    ProviderStatus,
    ReferenceFilters,
    ReferenceSummary,
)


class MarkerProvider(Protocol):
    def status(self) -> ProviderStatus: ...

    def search_cell_types(
        self,
        query: str,
        *,
        species: str | None = None,
        tissue: str | None = None,
        limit: int = 50,
    ) -> list[CellTypeSummary]: ...

    def get_markers(
        self,
        cell_type: str,
        *,
        species: str | None = None,
        tissue: str | None = None,
    ) -> MarkerSet: ...


class ReferenceProvider(Protocol):
    def status(self) -> ProviderStatus: ...

    def list_references(self, filters: ReferenceFilters) -> list[ReferenceSummary]: ...

    def load_reference(self, reference_id: str) -> LoadedReference: ...

