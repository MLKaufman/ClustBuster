"""Typed enrichment boundaries independent of any remote provider."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import pandas as pd


class EnrichmentError(ValueError):
    """Raised when enrichment input or provider output is invalid."""


@dataclass(slots=True)
class EnrichmentResult:
    genes: tuple[str, ...]
    library: str
    source: str
    values: pd.DataFrame


@dataclass(slots=True)
class AllClusterOraResult:
    library: str
    source: str
    markers: pd.DataFrame
    values: pd.DataFrame
    failures: tuple[str, ...]


class EnrichmentProvider(Protocol):
    def enrich(self, genes: tuple[str, ...], *, description: str) -> EnrichmentResult: ...
