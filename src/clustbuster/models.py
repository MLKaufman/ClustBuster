"""Shared typed models at ClustBuster's application boundaries."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from anndata import AnnData

    from clustbuster.core.annotations import AnnotationStore


class ObjectFormat(StrEnum):
    H5AD = "h5ad"
    SEURAT = "seurat"
    SINGLE_CELL_EXPERIMENT = "single_cell_experiment"
    UNKNOWN = "unknown"


class ExpressionKind(StrEnum):
    X = "X"
    RAW = "raw"
    LAYER = "layer"


@dataclass(frozen=True, slots=True)
class ExpressionSource:
    kind: ExpressionKind
    layer: str | None = None

    def __post_init__(self) -> None:
        if self.kind is ExpressionKind.LAYER and not self.layer:
            raise ValueError("A layer expression source requires a layer name")
        if self.kind is not ExpressionKind.LAYER and self.layer is not None:
            raise ValueError("Only a layer expression source may specify a layer name")

    @classmethod
    def x(cls) -> ExpressionSource:
        return cls(ExpressionKind.X)

    @classmethod
    def raw(cls) -> ExpressionSource:
        return cls(ExpressionKind.RAW)

    @classmethod
    def named_layer(cls, name: str) -> ExpressionSource:
        return cls(ExpressionKind.LAYER, name)

    @classmethod
    def from_label(cls, label: str) -> ExpressionSource:
        if label == ExpressionKind.X.value:
            return cls.x()
        if label == ExpressionKind.RAW.value:
            return cls.raw()
        prefix = "layer:"
        if label.startswith(prefix) and label[len(prefix) :]:
            return cls.named_layer(label[len(prefix) :])
        raise ValueError(f"Unknown expression source: {label}")

    @property
    def label(self) -> str:
        return f"layer:{self.layer}" if self.layer else self.kind.value


@dataclass(frozen=True, slots=True)
class ProbeResult:
    format: ObjectFormat
    confidence: float
    reason: str


@dataclass(frozen=True, slots=True)
class ImportOptions:
    expression_source: ExpressionSource | None = None


@dataclass(slots=True)
class ImportReport:
    source_format: ObjectFormat
    source_filename: str
    cell_count: int
    feature_count: int
    sparse: bool
    layers: tuple[str, ...] = ()
    embeddings: tuple[str, ...] = ()
    candidate_cluster_columns: tuple[str, ...] = ()
    expression_sources: tuple[ExpressionSource, ...] = ()
    warnings: list[str] = field(default_factory=list)
    mapping_decisions: dict[str, str] = field(default_factory=dict)
    unsupported_components: list[str] = field(default_factory=list)


@dataclass(slots=True)
class ImportResult:
    adata: AnnData
    report: ImportReport


@dataclass(frozen=True, slots=True)
class ExportArtifact:
    kind: str
    path: Path
    checksum: str
    media_type: str


@dataclass(frozen=True, slots=True)
class ExportResult:
    artifacts: tuple[ExportArtifact, ...]
    warnings: tuple[str, ...] = ()


@dataclass(slots=True)
class Workspace:
    adata: AnnData
    source_format: ObjectFormat
    source_filename: str
    expression_source: ExpressionSource
    annotations: AnnotationStore
    import_report: ImportReport
    cluster_column: str | None = None
    embedding_key: str | None = None
    workspace_id: str = ""


@dataclass(frozen=True, slots=True)
class ProviderStatus:
    available: bool
    message: str
    version: str | None = None


@dataclass(frozen=True, slots=True)
class CellTypeSummary:
    cell_type: str
    species: str | None = None
    tissue: str | None = None
    marker_count: int | None = None


@dataclass(frozen=True, slots=True)
class MarkerRecord:
    cell_type: str
    gene: str
    species: str | None = None
    tissue: str | None = None
    direction: str = "positive"
    evidence: str | None = None
    citation: str | None = None
    confidence: float | None = None
    resource_version: str | None = None


@dataclass(frozen=True, slots=True)
class MarkerSet:
    cell_type: str
    records: tuple[MarkerRecord, ...]


@dataclass(frozen=True, slots=True)
class ReferenceFilters:
    species: str | None = None
    tissue: str | None = None
    disease: str | None = None


@dataclass(frozen=True, slots=True)
class ReferenceSummary:
    reference_id: str
    name: str
    species: str
    tissue: str | None = None
    resource_version: str | None = None


@dataclass(frozen=True, slots=True)
class LoadedReference:
    summary: ReferenceSummary
    matrix_path: Path
    checksum: str
    metadata: dict[str, Any] = field(default_factory=dict)
