"""Typed, non-mutating adapter for reference annotation with pyclustifyr."""

from __future__ import annotations

from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from typing import Any, cast

import numpy as np
import pandas as pd
from scipy import sparse

from clustbuster.core.annotations import ClusterIdentifier, PreviewedAnnotation
from clustbuster.core.expression import expression_matrix
from clustbuster.models import LoadedReference, Workspace

PYCLUSTIFYR_REVISION = "ea5ee58ac3702d49ab873b65691513996b2d2cb0"
SUPPORTED_METHODS = ("spearman", "pearson", "cosine")


class ReferenceAnnotationError(ValueError):
    """Raised when query/reference inputs cannot be compared safely."""


@dataclass(frozen=True, slots=True)
class ReferenceAnnotationParameters:
    compute_method: str = "spearman"
    minimum_gene_overlap: int = 5
    threshold: float = 0.0
    query_is_log_normalized: bool = True

    def __post_init__(self) -> None:
        if self.compute_method not in SUPPORTED_METHODS:
            raise ValueError(f"Unsupported similarity method: {self.compute_method}")
        if self.minimum_gene_overlap < 2:
            raise ValueError("Minimum gene overlap must be at least 2")
        if not -1 <= self.threshold <= 1:
            raise ValueError("Prediction threshold must be between -1 and 1")


@dataclass(slots=True)
class ReferenceAnnotationResult:
    predictions: tuple[PreviewedAnnotation, ...]
    correlations: pd.DataFrame
    reference: LoadedReference
    parameters: ReferenceAnnotationParameters
    matched_genes: tuple[str, ...]
    missing_genes: tuple[str, ...]
    ambiguous_genes: tuple[str, ...]
    warnings: tuple[str, ...]
    package_version: str
    package_revision: str = PYCLUSTIFYR_REVISION


def _read_reference(reference: LoadedReference) -> pd.DataFrame:
    suffix = reference.matrix_path.suffix.casefold()
    if suffix == ".csv":
        table = pd.read_csv(reference.matrix_path, keep_default_na=False)
    elif suffix in {".tsv", ".txt"}:
        table = pd.read_csv(reference.matrix_path, sep="\t", keep_default_na=False)
    elif suffix == ".parquet":
        table = pd.read_parquet(reference.matrix_path)
    else:
        raise ReferenceAnnotationError(f"Unsupported reference format: {suffix}")
    if "gene" not in table:
        raise ReferenceAnnotationError("Reference matrix does not contain a gene column")
    return table.set_index("gene")


def _match_reference_genes(
    expression_names: pd.Index[str], reference_names: pd.Index[str]
) -> tuple[list[int], list[str], list[str], list[str]]:
    exact: dict[str, list[int]] = {}
    folded: dict[str, list[int]] = {}
    for index, name in enumerate(expression_names.astype(str)):
        exact.setdefault(name, []).append(index)
        folded.setdefault(name.casefold(), []).append(index)

    indices: list[int] = []
    matched: list[str] = []
    missing: list[str] = []
    ambiguous: list[str] = []
    for reference_gene in reference_names.astype(str):
        candidates = exact.get(reference_gene, [])
        if not candidates:
            candidates = folded.get(reference_gene.casefold(), [])
        if not candidates:
            missing.append(reference_gene)
        elif len(candidates) > 1:
            ambiguous.append(reference_gene)
        else:
            indices.append(candidates[0])
            matched.append(reference_gene)
    return indices, matched, missing, ambiguous


def _cluster_average(
    matrix: Any,
    feature_indices: list[int],
    clusters: list[Any],
    matched_genes: list[str],
    *,
    query_is_log_normalized: bool,
) -> tuple[pd.DataFrame, dict[str, ClusterIdentifier], dict[str, int]]:
    selected = matrix[:, feature_indices]
    if sparse.issparse(selected):
        transformed = selected.copy().astype(float)
        if query_is_log_normalized:
            transformed.data = np.expm1(transformed.data)
    else:
        transformed = np.asarray(selected, dtype=float)
        if query_is_log_normalized:
            transformed = np.expm1(transformed)
    values = transformed.data if sparse.issparse(transformed) else transformed
    if not np.isfinite(values).all():
        raise ReferenceAnnotationError("Query expression contains non-finite values")

    groups: dict[str, list[int]] = {}
    identifiers: dict[str, ClusterIdentifier] = {}
    for index, raw_cluster in enumerate(clusters):
        cluster_id = ClusterIdentifier.from_value(raw_cluster)
        identifiers[cluster_id.serialized] = cluster_id
        groups.setdefault(cluster_id.serialized, []).append(index)

    averages: dict[str, np.ndarray] = {}
    for serialized, cell_indices in groups.items():
        mean = np.asarray(transformed[cell_indices, :].mean(axis=0)).ravel()
        averages[serialized] = np.log1p(mean) if query_is_log_normalized else mean
    frame = pd.DataFrame(averages, index=matched_genes)
    return frame, identifiers, {key: len(value) for key, value in groups.items()}


class PyClustifyrAdapter:
    """Run pyclustifyr scoring while preserving ClustBuster annotation state."""

    def annotate(
        self,
        workspace: Workspace,
        reference: LoadedReference,
        parameters: ReferenceAnnotationParameters,
    ) -> ReferenceAnnotationResult:
        if workspace.cluster_column is None:
            raise ReferenceAnnotationError("Configure a cluster column before annotation")

        reference_matrix = _read_reference(reference)
        matrix, expression_names = expression_matrix(
            workspace.adata, workspace.expression_source
        )
        indices, matched, missing, ambiguous = _match_reference_genes(
            expression_names, pd.Index(reference_matrix.index.astype(str))
        )
        if len(matched) < parameters.minimum_gene_overlap:
            raise ReferenceAnnotationError(
                f"Only {len(matched)} reference genes overlap the query; "
                f"at least {parameters.minimum_gene_overlap} are required"
            )

        query_average, identifiers, cluster_sizes = _cluster_average(
            matrix,
            indices,
            workspace.adata.obs[workspace.cluster_column].tolist(),
            matched,
            query_is_log_normalized=parameters.query_is_log_normalized,
        )
        reference_subset = reference_matrix.loc[matched].astype(float)

        from pyclustifyr import clustify

        scored = clustify(
            query_average,
            reference_subset,
            metadata=None,
            per_cell=True,
            compute_method=parameters.compute_method,
            if_log=False,
            verbose=False,
        )
        correlations = cast(pd.DataFrame, scored).astype(float)
        non_finite = ~np.isfinite(correlations.to_numpy())
        warnings: list[str] = []
        if non_finite.any():
            correlations = correlations.replace([np.inf, -np.inf], np.nan).fillna(0.0)
            warnings.append("Non-finite similarity scores were replaced with zero")
        if missing:
            warnings.append(f"{len(missing)} reference gene(s) were absent from the query")
        if ambiguous:
            warnings.append(
                f"{len(ambiguous)} reference gene(s) had ambiguous query matches and were omitted"
            )
        small = [identifiers[key].display for key, count in cluster_sizes.items() if count < 10]
        if small:
            warnings.append("Clusters with fewer than 10 cells: " + ", ".join(small))

        predictions: list[PreviewedAnnotation] = []
        for serialized, scores in correlations.iterrows():
            ordered = scores.sort_values(ascending=False, kind="stable")
            best_score = float(ordered.iloc[0])
            second_score = float(ordered.iloc[1]) if len(ordered) > 1 else best_score
            tied = ordered.index[np.isclose(ordered.to_numpy(), best_score)].astype(str).tolist()
            if best_score < parameters.threshold:
                annotation = "unassigned"
            elif len(tied) > 1:
                annotation = "; ".join(tied) + "-CLASH!"
                warnings.append(
                    f"Cluster {identifiers[str(serialized)].display} has a tied best match"
                )
            else:
                annotation = str(ordered.index[0])
            cluster_id = identifiers[str(serialized)]
            predictions.append(
                PreviewedAnnotation(
                    cluster_id=cluster_id.serialized,
                    cluster_display=cluster_id.display,
                    annotation=annotation,
                    confidence=best_score,
                    margin=best_score - second_score,
                )
            )

        try:
            package_version = version("pyclustifyr")
        except PackageNotFoundError:
            package_version = "unknown"
        return ReferenceAnnotationResult(
            predictions=tuple(predictions),
            correlations=correlations,
            reference=reference,
            parameters=parameters,
            matched_genes=tuple(matched),
            missing_genes=tuple(missing),
            ambiguous_genes=tuple(ambiguous),
            warnings=tuple(dict.fromkeys(warnings)),
            package_version=package_version,
        )
