"""Sparse-safe differential marker ranking and heatmap aggregation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from scipy.stats import t as t_distribution

from clustbuster.core.annotations import ClusterIdentifier
from clustbuster.core.expression import expression_matrix
from clustbuster.models import ExpressionSource


class MarkerAnalysisError(ValueError):
    """Raised when a marker comparison cannot produce valid results."""


@dataclass(slots=True)
class MarkerResult:
    selected_cluster: str
    method: str
    values: pd.DataFrame
    heatmap: pd.DataFrame


@dataclass(slots=True)
class AllMarkerResult:
    method: str
    values: pd.DataFrame
    failures: tuple[str, ...]
    top_n_per_cluster: int


def _group_statistics(matrix: Any, mask: np.ndarray) -> tuple[np.ndarray, ...]:
    subset = matrix[mask, :]
    count = int(mask.sum())
    if sparse.issparse(subset):
        sums = np.asarray(subset.sum(axis=0)).ravel().astype(float)
        squared = subset.copy().astype(float)
        squared.data **= 2
        sum_squares = np.asarray(squared.sum(axis=0)).ravel()
        expressing = np.asarray((subset > 0).sum(axis=0)).ravel().astype(float)
    else:
        dense = np.asarray(subset, dtype=float)
        sums = dense.sum(axis=0)
        sum_squares = np.square(dense).sum(axis=0)
        expressing = (dense > 0).sum(axis=0).astype(float)
    means = sums / count
    variances = np.maximum((sum_squares - (sums**2 / count)) / (count - 1), 0)
    fractions = expressing / count
    return means, variances, fractions


def _adjust_benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted_ranked = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.clip(adjusted_ranked, 0, 1)
    return adjusted


def _cluster_means(
    matrix: Any,
    clusters: list[object],
    feature_indices: np.ndarray,
    feature_names: list[str],
) -> pd.DataFrame:
    groups: dict[str, tuple[ClusterIdentifier, list[int]]] = {}
    for index, raw_cluster in enumerate(clusters):
        cluster_id = ClusterIdentifier.from_value(raw_cluster)
        groups.setdefault(cluster_id.serialized, (cluster_id, []))[1].append(index)
    rows: list[np.ndarray] = []
    labels: list[str] = []
    for cluster_id, cell_indices in groups.values():
        subset = matrix[cell_indices, :][:, feature_indices]
        means = np.asarray(subset.mean(axis=0)).ravel()
        rows.append(means)
        labels.append(cluster_id.display)
    return pd.DataFrame(rows, index=labels, columns=feature_names)


def rank_markers(
    adata: AnnData,
    source: ExpressionSource,
    cluster_column: str,
    selected_cluster: str,
    *,
    top_n: int = 15,
    min_fraction: float = 0.1,
    min_log_fold_change: float = 0.0,
) -> MarkerResult:
    """Rank positive markers for one cluster versus all other cells using Welch's t-test."""

    if cluster_column not in adata.obs:
        raise MarkerAnalysisError(f"Unknown cluster column: {cluster_column}")
    if top_n < 1 or top_n > 100:
        raise MarkerAnalysisError("Top gene count must be between 1 and 100")
    if not 0 <= min_fraction <= 1:
        raise MarkerAnalysisError("Minimum expressing fraction must be between 0 and 1")
    if min_log_fold_change < 0:
        raise MarkerAnalysisError("Minimum log fold-change must be zero or greater")

    clusters = adata.obs[cluster_column].tolist()
    serialized = np.array(
        [ClusterIdentifier.from_value(value).serialized for value in clusters], dtype=object
    )
    selected_mask = serialized == selected_cluster
    rest_mask = ~selected_mask
    selected_count = int(selected_mask.sum())
    rest_count = int(rest_mask.sum())
    if selected_count < 2 or rest_count < 2:
        raise MarkerAnalysisError(
            "Marker ranking requires at least two cells in the selected cluster and comparison"
        )

    matrix, names = expression_matrix(adata, source)
    selected_mean, selected_var, selected_fraction = _group_statistics(matrix, selected_mask)
    rest_mean, rest_var, rest_fraction = _group_statistics(matrix, rest_mask)
    standard_error_squared = selected_var / selected_count + rest_var / rest_count
    standard_error = np.sqrt(standard_error_squared)
    mean_difference = selected_mean - rest_mean
    scores = np.divide(
        mean_difference,
        standard_error,
        out=np.zeros_like(mean_difference),
        where=standard_error > 0,
    )
    denominator = (selected_var / selected_count) ** 2 / (selected_count - 1)
    denominator += (rest_var / rest_count) ** 2 / (rest_count - 1)
    degrees_freedom = np.divide(
        standard_error_squared**2,
        denominator,
        out=np.ones_like(denominator),
        where=denominator > 0,
    )
    p_values = 2 * t_distribution.sf(np.abs(scores), degrees_freedom)
    p_values = np.where(np.isfinite(p_values), p_values, 1.0)
    adjusted = _adjust_benjamini_hochberg(p_values)

    table = pd.DataFrame(
        {
            "gene": names.astype(str),
            "score": scores,
            "p_value": p_values,
            "p_adjusted": adjusted,
            "mean_selected": selected_mean,
            "mean_rest": rest_mean,
            "mean_difference": mean_difference,
            "log_fold_change": mean_difference,
            "fraction_selected": selected_fraction,
            "fraction_rest": rest_fraction,
        }
    )
    table = table.loc[
        (table["fraction_selected"] >= min_fraction)
        & (table["log_fold_change"] >= min_log_fold_change)
        & (table["mean_difference"] > 0)
    ].sort_values(["p_adjusted", "score"], ascending=[True, False])
    table = table.head(top_n).reset_index(drop=True)
    if table.empty:
        raise MarkerAnalysisError("No positive markers passed the selected expression threshold")
    table.insert(0, "rank", np.arange(1, len(table) + 1))

    feature_indices = names.get_indexer(table["gene"])
    heatmap = _cluster_means(
        matrix,
        clusters,
        feature_indices,
        table["gene"].astype(str).tolist(),
    )
    selected_display = ClusterIdentifier.from_value(
        clusters[int(np.flatnonzero(selected_mask)[0])]
    ).display
    return MarkerResult(selected_display, "Welch t-test", table, heatmap)


def rank_all_markers(
    adata: AnnData,
    source: ExpressionSource,
    cluster_column: str,
    *,
    top_n_per_cluster: int = 25,
    min_fraction: float = 0.1,
    min_log_fold_change: float = 0.25,
) -> AllMarkerResult:
    """Rank positive one-vs-rest markers for every source cluster."""

    if cluster_column not in adata.obs:
        raise MarkerAnalysisError(f"Unknown cluster column: {cluster_column}")
    identifiers: dict[str, ClusterIdentifier] = {}
    for raw_cluster in adata.obs[cluster_column].tolist():
        identifier = ClusterIdentifier.from_value(raw_cluster)
        identifiers.setdefault(identifier.serialized, identifier)

    frames: list[pd.DataFrame] = []
    failures: list[str] = []
    for serialized, identifier in identifiers.items():
        try:
            result = rank_markers(
                adata,
                source,
                cluster_column,
                serialized,
                top_n=top_n_per_cluster,
                min_fraction=min_fraction,
                min_log_fold_change=min_log_fold_change,
            )
        except MarkerAnalysisError as exc:
            failures.append(f"Cluster {identifier.display}: {exc}")
            continue
        values = result.values.copy()
        values.insert(1, "cluster_id", serialized)
        values.insert(2, "cluster", identifier.display)
        frames.append(values)

    if not frames:
        details = "; ".join(failures)
        raise MarkerAnalysisError(f"No clusters produced marker genes. {details}")
    return AllMarkerResult(
        method="Welch t-test (one-vs-rest)",
        values=pd.concat(frames, ignore_index=True),
        failures=tuple(failures),
        top_n_per_cluster=top_n_per_cluster,
    )
