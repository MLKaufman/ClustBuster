"""Sparse-safe expression extraction, gene matching, and cluster aggregation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from clustbuster.core.annotations import ClusterIdentifier
from clustbuster.models import ExpressionKind, ExpressionSource


class GeneMatchError(ValueError):
    """Raised when no requested genes can be matched safely."""


@dataclass(frozen=True, slots=True)
class GeneMatchReport:
    matched: tuple[str, ...]
    missing: tuple[str, ...]
    ambiguous: tuple[str, ...]


@dataclass(slots=True)
class ExpressionResult:
    values: pd.DataFrame
    report: GeneMatchReport


@dataclass(slots=True)
class DotPlotResult:
    values: pd.DataFrame
    report: GeneMatchReport


def parse_gene_list(text: str, *, limit: int = 24) -> tuple[str, ...]:
    tokens = [token.strip() for token in text.replace(",", " ").split() if token.strip()]
    unique = tuple(dict.fromkeys(tokens))
    if not unique:
        raise GeneMatchError("Enter at least one gene symbol")
    if len(unique) > limit:
        raise GeneMatchError(f"Gene panels are limited to {limit} genes")
    return unique


def _matrix_and_names(
    adata: AnnData, source: ExpressionSource
) -> tuple[Any, pd.Index[str]]:
    if source.kind is ExpressionKind.X:
        return adata.X, pd.Index(adata.var_names.astype(str))
    if source.kind is ExpressionKind.RAW:
        if adata.raw is None:
            raise ValueError("The selected raw expression source is unavailable")
        return adata.raw.X, pd.Index(adata.raw.var_names.astype(str))
    assert source.layer is not None
    if source.layer not in adata.layers:
        raise ValueError(f"The selected expression layer is unavailable: {source.layer}")
    return adata.layers[source.layer], pd.Index(adata.var_names.astype(str))


def expression_matrix(
    adata: AnnData, source: ExpressionSource
) -> tuple[Any, pd.Index[str]]:
    """Return the selected matrix and its feature names without copying it."""

    return _matrix_and_names(adata, source)


def _match_indices(
    names: pd.Index[str], requested: tuple[str, ...]
) -> tuple[list[int], GeneMatchReport]:
    exact = {name: index for index, name in enumerate(names)}
    folded: dict[str, list[int]] = {}
    for index, name in enumerate(names):
        folded.setdefault(name.casefold(), []).append(index)

    indices: list[int] = []
    matched: list[str] = []
    missing: list[str] = []
    ambiguous: list[str] = []
    used_indices: set[int] = set()
    for query in requested:
        if query in exact:
            index = exact[query]
        else:
            candidates = folded.get(query.casefold(), [])
            if not candidates:
                missing.append(query)
                continue
            if len(candidates) > 1:
                ambiguous.append(query)
                continue
            index = candidates[0]
        if index in used_indices:
            continue
        used_indices.add(index)
        indices.append(index)
        matched.append(str(names[index]))
    report = GeneMatchReport(tuple(matched), tuple(missing), tuple(ambiguous))
    if not indices:
        details = ", ".join((*missing, *ambiguous))
        raise GeneMatchError(
            f"None of the requested genes matched the expression source: {details}"
        )
    return indices, report


def extract_expression(
    adata: AnnData, source: ExpressionSource, genes: tuple[str, ...]
) -> ExpressionResult:
    matrix, names = _matrix_and_names(adata, source)
    indices, report = _match_indices(names, genes)
    selected = matrix[:, indices]
    values = selected.toarray() if sparse.issparse(selected) else np.asarray(selected)
    frame = pd.DataFrame(values, index=adata.obs_names.astype(str), columns=report.matched)
    return ExpressionResult(frame, report)


def aggregate_dotplot(
    adata: AnnData,
    source: ExpressionSource,
    cluster_column: str,
    genes: tuple[str, ...],
) -> DotPlotResult:
    if cluster_column not in adata.obs:
        raise ValueError(f"Unknown cluster column: {cluster_column}")
    expression = extract_expression(adata, source, genes)
    clusters = adata.obs[cluster_column].tolist()
    groups: dict[str, tuple[ClusterIdentifier, list[int]]] = {}
    for index, value in enumerate(clusters):
        cluster_id = ClusterIdentifier.from_value(value)
        groups.setdefault(cluster_id.serialized, (cluster_id, []))[1].append(index)

    rows: list[dict[str, object]] = []
    data = expression.values.to_numpy()
    for cluster_id, indices in groups.values():
        subset = data[indices, :]
        means = np.asarray(subset.mean(axis=0)).ravel()
        fractions = np.asarray((subset > 0).mean(axis=0)).ravel()
        for gene_index, gene in enumerate(expression.report.matched):
            rows.append(
                {
                    "cluster_id": cluster_id.serialized,
                    "cluster_display": cluster_id.display,
                    "gene": gene,
                    "fraction_expressing": float(fractions[gene_index]),
                    "mean_expression": float(means[gene_index]),
                }
            )
    return DotPlotResult(pd.DataFrame(rows), expression.report)
