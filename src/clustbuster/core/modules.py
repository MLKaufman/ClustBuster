"""Deterministic, sparse-safe gene-set scoring."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from anndata import AnnData

from clustbuster.core.annotations import ClusterIdentifier
from clustbuster.core.expression import GeneMatchReport, extract_expression
from clustbuster.models import ExpressionSource


@dataclass(frozen=True, slots=True)
class ModulePreset:
    key: str
    label: str
    genes: tuple[str, ...]


MODULE_PRESETS = (
    ModulePreset("t_cell", "T cell", ("CD3D", "IL7R")),
    ModulePreset("myeloid", "Myeloid", ("LYZ", "S100A8", "FCGR3A")),
    ModulePreset("b_cell", "B cell", ("MS4A1", "CD79A")),
    ModulePreset("nk_cell", "NK cell", ("NKG7", "GNLY")),
    ModulePreset("platelet", "Platelet", ("PPBP",)),
)


@dataclass(slots=True)
class ModuleScoreResult:
    name: str
    values: pd.DataFrame
    cluster_summary: pd.DataFrame
    report: GeneMatchReport


def calculate_module_score(
    adata: AnnData,
    source: ExpressionSource,
    cluster_column: str,
    genes: tuple[str, ...],
    *,
    name: str = "Custom module",
) -> ModuleScoreResult:
    """Return the mean of per-gene standardized expression for each cell.

    Only the requested gene columns are materialized. Constant genes contribute
    zero, which keeps the score finite and deterministic for sparse and dense data.
    """

    if cluster_column not in adata.obs:
        raise ValueError(f"Unknown cluster column: {cluster_column}")
    expression = extract_expression(adata, source, genes)
    matrix = expression.values.to_numpy(dtype=float, copy=True)
    means = matrix.mean(axis=0)
    standard_deviations = matrix.std(axis=0)
    variable = standard_deviations > 0
    standardized = np.zeros_like(matrix, dtype=float)
    standardized[:, variable] = (matrix[:, variable] - means[variable]) / standard_deviations[
        variable
    ]
    scores = standardized.mean(axis=1)

    groups: dict[str, tuple[ClusterIdentifier, list[int]]] = {}
    serialized_clusters: list[str] = []
    display_clusters: list[str] = []
    for index, raw_cluster in enumerate(adata.obs[cluster_column].tolist()):
        cluster_id = ClusterIdentifier.from_value(raw_cluster)
        serialized_clusters.append(cluster_id.serialized)
        display_clusters.append(cluster_id.display)
        groups.setdefault(cluster_id.serialized, (cluster_id, []))[1].append(index)

    values = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "cluster_id": serialized_clusters,
            "cluster": display_clusters,
            "module_score": scores,
        },
        index=adata.obs_names.astype(str),
    )

    summary_rows: list[dict[str, object]] = []
    for cluster_id, indices in groups.values():
        cluster_scores = scores[indices]
        summary_rows.append(
            {
                "cluster_id": cluster_id.serialized,
                "cluster": cluster_id.display,
                "cells": len(indices),
                "mean_score": float(cluster_scores.mean()),
                "median_score": float(np.median(cluster_scores)),
            }
        )
    return ModuleScoreResult(
        name=name.strip() or "Custom module",
        values=values,
        cluster_summary=pd.DataFrame(summary_rows),
        report=expression.report,
    )
