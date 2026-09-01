"""Small hierarchical-clustering helpers shared by matrix-style plots."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage


@dataclass(frozen=True, slots=True)
class Hierarchy:
    order: tuple[int, ...]
    segments: tuple[tuple[tuple[float, ...], tuple[float, ...]], ...]


def cluster_hierarchy(observations: np.ndarray) -> Hierarchy:
    """Return optimal leaf order and Plotly-ready dendrogram line segments."""

    matrix = np.asarray(observations, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("Hierarchical clustering requires a two-dimensional matrix")
    if matrix.shape[0] < 2:
        return Hierarchy(tuple(range(matrix.shape[0])), ())
    tree = dendrogram(
        linkage(matrix, method="average", metric="euclidean", optimal_ordering=True),
        no_plot=True,
        color_threshold=0,
    )
    segments = tuple(
        (tuple(float(value) for value in leaves), tuple(float(value) for value in distances))
        for leaves, distances in zip(tree["icoord"], tree["dcoord"], strict=True)
    )
    return Hierarchy(tuple(int(value) for value in tree["leaves"]), segments)


def standardize_columns(values: np.ndarray) -> np.ndarray:
    """Z-score columns without producing NaNs for constant features."""

    matrix = np.asarray(values, dtype=float)
    means = matrix.mean(axis=0)
    deviations = matrix.std(axis=0)
    result = np.zeros_like(matrix)
    variable = deviations > 0
    result[:, variable] = (matrix[:, variable] - means[variable]) / deviations[variable]
    return result
