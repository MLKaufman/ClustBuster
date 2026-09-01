"""Static feature-plot builder."""

from __future__ import annotations

import numpy as np
from matplotlib.figure import Figure


def feature_gene_figure(
    coordinates: np.ndarray,
    gene: str,
    expression: np.ndarray,
    annotations: tuple[tuple[float, float, str], ...] = (),
) -> Figure:
    """Return a raster-friendly Matplotlib figure without browser-side point data."""

    order = np.argsort(expression, kind="stable")
    figure = Figure(figsize=(7.2, 5.2))
    axes = figure.subplots()
    figure.subplots_adjust(left=0.11, right=0.88, bottom=0.12, top=0.93)
    points = axes.scatter(
        coordinates[order, 0],
        coordinates[order, 1],
        c=expression[order],
        cmap="viridis",
        s=6,
        alpha=0.8,
        edgecolors="none",
        rasterized=True,
    )
    colorbar = figure.colorbar(points, ax=axes, pad=0.02)
    colorbar.set_label("Expression")
    for x_position, y_position, label in annotations:
        axes.annotate(
            label,
            (x_position, y_position),
            ha="center",
            va="center",
            color="#19324a",
            fontsize=11,
            fontweight="bold",
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "edgecolor": "#738695",
                "alpha": 0.82,
            },
        )
    axes.set_xlabel("Dimension 1")
    axes.set_ylabel("Dimension 2")
    axes.set_aspect("equal", adjustable="datalim")
    axes.set_title(gene)
    axes.spines[["top", "right"]].set_visible(False)
    return figure
