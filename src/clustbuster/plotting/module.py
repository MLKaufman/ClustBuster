"""Plotly/WebGL module-score embedding builder."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from matplotlib.figure import Figure

from clustbuster.core.modules import ModuleScoreResult


def module_score_figure(
    coordinates: np.ndarray, cell_ids: list[str], result: ModuleScoreResult
) -> go.Figure:
    scores = result.values["module_score"].to_numpy()
    figure = go.Figure(
        go.Scattergl(
            x=coordinates[:, 0],
            y=coordinates[:, 1],
            mode="markers",
            customdata=np.column_stack([cell_ids, scores]),
            hovertemplate=(
                "Cell: %{customdata[0]}<br>Module score: %{customdata[1]:.3f}<extra></extra>"
            ),
            marker={
                "size": 7,
                "opacity": 0.82,
                "color": scores,
                "colorscale": "RdBu_r",
                "cmid": 0,
                "showscale": True,
                "colorbar": {"title": "Score"},
            },
        )
    )
    figure.update_layout(
        title=result.name,
        template="plotly_white",
        margin={"l": 48, "r": 70, "t": 50, "b": 48},
        xaxis_title="Dimension 1",
        yaxis_title="Dimension 2",
        yaxis={"scaleanchor": "x", "scaleratio": 1},
        dragmode="lasso",
    )
    return figure


def module_score_violin_figure(
    result: ModuleScoreResult, labels: dict[str, str] | None = None,
) -> go.Figure:
    figure = go.Figure()
    for position, row in enumerate(result.cluster_summary.itertuples()):
        cluster = (labels or {}).get(str(row.cluster_id), str(row.cluster))
        values = result.values.loc[
            result.values["cluster_id"].astype(str) == str(row.cluster_id), "module_score"
        ].to_numpy()
        figure.add_trace(
            go.Violin(
                x=[position] * len(values),
                text=[cluster] * len(values),
                y=values,
                name=cluster,
                box_visible=True,
                meanline_visible=True,
                points="all",
                jitter=0.35,
                pointpos=0,
                marker={"size": 5, "opacity": 0.62, "color": "#2d8c88"},
                line={"color": "#19324a"},
                fillcolor="rgba(45, 140, 136, 0.28)",
                hovertemplate=("Cluster: %{text}<br>Module score: %{y:.3f}<extra></extra>"),
                showlegend=False,
            )
        )
    figure.update_layout(
        title=f"{result.name} by cluster",
        template="plotly_white",
        height=620,
        margin={"l": 60, "r": 30, "t": 55, "b": 70},
        xaxis_title="Cluster",
        xaxis={"tickmode": "array",
               "tickvals": list(range(len(result.cluster_summary))),
               "ticktext": [(labels or {}).get(str(row.cluster_id), str(row.cluster))
                            for row in result.cluster_summary.itertuples()]},
        yaxis_title="Module score",
        violinmode="overlay",
    )
    return figure


def module_score_static_figure(
    coordinates: np.ndarray, result: ModuleScoreResult,
    annotations: tuple[tuple[float, float, str], ...] = (),
) -> Figure:
    """Render a square, raster-friendly module-score embedding."""
    scores = result.values["module_score"].to_numpy(dtype=float)
    order = np.argsort(scores, kind="stable")
    limit = max(float(np.abs(scores).max()), 1e-8)
    figure = Figure(figsize=(10, 10))
    axes = figure.subplots()
    points = axes.scatter(
        coordinates[order, 0], coordinates[order, 1], c=scores[order],
        cmap="RdBu_r", vmin=-limit, vmax=limit, s=6, alpha=0.82,
        edgecolors="none", rasterized=True,
    )
    for x, y, label in annotations:
        axes.text(x, y, label, ha="center", va="center", fontsize=9,
                  bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"})
    # Equal ranges keep the embedding square without distorting its coordinates.
    centers = (coordinates[:, :2].min(axis=0) + coordinates[:, :2].max(axis=0)) / 2
    half_span = max(float(np.ptp(coordinates[:, :2], axis=0).max()) * 0.55, 0.5)
    axes.set_xlim(centers[0] - half_span, centers[0] + half_span)
    axes.set_ylim(centers[1] - half_span, centers[1] + half_span)
    axes.set_aspect("equal", adjustable="box")
    axes.set_box_aspect(1)
    # Resizing must adjust the axes box, never the data limits (which can crop cells).
    axes.set_adjustable("box")
    axes.set_title(result.name)
    axes.set_xlabel("Dimension 1")
    axes.set_ylabel("Dimension 2")
    axes.spines[["top", "right"]].set_visible(False)
    figure.colorbar(points, ax=axes, pad=0.02, shrink=0.8, label="Module score")
    figure.set_layout_engine("tight", pad=2.5)
    return figure
