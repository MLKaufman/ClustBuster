"""Plotly cluster-by-gene marker heatmap builder."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go

from clustbuster.core.markers import MarkerResult


def marker_heatmap_figure(result: MarkerResult) -> go.Figure:
    values = result.heatmap.to_numpy(dtype=float)
    means = values.mean(axis=0)
    standard_deviations = values.std(axis=0)
    variable = standard_deviations > 0
    standardized = np.zeros_like(values)
    standardized[:, variable] = (
        values[:, variable] - means[variable]
    ) / standard_deviations[variable]
    figure = go.Figure(
        go.Heatmap(
            z=standardized,
            x=result.heatmap.columns.tolist(),
            y=result.heatmap.index.tolist(),
            colorscale="RdBu_r",
            zmid=0,
            colorbar={"title": "Gene z-score"},
            hovertemplate=(
                "Cluster: %{y}<br>Gene: %{x}<br>Standardized mean: "
                "%{z:.3f}<extra></extra>"
            ),
        )
    )
    figure.update_layout(
        title=f"Top markers for cluster {result.selected_cluster}",
        template="plotly_white",
        margin={"l": 80, "r": 70, "t": 55, "b": 90},
        xaxis={"tickangle": -45},
        yaxis={"autorange": "reversed", "title": "Cluster"},
    )
    return figure
