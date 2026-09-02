"""Plotly/WebGL module-score embedding builder."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go

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


def module_score_violin_figure(result: ModuleScoreResult) -> go.Figure:
    figure = go.Figure()
    cluster_order = result.cluster_summary["cluster"].astype(str).tolist()
    for cluster in cluster_order:
        values = result.values.loc[
            result.values["cluster"].astype(str) == cluster, "module_score"
        ].to_numpy()
        figure.add_trace(
            go.Violin(
                x=[cluster] * len(values),
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
                hovertemplate=("Cluster: %{x}<br>Module score: %{y:.3f}<extra></extra>"),
                showlegend=False,
            )
        )
    figure.update_layout(
        title=f"{result.name} by source cluster",
        template="plotly_white",
        height=620,
        margin={"l": 60, "r": 30, "t": 55, "b": 70},
        xaxis_title="Source cluster",
        yaxis_title="Module score",
        violinmode="overlay",
    )
    return figure
