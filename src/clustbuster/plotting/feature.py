"""Plotly/WebGL feature-plot builder."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go


def feature_gene_figure(
    coordinates: np.ndarray,
    cell_ids: list[str],
    gene: str,
    expression: np.ndarray,
    annotations: tuple[tuple[float, float, str], ...] = (),
) -> go.Figure:
    figure = go.Figure(
        go.Scattergl(
            x=coordinates[:, 0],
            y=coordinates[:, 1],
            mode="markers",
            name=gene,
            showlegend=False,
            customdata=np.column_stack([cell_ids, expression]),
            hovertemplate=(
                "Cell: %{customdata[0]}<br>Expression: %{customdata[1]:.3g}<extra></extra>"
            ),
            marker={
                "size": 6,
                "opacity": 0.8,
                "color": expression,
                "colorscale": "Viridis",
                "showscale": True,
                "colorbar": {"title": "Expression"},
            },
        )
    )
    for x_position, y_position, label in annotations:
        figure.add_annotation(
            x=x_position,
            y=y_position,
            text=label,
            showarrow=False,
            font={"size": 15, "weight": 700, "color": "#19324a"},
            bgcolor="rgba(255, 255, 255, 0.82)",
            bordercolor="rgba(25, 50, 74, 0.35)",
            borderpad=3,
        )
    figure.update_layout(
        template="plotly_white",
        height=520,
        margin={"l": 50, "r": 80, "t": 25, "b": 50},
        dragmode="lasso",
        xaxis={"title": "Dimension 1"},
        yaxis={"title": "Dimension 2", "scaleanchor": "x", "scaleratio": 1},
    )
    return figure
