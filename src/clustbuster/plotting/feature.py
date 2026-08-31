"""Plotly/WebGL feature-plot builder."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from clustbuster.core.expression import ExpressionResult


def feature_figure(
    coordinates: np.ndarray, cell_ids: list[str], result: ExpressionResult
) -> go.Figure:
    genes = list(result.values.columns)
    figure = make_subplots(rows=1, cols=len(genes), subplot_titles=genes)
    for column, gene in enumerate(genes, start=1):
        expression = result.values[gene].to_numpy()
        figure.add_trace(
            go.Scattergl(
                x=coordinates[:, 0],
                y=coordinates[:, 1],
                mode="markers",
                name=gene,
                showlegend=False,
                customdata=np.column_stack([cell_ids, expression]),
                hovertemplate=(
                    "Cell: %{customdata[0]}<br>Expression: "
                    "%{customdata[1]:.3g}<extra></extra>"
                ),
                marker={
                    "size": 6,
                    "opacity": 0.8,
                    "color": expression,
                    "colorscale": "Viridis",
                    "showscale": column == len(genes),
                    "colorbar": {"title": "Expression"},
                },
            ),
            row=1,
            col=column,
        )
        figure.update_xaxes(title_text="Dimension 1", row=1, col=column)
        figure.update_yaxes(
            title_text="Dimension 2" if column == 1 else None,
            row=1,
            col=column,
            scaleanchor=f"x{column}" if column > 1 else "x",
            scaleratio=1,
        )
    figure.update_layout(
        template="plotly_white",
        height=560,
        margin={"l": 50, "r": 60, "t": 55, "b": 50},
        dragmode="lasso",
    )
    return figure
