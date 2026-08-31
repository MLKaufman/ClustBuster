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
                "Cell: %{customdata[0]}<br>Module score: "
                "%{customdata[1]:.3f}<extra></extra>"
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
        dragmode="lasso",
    )
    return figure
