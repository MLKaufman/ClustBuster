"""Cluster-by-gene dot-plot builder."""

import plotly.graph_objects as go

from clustbuster.core.expression import DotPlotResult


def dotplot_figure(result: DotPlotResult) -> go.Figure:
    values = result.values
    figure = go.Figure(
        go.Scatter(
            x=values["gene"],
            y=values["cluster_display"],
            mode="markers",
            customdata=values[["fraction_expressing", "mean_expression"]].to_numpy(),
            hovertemplate=(
                "Gene: %{x}<br>Cluster: %{y}<br>Fraction: %{customdata[0]:.1%}"
                "<br>Mean: %{customdata[1]:.3g}<extra></extra>"
            ),
            marker={
                "size": 6 + 28 * values["fraction_expressing"],
                "color": values["mean_expression"],
                "colorscale": "Viridis",
                "showscale": True,
                "colorbar": {"title": "Mean expression"},
                "line": {"width": 0.5, "color": "#34495e"},
            },
        )
    )
    figure.update_layout(
        template="plotly_white",
        height=max(360, 55 * values["cluster_display"].nunique()),
        margin={"l": 90, "r": 80, "t": 35, "b": 80},
        xaxis={"title": "Gene", "tickangle": -35},
        yaxis={"title": "Cluster", "autorange": "reversed"},
    )
    return figure

