"""Hierarchically clustered cluster-by-gene dot-plot builder."""

from __future__ import annotations

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from clustbuster.core.expression import DotPlotResult
from clustbuster.plotting.hierarchy import cluster_hierarchy, standardize_columns


def dotplot_figure(result: DotPlotResult) -> go.Figure:
    values = result.values.copy()
    cluster_ids = list(dict.fromkeys(values["cluster_id"].astype(str)))
    genes = list(dict.fromkeys(values["gene"].astype(str)))
    cluster_labels = (
        values.drop_duplicates("cluster_id")
        .set_index("cluster_id")["cluster_display"]
        .astype(str)
    )
    mean_matrix = (
        values.pivot(index="cluster_id", columns="gene", values="mean_expression")
        .reindex(index=cluster_ids, columns=genes)
        .to_numpy(dtype=float)
    )
    clustering_matrix = standardize_columns(mean_matrix)
    row_hierarchy = cluster_hierarchy(clustering_matrix)
    column_hierarchy = cluster_hierarchy(clustering_matrix.T)
    ordered_cluster_ids = [cluster_ids[index] for index in row_hierarchy.order]
    ordered_genes = [genes[index] for index in column_hierarchy.order]
    x_positions = {gene: 5 + 10 * index for index, gene in enumerate(ordered_genes)}
    y_positions = {
        cluster_id: 5 + 10 * index for index, cluster_id in enumerate(ordered_cluster_ids)
    }
    values["_x"] = values["gene"].astype(str).map(x_positions)
    values["_y"] = values["cluster_id"].astype(str).map(y_positions)

    figure = make_subplots(
        rows=2,
        cols=2,
        specs=[[None, {}], [{}, {}]],
        row_heights=[0.16, 0.84],
        column_widths=[0.16, 0.84],
        horizontal_spacing=0.01,
        vertical_spacing=0.01,
        shared_xaxes="columns",
        shared_yaxes="rows",
    )
    figure.add_trace(
        go.Scatter(
            x=values["_x"],
            y=values["_y"],
            mode="markers",
            customdata=values[
                ["gene", "cluster_display", "fraction_expressing", "mean_expression"]
            ].to_numpy(),
            hovertemplate=(
                "Gene: %{customdata[0]}<br>Cluster: %{customdata[1]}"
                "<br>Fraction: %{customdata[2]:.1%}"
                "<br>Mean: %{customdata[3]:.3g}<extra></extra>"
            ),
            marker={
                "size": 6 + 28 * values["fraction_expressing"],
                "color": values["mean_expression"],
                "colorscale": "Viridis",
                "showscale": True,
                "colorbar": {"title": "Mean expression"},
                "line": {"width": 0.5, "color": "#34495e"},
            },
            showlegend=False,
        ),
        row=2,
        col=2,
    )
    for leaves, distances in column_hierarchy.segments:
        figure.add_trace(
            go.Scatter(
                x=leaves,
                y=distances,
                mode="lines",
                line={"color": "#526b7a", "width": 1.4},
                hoverinfo="skip",
                showlegend=False,
            ),
            row=1,
            col=2,
        )
    for leaves, distances in row_hierarchy.segments:
        figure.add_trace(
            go.Scatter(
                x=distances,
                y=leaves,
                mode="lines",
                line={"color": "#526b7a", "width": 1.4},
                hoverinfo="skip",
                showlegend=False,
            ),
            row=2,
            col=1,
        )

    cluster_count = len(ordered_cluster_ids)
    figure.update_layout(
        template="plotly_white",
        height=min(1060, max(620, 28 * cluster_count + 240)),
        margin={"l": 120, "r": 95, "t": 30, "b": 120},
    )
    figure.update_xaxes(showticklabels=False, showgrid=False, zeroline=False, row=1, col=2)
    figure.update_yaxes(showticklabels=False, showgrid=False, zeroline=False, row=1, col=2)
    figure.update_xaxes(
        autorange="reversed",
        showticklabels=False,
        showgrid=False,
        zeroline=False,
        row=2,
        col=1,
    )
    figure.update_yaxes(
        title="Cluster",
        tickmode="array",
        tickvals=list(y_positions.values()),
        ticktext=[str(cluster_labels.loc[cluster_id]) for cluster_id in ordered_cluster_ids],
        autorange="reversed",
        showgrid=False,
        zeroline=False,
        row=2,
        col=1,
    )
    figure.update_xaxes(
        title="Gene",
        tickmode="array",
        tickvals=list(x_positions.values()),
        ticktext=ordered_genes,
        tickangle=-45,
        row=2,
        col=2,
    )
    figure.update_yaxes(
        showticklabels=False,
        autorange="reversed",
        row=2,
        col=2,
    )
    return figure
