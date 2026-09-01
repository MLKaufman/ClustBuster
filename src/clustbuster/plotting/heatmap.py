"""Plotly cluster-by-gene marker heatmap with row and column dendrograms."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from clustbuster.core.markers import MarkerResult
from clustbuster.plotting.hierarchy import cluster_hierarchy, standardize_columns


def marker_heatmap_figure(result: MarkerResult) -> go.Figure:
    standardized = standardize_columns(result.heatmap.to_numpy(dtype=float))
    row_hierarchy = cluster_hierarchy(standardized)
    column_hierarchy = cluster_hierarchy(standardized.T)
    row_order = list(row_hierarchy.order)
    column_order = list(column_hierarchy.order)
    ordered = standardized[np.ix_(row_order, column_order)]
    cluster_labels = [str(result.heatmap.index[index]) for index in row_order]
    gene_labels = [str(result.heatmap.columns[index]) for index in column_order]
    x_positions = [5 + 10 * index for index in range(len(gene_labels))]
    y_positions = [5 + 10 * index for index in range(len(cluster_labels))]
    customdata = np.empty((len(cluster_labels), len(gene_labels), 2), dtype=object)
    customdata[:, :, 0] = np.asarray(cluster_labels, dtype=object)[:, None]
    customdata[:, :, 1] = np.asarray(gene_labels, dtype=object)[None, :]

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
        go.Heatmap(
            z=ordered,
            x=x_positions,
            y=y_positions,
            customdata=customdata,
            colorscale="RdBu_r",
            zmid=0,
            colorbar={"title": "Gene z-score"},
            hovertemplate=(
                "Cluster: %{customdata[0]}<br>Gene: %{customdata[1]}"
                "<br>Standardized mean: %{z:.3f}<extra></extra>"
            ),
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

    figure.update_layout(
        title=f"Top markers for cluster {result.selected_cluster}",
        template="plotly_white",
        height=min(860, max(560, 24 * len(cluster_labels) + 220)),
        margin={"l": 35, "r": 90, "t": 65, "b": 115},
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
    figure.update_yaxes(showticklabels=False, showgrid=False, zeroline=False, row=2, col=1)
    figure.update_xaxes(
        tickmode="array",
        tickvals=x_positions,
        ticktext=gene_labels,
        tickangle=-45,
        row=2,
        col=2,
    )
    figure.update_yaxes(
        title="Cluster",
        tickmode="array",
        tickvals=y_positions,
        ticktext=cluster_labels,
        autorange="reversed",
        row=2,
        col=2,
    )
    return figure
