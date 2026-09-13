"""Clustered reference-correlation heatmap with row and column dendrograms."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from clustbuster.integrations.pyclustifyr import ReferenceAnnotationResult
from clustbuster.plotting.hierarchy import cluster_hierarchy


def reference_correlation_height(cluster_count: int) -> int:
    return min(1100, max(650, 28 * cluster_count + 280))


def reference_correlation_figure(
    result: ReferenceAnnotationResult, labels: dict[str, str] | None = None,
) -> go.Figure:
    display_by_id = {
        prediction.cluster_id: prediction.cluster_display for prediction in result.predictions
    }
    display_by_id.update(labels or {})
    correlations = result.correlations
    matrix = correlations.to_numpy(dtype=float)
    row_hierarchy = cluster_hierarchy(matrix)
    column_hierarchy = cluster_hierarchy(matrix.T)
    row_order = list(row_hierarchy.order)
    column_order = list(column_hierarchy.order)
    ordered = matrix[np.ix_(row_order, column_order)]
    cluster_labels = [
        display_by_id.get(str(correlations.index[index]), str(correlations.index[index]))
        for index in row_order
    ]
    reference_labels = [str(correlations.columns[index]) for index in column_order]
    x_positions = [5 + 10 * index for index in range(len(reference_labels))]
    y_positions = [5 + 10 * index for index in range(len(cluster_labels))]
    customdata = np.empty((len(cluster_labels), len(reference_labels), 2), dtype=object)
    customdata[:, :, 0] = np.asarray(cluster_labels, dtype=object)[:, None]
    customdata[:, :, 1] = np.asarray(reference_labels, dtype=object)[None, :]

    figure = make_subplots(
        rows=2,
        cols=2,
        specs=[[None, {}], [{}, {}]],
        row_heights=[0.16, 0.84],
        column_widths=[0.02, 0.98],
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
            colorscale="Viridis",
            zmin=0,
            zmax=1,
            colorbar={"title": result.parameters.compute_method.title(),
                      "tickvals": [0, 0.25, 0.5, 0.75, 1]},
            hovertemplate=(
                "Cluster: %{customdata[0]}<br>Reference type: %{customdata[1]}"
                "<br>Similarity: %{z:.3f}<extra></extra>"
            ),
        ),
        row=2,
        col=2,
    )
    # Apply consumes these predictions directly; do not independently choose an argmax.
    predictions = {prediction.cluster_id: prediction for prediction in result.predictions}
    star_x: list[int] = []
    star_y: list[int] = []
    star_data: list[list[str]] = []
    for row_position, original_row in enumerate(row_order):
        prediction = predictions.get(str(correlations.index[original_row]))
        if prediction is None or prediction.annotation == "unassigned":
            continue
        annotation = prediction.annotation
        selected = {annotation}
        if annotation not in reference_labels and annotation.endswith("-CLASH!"):
            selected = set(annotation.removesuffix("-CLASH!").split("; "))
        for column_position, label in enumerate(reference_labels):
            if label in selected:
                star_x.append(x_positions[column_position])
                star_y.append(y_positions[row_position])
                star_data.append([display_by_id[prediction.cluster_id], label, annotation])
    figure.add_trace(
        go.Scatter(
            x=star_x, y=star_y, mode="markers", customdata=star_data,
            marker={"symbol": "star", "size": 13, "color": "#ffffff",
                    "line": {"color": "#17212b", "width": 1.2}},
            name="Apply selection", showlegend=False,
            hovertemplate=("Cluster: %{customdata[0]}<br>Reference type: %{customdata[1]}"
                           "<br>Apply annotation: %{customdata[2]}<extra></extra>"),
        ), row=2, col=2,
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
    figure.update_layout(
        title=f"Similarity to {result.reference.summary.name}",
        template="plotly_white",
        height=reference_correlation_height(len(cluster_labels)),
        margin={"l": 35, "r": 90, "t": 70, "b": 130},
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
        title="Reference cell type",
        tickmode="array",
        tickvals=x_positions,
        ticktext=reference_labels,
        tickangle=-35,
        showticklabels=True,
        automargin=True,
        row=2,
        col=2,
    )
    figure.update_yaxes(
        title="Source cluster",
        tickmode="array",
        tickvals=y_positions,
        ticktext=cluster_labels,
        showticklabels=True,
        side="left",
        ticklabelposition="outside left",
        automargin=True,
        autorange="reversed",
        row=2,
        col=2,
    )
    return figure
