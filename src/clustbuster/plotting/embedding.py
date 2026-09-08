"""Interactive dimensional-reduction figure builder."""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import plotly.graph_objects as go

from clustbuster.models import Workspace


def embedding_figure(
    workspace: Workspace, *, color_by: str = "annotation", show_annotations: bool = True,
) -> go.Figure:
    if workspace.cluster_column is None or workspace.embedding_key is None:
        raise ValueError("The workspace must have a cluster column and embedding")
    coordinates = np.asarray(workspace.adata.obsm[workspace.embedding_key])
    clusters = workspace.adata.obs[workspace.cluster_column].tolist()
    records = [workspace.annotations.get(cluster) for cluster in clusters]
    source_labels = [record.cluster_id.display for record in records]
    current_annotations = [
        record.annotation.strip() or record.cluster_id.display for record in records
    ]
    if color_by == "cluster":
        labels = source_labels
        overlay_labels = source_labels
        legend_title = workspace.cluster_column
    elif color_by == "annotation":
        labels = current_annotations
        overlay_labels = current_annotations
        legend_title = "Annotation"
    else:
        raise ValueError(f"Unsupported embedding color mode: {color_by}")

    grouped_indices: dict[str, list[int]] = defaultdict(list)
    for index, label in enumerate(labels):
        grouped_indices[label].append(index)

    cluster_indices: dict[str, list[int]] = defaultdict(list)
    cluster_overlay_labels: dict[str, str] = {}
    for index, record in enumerate(records):
        cluster_key = record.cluster_id.serialized
        cluster_indices[cluster_key].append(index)
        cluster_overlay_labels[cluster_key] = overlay_labels[index]

    figure = go.Figure()
    cell_ids = workspace.adata.obs_names.astype(str).tolist()
    for label, indices in grouped_indices.items():
        figure.add_trace(
            go.Scattergl(
                x=coordinates[indices, 0],
                y=coordinates[indices, 1],
                mode="markers",
                name=label,
                customdata=[[cell_ids[index], labels[index]] for index in indices],
                hovertemplate="Cell: %{customdata[0]}<br>Cluster: %{customdata[1]}<extra></extra>",
                marker={"size": 6, "opacity": 0.78},
            )
        )
    for cluster_key, indices in cluster_indices.items():
        if not show_annotations:
            continue
        figure.add_annotation(
            x=float(np.median(coordinates[indices, 0])),
            y=float(np.median(coordinates[indices, 1])),
            text=cluster_overlay_labels[cluster_key],
            showarrow=False,
            font={"size": 15, "weight": 700, "color": "#19324a"},
            bgcolor="rgba(255, 255, 255, 0.82)",
            bordercolor="rgba(25, 50, 74, 0.35)",
            borderpad=3,
        )
    embedding_name = workspace.embedding_key.removeprefix("X_").upper()
    figure.update_layout(
        template="plotly_white",
        height=1160,
        autosize=True,
        margin={"l": 48, "r": 20, "t": 35, "b": 35},
        legend={"title": {"text": legend_title}, "itemsizing": "constant"},
        xaxis_title=f"{embedding_name} 1",
        yaxis_title=f"{embedding_name} 2",
        yaxis={"scaleanchor": "x", "scaleratio": 1},
        dragmode="lasso",
    )
    return figure
