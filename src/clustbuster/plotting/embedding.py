"""Interactive dimensional-reduction figure builder."""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import plotly.graph_objects as go

from clustbuster.models import Workspace


def embedding_figure(workspace: Workspace, *, color_by: str = "cluster") -> go.Figure:
    if workspace.cluster_column is None or workspace.embedding_key is None:
        raise ValueError("The workspace must have a cluster column and embedding")
    coordinates = np.asarray(workspace.adata.obsm[workspace.embedding_key])
    clusters = workspace.adata.obs[workspace.cluster_column].tolist()
    if color_by == "cluster":
        labels = ["<missing>" if value is None else str(value) for value in clusters]
        legend_title = workspace.cluster_column
    elif color_by == "annotation":
        labels = workspace.annotations.materialize(clusters)
        legend_title = "Annotation"
    else:
        raise ValueError(f"Unsupported embedding color mode: {color_by}")

    grouped_indices: dict[str, list[int]] = defaultdict(list)
    for index, label in enumerate(labels):
        grouped_indices[label].append(index)

    figure = go.Figure()
    cell_ids = workspace.adata.obs_names.astype(str).tolist()
    for label, indices in grouped_indices.items():
        figure.add_trace(
            go.Scattergl(
                x=coordinates[indices, 0],
                y=coordinates[indices, 1],
                mode="markers",
                name=label,
                customdata=[[cell_ids[index], str(clusters[index])] for index in indices],
                hovertemplate="Cell: %{customdata[0]}<br>Cluster: %{customdata[1]}<extra></extra>",
                marker={"size": 6, "opacity": 0.78},
            )
        )
    embedding_name = workspace.embedding_key.removeprefix("X_").upper()
    figure.update_layout(
        template="plotly_white",
        margin={"l": 48, "r": 20, "t": 35, "b": 48},
        legend={"title": {"text": legend_title}, "itemsizing": "constant"},
        xaxis_title=f"{embedding_name} 1",
        yaxis_title=f"{embedding_name} 2",
        dragmode="lasso",
    )
    return figure

