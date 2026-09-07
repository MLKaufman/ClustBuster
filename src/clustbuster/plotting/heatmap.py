"""Plotly cluster-by-gene marker heatmap with row and column dendrograms."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from anndata import AnnData
from matplotlib.figure import Figure
from plotly.subplots import make_subplots
from scipy import sparse

from clustbuster.core.annotations import ClusterIdentifier, _natural_cluster_key
from clustbuster.core.expression import expression_matrix
from clustbuster.core.markers import AllMarkerResult, MarkerResult
from clustbuster.models import ExpressionSource
from clustbuster.plotting.hierarchy import cluster_hierarchy, standardize_columns


def marker_heatmap_height(cluster_count: int) -> int:
    return min(860, max(560, 24 * cluster_count + 220))


def all_marker_heatmap_height(gene_count: int) -> int:
    return max(700, 24 * gene_count + 240)


def all_marker_heatmap_width(cluster_count: int) -> int:
    return max(1200, 48 * cluster_count + 240)


def all_marker_heatmap_figure(
    adata: AnnData,
    source: ExpressionSource,
    cluster_column: str,
    result: AllMarkerResult,
    labels: dict[str, str] | None = None,
) -> Figure:
    """Build a static Seurat-style gene-by-cell heatmap grouped by source cluster."""

    marker_rows = result.values.drop_duplicates("gene", keep="first")
    genes = marker_rows["gene"].astype(str).tolist()
    matrix, names = expression_matrix(adata, source)
    feature_indices = names.get_indexer(genes)
    if (feature_indices < 0).any():
        raise ValueError("All-marker heatmap genes are unavailable in the expression source")

    raw_clusters = adata.obs[cluster_column].tolist()
    serialized = np.asarray(
        [ClusterIdentifier.from_value(value).serialized for value in raw_clusters], dtype=object
    )
    cluster_rows = result.values.drop_duplicates("cluster_id", keep="first").copy()
    cluster_rows["_order"] = cluster_rows["cluster"].map(
        lambda value: _natural_cluster_key(ClusterIdentifier.from_value(value))
    )
    cluster_rows = cluster_rows.sort_values("_order", kind="stable")
    cluster_ids = cluster_rows["cluster_id"].astype(str).tolist()
    cluster_labels = [(labels or {}).get(str(row.cluster_id), str(row.cluster))
                      for row in cluster_rows.itertuples()]
    cell_groups = [np.flatnonzero(serialized == cluster_id) for cluster_id in cluster_ids]
    cell_order = np.concatenate(cell_groups)
    selected = matrix[cell_order, :][:, feature_indices]
    dense = selected.toarray() if sparse.issparse(selected) else np.asarray(selected)
    gene_by_cell = np.asarray(dense, dtype=float).T
    means = gene_by_cell.mean(axis=1, keepdims=True)
    deviations = gene_by_cell.std(axis=1, keepdims=True)
    scaled = np.divide(
        gene_by_cell - means,
        deviations,
        out=np.zeros_like(gene_by_cell),
        where=deviations > 0,
    )
    scaled = np.clip(scaled, -2.5, 2.5)

    figure_height = max(6.4, all_marker_heatmap_height(len(genes)) / 100)
    figure = Figure(figsize=(all_marker_heatmap_width(len(cluster_ids)) / 100, figure_height))
    heatmap_axes = figure.subplots()
    boundaries = np.cumsum([len(indices) for indices in cell_groups])
    starts = np.concatenate(([0], boundaries[:-1]))
    centers = (starts + boundaries - 1) / 2
    image = heatmap_axes.imshow(
        scaled,
        aspect="auto",
        interpolation="nearest",
        cmap="RdBu_r",
        vmin=-2.5,
        vmax=2.5,
        rasterized=True,
    )
    heatmap_axes.set_yticks(np.arange(len(genes)), labels=genes, fontsize=8)
    heatmap_axes.set_xticks(centers, labels=cluster_labels, rotation=90)
    heatmap_axes.xaxis.tick_top()
    heatmap_axes.tick_params(axis="x", length=0, pad=5)
    heatmap_axes.set_ylabel("Marker gene")
    heatmap_axes.set_title(
        f"Top {result.top_n_per_cluster} markers per source cluster", pad=28
    )
    for boundary in boundaries[:-1]:
        heatmap_axes.axvline(boundary - 0.5, color="#263746", linewidth=0.65, alpha=0.8)
    figure.colorbar(image, ax=heatmap_axes, pad=0.012, label="Scaled expression")
    figure.tight_layout(pad=1.5)
    return figure


def marker_heatmap_figure(result: MarkerResult, labels: dict[str, str] | None = None) -> go.Figure:
    standardized = standardize_columns(result.heatmap.to_numpy(dtype=float))
    row_hierarchy = cluster_hierarchy(standardized)
    column_hierarchy = cluster_hierarchy(standardized.T)
    row_order = list(row_hierarchy.order)
    column_order = list(column_hierarchy.order)
    ordered = standardized[np.ix_(row_order, column_order)]
    cluster_labels = [(labels or {}).get(str(result.heatmap.index[index]),
                                       str(result.heatmap.index[index])) for index in row_order]
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
    figure.update_layout(
        title=f"Top markers for cluster {result.selected_cluster}",
        template="plotly_white",
        height=marker_heatmap_height(len(cluster_labels)),
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
        showticklabels=True,
        automargin=True,
        side="bottom",
        row=2,
        col=2,
    )
    figure.update_yaxes(
        title="Cluster",
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
