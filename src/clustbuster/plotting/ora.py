"""Clustered pathway-by-cluster ORA heatmap."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go

from clustbuster.core.enrichment import AllClusterOraResult
from clustbuster.plotting.hierarchy import cluster_hierarchy


def ora_heatmap_height(pathway_count: int) -> int:
    return min(1400, max(700, 25 * pathway_count + 260))


def ora_heatmap_figure(result: AllClusterOraResult, *, top_n_pathways: int = 30) -> go.Figure:
    values = result.values.copy()
    values["significance"] = -np.log10(
        np.maximum(values["adjusted_p_value"].to_numpy(dtype=float), 1e-300)
    )
    selected_terms = (
        values.groupby("term", sort=False)["significance"]
        .max()
        .sort_values(ascending=False)
        .head(top_n_pathways)
        .index.tolist()
    )
    matrix = (
        values.loc[values["term"].isin(selected_terms)]
        .pivot_table(
            index="term",
            columns="cluster",
            values="significance",
            aggfunc="max",
            fill_value=0.0,
        )
        .reindex(index=selected_terms)
    )
    pathway_hierarchy = cluster_hierarchy(matrix.to_numpy(dtype=float))
    cluster_hierarchy_result = cluster_hierarchy(matrix.to_numpy(dtype=float).T)
    pathway_order = list(pathway_hierarchy.order)
    cluster_order = list(cluster_hierarchy_result.order)
    ordered = matrix.to_numpy(dtype=float)[np.ix_(pathway_order, cluster_order)]
    pathway_labels = [str(matrix.index[index]) for index in pathway_order]
    cluster_labels = [str(matrix.columns[index]) for index in cluster_order]
    customdata = np.empty((len(pathway_labels), len(cluster_labels), 2), dtype=object)
    customdata[:, :, 0] = np.asarray(pathway_labels, dtype=object)[:, None]
    customdata[:, :, 1] = np.asarray(cluster_labels, dtype=object)[None, :]

    figure = go.Figure(
        go.Heatmap(
            z=ordered,
            x=cluster_labels,
            y=pathway_labels,
            customdata=customdata,
            colorscale="YlGnBu",
            zmin=0,
            colorbar={"title": "-log10 adj. p"},
            hovertemplate=(
                "Pathway: %{customdata[0]}<br>Cluster: %{customdata[1]}"
                "<br>-log10 adjusted p: %{z:.3f}<extra></extra>"
            ),
        )
    )

    figure.update_layout(
        title=f"All-cluster ORA · {result.library.replace('_', ' ')}",
        template="plotly_white",
        height=ora_heatmap_height(len(pathway_labels)),
        margin={"l": 35, "r": 100, "t": 70, "b": 115},
    )
    figure.update_xaxes(
        title="Source cluster",
        tickangle=-45,
    )
    figure.update_yaxes(
        title="Pathway",
        showticklabels=True,
        side="left",
        ticklabelposition="outside left",
        automargin=True,
        autorange="reversed",
    )
    return figure
