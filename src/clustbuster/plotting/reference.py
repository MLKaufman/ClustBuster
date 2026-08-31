"""Plotly reference-correlation heatmap."""

from __future__ import annotations

import plotly.graph_objects as go

from clustbuster.integrations.pyclustifyr import ReferenceAnnotationResult


def reference_correlation_figure(result: ReferenceAnnotationResult) -> go.Figure:
    display_by_id = {
        prediction.cluster_id: prediction.cluster_display
        for prediction in result.predictions
    }
    correlations = result.correlations
    figure = go.Figure(
        go.Heatmap(
            z=correlations.to_numpy(dtype=float),
            x=correlations.columns.astype(str).tolist(),
            y=[display_by_id.get(str(value), str(value)) for value in correlations.index],
            colorscale="RdBu_r",
            zmid=0,
            zmin=-1,
            zmax=1,
            colorbar={"title": result.parameters.compute_method.title()},
            hovertemplate=(
                "Cluster: %{y}<br>Reference type: %{x}<br>Similarity: "
                "%{z:.3f}<extra></extra>"
            ),
        )
    )
    figure.update_layout(
        title=f"Similarity to {result.reference.summary.name}",
        template="plotly_white",
        margin={"l": 80, "r": 70, "t": 55, "b": 100},
        xaxis={"tickangle": -35, "title": "Reference cell type"},
        yaxis={"autorange": "reversed", "title": "Source cluster"},
    )
    return figure
