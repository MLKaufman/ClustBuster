"""Plotly enrichment summary builder."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go

from clustbuster.core.enrichment import EnrichmentResult


def enrichment_figure(result: EnrichmentResult, *, top_n: int = 12) -> go.Figure:
    table = result.values.head(top_n).iloc[::-1]
    significance = -np.log10(np.maximum(table["adjusted_p_value"].to_numpy(), 1e-300))
    figure = go.Figure(
        go.Bar(
            x=significance,
            y=table["term"],
            orientation="h",
            customdata=np.column_stack(
                [table["overlap_genes"], table["adjusted_p_value"]]
            ),
            hovertemplate=(
                "Term: %{y}<br>-log10 adjusted p: %{x:.3f}<br>Adjusted p: "
                "%{customdata[1]:.3g}<br>Genes: %{customdata[0]}<extra></extra>"
            ),
            marker={"color": significance, "colorscale": "Teal", "showscale": False},
        )
    )
    figure.update_layout(
        title=f"GO Biological Process enrichment · {result.source}",
        template="plotly_white",
        margin={"l": 260, "r": 30, "t": 55, "b": 50},
        xaxis_title="-log10 adjusted p-value",
        yaxis_title="",
    )
    return figure
