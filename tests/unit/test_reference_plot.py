from types import SimpleNamespace
from typing import cast

import pandas as pd

from clustbuster.core.annotations import PreviewedAnnotation
from clustbuster.integrations.pyclustifyr import ReferenceAnnotationResult
from clustbuster.plotting.reference import (
    reference_correlation_figure,
    reference_correlation_height,
)


def test_reference_heatmap_clusters_both_axes_and_labels_matrix_side() -> None:
    predictions = tuple(
        PreviewedAnnotation(
            cluster_id=f"cluster-{index}",
            cluster_display=str(index),
            annotation="Type A",
            confidence=0.8,
            margin=0.2,
        )
        for index in range(3)
    )
    result = cast(
        ReferenceAnnotationResult,
        SimpleNamespace(
            predictions=predictions,
            correlations=pd.DataFrame(
                [[0.9, 0.1, -0.2], [0.2, 0.8, 0.0], [-0.1, 0.3, 0.7]],
                index=["cluster-0", "cluster-1", "cluster-2"],
                columns=["Type A", "Type B", "Type C"],
            ),
            parameters=SimpleNamespace(compute_method="spearman"),
            reference=SimpleNamespace(summary=SimpleNamespace(name="Demo reference")),
        ),
    )

    figure = reference_correlation_figure(result)

    assert figure.data[0].type == "heatmap"
    assert sum(trace.mode == "lines" for trace in figure.data[1:]) == 2
    assert figure.layout.height == reference_correlation_height(3)
    assert figure.layout.yaxis2.showticklabels is False
    assert figure.layout.yaxis3.showticklabels is True
    assert figure.layout.yaxis3.side == "left"
    assert set(figure.layout.yaxis3.ticktext) == {"0", "1", "2"}


    stars = next(trace for trace in figure.data[1:] if trace.mode == "markers")
    # All previews deliberately say Type A despite different row maxima.
    assert len(stars.x) == 3
    assert {row[1] for row in stars.customdata} == {"Type A"}
    x_labels = dict(zip(figure.layout.xaxis3.tickvals, figure.layout.xaxis3.ticktext, strict=True))
    y_labels = dict(zip(figure.layout.yaxis3.tickvals, figure.layout.yaxis3.ticktext, strict=True))
    for x, y, row in zip(stars.x, stars.y, stars.customdata, strict=True):
        assert x_labels[x] == row[1]
        assert y_labels[y] == row[0]

    result.predictions = (
        PreviewedAnnotation("cluster-0", "0", "unassigned", 0.1, 0.0),
        PreviewedAnnotation("cluster-1", "1", "Type A; Type B-CLASH!", 0.8, 0.0),
        predictions[2],
    )
    figure = reference_correlation_figure(result)
    stars = next(trace for trace in figure.data[1:] if trace.mode == "markers")
    assert {(row[0], row[1]) for row in stars.customdata} == {
        ("1", "Type A"), ("1", "Type B"), ("2", "Type A"),
    }
