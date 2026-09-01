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
