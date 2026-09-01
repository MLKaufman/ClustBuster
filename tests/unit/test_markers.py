import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from clustbuster.core.annotations import ClusterIdentifier
from clustbuster.core.markers import MarkerAnalysisError, rank_markers
from clustbuster.models import ExpressionSource
from clustbuster.plotting.heatmap import marker_heatmap_figure


def _adata() -> ad.AnnData:
    matrix = sparse.csr_matrix(
        [
            [8, 0, 1],
            [9, 0, 1],
            [7, 1, 1],
            [0, 7, 1],
            [1, 8, 1],
            [0, 9, 1],
        ]
    )
    return ad.AnnData(
        X=matrix,
        obs=pd.DataFrame({"cluster": ["a", "a", "a", "b", "b", "b"]}, index=list("uvwxyz")),
        var=pd.DataFrame(index=["CD3D", "LYZ", "CONSTANT"]),
    )


def test_sparse_marker_ranking_finds_positive_cluster_marker() -> None:
    cluster = ClusterIdentifier.from_value("a").serialized
    result = rank_markers(
        _adata(), ExpressionSource.x(), "cluster", cluster, top_n=3, min_fraction=0
    )
    assert result.selected_cluster == "a"
    assert result.method == "Welch t-test"
    assert result.values.iloc[0]["gene"] == "CD3D"
    assert "LYZ" not in result.values["gene"].tolist()
    assert result.values.iloc[0]["p_adjusted"] < 0.05
    assert result.heatmap.index.tolist() == ["a", "b"]


def test_marker_ranking_validates_comparison_size_and_threshold() -> None:
    adata = _adata()
    singleton = ClusterIdentifier.from_value("singleton").serialized
    with pytest.raises(MarkerAnalysisError, match="at least two cells"):
        rank_markers(adata, ExpressionSource.x(), "cluster", singleton)
    cluster = ClusterIdentifier.from_value("a").serialized
    with pytest.raises(MarkerAnalysisError, match="between 0 and 1"):
        rank_markers(
            adata, ExpressionSource.x(), "cluster", cluster, min_fraction=1.1
        )


def test_marker_heatmap_is_standardized_and_labeled() -> None:
    cluster = ClusterIdentifier.from_value("a").serialized
    result = rank_markers(
        _adata(), ExpressionSource.x(), "cluster", cluster, top_n=2, min_fraction=0
    )
    figure = marker_heatmap_figure(result)
    assert figure.data[0].type == "heatmap"
    assert any(trace.mode == "lines" for trace in figure.data[1:])
    assert figure.layout.title.text == "Top markers for cluster a"
    assert np.isfinite(np.asarray(figure.data[0].z)).all()
