import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from clustbuster.core.expression import aggregate_dotplot, extract_expression
from clustbuster.models import ExpressionSource
from clustbuster.plotting.dotplot import dotplot_figure
from clustbuster.plotting.feature import feature_figure


def _adata() -> ad.AnnData:
    matrix = sparse.csr_matrix([[2, 0, 1], [4, 0, 0], [0, 3, 1], [0, 5, 0]])
    return ad.AnnData(
        X=matrix,
        obs=pd.DataFrame({"cluster": ["a", "a", "b", "b"]}, index=list("wxyz")),
        var=pd.DataFrame(index=["CD3D", "LYZ", "MS4A1"]),
    )


def test_feature_figure_uses_one_webgl_trace_per_gene() -> None:
    adata = _adata()
    result = extract_expression(adata, ExpressionSource.x(), ("CD3D", "LYZ"))
    coordinates = np.array([[0, 0], [1, 0], [2, 2], [3, 2]])
    figure = feature_figure(coordinates, adata.obs_names.tolist(), result)
    assert len(figure.data) == 2
    assert all(trace.type == "scattergl" for trace in figure.data)
    assert [trace.xaxis for trace in figure.data] == ["x", "x2"]
    assert figure.layout.height == 1000


def test_dotplot_figure_maps_fraction_to_marker_size() -> None:
    result = aggregate_dotplot(
        _adata(), ExpressionSource.x(), "cluster", ("CD3D", "LYZ")
    )
    figure = dotplot_figure(result)
    assert figure.data[0].type == "scatter"
    assert len(figure.data[0].marker.size) == 4
    assert [trace.mode for trace in figure.data[1:]] == ["lines", "lines"]
    assert figure.layout.margin.l == 120
    assert figure.layout.margin.b == 120
    assert figure.layout.height <= 1060
