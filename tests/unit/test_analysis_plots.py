import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from clustbuster.core.expression import aggregate_dotplot, extract_expression
from clustbuster.models import ExpressionSource
from clustbuster.plotting.dotplot import dotplot_figure
from clustbuster.plotting.feature import feature_gene_figure


def _adata() -> ad.AnnData:
    matrix = sparse.csr_matrix([[2, 0, 1], [4, 0, 0], [0, 3, 1], [0, 5, 0]])
    return ad.AnnData(
        X=matrix,
        obs=pd.DataFrame({"cluster": ["a", "a", "b", "b"]}, index=list("wxyz")),
        var=pd.DataFrame(index=["CD3D", "LYZ", "MS4A1"]),
    )


def test_feature_gene_figure_is_an_independent_webgl_plot() -> None:
    adata = _adata()
    result = extract_expression(adata, ExpressionSource.x(), ("CD3D", "LYZ"))
    coordinates = np.array([[0, 0], [1, 0], [2, 2], [3, 2]])
    figure = feature_gene_figure(
        coordinates,
        adata.obs_names.tolist(),
        "CD3D",
        result.values["CD3D"].to_numpy(),
    )
    assert len(figure.data) == 1
    assert figure.data[0].type == "scattergl"
    assert figure.data[0].name == "CD3D"
    assert figure.layout.height == 520

    labeled = feature_gene_figure(
        coordinates,
        adata.obs_names.tolist(),
        "CD3D",
        result.values["CD3D"].to_numpy(),
        ((0.5, 0.0, "T cell"), (2.5, 2.0, "Myeloid")),
    )
    assert [annotation.text for annotation in labeled.layout.annotations] == [
        "T cell",
        "Myeloid",
    ]


def test_dotplot_figure_maps_fraction_to_marker_size() -> None:
    result = aggregate_dotplot(_adata(), ExpressionSource.x(), "cluster", ("CD3D", "LYZ"))
    figure = dotplot_figure(result)
    assert figure.data[0].type == "scatter"
    assert len(figure.data[0].marker.size) == 4
    assert [trace.mode for trace in figure.data[1:]] == ["lines"]
    assert figure.layout.margin.l == 120
    assert figure.layout.margin.b == 120
    assert figure.layout.height <= 1060
    assert figure.layout.yaxis2.showticklabels is False
    assert figure.layout.yaxis3.showticklabels is True
    assert figure.layout.yaxis3.side == "left"
