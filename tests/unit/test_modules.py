import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from clustbuster.core.modules import calculate_module_score
from clustbuster.models import ExpressionSource
from clustbuster.plotting.module import module_score_figure, module_score_violin_figure


def _adata() -> ad.AnnData:
    return ad.AnnData(
        X=sparse.csr_matrix([[1, 0, 2], [3, 0, 2], [0, 4, 2], [0, 6, 2]]),
        obs=pd.DataFrame({"cluster": ["a", "a", "b", "b"]}, index=list("wxyz")),
        var=pd.DataFrame(index=["CD3D", "LYZ", "CONSTANT"]),
        obsm={"X_umap": np.array([[0, 0], [1, 0], [2, 2], [3, 2]])},
    )


def test_module_score_is_standardized_sparse_safe_and_clustered() -> None:
    result = calculate_module_score(
        _adata(),
        ExpressionSource.x(),
        "cluster",
        ("CD3D", "CONSTANT", "missing"),
        name="T cell",
    )
    assert result.name == "T cell"
    assert result.report.matched == ("CD3D", "CONSTANT")
    assert result.report.missing == ("missing",)
    assert result.values["module_score"].mean() == pytest.approx(0.0)
    assert np.isfinite(result.values["module_score"]).all()
    assert result.values["cluster"].tolist() == ["a", "a", "b", "b"]
    assert result.cluster_summary["cluster"].tolist() == ["a", "b"]
    assert result.cluster_summary["cells"].tolist() == [2, 2]
    assert result.cluster_summary.loc[0, "mean_score"] > result.cluster_summary.loc[1, "mean_score"]


def test_module_score_figure_uses_webgl_and_cell_ids() -> None:
    adata = _adata()
    result = calculate_module_score(adata, ExpressionSource.x(), "cluster", ("CD3D", "LYZ"))
    figure = module_score_figure(
        np.asarray(adata.obsm["X_umap"]), adata.obs_names.astype(str).tolist(), result
    )
    assert figure.data[0].type == "scattergl"
    assert figure.data[0].customdata[0][0] == "w"


def test_module_score_violin_shows_every_cell_as_a_point() -> None:
    result = calculate_module_score(_adata(), ExpressionSource.x(), "cluster", ("CD3D", "LYZ"))
    figure = module_score_violin_figure(result)

    assert [trace.type for trace in figure.data] == ["violin", "violin"]
    assert all(trace.points == "all" for trace in figure.data)
    assert sum(len(trace.y) for trace in figure.data) == 4
    assert figure.layout.yaxis.title.text == "Module score"
