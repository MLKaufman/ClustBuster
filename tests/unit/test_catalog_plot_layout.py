import anndata as ad
import numpy as np
import pandas as pd
import pytest

from clustbuster.app import _catalog_plot_output
from clustbuster.core.expression import aggregate_dotplot
from clustbuster.core.markers import MarkerResult
from clustbuster.models import ExpressionSource
from clustbuster.plotting.dotplot import dotplot_figure
from clustbuster.plotting.heatmap import marker_heatmap_figure


@pytest.mark.parametrize("cluster_count", [5, 35, 80])
@pytest.mark.parametrize("kind", ["dot", "heatmap"])
def test_catalog_container_fits_cluster_dependent_figure(cluster_count: int, kind: str) -> None:
    adata = ad.AnnData(
        X=np.arange(cluster_count * 3).reshape(cluster_count, 3).astype(float),
        obs=pd.DataFrame({"cluster": [str(i) for i in range(cluster_count)]},
                         index=[f"cell-{i}" for i in range(cluster_count)]),
        var=pd.DataFrame(index=["CD3D", "IL7R", "NKG7"]),
    )
    result = aggregate_dotplot(adata, ExpressionSource.x(), "cluster", tuple(adata.var_names))
    if kind == "dot":
        figure = dotplot_figure(result)
    else:
        means = result.values.pivot(index="cluster_id", columns="gene", values="mean_expression")
        figure = marker_heatmap_figure(MarkerResult("T cell", "catalog", pd.DataFrame(), means))
    html = str(_catalog_plot_output(0, figure, kind))
    assert f"height:{figure.layout.height}px" in html
    assert f"min-height:{figure.layout.height}px" in html
    assert figure.layout.xaxis3.showticklabels is True
    assert set(figure.layout.xaxis3.ticktext) == set(adata.var_names)
    if cluster_count >= 35:
        assert figure.layout.height > 700
