import anndata as ad
import numpy as np
import pandas as pd

from clustbuster.core.annotations import AnnotationStore
from clustbuster.models import (
    ExpressionSource,
    ImportReport,
    ObjectFormat,
    Workspace,
)
from clustbuster.plotting.embedding import embedding_figure


def _workspace() -> Workspace:
    adata = ad.AnnData(
        X=np.ones((3, 2)),
        obs=pd.DataFrame({"cluster": ["a", "a", "b"]}, index=["c1", "c2", "c3"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    adata.obsm["X_umap"] = np.array([[0, 0], [1, 0], [2, 3]])
    report = ImportReport(
        source_format=ObjectFormat.H5AD,
        source_filename="fixture.h5ad",
        cell_count=3,
        feature_count=2,
        sparse=False,
        candidate_cluster_columns=("cluster",),
        embeddings=("X_umap",),
        expression_sources=(ExpressionSource.x(),),
    )
    return Workspace(
        adata=adata,
        source_format=ObjectFormat.H5AD,
        source_filename="fixture.h5ad",
        expression_source=ExpressionSource.x(),
        annotations=AnnotationStore.from_clusters(adata.obs["cluster"].tolist()),
        import_report=report,
        cluster_column="cluster",
        embedding_key="X_umap",
    )


def test_embedding_uses_webgl_and_stable_cell_customdata() -> None:
    figure = embedding_figure(_workspace())
    assert len(figure.data) == 2
    assert all(trace.type == "scattergl" for trace in figure.data)
    assert figure.data[0].customdata[0][0] == "c1"


def test_embedding_can_color_by_current_annotations() -> None:
    workspace = _workspace()
    workspace.annotations.assign("a", "T cell")
    figure = embedding_figure(workspace, color_by="annotation")
    assert {trace.name for trace in figure.data} == {"T cell", "b"}
