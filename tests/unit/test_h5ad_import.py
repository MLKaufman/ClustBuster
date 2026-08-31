from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from clustbuster.io.h5ad import H5adImporter, H5adImportError
from clustbuster.models import ObjectFormat


def _fixture_adata(matrix: object) -> ad.AnnData:
    adata = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(
            {"leiden": pd.Categorical(["0", "0", "1"]), "sample": ["a", "b", "c"]},
            index=["cell-1", "cell-2", "cell-3"],
        ),
        var=pd.DataFrame(index=["CD3D", "LYZ"]),
    )
    adata.layers["counts"] = matrix.copy()  # type: ignore[attr-defined]
    adata.obsm["X_pca"] = np.ones((3, 4))
    adata.obsm["X_umap"] = np.ones((3, 2))
    return adata


@pytest.mark.parametrize(
    "matrix",
    [
        np.array([[1, 0], [2, 1], [0, 3]], dtype=float),
        sparse.csr_matrix([[1, 0], [2, 1], [0, 3]]),
    ],
)
def test_load_reports_semantic_mapping(tmp_path: Path, matrix: object) -> None:
    path = tmp_path / "fixture.h5ad"
    _fixture_adata(matrix).write_h5ad(path)
    result = H5adImporter().load(path)
    report = result.report
    assert report.source_format is ObjectFormat.H5AD
    assert (report.cell_count, report.feature_count) == (3, 2)
    assert report.embeddings[0] == "X_umap"
    assert "leiden" in report.candidate_cluster_columns
    assert [source.label for source in report.expression_sources] == ["X", "layer:counts"]


def test_probe_uses_signature_not_only_extension(tmp_path: Path) -> None:
    path = tmp_path / "fixture.data"
    _fixture_adata(np.ones((3, 2))).write_h5ad(path)
    result = H5adImporter().probe(path)
    assert result.format is ObjectFormat.H5AD
    assert result.confidence == 0.6


def test_invalid_file_has_actionable_error(tmp_path: Path) -> None:
    path = tmp_path / "bad.h5ad"
    path.write_text("not hdf5")
    with pytest.raises(H5adImportError, match="truncated, corrupt"):
        H5adImporter().load(path)
