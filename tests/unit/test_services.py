from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from clustbuster.core.workspace import workspace_from_import
from clustbuster.models import ExpressionSource
from clustbuster.services.imports import ImportService, SessionFiles, UploadError
from clustbuster.services.workspaces import configure_workspace


def _write_fixture(path: Path) -> None:
    adata = ad.AnnData(
        X=np.array([[1, 0], [2, 1], [0, 3]], dtype=float),
        obs=pd.DataFrame(
            {"leiden": pd.Categorical(["0", "0", "1"])},
            index=["cell-1", "cell-2", "cell-3"],
        ),
        var=pd.DataFrame(index=["CD3D", "LYZ"]),
    )
    adata.layers["counts"] = adata.X.copy()
    adata.obsm["X_umap"] = np.array([[0, 0], [1, 0], [3, 2]], dtype=float)
    adata.write_h5ad(path)


def test_import_service_copies_upload_into_isolated_session(tmp_path: Path) -> None:
    upload_path = tmp_path / "shiny-upload"
    _write_fixture(upload_path)
    sessions_root = tmp_path / "sessions"
    session_files = SessionFiles.create(sessions_root)
    result = ImportService(max_upload_mb=10).import_upload(
        {"name": "study.h5ad", "datapath": str(upload_path)}, session_files
    )
    assert result.report.source_filename == "study.h5ad"
    assert len(list(session_files.uploads.glob("*.h5ad"))) == 1
    assert upload_path.is_file()
    session_files.cleanup()
    assert not session_files.root.exists()


def test_import_service_rejects_empty_upload(tmp_path: Path) -> None:
    upload_path = tmp_path / "empty"
    upload_path.touch()
    session_files = SessionFiles.create(tmp_path / "sessions")
    with pytest.raises(UploadError, match="empty"):
        ImportService(max_upload_mb=10).import_upload(
            {"name": "empty.h5ad", "datapath": str(upload_path)}, session_files
        )


def test_workspace_configuration_initializes_annotation_state(tmp_path: Path) -> None:
    upload_path = tmp_path / "fixture.h5ad"
    _write_fixture(upload_path)
    session_files = SessionFiles.create(tmp_path / "sessions")
    imported = ImportService(max_upload_mb=10).import_upload(
        {"name": "fixture.h5ad", "datapath": str(upload_path)}, session_files
    )
    workspace = workspace_from_import(imported)
    configure_workspace(
        workspace,
        cluster_column="leiden",
        embedding_key="X_umap",
        expression_source=ExpressionSource.named_layer("counts"),
    )
    assert workspace.expression_source.label == "layer:counts"
    assert [record.annotation for record in workspace.annotations.records()] == ["0", "1"]

