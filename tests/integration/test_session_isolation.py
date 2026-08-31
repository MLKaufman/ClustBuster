from __future__ import annotations

import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from clustbuster.core.workspace import workspace_from_import
from clustbuster.models import ExpressionSource
from clustbuster.services.exports import WorkspaceExportService
from clustbuster.services.imports import ImportService, SessionFiles
from clustbuster.services.workspaces import configure_workspace


def _write_fixture(path: Path) -> None:
    adata = ad.AnnData(
        X=np.array([[1.0, 0.0], [2.0, 1.0], [0.0, 3.0]]),
        obs=pd.DataFrame(
            {"leiden": pd.Categorical(["0", "0", "1"])},
            index=["cell-1", "cell-2", "cell-3"],
        ),
        var=pd.DataFrame(index=["CD3D", "LYZ"]),
    )
    adata.obsm["X_umap"] = np.array([[0.0, 0.0], [1.0, 0.0], [3.0, 2.0]])
    adata.write_h5ad(path)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_concurrent_sessions_keep_files_state_and_exports_independent(tmp_path: Path) -> None:
    source = tmp_path / "source.h5ad"
    _write_fixture(source)
    source_checksum = _sha256(source)
    sessions_root = tmp_path / "sessions"
    session_a = SessionFiles.create(sessions_root)
    session_b = SessionFiles.create(sessions_root)

    assert session_a.session_id != session_b.session_id
    assert session_a.root.parent == session_b.root.parent == sessions_root.resolve()
    assert session_a.root != session_b.root

    importer = ImportService(max_upload_mb=10)
    upload = {"name": "source.h5ad", "datapath": str(source)}
    workspace_a = workspace_from_import(importer.import_upload(upload, session_a))
    workspace_b = workspace_from_import(importer.import_upload(upload, session_b))
    for workspace in (workspace_a, workspace_b):
        configure_workspace(
            workspace,
            cluster_column="leiden",
            embedding_key="X_umap",
            expression_source=ExpressionSource.x(),
        )

    workspace_a.annotations.assign("0", "T cell")
    assert workspace_a.annotations.annotation_for("0") == "T cell"
    assert workspace_b.annotations.annotation_for("0") == "0"

    export = WorkspaceExportService().export_h5ad(workspace_a, session_a.exports)
    assert export.artifacts[0].path.is_relative_to(session_a.root)
    assert not export.artifacts[0].path.is_relative_to(session_b.root)
    assert not any(session_b.exports.iterdir())
    assert _sha256(source) == source_checksum

    session_a.cleanup()
    assert not session_a.root.exists()
    assert session_b.root.exists()
    assert len(list(session_b.uploads.glob("*.h5ad"))) == 1
    session_b.cleanup()


def test_cleanup_refuses_a_tampered_session_boundary(tmp_path: Path) -> None:
    session = SessionFiles.create(tmp_path / "sessions")
    tampered = SessionFiles(
        workspace_root=session.workspace_root,
        session_id=session.session_id,
        root=tmp_path,
        uploads=session.uploads,
        cache=session.cache,
        state=session.state,
        exports=session.exports,
    )

    try:
        tampered.cleanup()
    except RuntimeError as exc:
        assert "outside its session boundary" in str(exc)
    else:
        raise AssertionError("Tampered cleanup should have been rejected")
    finally:
        session.cleanup()
