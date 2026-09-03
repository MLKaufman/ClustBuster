from __future__ import annotations

import io
import zipfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from clustbuster.core.annotations import AnnotationStore
from clustbuster.io.export import ExportCollisionError
from clustbuster.models import ExpressionSource, ImportReport, ObjectFormat, Workspace
from clustbuster.services.exports import (
    ANNOTATION_COLUMN,
    CLUSTER_COLUMN,
    PROVENANCE_KEY,
    WorkspaceExportService,
)


def _workspace(*, mixed_clusters: bool = False) -> Workspace:
    clusters: list[object] = [1, "1"] if mixed_clusters else ["a", "a", "b"]
    cell_count = len(clusters)
    matrix = sparse.csr_matrix(np.arange(cell_count * 2).reshape(cell_count, 2))
    adata = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(
            {"source_cluster": clusters},
            index=[f"cell-{index}" for index in range(cell_count)],
        ),
        var=pd.DataFrame(index=["gene-1", "gene-2"]),
    )
    adata.layers["counts"] = matrix.copy()
    adata.obsm["X_umap"] = np.arange(cell_count * 2).reshape(cell_count, 2)
    report = ImportReport(
        source_format=ObjectFormat.H5AD,
        source_filename="study, one.h5ad",
        cell_count=cell_count,
        feature_count=2,
        sparse=True,
        candidate_cluster_columns=("source_cluster",),
        embeddings=("X_umap",),
        expression_sources=(ExpressionSource.x(), ExpressionSource.named_layer("counts")),
    )
    return Workspace(
        adata=adata,
        source_format=ObjectFormat.H5AD,
        source_filename="study, one.h5ad",
        expression_source=ExpressionSource.x(),
        annotations=AnnotationStore.from_clusters(clusters),
        import_report=report,
        cluster_column="source_cluster",
        embedding_key="X_umap",
    )


def test_annotation_zip_round_trips_unicode_quotes_and_typed_ids(tmp_path: Path) -> None:
    workspace = _workspace(mixed_clusters=True)
    workspace.annotations.assign(1, 'Myeloid, "reviewed"', notes="naïve ✓")
    workspace.annotations.assign("1", "T cell")
    result = WorkspaceExportService().export_annotation_zip(workspace, tmp_path)
    artifact = result.artifacts[0]
    assert artifact.path.is_file()
    assert len(artifact.checksum) == 64

    with zipfile.ZipFile(artifact.path) as archive:
        names = archive.namelist()
        assert len(names) == 2
        cluster_csv = archive.read(next(name for name in names if "cluster" in name)).decode()
        cell_csv = archive.read(next(name for name in names if "cell" in name)).decode()
    clusters = pd.read_csv(io.StringIO(cluster_csv))
    cells = pd.read_csv(io.StringIO(cell_csv))
    assert clusters["annotation"].tolist() == ['Myeloid, "reviewed"', "T cell"]
    assert clusters["notes"].iloc[0] == "naïve ✓"
    assert len(set(cells["cluster_id"])) == 2


def test_cluster_annotation_csv_exports_current_table(tmp_path: Path) -> None:
    workspace = _workspace()
    workspace.annotations.assign("a", "T cell", notes="reviewed")

    result = WorkspaceExportService().export_cluster_annotations_csv(workspace, tmp_path)
    artifact = result.artifacts[0]
    exported = pd.read_csv(artifact.path)

    assert artifact.kind == "cluster_annotation_csv"
    assert artifact.media_type == "text/csv"
    assert exported["cluster_display"].tolist() == ["a", "b"]
    assert exported["annotation"].tolist() == ["T cell", "b"]
    assert exported["notes"].iloc[0] == "reviewed"


def test_reference_matrix_csv_groups_by_current_annotations(tmp_path: Path) -> None:
    workspace = _workspace()
    workspace.annotations.assign("a", "T cell")

    result = WorkspaceExportService().export_reference_matrix_csv(workspace, tmp_path)
    artifact = result.artifacts[0]
    exported = pd.read_csv(artifact.path)

    assert artifact.kind == "reference_matrix_csv"
    assert artifact.media_type == "text/csv"
    assert exported.columns.tolist() == ["gene", "T cell", "b"]
    assert exported["gene"].tolist() == ["gene-1", "gene-2"]
    np.testing.assert_allclose(exported["T cell"], [1.0, 2.0])
    np.testing.assert_allclose(exported["b"], [4.0, 5.0])


def test_reference_matrix_csv_groups_by_metadata_and_omits_missing_cells(
    tmp_path: Path,
) -> None:
    workspace = _workspace()
    workspace.adata.obs["sample"] = ["control", "treated", None]

    result = WorkspaceExportService().export_reference_matrix_csv(
        workspace, tmp_path, metadata_column="sample"
    )
    exported = pd.read_csv(result.artifacts[0].path)

    assert exported.columns.tolist() == ["gene", "control", "treated"]
    np.testing.assert_allclose(exported["control"], [0.0, 1.0])
    np.testing.assert_allclose(exported["treated"], [2.0, 3.0])


def test_reference_matrix_csv_keeps_gene_column_name_unique(tmp_path: Path) -> None:
    workspace = _workspace()
    workspace.adata.obs["label"] = ["gene", "gene [group 2]", "other"]

    result = WorkspaceExportService().export_reference_matrix_csv(
        workspace, tmp_path, metadata_column="label"
    )
    exported = pd.read_csv(result.artifacts[0].path)

    assert exported.columns.tolist() == ["gene", "gene [group 1]", "gene [group 2]", "other"]


def test_h5ad_export_preserves_source_and_validates_round_trip(tmp_path: Path) -> None:
    workspace = _workspace()
    source_path = tmp_path / "source.h5ad"
    workspace.adata.write_h5ad(source_path)
    source_bytes = source_path.read_bytes()
    workspace.annotations.assign("a", "Astrocyte")

    result = WorkspaceExportService().export_h5ad(workspace, tmp_path / "exports")
    exported = ad.read_h5ad(result.artifacts[0].path)
    assert exported.obs[ANNOTATION_COLUMN].tolist() == ["Astrocyte", "Astrocyte", "b"]
    assert CLUSTER_COLUMN in exported.obs
    assert PROVENANCE_KEY in exported.uns
    assert set(exported.layers) == {"counts"}
    assert set(exported.obsm) == {"X_umap"}
    assert source_path.read_bytes() == source_bytes
    assert CLUSTER_COLUMN not in workspace.adata.obs
    assert PROVENANCE_KEY not in workspace.adata.uns


@pytest.mark.parametrize("collision", [CLUSTER_COLUMN, ANNOTATION_COLUMN, PROVENANCE_KEY])
def test_h5ad_export_refuses_namespace_collisions(tmp_path: Path, collision: str) -> None:
    workspace = _workspace()
    if collision == PROVENANCE_KEY:
        workspace.adata.uns[collision] = {"existing": True}
    else:
        workspace.adata.obs[collision] = "existing"
    with pytest.raises(ExportCollisionError, match="would be overwritten"):
        WorkspaceExportService().export_h5ad(workspace, tmp_path)
