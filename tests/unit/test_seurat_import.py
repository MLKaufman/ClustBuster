from pathlib import Path

import anndata as ad
import pytest
from scipy import sparse

from clustbuster.core.workspace import workspace_from_import
from clustbuster.io.seurat import SeuratImporter, SeuratImportError
from clustbuster.models import ObjectFormat
from clustbuster.services.exports import ANNOTATION_COLUMN, PROVENANCE_KEY, WorkspaceExportService
from clustbuster.services.workspaces import configure_workspace

FIXTURE = Path(__file__).parents[2] / "testdata" / "so.rds"
H5SEURAT_FIXTURE = Path(__file__).parents[2] / "testdata" / "so.h5seurat"


def test_load_genuine_seurat_v5_fixture() -> None:
    result = SeuratImporter().load(FIXTURE)

    assert result.report.source_format is ObjectFormat.SEURAT
    assert result.adata.shape == (40, 8)
    assert list(result.adata.var_names) == [
        "CD3D",
        "IL7R",
        "LYZ",
        "S100A8",
        "MS4A1",
        "CD79A",
        "NKG7",
        "GNLY",
    ]
    assert result.report.candidate_cluster_columns[0] == "seurat_clusters"
    assert result.report.embeddings == ("X_umap",)
    assert [source.label for source in result.report.expression_sources] == [
        "X",
        "layer:assay:ALT:counts",
        "layer:counts",
    ]
    assert sparse.issparse(result.adata.X)
    assert (result.adata.layers["assay:ALT:counts"] != result.adata.layers["counts"] * 2).nnz == 0
    assert any("Recovered Seurat v5 feature identifiers" in item for item in result.report.warnings)
    assert any("Mapped 1 compatible secondary-assay" in item for item in result.report.warnings)
    assert result.report.unsupported_components == [
        "Non-active assay 'ADT' has a different cell or feature space"
    ]


def test_probe_recognizes_compressed_rds_fixture() -> None:
    result = SeuratImporter().probe(FIXTURE)
    assert result.format is ObjectFormat.SEURAT
    assert result.confidence == 0.8


def test_load_h5seurat_fixture() -> None:
    result = SeuratImporter().load(H5SEURAT_FIXTURE)

    assert result.report.source_format is ObjectFormat.SEURAT
    assert result.adata.shape == (12, 5)
    assert list(result.adata.obs_names) == [f"h5-cell-{index:03d}" for index in range(12)]
    assert list(result.adata.var_names) == ["CD3D", "IL7R", "LYZ", "MS4A1", "NKG7"]
    assert result.adata.obs["seurat_clusters"].tolist() == ["0"] * 4 + ["1"] * 4 + ["2"] * 4
    assert result.report.candidate_cluster_columns[0] == "seurat_clusters"
    assert result.report.embeddings == ("X_umap",)
    assert list(result.adata.layers) == ["counts"]
    assert sparse.issparse(result.adata.X)
    assert result.report.mapping_decisions["X"] == "Seurat RNA primary matrix -> adata.X"
    assert result.report.unsupported_components == []


def test_probe_recognizes_h5seurat_fixture() -> None:
    result = SeuratImporter().probe(H5SEURAT_FIXTURE)
    assert result.format is ObjectFormat.SEURAT
    assert result.confidence == 1.0


def test_non_seurat_rds_has_actionable_error(tmp_path: Path) -> None:
    path = tmp_path / "not-seurat.rds"
    path.write_bytes(b"not an R object")
    with pytest.raises(SeuratImportError, match="supported in-memory Seurat v5"):
        SeuratImporter().load(path)


def test_seurat_annotation_and_h5ad_export_round_trip(tmp_path: Path) -> None:
    source_bytes = FIXTURE.read_bytes()
    workspace = workspace_from_import(SeuratImporter().load(FIXTURE))
    configure_workspace(
        workspace,
        cluster_column="seurat_clusters",
        embedding_key="X_umap",
        expression_source=workspace.expression_source,
    )
    workspace.annotations.assign("0", "T cell")

    artifact = WorkspaceExportService().export_h5ad(workspace, tmp_path).artifacts[0]

    exported = ad.read_h5ad(artifact.path)
    assert exported.obs[ANNOTATION_COLUMN].iloc[:10].tolist() == ["T cell"] * 10
    assert exported.uns[PROVENANCE_KEY]["source_format"] == "seurat"
    assert list(exported.var_names) == list(workspace.adata.var_names)
    assert set(exported.layers) == {"counts", "assay:ALT:counts"}
    assert FIXTURE.read_bytes() == source_bytes
