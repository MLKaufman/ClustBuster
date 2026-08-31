from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pytest
from scipy import sparse

import clustbuster.io.sce as sce_module
from clustbuster.core.workspace import workspace_from_import
from clustbuster.io.sce import SceImporter, SceImportError
from clustbuster.models import ObjectFormat
from clustbuster.services.exports import ANNOTATION_COLUMN, PROVENANCE_KEY, WorkspaceExportService
from clustbuster.services.workspaces import configure_workspace

FIXTURE = Path(__file__).parents[2] / "testdata" / "sce.rds"


def test_load_genuine_single_cell_experiment_fixture() -> None:
    result = SceImporter().load(FIXTURE)

    assert result.report.source_format is ObjectFormat.SINGLE_CELL_EXPERIMENT
    assert result.adata.shape == (40, 8)
    assert list(result.adata.obs_names) == [f"demo-cell-{index:03d}" for index in range(1, 41)]
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
    assert result.adata.obs["sce_clusters"].tolist() == [
        str(cluster) for cluster in range(4) for _ in range(10)
    ]
    assert result.report.candidate_cluster_columns[0] == "sce_clusters"
    assert result.report.embeddings == ("X_umap",)
    assert set(result.adata.layers) == {"counts", "logcounts"}
    assert sparse.issparse(result.adata.layers["counts"])
    assert not sparse.issparse(result.adata.layers["logcounts"])
    assert np.allclose(result.adata.X, result.adata.layers["logcounts"])
    assert result.report.mapping_decisions["X"] == "SCE assay 'logcounts' -> adata.X"
    assert result.report.unsupported_components == ["Alternative experiments are not imported"]


def test_probe_recognizes_single_cell_experiment_fixture() -> None:
    result = SceImporter().probe(FIXTURE)
    assert result.format is ObjectFormat.SINGLE_CELL_EXPERIMENT
    assert result.confidence == 1.0


def test_non_sce_rds_has_actionable_error(tmp_path: Path) -> None:
    path = tmp_path / "not-sce.rds"
    path.write_bytes(b"not an R object")
    with pytest.raises(SceImportError, match="supported in-memory SingleCellExperiment"):
        SceImporter().load(path)


def test_delayed_or_custom_assay_fails_safely(monkeypatch: pytest.MonkeyPatch) -> None:
    source = sce_module._decode_source(FIXTURE)
    source.assays.data.listData["counts"] = SimpleNamespace(shape=(8, 40))
    monkeypatch.setattr(sce_module, "_decode_source", lambda path: source)

    with pytest.raises(SceImportError, match="unsupported delayed or custom representation"):
        SceImporter().load(FIXTURE)


def test_sce_annotation_and_h5ad_export_round_trip(tmp_path: Path) -> None:
    source_bytes = FIXTURE.read_bytes()
    workspace = workspace_from_import(SceImporter().load(FIXTURE))
    configure_workspace(
        workspace,
        cluster_column="sce_clusters",
        embedding_key="X_umap",
        expression_source=workspace.expression_source,
    )
    workspace.annotations.assign("0", "T cell")

    artifact = WorkspaceExportService().export_h5ad(workspace, tmp_path).artifacts[0]

    exported = ad.read_h5ad(artifact.path)
    assert exported.obs[ANNOTATION_COLUMN].iloc[:10].tolist() == ["T cell"] * 10
    assert exported.uns[PROVENANCE_KEY]["source_format"] == "single_cell_experiment"
    assert list(exported.var_names) == list(workspace.adata.var_names)
    assert set(exported.layers) == {"counts", "logcounts"}
    assert FIXTURE.read_bytes() == source_bytes
