from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from clustbuster.core.workspace import workspace_from_import
from clustbuster.integrations.pyclustifyr import (
    PYCLUSTIFYR_REVISION,
    PyClustifyrAdapter,
    ReferenceAnnotationError,
    ReferenceAnnotationParameters,
)
from clustbuster.io.h5ad import H5adImporter
from clustbuster.models import ExpressionSource
from clustbuster.resources.references.local import LocalReferenceProvider

ROOT = Path(__file__).resolve().parents[2]


def _workspace():
    result = H5adImporter().load(ROOT / "testdata" / "clustbuster-demo.h5ad")
    workspace = workspace_from_import(result)
    workspace.cluster_column = "leiden"
    workspace.expression_source = ExpressionSource.x()
    return workspace


def test_pyclustifyr_contract_recovers_demo_clusters_without_mutating_obs() -> None:
    workspace = _workspace()
    reference = LocalReferenceProvider(
        ROOT / "resources" / "reference_matrices"
    ).load_reference("pbmc-demo-v1")
    before = workspace.adata.obs.copy(deep=True)

    result = PyClustifyrAdapter().annotate(
        workspace,
        reference,
        ReferenceAnnotationParameters(),
    )

    assert [item.annotation for item in result.predictions] == [
        "T cell",
        "Myeloid cell",
        "B cell",
        "NK cell",
        "Platelet",
    ]
    assert result.correlations.shape == (5, 5)
    assert len(result.matched_genes) == 10
    assert result.package_version == "1.0.0"
    assert result.package_revision == PYCLUSTIFYR_REVISION
    pd.testing.assert_frame_equal(workspace.adata.obs, before)


def test_reference_annotation_enforces_minimum_overlap() -> None:
    workspace = _workspace()
    reference = LocalReferenceProvider(
        ROOT / "resources" / "reference_matrices"
    ).load_reference("pbmc-demo-v1")
    with pytest.raises(ReferenceAnnotationError, match="Only 10 reference genes"):
        PyClustifyrAdapter().annotate(
            workspace,
            reference,
            ReferenceAnnotationParameters(minimum_gene_overlap=11),
        )
