from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from clustbuster.models import ReferenceFilters
from clustbuster.resources.errors import ProviderSchemaError, ResourceNotFoundError
from clustbuster.resources.markers.csv import CsvMarkerProvider
from clustbuster.resources.references.local import LocalReferenceProvider


def _marker_catalog(path: Path) -> None:
    pd.DataFrame(
        [
            {
                "cell_type": "T cell",
                "gene": "CD3D",
                "species": "human",
                "tissue": "blood",
                "direction": "positive",
                "evidence": "canonical",
                "citation": "fixture",
                "confidence": 0.95,
                "resource_version": "1",
            },
            {
                "cell_type": "T cell",
                "gene": "IL7R",
                "species": "human",
                "tissue": "blood",
                "direction": "positive",
                "evidence": "canonical",
                "citation": "fixture",
                "confidence": 0.9,
                "resource_version": "1",
            },
            {
                "cell_type": "B cell",
                "gene": "MS4A1",
                "species": "human",
                "tissue": "blood",
                "direction": "positive",
                "evidence": "canonical",
                "citation": "fixture",
                "confidence": 0.97,
                "resource_version": "1",
            },
        ]
    ).to_csv(path, index=False)


def test_csv_marker_provider_searches_filters_and_loads_sets(tmp_path: Path) -> None:
    path = tmp_path / "markers.csv"
    _marker_catalog(path)
    provider = CsvMarkerProvider(path)
    assert provider.status().available
    summaries = provider.search_cell_types("cell", species="HUMAN", tissue="blood")
    assert [item.cell_type for item in summaries] == ["B cell", "T cell"]
    assert [item.marker_count for item in summaries] == [1, 2]
    marker_set = provider.get_markers("t CELL", species="human")
    assert marker_set.cell_type == "T cell"
    assert [record.gene for record in marker_set.records] == ["CD3D", "IL7R"]
    with pytest.raises(ResourceNotFoundError, match="No marker set"):
        provider.get_markers("NK cell")


def test_csv_marker_provider_reports_invalid_schema(tmp_path: Path) -> None:
    path = tmp_path / "markers.csv"
    path.write_text("cell_type,gene\nT cell,CD3D\n")
    provider = CsvMarkerProvider(path)
    assert not provider.status().available
    with pytest.raises(ProviderSchemaError, match="missing required columns"):
        provider.search_cell_types("")


def _reference_fixture(root: Path, *, bad_checksum: bool = False) -> None:
    matrix = root / "reference.csv"
    matrix.write_text("gene,T cell,B cell\nCD3D,4.0,0.1\nMS4A1,0.0,5.0\n")
    checksum = hashlib.sha256(matrix.read_bytes()).hexdigest()
    metadata = {
        "schema_version": "1",
        "reference_id": "demo",
        "name": "Demo reference",
        "species": "human",
        "tissue": "blood",
        "disease": "healthy",
        "assay": "RNA",
        "platform": "synthetic",
        "resource_version": "1",
        "matrix": "reference.csv",
        "checksum": "wrong" if bad_checksum else checksum,
    }
    (root / "reference.json").write_text(json.dumps(metadata))


def test_local_reference_provider_filters_validates_and_loads(tmp_path: Path) -> None:
    _reference_fixture(tmp_path)
    provider = LocalReferenceProvider(tmp_path)
    assert provider.status().available
    summaries = provider.list_references(
        ReferenceFilters(species="human", tissue="BLOOD", disease="healthy")
    )
    assert [item.reference_id for item in summaries] == ["demo"]
    loaded = provider.load_reference("demo")
    assert len(loaded.checksum) == 64
    assert loaded.metadata["gene_count"] == 2
    assert loaded.metadata["cell_types"] == ["T cell", "B cell"]
    with pytest.raises(ResourceNotFoundError, match="Unknown reference"):
        provider.load_reference("missing")


def test_local_reference_provider_enforces_checksum(tmp_path: Path) -> None:
    _reference_fixture(tmp_path, bad_checksum=True)
    with pytest.raises(ProviderSchemaError, match="checksum"):
        LocalReferenceProvider(tmp_path).load_reference("demo")
