from __future__ import annotations

import hashlib
from pathlib import Path

import duckdb
import pytest

from clustbuster.models import ReferenceFilters
from clustbuster.resources.errors import ProviderSchemaError
from clustbuster.resources.markers.duckdb import DuckDbMarkerProvider
from clustbuster.resources.references.duckdb import DuckDbReferenceProvider


def _marker_database(path: Path) -> None:
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            """CREATE TABLE marker_records (
                gene_symbol VARCHAR,
                species VARCHAR,
                cell_type VARCHAR,
                marker_direction VARCHAR,
                tissue VARCHAR,
                confidence VARCHAR,
                bbsr_verified BOOLEAN,
                source_titles VARCHAR,
                source_identifiers VARCHAR
            )"""
        )
        connection.executemany(
            "INSERT INTO marker_records VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                (
                    "CD3D",
                    "Homo sapiens",
                    "T cell",
                    "positive",
                    "Blood",
                    "high",
                    True,
                    "Paper A",
                    "PMID:1",
                ),
                (
                    "IL7R",
                    "Homo sapiens",
                    "T cell",
                    "positive",
                    "Blood",
                    "moderate",
                    False,
                    "Paper B",
                    "10.1/example",
                ),
                (
                    "Cd3d",
                    "Mus musculus",
                    "T cell",
                    "positive",
                    "Spleen",
                    "high",
                    False,
                    "Paper C",
                    "PMID:2",
                ),
            ],
        )
        connection.execute("CREATE VIEW marker_atlas AS SELECT * FROM marker_records")


def test_duckdb_marker_provider_queries_consumer_view(tmp_path: Path) -> None:
    path = tmp_path / "markers.duckdb"
    _marker_database(path)
    provider = DuckDbMarkerProvider(path, resource_version="preview-1")

    assert provider.status().available
    assert provider.status().version == "preview-1"
    assert provider.list_facets().species == ("Homo sapiens", "Mus musculus")
    results = provider.search_cell_types("t cell", species="homo sapiens", tissue="blood")
    assert [(result.cell_type, result.marker_count) for result in results] == [("T cell", 2)]
    marker_set = provider.get_markers("T CELL", species="Homo sapiens", tissue="Blood")
    assert [record.gene for record in marker_set.records] == ["CD3D", "IL7R"]
    assert marker_set.records[0].confidence is None
    assert marker_set.records[0].confidence_label == "high"
    assert marker_set.records[0].verified is True


def test_duckdb_marker_provider_reports_incompatible_view(tmp_path: Path) -> None:
    path = tmp_path / "bad.duckdb"
    with duckdb.connect(str(path)) as connection:
        connection.execute("CREATE VIEW marker_atlas AS SELECT 1 AS unexpected")
    provider = DuckDbMarkerProvider(path)
    assert not provider.status().available
    with pytest.raises(ProviderSchemaError, match="missing consumer fields"):
        provider.list_facets()


def _reference_catalog(path: Path, files_root: Path) -> tuple[str, Path]:
    matrix_id = "fad84409-9f3d-4005-a5c9-64719d8100e4"
    matrix = files_root / "reference.tsv"
    matrix.write_text("gene\tT cell\tB cell\nCD3D\t4\t0\nMS4A1\t0\t5\n")
    checksum = hashlib.sha256(matrix.read_bytes()).hexdigest()
    with duckdb.connect(str(path)) as connection:
        connection.execute(
            """CREATE TABLE reference_matrices (
                matrix_id UUID,
                title VARCHAR,
                species_scientific_name VARCHAR,
                tissue VARCHAR,
                condition VARCHAR,
                assay VARCHAR,
                platform VARCHAR,
                stored_filename VARCHAR,
                sha256 VARCHAR,
                row_count BIGINT,
                column_count BIGINT,
                column_names VARCHAR[],
                normalization VARCHAR,
                value_type VARCHAR,
                source_title VARCHAR,
                citation VARCHAR,
                doi VARCHAR,
                pmid VARCHAR,
                source_url VARCHAR,
                data_license VARCHAR
            )"""
        )
        connection.execute(
            """INSERT INTO reference_matrices VALUES (
                ?, 'PBMC reference', 'Homo sapiens', 'Blood', 'healthy', 'scRNA-seq',
                '10x', 'reference.tsv', ?, 2, 3, ['gene', 'T cell', 'B cell'],
                'log1p', 'mean expression', 'Fixture', NULL, NULL, NULL, NULL, 'CC0'
            )""",
            [matrix_id, checksum],
        )
    return matrix_id, matrix


def test_duckdb_reference_provider_resolves_and_validates_matrix(tmp_path: Path) -> None:
    files_root = tmp_path / "files"
    files_root.mkdir()
    catalog = tmp_path / "references.duckdb"
    matrix_id, matrix = _reference_catalog(catalog, files_root)
    provider = DuckDbReferenceProvider(catalog, files_root, resource_version="preview-1")

    assert provider.status().available
    summaries = provider.list_references(
        ReferenceFilters(species="homo sapiens", tissue="blood", disease="HEALTHY")
    )
    assert [summary.reference_id for summary in summaries] == [matrix_id]
    loaded = provider.load_reference(matrix_id)
    assert loaded.matrix_path == matrix
    assert loaded.metadata["gene_count"] == 2
    assert loaded.metadata["cell_types"] == ["T cell", "B cell"]
    assert loaded.summary.resource_version == "preview-1"


def test_duckdb_reference_provider_rejects_changed_bytes(tmp_path: Path) -> None:
    files_root = tmp_path / "files"
    files_root.mkdir()
    catalog = tmp_path / "references.duckdb"
    matrix_id, matrix = _reference_catalog(catalog, files_root)
    matrix.write_text(matrix.read_text() + "NKG7\t3\t0\n")

    with pytest.raises(ProviderSchemaError, match="checksum"):
        DuckDbReferenceProvider(catalog, files_root).load_reference(matrix_id)


@pytest.mark.parametrize("query", ["cd3d", "sapiens", "blood", "high", "true", "paper a", "PMID:1"])
def test_duckdb_marker_search_matches_all_fields(tmp_path: Path, query: str) -> None:
    path = tmp_path / "markers.duckdb"
    _marker_database(path)
    provider = DuckDbMarkerProvider(path)
    results = provider.search_cell_types(query, species="Homo sapiens")
    assert [item.cell_type for item in results] == ["T cell"]
    assert provider.search_cell_types(query, tissue="missing") == []
    assert provider.search_cell_types("' OR 1=1 --") == []
