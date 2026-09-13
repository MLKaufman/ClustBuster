from pathlib import Path

from clustbuster.models import ReferenceFilters, ReferenceSummary
from clustbuster.resources.references.details import (
    cell_type_table,
    reference_compatibility,
    reference_label,
    search_references,
)
from clustbuster.resources.references.local import LocalReferenceProvider


def test_search_covers_cell_types_metadata_and_punctuation() -> None:
    ref = ReferenceSummary(
        'id-12345', 'Atlas', 'Mus musculus', tissue='mammary', assay='scRNA-seq',
        cell_types=('Prolif.Fibroblasts', 'Luminal_progenitor'),
        metadata={'submitter': 'Alice', 'description': 'Healthy adult tissue'},
    )
    assert search_references([ref], 'prolif fibroblast') == [ref]
    assert search_references([ref], 'luminal progenitor') == [ref]
    assert search_references([ref], 'alice adult') == [ref]
    assert search_references([ref], 'missing') == []
    assert search_references([ref], '') == [ref]
    assert 'mammary | scRNA-seq | 2 cell types | id-12345' in reference_label(ref)


def test_cell_coverage_includes_optional_metadata_and_synonyms() -> None:
    ref = ReferenceSummary(
        'r', 'Atlas', 'human', cell_types=('T cell', 'B cell'),
        metadata={
            'cell_type_metadata': '{"T cell": {"cell_count": 120, "aliases": "T lymphocyte"}}',
        },
    )
    table = cell_type_table(ref, 'lymphocyte')
    assert table['Cell type'].tolist() == ['T cell']
    assert table['Cells'].tolist() == [120]
    assert 'Donors' not in table


def test_local_header_discovery_and_compatibility_use_actual_genes(tmp_path: Path) -> None:
    import json
    (tmp_path / 'matrix.tsv').write_text('gene\tT cell\nCD3D\t1\nLYZ\t2\nAPOE\t3\nMISSING\t4\n')
    (tmp_path / 'ref.json').write_text(json.dumps({
        'schema_version': 1, 'reference_id': 'r', 'name': 'Atlas', 'species': 'human',
        'matrix': 'matrix.tsv', 'resource_version': 'v1', 'description': 'Study notes',
    }))
    provider = LocalReferenceProvider(tmp_path)
    summary = provider.list_references(ReferenceFilters())[0]
    assert summary.cell_types == ('T cell',)
    assert summary.metadata['description'] == 'Study notes'
    table = reference_compatibility(provider.load_reference('r'), ('cd3d', 'LYZ', 'Apoe', 'apoe'))
    assert dict(zip(table['Reference gene'], table['Match'], strict=True)) == {
        'CD3D': 'Shared', 'LYZ': 'Shared', 'APOE': 'Ambiguous', 'MISSING': 'Missing from dataset',
    }
