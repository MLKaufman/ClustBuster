import anndata as ad
import pandas as pd
import pytest
from scipy import sparse

from clustbuster.core.expression import (
    GeneMatchError,
    aggregate_dotplot,
    extract_expression,
    parse_gene_list,
)
from clustbuster.models import ExpressionSource


def _adata() -> ad.AnnData:
    matrix = sparse.csr_matrix([[2, 0, 1], [4, 0, 0], [0, 3, 1], [0, 5, 0]])
    adata = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame({"cluster": ["a", "a", "b", "b"]}, index=list("wxyz")),
        var=pd.DataFrame(index=["CD3D", "LYZ", "MS4A1"]),
    )
    adata.raw = adata.copy()
    adata.layers["double"] = matrix * 2
    return adata


def test_parse_gene_list_deduplicates_and_limits() -> None:
    assert parse_gene_list("CD3D, LYZ\nCD3D") == ("CD3D", "LYZ")
    assert len(parse_gene_list(" ".join(f"G{i}" for i in range(30)), limit=None)) == 30
    with pytest.raises(GeneMatchError, match="at least one"):
        parse_gene_list("  , ")
    with pytest.raises(GeneMatchError, match="limited"):
        parse_gene_list("a b c", limit=2)


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        (ExpressionSource.x(), [2, 4, 0, 0]),
        (ExpressionSource.raw(), [2, 4, 0, 0]),
        (ExpressionSource.named_layer("double"), [4, 8, 0, 0]),
    ],
)
def test_sparse_expression_extraction_supports_every_source(
    source: ExpressionSource, expected: list[int]
) -> None:
    result = extract_expression(_adata(), source, ("cd3d", "missing"))
    assert result.report.matched == ("CD3D",)
    assert result.report.missing == ("missing",)
    assert result.values["CD3D"].tolist() == expected


def test_no_matching_genes_has_actionable_error() -> None:
    with pytest.raises(GeneMatchError, match="None of the requested genes"):
        extract_expression(_adata(), ExpressionSource.x(), ("NOT_A_GENE",))


def test_dotplot_aggregation_reports_fraction_and_mean() -> None:
    result = aggregate_dotplot(_adata(), ExpressionSource.x(), "cluster", ("CD3D", "LYZ"))
    cluster_a_cd3d = result.values.query("cluster_display == 'a' and gene == 'CD3D'").iloc[0]
    cluster_b_lyz = result.values.query("cluster_display == 'b' and gene == 'LYZ'").iloc[0]
    assert cluster_a_cd3d["fraction_expressing"] == 1.0
    assert cluster_a_cd3d["mean_expression"] == 3.0
    assert cluster_b_lyz["fraction_expressing"] == 1.0
    assert cluster_b_lyz["mean_expression"] == 4.0
