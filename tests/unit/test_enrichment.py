from urllib.error import URLError
from urllib.request import Request

import pytest

from clustbuster.core.enrichment import EnrichmentError
from clustbuster.integrations.enrichr import DEFAULT_LIBRARY, EnrichrClient
from clustbuster.plotting.enrichment import enrichment_figure


def _fake_reader(request: Request, timeout: float) -> bytes:
    assert timeout == 5
    if request.full_url.endswith("/addList"):
        assert request.method == "POST"
        assert request.data is not None and b"CD3D\nIL7R" in request.data
        return b'{"userListId": 123, "shortId": "demo"}'
    assert "backgroundType=GO_Biological_Process_2025" in request.full_url
    return (
        b'{"GO_Biological_Process_2025": [[1, "T cell activation (GO:0042110)", '
        b'0.001, 9.5, 65.6, ["CD3D", "IL7R"], 0.01, 0.0, 0.0]]}'
    )


def test_enrichr_adapter_normalizes_results_and_plot() -> None:
    result = EnrichrClient(timeout_seconds=5, reader=_fake_reader).enrich(
        ("CD3D", "IL7R"), description="cluster 0 markers"
    )
    assert result.library == DEFAULT_LIBRARY
    assert result.values.loc[0, "term"] == "T cell activation (GO:0042110)"
    assert result.values.loc[0, "overlap_genes"] == "CD3D, IL7R"
    figure = enrichment_figure(result)
    assert figure.data[0].type == "bar"
    assert figure.layout.xaxis.title.text == "-log10 adjusted p-value"


def test_enrichr_adapter_requires_multiple_genes() -> None:
    with pytest.raises(EnrichmentError, match="at least two"):
        EnrichrClient(reader=_fake_reader).enrich(("CD3D",), description="one gene")


def test_enrichr_network_error_is_actionable() -> None:
    def unavailable(request: Request, timeout: float) -> bytes:
        raise URLError("offline")

    with pytest.raises(EnrichmentError, match="service is unavailable"):
        EnrichrClient(reader=unavailable).enrich(
            ("CD3D", "IL7R"), description="offline test"
        )
