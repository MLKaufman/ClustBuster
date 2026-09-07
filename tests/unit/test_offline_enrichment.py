from pathlib import Path
from urllib.request import Request

import pytest
from scipy.stats import fisher_exact

from clustbuster.core.enrichment import EnrichmentError
from clustbuster.integrations.offline_enrichment import (
    GeneSetCache,
    OfflineEnrichmentClient,
    parse_gmt,
)


@pytest.fixture
def cache(tmp_path: Path) -> GeneSetCache:
    cache = GeneSetCache(tmp_path)
    cache.path("test").write_text("A term\tdesc\tA\tB\nB term\t\tC\tD\tE\nC term\t\tA\tC\n")
    return cache


def test_offline_fisher_and_bh_include_zero_overlap_terms(cache: GeneSetCache) -> None:
    result = OfflineEnrichmentClient(cache, "test", tuple("ABCDEFGHIJ")).enrich(
        ("A", "B", "A", "outside"), description="test",
    )
    values = result.values.set_index("term")
    assert values.loc["A term", "p_value"] == pytest.approx(
        fisher_exact([[2, 0], [0, 8]], alternative="greater").pvalue
    )
    assert values.loc["A term", "adjusted_p_value"] == pytest.approx(3 / 45)
    assert "B term" not in values.index
    assert result.genes == ("A", "B")
    assert set(values["background_size"]) == {10}
    assert values["combined_score"].isna().all()


def test_offline_never_accesses_network(
    cache: GeneSetCache, monkeypatch: pytest.MonkeyPatch,
) -> None:
    def forbidden(*args: object, **kwargs: object) -> None:
        raise AssertionError("Offline enrichment attempted a network request")
    monkeypatch.setattr("clustbuster.integrations.enrichr.urlopen", forbidden)
    monkeypatch.setattr(GeneSetCache, "download", forbidden)
    result = OfflineEnrichmentClient(cache, "test", tuple("ABCDEFGHIJ")).enrich(
        ("A", "B"), description="offline",
    )
    assert result.source.startswith("Offline")
    with pytest.raises(EnrichmentError, match="Download it in Settings"):
        OfflineEnrichmentClient(cache, "missing", tuple("ABCDEFGHIJ"))


def test_offline_uses_only_measured_genes(cache: GeneSetCache) -> None:
    client = OfflineEnrichmentClient(cache, "test", ("A", "B", "J", "K"))
    assert client.terms["C term"] == {"A"}
    assert "B term" not in client.terms
    result = client.enrich(("A", "B"), description="measured")
    assert result.values.iloc[0]["p_value"] == pytest.approx(1 / 6)
    with pytest.raises(EnrichmentError, match="at least two"):
        client.enrich(("A", "missing"), description="invalid")
    with pytest.raises(EnrichmentError, match="No gene-set terms overlap"):
        client.enrich(("J", "K"), description="no overlap")


def test_cache_validates_download_and_preserves_existing_library(cache: GeneSetCache) -> None:
    original = cache.path("test").read_bytes()
    def invalid(request: Request, timeout: float) -> bytes:
        return b"<html>service unavailable</html>"
    with pytest.raises(EnrichmentError, match="Invalid gene-set library"):
        cache.download("test", reader=invalid)
    assert cache.path("test").read_bytes() == original
    def valid(request: Request, timeout: float) -> bytes:
        assert "mode=text" in request.full_url and "libraryName=test" in request.full_url
        return b"New term\t\tA\tB\n"
    assert cache.download("test", reader=valid) == 1
    assert cache.load("test") == {"New term": {"A", "B"}}
    assert cache.path("test").with_suffix(".json").is_file()
    assert not list(cache.root.glob("*.tmp"))


def test_gmt_merges_duplicate_terms_and_rejects_unsafe_names(cache: GeneSetCache) -> None:
    assert parse_gmt("Term\tdesc\tA\tA\nTerm\t\tB,1\t\n") == {"Term": {"A", "B"}}
    with pytest.raises(EnrichmentError, match="Invalid gene-set library name"):
        cache.path("../escape")


def test_species_specific_libraries_preserve_mouse_biology(tmp_path: Path) -> None:
    from clustbuster.integrations.gene_species import detect_species
    library = "GO_Biological_Process_2025"
    cache = GeneSetCache(tmp_path)
    cache.path(library).write_text("Human term\t\tAPOE\tLYZ\tCD3D\n")
    mouse = GeneSetCache(tmp_path / "mouse")
    mouse.root.mkdir()
    mouse.path(library).write_text("Mouse specific term\t\tApoe\tLyz2\tCd3d\n")
    assert detect_species(tmp_path, ("APOE", "LYZ", "CD3D")) == "human"
    assert detect_species(tmp_path, ("Apoe", "Lyz2", "Cd3d")) == "mouse"
    client = OfflineEnrichmentClient(cache, library, ("Apoe", "Lyz2", "Cd3d", "Other"), "auto")
    result = client.enrich(("Apoe", "Lyz2"), description="mouse")
    assert result.genes == ("Apoe", "Lyz2")
    assert result.values.iloc[0]["term"] == "Mouse specific term"
    assert "2025.1.Mm" in result.source
    with pytest.raises(EnrichmentError, match="ambiguous"):
        detect_species(tmp_path, ("APOE", "LYZ", "Apoe", "Lyz2"))
    with pytest.raises(EnrichmentError, match="No native mouse"):
        OfflineEnrichmentClient(cache, "KEGG_2021_Human", ("Apoe", "Lyz2"), "mouse")
    with pytest.raises(EnrichmentError, match="Ensembl"):
        detect_species(tmp_path, ("ENSMUSG000001",))


def test_mouse_download_is_separate_and_network_free_enrichment(tmp_path: Path) -> None:
    cache = GeneSetCache(tmp_path)
    library = "GO_Biological_Process_2025"
    cache.path(library).write_text("Human\t\tAPOE\tLYZ\n")
    def reader(request: Request, timeout: float) -> bytes:
        assert "2025.1.Mm/m5.go.bp" in request.full_url
        return b"Mouse\t\tApoe\tLyz2\n"
    cache.download(library, reader=reader, species="mouse")
    assert "Human" in cache.load(library)
    result = OfflineEnrichmentClient(cache, library, ("Apoe", "Lyz2", "Cd3d"), "mouse").enrich(
        ("Apoe", "Lyz2"), description="mouse",
    )
    assert result.genes == ("Apoe", "Lyz2")
