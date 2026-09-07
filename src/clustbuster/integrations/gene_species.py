"""Detect species from exact symbol overlap with local species-specific libraries."""
from pathlib import Path

from clustbuster.core.enrichment import EnrichmentError

MOUSE_LIBRARIES = {
    "GO_Biological_Process_2025": "m5.go.bp",
    "GO_Molecular_Function_2025": "m5.go.mf",
    "GO_Cellular_Component_2025": "m5.go.cc",
    "Reactome_Pathways_2024": "m2.cp.reactome",
}
MOUSE_RELEASE = "2025.1.Mm"


def detect_species(root: Path, genes: tuple[str, ...]) -> str:
    symbols = set(genes)
    if any(g.startswith(("ENSG", "ENSMUSG")) for g in symbols):
        raise EnrichmentError(
            "Ensembl IDs detected. Use gene symbols as the expression source's gene names."
        )
    references: dict[str, set[str]] = {"human": set(), "mouse": set()}
    for species in references:
        path = root / ("mouse" if species == "mouse" else "")
        for gmt in path.glob("*.gmt"):
            for line in gmt.read_text(encoding="utf-8").splitlines():
                references[species].update(g.split(",", 1)[0] for g in line.split("\t")[2:])
    if not all(references.values()):
        raise EnrichmentError(
            "Download human and mouse libraries in Settings to enable species detection, "
            "or select the dataset species manually."
        )
    hits = {s: len(symbols & (ref - references["mouse" if s == "human" else "human"]))
            for s, ref in references.items()}
    total = sum(hits.values())
    best = max(hits, key=lambda s: hits[s])
    if total < 3 or hits[best] / total < 0.8:
        raise EnrichmentError(
            "Species detection is ambiguous (human: " + str(hits["human"])
            + ", mouse: " + str(hits["mouse"]) + " matching symbols). "
            "Select Human or Mouse in Settings."
        )
    return best
