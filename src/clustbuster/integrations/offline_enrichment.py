"""Downloaded GMT libraries and network-free, measured-background ORA."""

from __future__ import annotations

import hashlib
import json
import re
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request
from uuid import uuid4

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from clustbuster.core.enrichment import EnrichmentError, EnrichmentResult
from clustbuster.integrations.enrichr import BASE_URL, UrlReader, _read_url
from clustbuster.integrations.gene_species import MOUSE_LIBRARIES, MOUSE_RELEASE, detect_species


def parse_gmt(text: str) -> dict[str, frozenset[str]]:
    terms: dict[str, set[str]] = {}
    for line_number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 3 or not fields[0].strip():
            raise EnrichmentError(f"Invalid gene-set library at line {line_number}")
        genes = {gene.split(",", 1)[0].strip() for gene in fields[2:] if gene.strip()}
        genes.discard("")
        if genes:
            terms.setdefault(fields[0].strip(), set()).update(genes)
    if not terms:
        raise EnrichmentError("The gene-set library contains no usable gene sets")
    return {term: frozenset(genes) for term, genes in terms.items()}


class GeneSetCache:
    def __init__(self, root: Path) -> None:
        self.root = root

    def path(self, library: str) -> Path:
        if not re.fullmatch(r"[A-Za-z0-9_-]+", library):
            raise EnrichmentError("Invalid gene-set library name")
        return self.root / f"{library}.gmt"

    def load(self, library: str) -> dict[str, frozenset[str]]:
        try:
            return parse_gmt(self.path(library).read_text(encoding="utf-8"))
        except FileNotFoundError as exc:
            raise EnrichmentError(
                f"{library} is not available offline. Download it in Settings while online."
            ) from exc
        except (OSError, UnicodeError) as exc:
            raise EnrichmentError(f"Cannot read offline library {library}: {exc}") from exc

    def download(
        self, library: str, *, reader: UrlReader = _read_url, timeout: float = 120,
        species: str = "human",
    ) -> int:
        path = self.path(library)
        url = f"{BASE_URL}/geneSetLibrary?" + urlencode({"mode": "text", "libraryName": library})
        if species == "mouse":
            if library not in MOUSE_LIBRARIES:
                raise EnrichmentError("No native mouse version of this library is configured. "
                                      "Choose GO or Reactome.")
            path = GeneSetCache(self.root / "mouse").path(library)
            collection = MOUSE_LIBRARIES[library]
            url = ("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
                   f"{MOUSE_RELEASE}/{collection}.v{MOUSE_RELEASE}.symbols.gmt")
        elif species != "human":
            raise EnrichmentError("Select Human or Mouse for downloading libraries")
        try:
            payload = reader(Request(url), timeout)
            terms = parse_gmt(payload.decode("utf-8"))
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_suffix(f".{uuid4().hex}.tmp")
            try:
                temporary.write_bytes(payload)
                temporary.replace(path)
            finally:
                temporary.unlink(missing_ok=True)
            metadata = {
                "library": library, "species": species, "source_url": url,
                "downloaded_at": datetime.now(UTC).isoformat(),
                "sha256": hashlib.sha256(payload).hexdigest(), "terms": len(terms),
            }
            path.with_suffix(".json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
            return len(terms)
        except (OSError, UnicodeError) as exc:
            raise EnrichmentError(f"Could not download {library}: {exc}") from exc


class OfflineEnrichmentClient:
    """One-sided Fisher equivalent with BH across all testable terms in one library."""

    def __init__(
        self, cache: GeneSetCache, library: str, background: tuple[str, ...],
        species: str = "human",
    ) -> None:
        self.library = library
        self.species = detect_species(cache.root, background) if species == "auto" else species
        if self.species not in {"human", "mouse"}:
            raise EnrichmentError("Select Human or Mouse in Settings")
        if self.species == "mouse":
            if library not in MOUSE_LIBRARIES:
                raise EnrichmentError("No native mouse version of this library is configured. "
                                      "Choose GO or Reactome; no homolog conversion is used.")
            cache = GeneSetCache(cache.root / "mouse")
        self.background = frozenset(gene.strip() for gene in background if gene.strip())
        if len(self.background) < 2:
            raise EnrichmentError("Offline enrichment requires at least two background genes")
        self.terms = {
            term: overlap for term, genes in cache.load(library).items()
            if (overlap := genes & self.background)
        }
        if not self.terms:
            raise EnrichmentError("No library genes match the expression source's gene symbols")

    def enrich(self, genes: tuple[str, ...], *, description: str) -> EnrichmentResult:
        query = frozenset(gene.strip() for gene in genes if gene.strip()) & self.background
        if len(query) < 2:
            raise EnrichmentError("Enrichment requires at least two marker genes in the background")
        population = len(self.background)
        query_size = len(query)
        rows = []
        for term, members in self.terms.items():
            overlap = query & members
            a = len(overlap)
            b, c = query_size - a, len(members) - a
            d = population - a - b - c
            denominator = b * c
            odds = a * d / denominator if denominator else (np.inf if a * d else np.nan)
            rows.append({
                "term": term,
                "odds_ratio": odds,
                # Enrichr's combined score is not calculated by this local method.
                "combined_score": np.nan,
                "overlap_genes": ", ".join(sorted(overlap)),
                "overlap_count": a, "term_size": len(members),
                "query_size": query_size, "background_size": population,
            })
        table = pd.DataFrame(rows)
        table["p_value"] = hypergeom.sf(
            table["overlap_count"].to_numpy() - 1, population,
            table["term_size"].to_numpy(), query_size,
        )
        table = table.sort_values(["p_value", "term"]).reset_index(drop=True)
        count = len(table)
        adjusted = table["p_value"].to_numpy() * count / np.arange(1, count + 1)
        table["adjusted_p_value"] = np.minimum(np.minimum.accumulate(adjusted[::-1])[::-1], 1)
        table = table.loc[table["overlap_count"] > 0].reset_index(drop=True)
        if table.empty:
            raise EnrichmentError("No gene-set terms overlap these marker genes")
        table.insert(0, "rank", np.arange(1, len(table) + 1))
        return EnrichmentResult(
            genes=tuple(sorted(query)), library=self.library,
            source=(f"Offline ORA ({self.species}; measured-gene background"
                    + (f"; MSigDB {MOUSE_RELEASE}" if self.species == "mouse" else "") + ")"),
            values=table,
        )
