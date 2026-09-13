"""Reference discovery and compatibility using the annotation engine's gene matcher."""
from __future__ import annotations

import json
import re
from typing import Any

import pandas as pd

from clustbuster.integrations.pyclustifyr import _match_reference_genes, _read_reference
from clustbuster.models import LoadedReference, ReferenceSummary


def reference_label(reference: ReferenceSummary) -> str:
    fields = [reference.name, reference.species, reference.tissue, reference.disease,
              reference.assay]
    parts = list(dict.fromkeys(str(value) for value in fields if value))
    if reference.cell_types:
        parts.append(f"{len(reference.cell_types)} cell types")
    parts.append(reference.reference_id.removeprefix("upload:refmat-")[:8])
    return " | ".join(parts)


def search_references(references: list[ReferenceSummary], query: str) -> list[ReferenceSummary]:
    def normalize(value: Any) -> str:
        return re.sub(r"[\W_]+", " ", str(value).casefold())
    tokens = normalize(query).split()
    return [reference for reference in references if all(
        token in normalize(" ".join([
            reference_label(reference), *reference.cell_types,
            *(str(value) for value in reference.metadata.values() if value is not None),
        ])) for token in tokens
    )]


def reference_compatibility(
    reference: LoadedReference, query_genes: tuple[str, ...],
) -> pd.DataFrame:
    matrix = _read_reference(reference)
    _, matched, missing, ambiguous = _match_reference_genes(
        pd.Index(query_genes), pd.Index(matrix.index.astype(str)),
    )
    return pd.DataFrame(
        [(gene, status) for status, genes in
         (("Shared", matched), ("Missing from dataset", missing), ("Ambiguous", ambiguous))
         for gene in genes], columns=["Reference gene", "Match"],
    )


def cell_type_table(reference: ReferenceSummary, query: str = "") -> pd.DataFrame:
    metadata = reference.metadata.get("cell_type_metadata") or {}
    if isinstance(metadata, str):
        try:
            metadata = json.loads(metadata)
        except ValueError:
            metadata = {}
    if not isinstance(metadata, dict):
        metadata = {}
    rows = []
    for cell in reference.cell_types:
        details = metadata.get(cell, {})
        if not isinstance(details, dict):
            details = {}
        row = {"Cell type": cell, **{label: details.get(key, "Not provided")
               for label, key in (("Cells", "cell_count"), ("Samples", "sample_count"),
                                  ("Donors", "donor_count"), ("Description", "description"),
                                  ("Ontology ID", "ontology_id"), ("Aliases", "aliases"))}}
        normalized = re.sub(r"[\W_]+", " ", str(row).casefold())
        if all(token in normalized for token in re.sub(r"[\W_]+", " ", query.casefold()).split()):
            rows.append(row)
    table = pd.DataFrame(rows, columns=["Cell type", "Cells", "Samples", "Donors", "Description",
                                       "Ontology ID", "Aliases"])
    return table[[column for column in table if column == "Cell type" or
                  any(str(value) not in {"Not provided", "None", ""} for value in table[column])]]
