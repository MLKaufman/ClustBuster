"""Narrow Enrichr adapter for GO Biological Process enrichment."""

from __future__ import annotations

import json
from collections.abc import Callable
from typing import Any, cast
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from uuid import uuid4

import pandas as pd

from clustbuster.core.enrichment import EnrichmentError, EnrichmentResult

DEFAULT_LIBRARY = "GO_Biological_Process_2025"
BASE_URL = "https://maayanlab.cloud/Enrichr"
UrlReader = Callable[[Request, float], bytes]


def _read_url(request: Request, timeout: float) -> bytes:
    with urlopen(request, timeout=timeout) as response:
        return cast(bytes, response.read())


def _multipart_body(fields: dict[str, str]) -> tuple[bytes, str]:
    boundary = f"clustbuster-{uuid4().hex}"
    chunks: list[bytes] = []
    for name, value in fields.items():
        chunks.extend(
            [
                f"--{boundary}\r\n".encode(),
                f'Content-Disposition: form-data; name="{name}"\r\n\r\n'.encode(),
                value.encode("utf-8"),
                b"\r\n",
            ]
        )
    chunks.append(f"--{boundary}--\r\n".encode())
    return b"".join(chunks), boundary


class EnrichrClient:
    """Submit gene symbols only and normalize Enrichr's JSON response."""

    def __init__(
        self,
        *,
        library: str = DEFAULT_LIBRARY,
        timeout_seconds: float = 20,
        reader: UrlReader = _read_url,
    ) -> None:
        self.library = library
        self.timeout_seconds = timeout_seconds
        self._reader = reader

    def _json(self, request: Request) -> Any:
        try:
            payload = self._reader(request, self.timeout_seconds)
            return json.loads(payload)
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            raise EnrichmentError(
                "The Enrichr service is unavailable; check network access and try again"
            ) from exc
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            raise EnrichmentError("Enrichr returned an unreadable response") from exc

    def enrich(self, genes: tuple[str, ...], *, description: str) -> EnrichmentResult:
        unique_genes = tuple(dict.fromkeys(gene.strip() for gene in genes if gene.strip()))
        if len(unique_genes) < 2:
            raise EnrichmentError("Enrichment requires at least two marker genes")
        body, boundary = _multipart_body(
            {"list": "\n".join(unique_genes), "description": description}
        )
        add_request = Request(
            f"{BASE_URL}/addList",
            data=body,
            headers={"Content-Type": f"multipart/form-data; boundary={boundary}"},
            method="POST",
        )
        added = self._json(add_request)
        try:
            user_list_id = int(added["userListId"])
        except (KeyError, TypeError, ValueError) as exc:
            raise EnrichmentError("Enrichr did not return a valid gene-list identifier") from exc

        query = urlencode(
            {"userListId": user_list_id, "backgroundType": self.library}
        )
        result_request = Request(f"{BASE_URL}/enrich?{query}", method="GET")
        response = self._json(result_request)
        try:
            rows = response[self.library]
        except (KeyError, TypeError) as exc:
            raise EnrichmentError(
                f"Enrichr did not return results for {self.library}"
            ) from exc
        normalized: list[dict[str, object]] = []
        for row in rows:
            if not isinstance(row, list) or len(row) < 7:
                continue
            overlap = row[5] if isinstance(row[5], list) else []
            normalized.append(
                {
                    "rank": int(row[0]),
                    "term": str(row[1]),
                    "p_value": float(row[2]),
                    "odds_ratio": float(row[3]),
                    "combined_score": float(row[4]),
                    "overlap_genes": ", ".join(str(gene) for gene in overlap),
                    "adjusted_p_value": float(row[6]),
                }
            )
        if not normalized:
            raise EnrichmentError("No GO Biological Process terms were returned for these genes")
        values = pd.DataFrame(normalized).sort_values(
            ["adjusted_p_value", "rank"], ascending=True
        )
        return EnrichmentResult(
            genes=unique_genes,
            library=self.library,
            source="Enrichr",
            values=values.reset_index(drop=True),
        )
