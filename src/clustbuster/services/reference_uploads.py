"""Validate session-local user reference matrices without changing catalog resources."""
from __future__ import annotations

import csv
import hashlib
import shutil
from collections.abc import Mapping
from pathlib import Path
from typing import Any
from uuid import uuid4

import numpy as np
import pandas as pd

from clustbuster.models import LoadedReference, ReferenceSummary
from clustbuster.services.imports import SessionFiles, UploadError


def import_reference_upload(
    upload: Mapping[str, Any], session_files: SessionFiles, *, max_upload_mb: int,
    name: str = "", species: str = "Unknown", tissue: str = "", normalization: str = "",
) -> LoadedReference:
    original_name = Path(str(upload.get("name", "reference.tsv"))).name
    suffix = Path(original_name).suffix.casefold()
    if suffix not in {".csv", ".tsv", ".txt"}:
        raise UploadError("Upload a CSV or tab-delimited TSV/TXT reference matrix.")
    source = Path(str(upload.get("datapath", "")))
    if not source.is_file() or not source.stat().st_size:
        raise UploadError("The uploaded reference is empty or unavailable. Upload it again.")
    size = source.stat().st_size
    if size > max_upload_mb * 1024 * 1024:
        raise UploadError(f"The reference exceeds the {max_upload_mb} MB upload limit.")
    if shutil.disk_usage(session_files.root).free < size * 2:
        raise UploadError("Not enough session disk space for this reference.")
    folder = session_files.uploads / f"refmat-{uuid4().hex}"
    folder.mkdir(mode=0o700)
    try:
        delimiter = "," if suffix == ".csv" else "\t"
        with source.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.reader(stream, delimiter=delimiter)
            header = next(reader)
            for line_number, row in enumerate(reader, 2):
                if not row:
                    continue
                if len(row) != len(header):
                    raise UploadError(
                        f"Row {line_number} has a different number of columns than the header."
                    )
        if len(header) < 2:
            raise UploadError("Use genes in the first column and one or more cell-type columns.")
        cell_types = [column.strip() for column in header[1:]]
        if any(not column for column in cell_types) or len(set(cell_types)) != len(cell_types):
            raise UploadError("Cell-type column names must be nonempty and unique.")
        if "gene" in cell_types:
            raise UploadError("The name 'gene' is reserved for the first column.")
        table = pd.read_csv(source, sep=delimiter, dtype=str, keep_default_na=False,
                            encoding="utf-8-sig", index_col=False)
        table.columns = ["gene", *cell_types]
        table["gene"] = table["gene"].str.strip()
        if table.empty or (table["gene"] == "").any() or not table["gene"].is_unique:
            raise UploadError("Reference gene identifiers must be nonempty and unique.")
        numeric = table[cell_types].apply(pd.to_numeric, errors="coerce")
        if not np.isfinite(numeric.to_numpy(dtype=float)).all():
            raise UploadError(
                "All cell-type expression values must be numeric and finite (no missing values)."
            )
        table[cell_types] = numeric
        path = folder / "reference.tsv"
        table.to_csv(path, sep="\t", index=False)
        checksum = hashlib.sha256(path.read_bytes()).hexdigest()
        metadata = {
            "description": "User-uploaded reference; available only in this session.",
            "source_title": original_name, "submitter": "Session upload",
            "gene_count": len(table), "row_count": len(table),
            "cell_types": cell_types, "normalization": normalization.strip() or None,
            "original_filename": original_name,
        }
        summary = ReferenceSummary(
            reference_id=f"upload:{folder.name}", name=name.strip() or Path(original_name).stem,
            species=species, tissue=tissue.strip() or None,
            resource_version=f"sha256:{checksum[:12]}",
            cell_types=tuple(cell_types), metadata=metadata,
        )
        return LoadedReference(summary, path, checksum, metadata)
    except Exception:
        shutil.rmtree(folder, ignore_errors=True)
        raise
