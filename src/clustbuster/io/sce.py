"""Pure-Python importer for the tested in-memory SingleCellExperiment RDS contract."""

from __future__ import annotations

import warnings as python_warnings
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from clustbuster.io.validation import (
    AnnDataValidationError,
    candidate_cluster_columns,
    compatible_embeddings,
    expression_sources,
    make_identifiers_unique,
    validate_adata,
)
from clustbuster.models import ImportOptions, ImportReport, ImportResult, ObjectFormat, ProbeResult


class SceImportError(ValueError):
    """An actionable error raised at the SCE compatibility boundary."""


def _first_string(value: Any) -> str | None:
    if isinstance(value, str):
        return value
    try:
        values = list(value)
    except TypeError:
        return None
    return str(values[0]) if values else None


def _coordinates(value: Any, coordinate: str) -> tuple[str, ...]:
    coordinates = getattr(value, "coords", {})
    selected = coordinates.get(coordinate)
    values = getattr(selected, "values", selected)
    if values is None:
        return ()
    return tuple(str(item) for item in values)


def _list_data(value: Any) -> dict[Any, Any]:
    return value.__dict__.get("listData", {}) if value is not None else {}


def _dframe(value: Any, names: tuple[str, ...], label: str) -> pd.DataFrame:
    if value is None:
        return pd.DataFrame(index=names)
    raw_rows = value.__dict__.get("rownames")
    if not isinstance(raw_rows, str):
        rows = tuple(str(item) for item in raw_rows)
        if rows and rows != names:
            raise SceImportError(f"{label} row names do not align with object identifiers")
    columns: dict[str, Any] = {}
    for raw_name, raw_values in _list_data(value).items():
        values = np.asarray(raw_values)
        if values.ndim != 1 or len(values) != len(names):
            raise SceImportError(f"{label} column {str(raw_name)!r} has incompatible dimensions")
        columns[str(raw_name)] = raw_values
    return pd.DataFrame(columns, index=names)


def _decode_source(path: Path) -> Any:
    try:
        from readseurat.rdata import read_rds

        with python_warnings.catch_warnings():
            python_warnings.simplefilter("ignore")
            source = read_rds(str(path))
    except Exception as exc:
        raise SceImportError(
            "The RDS could not be decoded as a supported in-memory SingleCellExperiment "
            f"object ({type(exc).__name__})."
        ) from exc
    source_class = _first_string(source.__dict__.get("class"))
    if source_class != "SingleCellExperiment":
        raise SceImportError(
            f"The RDS contains {source_class or 'an unknown R object'}, not a "
            "SingleCellExperiment object"
        )
    return source


class SceImporter:
    """Import a tested, in-memory SingleCellExperiment RDS into AnnData."""

    def probe(self, path: Path) -> ProbeResult:
        if path.suffix.lower() != ".rds":
            return ProbeResult(ObjectFormat.UNKNOWN, 0.0, "Not an RDS filename")
        try:
            source = _decode_source(path)
        except SceImportError as exc:
            return ProbeResult(ObjectFormat.UNKNOWN, 0.0, str(exc))
        if _first_string(source.__dict__.get("class")) == "SingleCellExperiment":
            return ProbeResult(
                ObjectFormat.SINGLE_CELL_EXPERIMENT,
                1.0,
                "RDS containing a SingleCellExperiment S4 object",
            )
        return ProbeResult(ObjectFormat.UNKNOWN, 0.0, "Not a SingleCellExperiment object")

    def load(self, path: Path, options: ImportOptions | None = None) -> ImportResult:
        del options
        if not path.is_file():
            raise SceImportError(f"Input file does not exist: {path.name}")
        unsupported: list[str] = []
        warnings: list[str] = []
        try:
            source = _decode_source(path)
            source_data = source.__dict__
            assays_container = source_data.get("assays")
            assays_data = getattr(assays_container, "data", None)
            assays = _list_data(assays_data)
            if not assays:
                raise SceImportError("The SingleCellExperiment contains no readable assays")

            row_ranges = source_data.get("rowRanges")
            partitioning = getattr(row_ranges, "partitioning", None)
            raw_features = getattr(partitioning, "NAMES", ())
            features = tuple(str(value) for value in raw_features)
            col_data = source_data.get("colData")
            raw_cells = getattr(col_data, "rownames", ())
            cells = (
                tuple(str(value) for value in raw_cells)
                if not isinstance(raw_cells, str)
                else ()
            )
            if not features or not cells:
                raise SceImportError("The SCE does not expose feature and cell identifiers")
            if len(features) != len(set(features)) or len(cells) != len(set(cells)):
                raise SceImportError("The SCE contains missing or duplicate identifiers")

            aligned: dict[str, Any] = {}
            for raw_name, matrix in assays.items():
                name = str(raw_name)
                assay_features = _coordinates(matrix, "dim_0")
                assay_cells = _coordinates(matrix, "dim_1")
                if assay_features and assay_features != features:
                    raise SceImportError(f"Assay {name!r} does not align to SCE identifiers")
                if assay_cells and assay_cells != cells:
                    raise SceImportError(f"Assay {name!r} does not align to SCE identifiers")
                values = getattr(matrix, "values", matrix)
                if getattr(values, "shape", None) != (len(features), len(cells)):
                    raise SceImportError(f"Assay {name!r} has incompatible dimensions")
                if not sparse.issparse(values) and not isinstance(values, np.ndarray):
                    raise SceImportError(
                        f"Assay {name!r} uses an unsupported delayed or custom representation"
                    )
                aligned[name] = values.T

            primary_name = "logcounts" if "logcounts" in aligned else (
                "counts" if "counts" in aligned else next(iter(aligned))
            )
            obs = _dframe(col_data, cells, "colData")
            var = _dframe(source_data.get("elementMetadata"), features, "rowData")
            adata = ad.AnnData(X=aligned[primary_name], obs=obs, var=var)
            for name, matrix in aligned.items():
                adata.layers[name] = matrix

            internal_col_data = source_data.get("int_colData")
            internal_items = _list_data(internal_col_data)
            reductions = _list_data(internal_items.get("reducedDims"))
            for raw_name, values in reductions.items():
                reduction_cells = _coordinates(values, "dim_0")
                if reduction_cells != cells:
                    raise SceImportError(
                        f"Reduced dimension {str(raw_name)!r} does not align to cell identifiers"
                    )
                array = np.asarray(values)
                if array.ndim != 2 or array.shape[0] != len(cells):
                    raise SceImportError(
                        f"Reduced dimension {str(raw_name)!r} has incompatible dimensions"
                    )
                adata.obsm[f"X_{str(raw_name).lower()}"] = array

            for component, description in (
                ("altExps", "Alternative experiments"),
                ("colPairs", "Column pair data"),
            ):
                if _list_data(internal_items.get(component)):
                    unsupported.append(f"{description} are not imported")
            internal_rows = _list_data(source_data.get("int_elementMetadata"))
            if _list_data(internal_rows.get("rowPairs")):
                unsupported.append("Row pair data are not imported")

            validate_adata(adata)
            warnings.extend(make_identifiers_unique(adata))
            validate_adata(adata)
        except (AnnDataValidationError, SceImportError):
            raise
        except Exception as exc:
            raise SceImportError(
                "The file could not be imported as a supported SingleCellExperiment. "
                "ClustBuster currently supports tested in-memory assays, DFrame metadata, "
                f"and reduced dimensions ({type(exc).__name__})."
            ) from exc

        report = ImportReport(
            source_format=ObjectFormat.SINGLE_CELL_EXPERIMENT,
            source_filename=path.name,
            cell_count=adata.n_obs,
            feature_count=adata.n_vars,
            sparse=sparse.issparse(adata.X),
            layers=tuple(str(name) for name in adata.layers),
            embeddings=compatible_embeddings(adata),
            candidate_cluster_columns=candidate_cluster_columns(adata.obs),
            expression_sources=expression_sources(adata),
            warnings=warnings,
            mapping_decisions={
                "X": f"SCE assay {primary_name!r} -> adata.X",
                "obs": "SCE colData -> adata.obs",
                "var": "SCE rowData -> adata.var",
                "reductions": "SCE reducedDims -> adata.obsm['X_<name>']",
                "layers": "All aligned SCE assays -> adata.layers",
            },
            unsupported_components=unsupported,
        )
        return ImportResult(adata=adata, report=report)
