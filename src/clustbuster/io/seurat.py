"""Pure-Python Seurat import boundary backed by a pinned readseurat release."""

from __future__ import annotations

import warnings as python_warnings
from pathlib import Path
from typing import Any

import numpy as np
import readseurat
from scipy import sparse

from clustbuster.io.validation import (
    AnnDataValidationError,
    candidate_cluster_columns,
    compatible_embeddings,
    expression_sources,
    make_identifiers_unique,
    validate_adata,
)
from clustbuster.models import (
    ImportOptions,
    ImportReport,
    ImportResult,
    ObjectFormat,
    ProbeResult,
)


class SeuratImportError(ValueError):
    """An actionable error raised at the Seurat compatibility boundary."""


def _first_string(value: Any) -> str | None:
    if isinstance(value, str):
        return value
    try:
        items = list(value)
    except TypeError:
        return None
    return str(items[0]) if items else None


def _coordinate_names(value: Any, coordinate: str = "dim_0") -> tuple[str, ...]:
    coordinates = getattr(value, "coords", {})
    selected = coordinates.get(coordinate)
    values = getattr(selected, "values", selected)
    if values is None:
        return ()
    return tuple(str(item) for item in values)


def _recover_rds_metadata(path: Path) -> tuple[tuple[str, ...], str, Any]:
    """Recover Assay5 coordinates that readseurat 0.1.0 converts as booleans."""

    try:
        from readseurat.rdata import read_rds

        with python_warnings.catch_warnings():
            python_warnings.simplefilter("ignore")
            source = read_rds(str(path))
        source_data = source.__dict__
        source_class = _first_string(source_data.get("class"))
        if source_class != "Seurat":
            raise SeuratImportError(
                f"The RDS contains {source_class or 'an unknown R object'}, not a Seurat object. "
                "Use the SingleCellExperiment importer for SCE files or convert the object to "
                "H5AD in R."
            )
        active_assay = _first_string(source_data.get("active.assay"))
        assays = source_data.get("assays")
        if not active_assay or not isinstance(assays, dict) or active_assay not in assays:
            raise SeuratImportError("The Seurat object has no readable active assay")
        features = _coordinate_names(getattr(assays[active_assay], "features", None))
        if not features:
            raise SeuratImportError("The active assay has no feature-name coordinate")
        return features, active_assay, source
    except SeuratImportError:
        raise
    except Exception as exc:
        raise SeuratImportError(
            "The RDS could not be decoded as a supported in-memory Seurat v5 object. "
            "On-disk layers, custom S4 extensions, and non-Seurat RDS files are not supported "
            f"({type(exc).__name__})."
        ) from exc


def _map_secondary_assays(
    source: Any,
    active_assay: str,
    adata: Any,
) -> tuple[list[str], list[str]]:
    """Map only secondary layers whose identifiers can be aligned exactly."""

    mapped: list[str] = []
    unsupported: list[str] = []
    assays = source.__dict__.get("assays", {})
    target_features = tuple(str(value) for value in adata.var_names)
    target_cells = tuple(str(value) for value in adata.obs_names)

    for raw_name, assay in assays.items():
        assay_name = str(raw_name)
        if assay_name == active_assay:
            continue
        features = _coordinate_names(getattr(assay, "features", None))
        cells = _coordinate_names(getattr(assay, "cells", None))
        layers = getattr(assay, "layers", {})
        if (
            len(features) != len(set(features))
            or len(cells) != len(set(cells))
            or set(features) != set(target_features)
            or set(cells) != set(target_cells)
        ):
            unsupported.append(
                f"Non-active assay {assay_name!r} has a different cell or feature space"
            )
            continue

        feature_positions = {name: index for index, name in enumerate(features)}
        cell_positions = {name: index for index, name in enumerate(cells)}
        feature_order = [feature_positions[name] for name in target_features]
        cell_order = [cell_positions[name] for name in target_cells]
        mapped_for_assay = 0
        for raw_layer_name, raw_matrix in layers.items():
            matrix = raw_matrix
            if matrix.shape == (len(features), len(cells)):
                aligned = matrix[feature_order, :][:, cell_order].T
            elif matrix.shape == (len(cells), len(features)):
                aligned = matrix[cell_order, :][:, feature_order]
            else:
                unsupported.append(
                    f"Non-active assay {assay_name!r} layer {str(raw_layer_name)!r} has "
                    f"incompatible shape {matrix.shape}"
                )
                continue
            layer_name = f"assay:{assay_name}:{raw_layer_name!s}"
            adata.layers[layer_name] = (
                aligned if sparse.issparse(aligned) else np.asarray(aligned)
            )
            mapped.append(layer_name)
            mapped_for_assay += 1
        if not mapped_for_assay and not layers:
            unsupported.append(f"Non-active assay {assay_name!r} has no readable layers")
    return mapped, unsupported


class SeuratImporter:
    """Import tested Seurat RDS/H5Seurat objects into canonical AnnData."""

    _RDS_SUFFIX = ".rds"
    _H5SEURAT_SUFFIX = ".h5seurat"

    def probe(self, path: Path) -> ProbeResult:
        suffix = path.suffix.lower()
        try:
            with path.open("rb") as handle:
                signature = handle.read(8)
        except OSError:
            return ProbeResult(ObjectFormat.UNKNOWN, 0.0, "The file cannot be read")
        if suffix == self._H5SEURAT_SUFFIX and signature == b"\x89HDF\r\n\x1a\n":
            return ProbeResult(ObjectFormat.SEURAT, 1.0, "H5Seurat extension and HDF5 signature")
        if suffix == self._RDS_SUFFIX:
            compressed = signature.startswith((b"\x1f\x8b", b"BZh", b"\xfd7zXZ"))
            reason = "RDS extension" + (" and recognized compression" if compressed else "")
            return ProbeResult(ObjectFormat.SEURAT, 0.8 if compressed else 0.7, reason)
        return ProbeResult(ObjectFormat.UNKNOWN, 0.0, "Not a recognized Seurat filename")

    def load(self, path: Path, options: ImportOptions | None = None) -> ImportResult:
        del options
        if not path.is_file():
            raise SeuratImportError(f"Input file does not exist: {path.name}")

        report_warnings: list[str] = []
        unsupported: list[str] = []
        active_assay = "active assay"
        source_object: Any | None = None
        try:
            if path.suffix.lower() == self._H5SEURAT_SUFFIX:
                adata = readseurat.read_h5seurat(str(path))
            else:
                features, active_assay, source_object = _recover_rds_metadata(path)
                with python_warnings.catch_warnings():
                    python_warnings.simplefilter("ignore")
                    adata = readseurat.read_seurat(str(path))
                if len(features) != adata.n_vars or not all(features):
                    raise SeuratImportError(
                        "The active assay feature map is missing or does not match the matrix"
                    )
                if tuple(str(name) for name in adata.var_names) != features:
                    adata.var_names = list(features)
                    report_warnings.append(
                        "Recovered Seurat v5 feature identifiers from the Assay5 coordinate map"
                    )

                mapped_layers, unsupported_assays = _map_secondary_assays(
                    source_object, active_assay, adata
                )
                unsupported.extend(unsupported_assays)
                if mapped_layers:
                    report_warnings.append(
                        f"Mapped {len(mapped_layers)} compatible secondary-assay layer(s)"
                    )

            validate_adata(adata)
            report_warnings.extend(make_identifiers_unique(adata))
            validate_adata(adata)
        except (AnnDataValidationError, SeuratImportError):
            raise
        except Exception as exc:
            raise SeuratImportError(
                "The file could not be imported as a supported Seurat object. ClustBuster "
                "currently supports in-memory Seurat v5 RDS objects and standard H5Seurat "
                f"files ({type(exc).__name__})."
            ) from exc

        report = ImportReport(
            source_format=ObjectFormat.SEURAT,
            source_filename=path.name,
            cell_count=adata.n_obs,
            feature_count=adata.n_vars,
            sparse=sparse.issparse(adata.X),
            layers=tuple(str(name) for name in adata.layers),
            embeddings=compatible_embeddings(adata),
            candidate_cluster_columns=candidate_cluster_columns(adata.obs),
            expression_sources=expression_sources(adata),
            warnings=report_warnings,
            mapping_decisions={
                "X": f"Seurat {active_assay} primary matrix -> adata.X",
                "obs": "Seurat meta.data -> adata.obs",
                "var": f"Seurat {active_assay} feature metadata -> adata.var",
                "reductions": "Seurat reductions -> adata.obsm['X_<name>']",
                "layers": "Compatible active and namespaced secondary-assay layers -> adata.layers",
            },
            unsupported_components=unsupported,
        )
        return ImportResult(adata=adata, report=report)
