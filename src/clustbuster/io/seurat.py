"""Pure-Python Seurat import boundary backed by a pinned readseurat release."""

from __future__ import annotations

import gc
import warnings as python_warnings
from pathlib import Path
from typing import Any

import anndata as ad
import h5py  # type: ignore[import-untyped]
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


def _inspect_h5seurat(path: Path) -> tuple[str, list[str]]:
    """Read lightweight assay metadata without materializing expression matrices."""

    try:
        with h5py.File(path, "r") as source:
            active_value = source.attrs.get("active.assay")
            if isinstance(active_value, bytes):
                active_assay = active_value.decode("utf-8")
            else:
                active_assay = str(active_value) if active_value is not None else ""
            assays = source.get("assays")
            assay_names = [str(name) for name in assays] if assays is not None else []
    except (OSError, UnicodeDecodeError) as exc:
        raise SeuratImportError("The H5Seurat header could not be read") from exc

    if not active_assay:
        raise SeuratImportError("The H5Seurat file does not declare an active assay")
    if active_assay not in assay_names:
        raise SeuratImportError(
            f"The H5Seurat active assay {active_assay!r} is not present in the assay table"
        )
    return active_assay, [name for name in assay_names if name != active_assay]


def _has_canonical_seurat_slots(source_data: dict[str, Any]) -> bool:
    """Recognize a Seurat S4 payload when rdata omits its top-level class attribute."""

    required_slots = {
        "active.assay",
        "active.ident",
        "assays",
        "meta.data",
        "project.name",
        "reductions",
        "version",
    }
    if not required_slots.issubset(source_data):
        return False
    active_assay = _first_string(source_data.get("active.assay"))
    assays = source_data.get("assays")
    if not active_assay or not isinstance(assays, dict) or active_assay not in assays:
        return False
    assay_data = getattr(assays[active_assay], "__dict__", {})
    assay_class = _first_string(assay_data.get("class"))
    return assay_class in {"Assay", "Assay5"}


def _read_rds_source(path: Path) -> tuple[str, Any, bool]:
    """Decode the source once for strict class and assay inspection."""

    try:
        from readseurat.rdata import read_rds

        with python_warnings.catch_warnings():
            python_warnings.simplefilter("ignore")
        source = read_rds(str(path))
        source_data = source.__dict__
        source_class = _first_string(source_data.get("class"))
        inferred_class = source_class is None and _has_canonical_seurat_slots(source_data)
        if source_class != "Seurat" and not inferred_class:
            raise SeuratImportError(
                f"The RDS contains {source_class or 'an unknown R object'}, not a Seurat object. "
                "Use the SingleCellExperiment importer for SCE files or convert the object to "
                "H5AD in R."
            )
        active_assay = _first_string(source_data.get("active.assay"))
        assays = source_data.get("assays")
        if not active_assay or not isinstance(assays, dict) or active_assay not in assays:
            raise SeuratImportError("The Seurat object has no readable active assay")
        return active_assay, source, inferred_class
    except SeuratImportError:
        raise
    except Exception as exc:
        raise SeuratImportError(
            "The RDS could not be decoded as a supported in-memory Seurat v4/v5 object. "
            "On-disk layers, custom S4 extensions, and non-Seurat RDS files are not supported "
            f"({type(exc).__name__})."
        ) from exc


def _legacy_assay_features(assay: Any) -> tuple[str, ...]:
    metadata = assay.__dict__.get("meta.features")
    index = getattr(metadata, "index", None)
    return tuple(str(value) for value in index) if index is not None else ()


def _convert_legacy_seurat(source: Any, active_assay: str) -> tuple[ad.AnnData, list[str]]:
    """Convert the tested Seurat v4 Assay representation without guessing identifiers."""

    source_data = source.__dict__
    assays = source_data.get("assays", {})
    assay = assays[active_assay]
    assay_data = assay.__dict__
    obs = source_data.get("meta.data")
    var = assay_data.get("meta.features")
    if obs is None or var is None:
        raise SeuratImportError("The legacy Seurat object is missing meta.data or meta.features")

    cells = tuple(str(value) for value in obs.index)
    features = _legacy_assay_features(assay)
    if not cells or len(cells) != len(set(cells)):
        raise SeuratImportError("The legacy Seurat object has missing or duplicate cell IDs")
    if not features or len(features) != len(set(features)):
        raise SeuratImportError("The legacy active assay has missing or duplicate feature IDs")

    counts = assay_data.get("counts")
    normalized = assay_data.get("data")
    primary = (
        normalized
        if getattr(normalized, "shape", None) == (len(features), len(cells))
        else counts
    )
    if getattr(primary, "shape", None) != (len(features), len(cells)):
        raise SeuratImportError("The legacy active assay matrix does not match its identifiers")

    adata = ad.AnnData(X=primary.T, obs=obs.copy(), var=var.copy())
    if getattr(counts, "shape", None) == (len(features), len(cells)):
        adata.layers["counts"] = counts.T

    unsupported: list[str] = []
    target_features = set(features)
    for raw_name, secondary in assays.items():
        assay_name = str(raw_name)
        if assay_name == active_assay:
            continue
        secondary_features = _legacy_assay_features(secondary)
        secondary_data = secondary.__dict__
        if set(secondary_features) != target_features or len(secondary_features) != len(features):
            unsupported.append(
                f"Non-active assay {assay_name!r} has a different cell or feature space"
            )
            continue
        positions = {name: index for index, name in enumerate(secondary_features)}
        feature_order = [positions[name] for name in features]
        for layer_name in ("counts", "data"):
            matrix = secondary_data.get(layer_name)
            if getattr(matrix, "shape", None) != (len(features), len(cells)):
                unsupported.append(
                    f"Non-active assay {assay_name!r} layer {layer_name!r} has an "
                    "incompatible shape"
                )
                continue
            adata.layers[f"assay:{assay_name}:{layer_name}"] = matrix[feature_order, :].T

    reductions = source_data.get("reductions", {})
    cell_positions = {name: index for index, name in enumerate(cells)}
    for raw_name, reduction in reductions.items():
        values = reduction.__dict__.get("cell.embeddings")
        reduction_cells = _coordinate_names(values)
        if set(reduction_cells) != set(cells) or len(reduction_cells) != len(cells):
            unsupported.append(f"Reduction {str(raw_name)!r} does not align to cell IDs")
            continue
        reduction_positions = {name: index for index, name in enumerate(reduction_cells)}
        order = [reduction_positions[name] for name in cell_positions]
        adata.obsm[f"X_{str(raw_name).lower()}"] = np.asarray(values)[order, :]
    return adata, unsupported


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
                active_assay, non_active_assays = _inspect_h5seurat(path)
                adata = readseurat.read_h5seurat(str(path))
                unsupported.extend(
                    f"Non-active H5Seurat assay {name!r} is not imported by the pinned reader"
                    for name in non_active_assays
                )
            else:
                active_assay, source_object, inferred_class = _read_rds_source(path)
                if inferred_class:
                    report_warnings.append(
                        "Recovered the Seurat class from its canonical S4 slot structure"
                    )
                source_assay = source_object.__dict__["assays"][active_assay]
                assay_class = _first_string(source_assay.__dict__.get("class"))
                if assay_class == "Assay":
                    adata, unsupported_assays = _convert_legacy_seurat(
                        source_object, active_assay
                    )
                    unsupported.extend(unsupported_assays)
                    report_warnings.append(
                        "Imported legacy Seurat v4 Assay through the compatibility converter"
                    )
                elif assay_class == "Assay5":
                    features = _coordinate_names(getattr(source_assay, "features", None))
                    if not features:
                        raise SeuratImportError(
                            "The active assay has no feature-name coordinate"
                        )
                    has_secondary_assays = len(source_object.__dict__["assays"]) > 1
                    if not has_secondary_assays:
                        # Avoid retaining one complete decoded R object while readseurat
                        # decodes it again. This is material for real-world multi-GB RDS files.
                        source_assay = None
                        source_object = None
                        gc.collect()
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

                    if source_object is not None:
                        mapped_layers, unsupported_assays = _map_secondary_assays(
                            source_object, active_assay, adata
                        )
                        unsupported.extend(unsupported_assays)
                        if mapped_layers:
                            report_warnings.append(
                                f"Mapped {len(mapped_layers)} compatible secondary-assay layer(s)"
                            )
                else:
                    raise SeuratImportError(
                        f"The active assay class {assay_class or 'unknown'!r} is not supported"
                    )

            validate_adata(adata)
            report_warnings.extend(make_identifiers_unique(adata))
            validate_adata(adata)
        except (AnnDataValidationError, SeuratImportError):
            raise
        except Exception as exc:
            raise SeuratImportError(
                "The file could not be imported as a supported Seurat object. ClustBuster "
                "currently supports tested in-memory Seurat v4/v5 RDS objects and standard "
                "H5Seurat "
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
