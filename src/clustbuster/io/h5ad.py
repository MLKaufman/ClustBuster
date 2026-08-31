"""Reference H5AD importer."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
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


class H5adImportError(ValueError):
    """An actionable error raised at the H5AD compatibility boundary."""


class H5adImporter:
    def probe(self, path: Path) -> ProbeResult:
        suffix_match = path.suffix.lower() == ".h5ad"
        try:
            with path.open("rb") as handle:
                hdf5_signature = handle.read(8) == b"\x89HDF\r\n\x1a\n"
        except OSError:
            return ProbeResult(ObjectFormat.UNKNOWN, 0.0, "The file cannot be read")
        if suffix_match and hdf5_signature:
            return ProbeResult(ObjectFormat.H5AD, 1.0, "H5AD extension and HDF5 signature")
        if hdf5_signature:
            return ProbeResult(ObjectFormat.H5AD, 0.6, "HDF5 signature; extension is advisory")
        return ProbeResult(ObjectFormat.UNKNOWN, 0.0, "Not an HDF5 file")

    def load(self, path: Path, options: ImportOptions | None = None) -> ImportResult:
        del options  # Selection is applied after the import summary is reviewed.
        if not path.is_file():
            raise H5adImportError(f"Input file does not exist: {path.name}")
        try:
            adata = ad.read_h5ad(path)
            validate_adata(adata)
            warnings = make_identifiers_unique(adata)
            validate_adata(adata)
        except AnnDataValidationError:
            raise
        except Exception as exc:
            raise H5adImportError(
                "The file could not be parsed as H5AD. It may be truncated, corrupt, or use an "
                f"unsupported AnnData encoding ({type(exc).__name__})."
            ) from exc

        report = ImportReport(
            source_format=ObjectFormat.H5AD,
            source_filename=path.name,
            cell_count=adata.n_obs,
            feature_count=adata.n_vars,
            sparse=sparse.issparse(adata.X),
            layers=tuple(str(name) for name in adata.layers),
            embeddings=compatible_embeddings(adata),
            candidate_cluster_columns=candidate_cluster_columns(adata.obs),
            expression_sources=expression_sources(adata),
            warnings=warnings,
            mapping_decisions={"X": "adata.X", "obs": "adata.obs", "var": "adata.var"},
        )
        return ImportResult(adata=adata, report=report)

