"""Structural validation and discovery for canonical AnnData workspaces."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd
from anndata import AnnData

from clustbuster.models import ExpressionSource


class AnnDataValidationError(ValueError):
    """Raised when an AnnData object cannot safely become a workspace."""


def validate_adata(adata: AnnData) -> None:
    """Validate invariants required by all downstream application layers."""
    if adata.n_obs == 0:
        raise AnnDataValidationError("The object contains no cells")
    if adata.n_vars == 0:
        raise AnnDataValidationError("The object contains no features")
    if adata.obs_names.hasnans or adata.var_names.hasnans:
        raise AnnDataValidationError("Cell and feature identifiers must not be missing")
    for name, matrix in adata.layers.items():
        if matrix.shape != adata.shape:
            raise AnnDataValidationError(
                f"Layer {name!r} has shape {matrix.shape}; expected {adata.shape}"
            )
    for name, embedding in adata.obsm.items():
        if getattr(embedding, "ndim", 0) != 2 or embedding.shape[0] != adata.n_obs:
            raise AnnDataValidationError(
                f"Embedding {name!r} is incompatible with {adata.n_obs} cells"
            )


def make_identifiers_unique(adata: AnnData) -> list[str]:
    """Make identifiers unique and return explicit mapping warnings."""
    warnings: list[str] = []
    if not adata.obs_names.is_unique:
        duplicate_count = int(adata.obs_names.duplicated().sum())
        adata.obs_names_make_unique()
        warnings.append(f"Made {duplicate_count} duplicate cell identifier(s) unique")
    if not adata.var_names.is_unique:
        duplicate_count = int(adata.var_names.duplicated().sum())
        adata.var_names_make_unique()
        warnings.append(f"Made {duplicate_count} duplicate feature identifier(s) unique")
    return warnings


def candidate_cluster_columns(obs: pd.DataFrame) -> tuple[str, ...]:
    """Return plausible clustering columns without selecting one silently."""
    candidates: list[str] = []
    cardinality_limit = max(50, int(len(obs) * 0.05))
    for column in obs.columns:
        series = obs[column]
        unique_count = int(series.nunique(dropna=True))
        if unique_count == 0 or unique_count >= len(obs):
            continue
        is_candidate = (
            isinstance(series.dtype, pd.CategoricalDtype)
            or pd.api.types.is_integer_dtype(series.dtype)
            or unique_count <= cardinality_limit
        )
        if is_candidate:
            candidates.append(str(column))
    return tuple(candidates)


def compatible_embeddings(adata: AnnData) -> tuple[str, ...]:
    compatible = [
        str(name)
        for name, value in adata.obsm.items()
        if getattr(value, "ndim", 0) == 2
        and value.shape[0] == adata.n_obs
        and value.shape[1] >= 2
    ]
    priorities = {"X_umap": 0, "X_tsne": 1}
    return tuple(sorted(compatible, key=lambda name: (priorities.get(name, 2), name)))


def expression_sources(adata: AnnData) -> tuple[ExpressionSource, ...]:
    sources = [ExpressionSource.x()]
    if adata.raw is not None:
        sources.append(ExpressionSource.raw())
    sources.extend(ExpressionSource.named_layer(str(name)) for name in sorted(adata.layers))
    return tuple(sources)


def validate_selected_columns(
    adata: AnnData, cluster_column: str | None, embedding_key: str | None
) -> Sequence[str]:
    errors: list[str] = []
    if cluster_column is not None and cluster_column not in adata.obs:
        errors.append(f"Unknown cluster column: {cluster_column}")
    if embedding_key is not None and embedding_key not in compatible_embeddings(adata):
        errors.append(f"Incompatible embedding: {embedding_key}")
    return errors

