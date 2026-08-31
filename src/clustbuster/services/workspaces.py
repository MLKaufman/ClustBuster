"""Workspace configuration use cases."""

from clustbuster.core.annotations import AnnotationStore
from clustbuster.io.validation import compatible_embeddings
from clustbuster.models import ExpressionSource, Workspace


def configure_workspace(
    workspace: Workspace,
    *,
    cluster_column: str,
    embedding_key: str,
    expression_source: ExpressionSource,
) -> None:
    if cluster_column not in workspace.adata.obs:
        raise ValueError(f"Unknown cluster column: {cluster_column}")
    if embedding_key not in compatible_embeddings(workspace.adata):
        raise ValueError(f"Incompatible embedding: {embedding_key}")
    if expression_source not in workspace.import_report.expression_sources:
        raise ValueError(f"Unavailable expression source: {expression_source.label}")
    workspace.cluster_column = cluster_column
    workspace.embedding_key = embedding_key
    workspace.expression_source = expression_source
    workspace.annotations = AnnotationStore.from_clusters(
        workspace.adata.obs[cluster_column].tolist()
    )

