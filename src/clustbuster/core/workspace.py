"""Workspace construction helpers."""

from __future__ import annotations

from uuid import uuid4

from clustbuster.core.annotations import AnnotationStore
from clustbuster.models import ExpressionSource, ImportResult, Workspace


def workspace_from_import(result: ImportResult) -> Workspace:
    report = result.report
    source = report.expression_sources[0] if report.expression_sources else ExpressionSource.x()
    cluster_column = (
        report.candidate_cluster_columns[0] if report.candidate_cluster_columns else None
    )
    clusters = result.adata.obs[cluster_column].tolist() if cluster_column else []
    return Workspace(
        adata=result.adata,
        source_format=report.source_format,
        source_filename=report.source_filename,
        expression_source=source,
        annotations=AnnotationStore.from_clusters(clusters),
        import_report=report,
        cluster_column=cluster_column,
        embedding_key=report.embeddings[0] if report.embeddings else None,
        workspace_id=uuid4().hex,
    )
