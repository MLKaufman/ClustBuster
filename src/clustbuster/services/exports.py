"""CSV and annotated-H5AD export use cases."""

from __future__ import annotations

import hashlib
import re
import zipfile
from datetime import UTC, datetime
from pathlib import Path
from uuid import uuid4

import anndata as ad
import pandas as pd

from clustbuster import __version__
from clustbuster.core.annotations import ClusterIdentifier
from clustbuster.io.export import ExportCollisionError, ExportError
from clustbuster.models import ExportArtifact, ExportResult, Workspace

CLUSTER_COLUMN = "clustbuster_cluster"
ANNOTATION_COLUMN = "clustbuster_annotation"
PROVENANCE_KEY = "clustbuster"


def _checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_stem(filename: str) -> str:
    stem = Path(filename).stem
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "-", stem).strip("-.")
    return cleaned[:80] or "workspace"


class WorkspaceExportService:
    """Generate download artifacts without mutating the source AnnData object."""

    def annotation_tables(self, workspace: Workspace) -> tuple[pd.DataFrame, pd.DataFrame]:
        if workspace.cluster_column is None:
            raise ExportError("Configure a cluster column before exporting")
        records = workspace.annotations.records()
        cluster_rows = [record.to_dict() for record in records]
        cluster_table = pd.DataFrame(cluster_rows)[
            [
                "cluster_id",
                "cluster_display",
                "annotation",
                "notes",
                "source",
                "confidence",
                "reference_id",
                "updated_at",
            ]
        ]

        cell_rows: list[dict[str, object]] = []
        clusters = workspace.adata.obs[workspace.cluster_column].tolist()
        for cell_id, cluster in zip(workspace.adata.obs_names.astype(str), clusters, strict=True):
            record = workspace.annotations.get(cluster)
            cell_rows.append(
                {
                    "cell_id": cell_id,
                    "cluster_id": record.cluster_id.serialized,
                    "cluster_display": record.cluster_id.display,
                    "annotation": record.annotation,
                    "source": record.source,
                    "confidence": record.confidence,
                    "reference_id": record.reference_id,
                }
            )
        cell_table = pd.DataFrame(cell_rows)
        return cluster_table, cell_table

    def export_annotation_zip(self, workspace: Workspace, destination: Path) -> ExportResult:
        destination.mkdir(parents=True, exist_ok=True)
        cluster_table, cell_table = self.annotation_tables(workspace)
        stem = _safe_stem(workspace.source_filename)
        artifact_path = destination / f"{uuid4().hex}-{stem}-annotations.zip"
        with zipfile.ZipFile(artifact_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr(
                f"{stem}-cluster-annotations.csv",
                cluster_table.to_csv(index=False, lineterminator="\n"),
            )
            archive.writestr(
                f"{stem}-cell-annotations.csv",
                cell_table.to_csv(index=False, lineterminator="\n"),
            )
        artifact = ExportArtifact(
            kind="annotation_zip",
            path=artifact_path,
            checksum=_checksum(artifact_path),
            media_type="application/zip",
        )
        return ExportResult((artifact,))

    def export_h5ad(self, workspace: Workspace, destination: Path) -> ExportResult:
        if workspace.cluster_column is None:
            raise ExportError("Configure a cluster column before exporting")
        collisions = {
            name for name in (CLUSTER_COLUMN, ANNOTATION_COLUMN) if name in workspace.adata.obs
        }
        if PROVENANCE_KEY in workspace.adata.uns:
            collisions.add(f"uns:{PROVENANCE_KEY}")
        if collisions:
            raise ExportCollisionError(
                "Existing ClustBuster export fields would be overwritten: "
                + ", ".join(sorted(collisions))
            )

        destination.mkdir(parents=True, exist_ok=True)
        exported = workspace.adata.copy()
        clusters = workspace.adata.obs[workspace.cluster_column].tolist()
        serialized_clusters = [ClusterIdentifier.from_value(value).serialized for value in clusters]
        annotations = workspace.annotations.materialize(clusters)
        exported.obs[CLUSTER_COLUMN] = serialized_clusters
        exported.obs[ANNOTATION_COLUMN] = annotations
        exported.uns[PROVENANCE_KEY] = {
            "clustbuster_version": __version__,
            "annotation_schema_version": "1",
            "source_filename": workspace.source_filename,
            "source_format": workspace.source_format.value,
            "source_cluster_column": workspace.cluster_column,
            "expression_source": workspace.expression_source.label,
            "exported_at": datetime.now(UTC).isoformat(),
        }

        stem = _safe_stem(workspace.source_filename)
        artifact_path = destination / f"{uuid4().hex}-{stem}-annotated.h5ad"
        exported.write_h5ad(artifact_path)
        reopened = ad.read_h5ad(artifact_path)
        self._validate_round_trip(workspace, reopened, serialized_clusters, annotations)
        artifact = ExportArtifact(
            kind="annotated_h5ad",
            path=artifact_path,
            checksum=_checksum(artifact_path),
            media_type="application/x-hdf5",
        )
        return ExportResult((artifact,))

    @staticmethod
    def _validate_round_trip(
        source: Workspace,
        reopened: ad.AnnData,
        serialized_clusters: list[str],
        annotations: list[str],
    ) -> None:
        if reopened.shape != source.adata.shape:
            raise ExportError("Annotated H5AD failed shape validation after writing")
        if set(reopened.layers) != set(source.adata.layers):
            raise ExportError("Annotated H5AD did not preserve expression layers")
        if set(reopened.obsm) != set(source.adata.obsm):
            raise ExportError("Annotated H5AD did not preserve embeddings")
        if reopened.obs[CLUSTER_COLUMN].astype(str).tolist() != serialized_clusters:
            raise ExportError("Annotated H5AD failed cluster identifier validation")
        if reopened.obs[ANNOTATION_COLUMN].astype(str).tolist() != annotations:
            raise ExportError("Annotated H5AD failed annotation validation")
        if PROVENANCE_KEY not in reopened.uns:
            raise ExportError("Annotated H5AD is missing export provenance")

