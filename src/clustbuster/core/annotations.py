"""Independent, session-scoped cluster annotation state."""

from __future__ import annotations

import json
import math
from collections.abc import Iterable
from dataclasses import asdict, dataclass, replace
from datetime import UTC, datetime
from typing import Any


def _utc_now() -> str:
    return datetime.now(UTC).isoformat()


@dataclass(frozen=True, slots=True)
class ClusterIdentifier:
    """A type-tagged cluster ID that keeps numeric and text labels distinct."""

    type_name: str
    value: str

    @classmethod
    def from_value(cls, value: Any) -> ClusterIdentifier:
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return cls("missing", "null")
        if isinstance(value, bool):
            return cls("bool", json.dumps(value))
        if isinstance(value, int):
            return cls("int", str(value))
        if isinstance(value, float):
            return cls("float", repr(value))
        if isinstance(value, str):
            return cls("str", value)
        return cls(type(value).__name__, str(value))

    @property
    def serialized(self) -> str:
        return json.dumps({"type": self.type_name, "value": self.value}, sort_keys=True)

    @property
    def display(self) -> str:
        return "<missing>" if self.type_name == "missing" else self.value


@dataclass(frozen=True, slots=True)
class AnnotationRecord:
    cluster_id: ClusterIdentifier
    annotation: str
    source: str
    updated_at: str
    notes: str | None = None
    confidence: str | float | None = None
    reference_id: str | None = None

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["cluster_id"] = self.cluster_id.serialized
        result["cluster_display"] = self.cluster_id.display
        return result


class AnnotationStore:
    """Mutable annotation state that never mutates source observations."""

    def __init__(self, records: Iterable[AnnotationRecord]) -> None:
        self._records = {record.cluster_id.serialized: record for record in records}

    @classmethod
    def from_clusters(cls, clusters: Iterable[Any]) -> AnnotationStore:
        records: list[AnnotationRecord] = []
        seen: set[str] = set()
        timestamp = _utc_now()
        for raw_cluster in clusters:
            cluster_id = ClusterIdentifier.from_value(raw_cluster)
            if cluster_id.serialized in seen:
                continue
            seen.add(cluster_id.serialized)
            records.append(
                AnnotationRecord(
                    cluster_id=cluster_id,
                    annotation=cluster_id.display,
                    source="imported",
                    updated_at=timestamp,
                )
            )
        return cls(records)

    def __len__(self) -> int:
        return len(self._records)

    def records(self) -> tuple[AnnotationRecord, ...]:
        return tuple(self._records.values())

    def get(self, cluster: Any) -> AnnotationRecord:
        key = ClusterIdentifier.from_value(cluster).serialized
        try:
            return self._records[key]
        except KeyError as exc:
            raise KeyError(f"Unknown cluster identifier: {cluster!r}") from exc

    def get_serialized(self, serialized_cluster: str) -> AnnotationRecord:
        try:
            return self._records[serialized_cluster]
        except KeyError as exc:
            raise KeyError(f"Unknown serialized cluster identifier: {serialized_cluster}") from exc

    def assign(
        self,
        cluster: Any,
        annotation: str,
        *,
        source: str = "manual",
        notes: str | None = None,
        confidence: str | float | None = None,
        reference_id: str | None = None,
    ) -> AnnotationRecord:
        if not annotation.strip():
            raise ValueError("Annotation must not be empty")
        existing = self.get(cluster)
        updated = replace(
            existing,
            annotation=annotation.strip(),
            source=source,
            notes=notes,
            confidence=confidence,
            reference_id=reference_id,
            updated_at=_utc_now(),
        )
        self._records[existing.cluster_id.serialized] = updated
        return updated

    def assign_serialized(
        self, serialized_cluster: str, annotation: str, **metadata: Any
    ) -> AnnotationRecord:
        existing = self.get_serialized(serialized_cluster)
        if not annotation.strip():
            raise ValueError("Annotation must not be empty")
        updated = replace(
            existing,
            annotation=annotation.strip(),
            source=str(metadata.get("source", "manual")),
            notes=metadata.get("notes"),
            confidence=metadata.get("confidence"),
            reference_id=metadata.get("reference_id"),
            updated_at=_utc_now(),
        )
        self._records[serialized_cluster] = updated
        return updated

    def assign_many(self, clusters: Iterable[Any], annotation: str, **metadata: Any) -> None:
        cluster_values = list(clusters)
        keys = [ClusterIdentifier.from_value(cluster).serialized for cluster in cluster_values]
        missing = [key for key in keys if key not in self._records]
        if missing:
            raise KeyError("One or more cluster identifiers are unknown")
        for cluster in cluster_values:
            self.assign(cluster, annotation, **metadata)

    def reset(self, clusters: Iterable[Any] | None = None) -> None:
        targets = (
            list(clusters)
            if clusters is not None
            else [record.cluster_id for record in self.records()]
        )
        for target in targets:
            cluster_id = (
                target
                if isinstance(target, ClusterIdentifier)
                else ClusterIdentifier.from_value(target)
            )
            existing = self._records.get(cluster_id.serialized)
            if existing is None:
                raise KeyError(f"Unknown cluster identifier: {target!r}")
            self._records[cluster_id.serialized] = AnnotationRecord(
                cluster_id=cluster_id,
                annotation=cluster_id.display,
                source="imported",
                updated_at=_utc_now(),
            )

    def annotation_for(self, cluster: Any) -> str:
        return self.get(cluster).annotation

    def materialize(self, clusters: Iterable[Any]) -> list[str]:
        return [self.annotation_for(cluster) for cluster in clusters]
