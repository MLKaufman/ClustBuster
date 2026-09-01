"""Independent, session-scoped cluster annotation state."""

from __future__ import annotations

import json
import math
import re
from collections.abc import Iterable
from dataclasses import asdict, dataclass, replace
from datetime import UTC, datetime
from typing import Any


def _utc_now() -> str:
    return datetime.now(UTC).isoformat()


def _natural_cluster_key(identifier: ClusterIdentifier) -> tuple[Any, ...]:
    """Sort cluster displays naturally while keeping typed IDs deterministic."""

    if identifier.type_name == "missing":
        return (1, (), identifier.type_name, identifier.value)
    parts = tuple(
        (0, int(part)) if part.isdigit() else (1, part.casefold())
        for part in re.split(r"(\d+)", identifier.display)
        if part
    )
    return (0, parts, identifier.type_name, identifier.value)


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


@dataclass(frozen=True, slots=True)
class PreviewedAnnotation:
    """A non-mutating cluster-label prediction ready for explicit application."""

    cluster_id: str
    cluster_display: str
    annotation: str
    confidence: float
    margin: float


class AnnotationStore:
    """Mutable annotation state that never mutates source observations."""

    def __init__(
        self, records: Iterable[AnnotationRecord], *, history_limit: int = 50
    ) -> None:
        if history_limit < 1:
            raise ValueError("Annotation history limit must be positive")
        self._records = {record.cluster_id.serialized: record for record in records}
        self._history_limit = history_limit
        self._undo_history: list[dict[str, AnnotationRecord]] = []
        self._redo_history: list[dict[str, AnnotationRecord]] = []

    def _checkpoint(self) -> None:
        self._undo_history.append(self._records.copy())
        if len(self._undo_history) > self._history_limit:
            self._undo_history.pop(0)
        self._redo_history.clear()

    @property
    def can_undo(self) -> bool:
        return bool(self._undo_history)

    @property
    def can_redo(self) -> bool:
        return bool(self._redo_history)

    @classmethod
    def from_clusters(
        cls, clusters: Iterable[Any], *, history_limit: int = 50
    ) -> AnnotationStore:
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
        return cls(records, history_limit=history_limit)

    def __len__(self) -> int:
        return len(self._records)

    def records(self) -> tuple[AnnotationRecord, ...]:
        return tuple(
            sorted(
                self._records.values(),
                key=lambda record: _natural_cluster_key(record.cluster_id),
            )
        )

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
        self._checkpoint()
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
        self._checkpoint()
        self._records[serialized_cluster] = updated
        return updated

    def assign_many(self, clusters: Iterable[Any], annotation: str, **metadata: Any) -> None:
        if not annotation.strip():
            raise ValueError("Annotation must not be empty")
        cluster_values = list(clusters)
        keys = [ClusterIdentifier.from_value(cluster).serialized for cluster in cluster_values]
        missing = [key for key in keys if key not in self._records]
        if missing:
            raise KeyError("One or more cluster identifiers are unknown")
        if not keys:
            return
        timestamp = _utc_now()
        self._checkpoint()
        for key in dict.fromkeys(keys):
            existing = self._records[key]
            self._records[key] = replace(
                existing,
                annotation=annotation.strip(),
                source=str(metadata.get("source", "manual")),
                notes=metadata.get("notes", existing.notes),
                confidence=metadata.get("confidence"),
                reference_id=metadata.get("reference_id"),
                updated_at=timestamp,
            )

    def set_notes_serialized(
        self, serialized_cluster: str, notes: str | None
    ) -> AnnotationRecord:
        """Update notes without changing label or prediction provenance."""

        existing = self.get_serialized(serialized_cluster)
        cleaned = notes.strip() if notes and notes.strip() else None
        if cleaned == existing.notes:
            return existing
        updated = replace(existing, notes=cleaned, updated_at=_utc_now())
        self._checkpoint()
        self._records[serialized_cluster] = updated
        return updated

    def apply_table_edits(
        self, edits: Iterable[tuple[str, str, str | None]]
    ) -> tuple[AnnotationRecord, ...]:
        """Atomically apply editable-table labels and notes as one history step."""

        items = tuple(edits)
        keys = [serialized for serialized, _, _ in items]
        if len(keys) != len(set(keys)):
            raise ValueError("Annotation table edits must contain each cluster once")
        if any(not annotation.strip() for _, annotation, _ in items):
            raise ValueError("Annotation must not be empty")
        if any(key not in self._records for key in keys):
            raise KeyError("One or more annotation table clusters are unknown")

        changed: list[tuple[AnnotationRecord, str, str | None]] = []
        for serialized, annotation, notes in items:
            existing = self._records[serialized]
            cleaned_annotation = annotation.strip()
            cleaned_notes = notes.strip() if notes and notes.strip() else None
            if cleaned_annotation != existing.annotation or cleaned_notes != existing.notes:
                changed.append((existing, cleaned_annotation, cleaned_notes))
        if not changed:
            return ()

        timestamp = _utc_now()
        self._checkpoint()
        updated: list[AnnotationRecord] = []
        for existing, annotation, notes in changed:
            annotation_changed = annotation != existing.annotation
            record = replace(
                existing,
                annotation=annotation,
                notes=notes,
                source="manual" if annotation_changed else existing.source,
                confidence=None if annotation_changed else existing.confidence,
                reference_id=None if annotation_changed else existing.reference_id,
                updated_at=timestamp,
            )
            self._records[record.cluster_id.serialized] = record
            updated.append(record)
        return tuple(updated)

    def apply_previewed(
        self,
        predictions: Iterable[PreviewedAnnotation],
        *,
        source: str,
        reference_id: str,
    ) -> tuple[AnnotationRecord, ...]:
        """Atomically apply a validated prediction set to annotation state."""

        items = tuple(predictions)
        keys = [item.cluster_id for item in items]
        if len(keys) != len(set(keys)):
            raise ValueError("A prediction set must contain each cluster at most once")
        if any(not item.annotation.strip() for item in items):
            raise ValueError("Predicted annotations must not be empty")
        if any(key not in self._records for key in keys):
            raise KeyError("One or more predicted cluster identifiers are unknown")
        if not items:
            return ()

        timestamp = _utc_now()
        updated = tuple(
            replace(
                self._records[item.cluster_id],
                annotation=item.annotation.strip(),
                source=source,
                confidence=item.confidence,
                reference_id=reference_id,
                updated_at=timestamp,
            )
            for item in items
        )
        self._checkpoint()
        for record in updated:
            self._records[record.cluster_id.serialized] = record
        return updated

    def reset(self, clusters: Iterable[Any] | None = None) -> None:
        targets = (
            list(clusters)
            if clusters is not None
            else [record.cluster_id for record in self.records()]
        )
        resolved: list[ClusterIdentifier] = []
        for target in targets:
            cluster_id = (
                target
                if isinstance(target, ClusterIdentifier)
                else ClusterIdentifier.from_value(target)
            )
            existing = self._records.get(cluster_id.serialized)
            if existing is None:
                raise KeyError(f"Unknown cluster identifier: {target!r}")
            resolved.append(cluster_id)
        if not resolved:
            return
        self._checkpoint()
        timestamp = _utc_now()
        for cluster_id in resolved:
            self._records[cluster_id.serialized] = AnnotationRecord(
                cluster_id=cluster_id,
                annotation=cluster_id.display,
                source="imported",
                updated_at=timestamp,
            )

    def undo(self) -> tuple[AnnotationRecord, ...]:
        if not self._undo_history:
            raise ValueError("There are no annotation changes to undo")
        self._redo_history.append(self._records.copy())
        self._records = self._undo_history.pop()
        return self.records()

    def redo(self) -> tuple[AnnotationRecord, ...]:
        if not self._redo_history:
            raise ValueError("There are no annotation changes to redo")
        self._undo_history.append(self._records.copy())
        self._records = self._redo_history.pop()
        return self.records()

    def annotation_for(self, cluster: Any) -> str:
        return self.get(cluster).annotation

    def materialize(self, clusters: Iterable[Any]) -> list[str]:
        return [self.annotation_for(cluster) for cluster in clusters]
