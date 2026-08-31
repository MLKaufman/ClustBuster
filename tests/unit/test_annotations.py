import pytest

from clustbuster.core.annotations import (
    AnnotationStore,
    ClusterIdentifier,
    PreviewedAnnotation,
)


def test_cluster_identifiers_preserve_python_scalar_type() -> None:
    numeric = ClusterIdentifier.from_value(1)
    text = ClusterIdentifier.from_value("1")
    assert numeric.serialized != text.serialized


def test_assignment_and_materialization_do_not_collapse_ids() -> None:
    store = AnnotationStore.from_clusters([1, "1", "T cells", 1])
    assert len(store) == 3
    store.assign(1, "Myeloid", notes="review", confidence=0.8)
    store.assign("1", "Lymphoid")
    assert store.materialize([1, "1", "T cells"]) == ["Myeloid", "Lymphoid", "T cells"]
    assert store.get(1).notes == "review"
    assert store.get(1).source == "manual"


def test_empty_annotations_and_unknown_clusters_fail() -> None:
    store = AnnotationStore.from_clusters(["a"])
    with pytest.raises(ValueError, match="must not be empty"):
        store.assign("a", " ")
    with pytest.raises(KeyError, match="Unknown cluster"):
        store.assign("b", "B cell")


def test_reset_restores_imported_label_and_metadata() -> None:
    store = AnnotationStore.from_clusters(["a", "b"])
    store.assign("a", "Astrocyte", reference_id="ref-1")
    store.reset(["a"])
    record = store.get("a")
    assert record.annotation == "a"
    assert record.source == "imported"
    assert record.reference_id is None


def test_assign_many_is_atomic_when_a_cluster_is_unknown() -> None:
    store = AnnotationStore.from_clusters(["a", "b"])
    with pytest.raises(KeyError, match="One or more"):
        store.assign_many(["a", "unknown"], "Changed")
    assert store.annotation_for("a") == "a"


def test_missing_cluster_is_supported() -> None:
    store = AnnotationStore.from_clusters([None, float("nan")])
    assert len(store) == 1
    assert store.annotation_for(None) == "<missing>"


def test_previewed_predictions_apply_atomically_with_provenance() -> None:
    store = AnnotationStore.from_clusters(["a", "b"])
    predictions = (
        PreviewedAnnotation(
            cluster_id=ClusterIdentifier.from_value("a").serialized,
            cluster_display="a",
            annotation="T cell",
            confidence=0.91,
            margin=0.4,
        ),
        PreviewedAnnotation(
            cluster_id=ClusterIdentifier.from_value("missing").serialized,
            cluster_display="missing",
            annotation="B cell",
            confidence=0.8,
            margin=0.3,
        ),
    )
    with pytest.raises(KeyError, match="predicted cluster"):
        store.apply_previewed(predictions, source="pyclustifyr", reference_id="ref-v1")
    assert store.annotation_for("a") == "a"

    applied = store.apply_previewed(
        predictions[:1], source="pyclustifyr", reference_id="ref-v1"
    )
    assert applied[0].annotation == "T cell"
    assert applied[0].source == "pyclustifyr"
    assert applied[0].reference_id == "ref-v1"
    assert applied[0].confidence == pytest.approx(0.91)


def test_note_updates_preserve_annotation_provenance() -> None:
    store = AnnotationStore.from_clusters(["a"])
    store.assign(
        "a",
        "T cell",
        source="pyclustifyr",
        confidence=0.92,
        reference_id="pbmc@1.0",
    )
    cluster_id = ClusterIdentifier.from_value("a").serialized

    updated = store.set_notes_serialized(cluster_id, "  Review CD4/CD8 state  ")

    assert updated.notes == "Review CD4/CD8 state"
    assert updated.annotation == "T cell"
    assert updated.source == "pyclustifyr"
    assert updated.confidence == pytest.approx(0.92)
    assert updated.reference_id == "pbmc@1.0"

    assert store.set_notes_serialized(cluster_id, "  ").notes is None


def test_undo_redo_tracks_labels_notes_and_resets_as_atomic_changes() -> None:
    store = AnnotationStore.from_clusters(["a", "b"])
    cluster_a = ClusterIdentifier.from_value("a").serialized
    store.assign("a", "T cell")
    store.set_notes_serialized(cluster_a, "reviewed")
    store.assign_many(["a", "b"], "Lymphoid")
    store.reset()

    assert store.annotation_for("a") == "a"
    assert store.annotation_for("b") == "b"
    assert store.can_undo
    store.undo()
    assert store.materialize(["a", "b"]) == ["Lymphoid", "Lymphoid"]
    assert store.get("a").notes == "reviewed"
    store.undo()
    assert store.annotation_for("a") == "T cell"
    assert store.get("a").notes == "reviewed"
    assert store.can_redo
    store.redo()
    assert store.materialize(["a", "b"]) == ["Lymphoid", "Lymphoid"]

    store.assign("b", "B cell")
    assert not store.can_redo


def test_annotation_history_is_bounded() -> None:
    store = AnnotationStore.from_clusters(["a"], history_limit=2)
    store.assign("a", "first")
    store.assign("a", "second")
    store.assign("a", "third")
    store.undo()
    store.undo()
    assert store.annotation_for("a") == "first"
    with pytest.raises(ValueError, match="no annotation changes"):
        store.undo()


def test_undo_and_redo_fail_cleanly_without_history() -> None:
    store = AnnotationStore.from_clusters(["a"])
    with pytest.raises(ValueError, match="no annotation changes"):
        store.undo()
    with pytest.raises(ValueError, match="no annotation changes"):
        store.redo()


def test_table_edits_are_atomic_and_note_only_edits_keep_provenance() -> None:
    store = AnnotationStore.from_clusters(["a", "b"])
    store.assign(
        "a",
        "T cell",
        source="pyclustifyr",
        confidence=0.9,
        reference_id="pbmc@1",
    )
    records = store.records()
    updated = store.apply_table_edits(
        [
            (records[0].cluster_id.serialized, "T cell", "review subtype"),
            (records[1].cluster_id.serialized, "B cell", "manual call"),
        ]
    )

    assert len(updated) == 2
    assert store.get("a").source == "pyclustifyr"
    assert store.get("a").confidence == pytest.approx(0.9)
    assert store.get("a").notes == "review subtype"
    assert store.get("b").source == "manual"
    assert store.get("b").reference_id is None

    store.undo()
    assert store.get("a").notes is None
    assert store.annotation_for("b") == "b"

    with pytest.raises(ValueError, match="must not be empty"):
        store.apply_table_edits(
            [(records[0].cluster_id.serialized, " ", None)]
        )
