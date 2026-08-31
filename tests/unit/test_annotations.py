import pytest

from clustbuster.core.annotations import AnnotationStore, ClusterIdentifier


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
