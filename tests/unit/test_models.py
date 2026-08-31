import pytest

from clustbuster.models import ExpressionKind, ExpressionSource


def test_expression_sources_are_unambiguous() -> None:
    assert ExpressionSource.x().label == "X"
    assert ExpressionSource.raw().kind is ExpressionKind.RAW
    assert ExpressionSource.named_layer("counts").label == "layer:counts"


def test_layer_name_is_required_only_for_layer_source() -> None:
    with pytest.raises(ValueError, match="requires a layer name"):
        ExpressionSource(ExpressionKind.LAYER)
    with pytest.raises(ValueError, match="Only a layer"):
        ExpressionSource(ExpressionKind.X, "counts")


def test_expression_source_labels_round_trip() -> None:
    sources = (
        ExpressionSource.x(),
        ExpressionSource.raw(),
        ExpressionSource.named_layer("log1p"),
    )
    for source in sources:
        assert ExpressionSource.from_label(source.label) == source
    with pytest.raises(ValueError, match="Unknown expression source"):
        ExpressionSource.from_label("counts")
