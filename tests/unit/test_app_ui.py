from clustbuster import __version__
from clustbuster.app import app_ui


def _rendered_html() -> str:
    return str(app_ui.render()["html"])


def test_workspace_sidebar_owns_branding_version_and_embedding_color_control() -> None:
    html = _rendered_html()

    assert html.index("cb-workspace-brand") < html.index('id="dataset"')
    assert html.index("ClustBuster logo") < html.index('class="cb-title"')
    assert "ClustBuster logo" in html
    assert f"Version {__version__}" in html
    assert 'id="color_by"' in html
    assert "AnnData MVP" not in html
    assert 'data-value="Import report"' not in html
    assert "Import report" in html
    assert 'class="cb-overview-plot"' in html
    assert 'class="cb-overview-stack"' in html
    assert "height:1160px; min-height:1160px" in html
    assert 'id="embedding_plot"' in html and "height:1160px" in html
    assert "cb-import-report html-fill-container" in html


def test_annotation_sidebar_uses_only_the_autosaving_table_editor() -> None:
    html = _rendered_html()

    assert 'id="annotation_cluster"' not in html
    assert 'id="annotation_label"' not in html
    assert 'id="save_annotation"' not in html
    assert 'id="save_annotation_table"' not in html
    assert 'id="reset_selected_annotation"' not in html
    assert "Changes save automatically when you leave a field" in html
