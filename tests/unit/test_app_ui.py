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
    assert 'id="dot_plot"' in html and "height:1100px" in html
    assert 'id="marker_heatmap"' in html and "height:900px" in html


def test_annotation_sidebar_uses_only_the_autosaving_table_editor() -> None:
    html = _rendered_html()

    assert 'id="annotation_cluster"' not in html
    assert 'id="annotation_label"' not in html
    assert 'id="save_annotation"' not in html
    assert 'id="save_annotation_table"' not in html
    assert 'id="reset_selected_annotation"' not in html
    assert 'id="undo_annotation"' not in html
    assert 'id="redo_annotation"' not in html
    assert 'id="request_reset_annotations"' in html
    assert 'id="download_annotation_csv"' in html
    assert "Changes save automatically when you leave a field" in html
    assert html.index('id="download_annotation_csv"') < html.index(
        "Changes save automatically when you leave a field"
    )
    assert ".cb-cluster-summary { background:#fff" in html


def test_workspace_and_tab_labels_match_product_language() -> None:
    html = _rendered_html()

    assert 'id="initialize_workspace_control"' in html
    assert 'id="feature_plot_container"' in html
    assert 'data-value="MarkerCodex"' in html
    assert 'data-value="Refmats"' in html
    assert 'data-value="Top Markers"' in html
    assert 'data-value="Resources"' not in html
    assert 'data-value="Reference annotation"' not in html
    assert 'data-value="Markers &amp; heatmap"' not in html
    assert html.index('data-value="Module scores"') < html.index('data-value="MarkerCodex"')
    assert html.index('data-value="Top Markers"') < html.index('data-value="MarkerCodex"')
