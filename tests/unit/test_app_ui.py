from clustbuster import __version__
from clustbuster.app import app_ui


def _rendered_html() -> str:
    return str(app_ui.render()["html"])


def test_workspace_sidebar_owns_branding_version_and_embedding_color_control() -> None:
    html = _rendered_html()

    assert '<span class="sidebar-title">Workspace</span>' not in html
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
    assert "#dataset_progress.shiny-file-input-progress { height:1.5rem" in html
    assert 'setInputValue("plot_refresh"' in html
    assert 'addEventListener("visibilitychange"' in html
    assert 'id="marker_heatmap_container"' in html
    assert "height:900px" not in html


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
    assert 'data-value="ORA"' in html
    assert 'data-value="Enrichment"' not in html
    assert 'data-value="Resources"' not in html
    assert 'data-value="Reference annotation"' not in html
    assert 'data-value="Markers &amp; heatmap"' not in html
    assert html.index('data-value="Module scores"') < html.index('data-value="MarkerCodex"')
    assert html.index('data-value="Top Markers"') < html.index('data-value="MarkerCodex"')


def test_feature_and_dot_plots_update_without_run_buttons() -> None:
    html = _rendered_html()
    default_genes = (
        "PTPRC, CD3E, CD4, CD8A, MS4A1, CD14, NKG7, EPCAM, PECAM1, "
        "COL1A1, ACTA2, MKI67"
    )

    assert 'id="run_feature"' not in html
    assert 'id="run_dotplot"' not in html
    assert html.count(default_genes) == 2
    assert "Plots update automatically when the gene list changes." in html
    assert "The dot plot updates automatically when the gene list changes." in html
    assert 'id="feature_show_annotations"' in html


def test_module_scores_include_per_cell_violin_plot() -> None:
    html = _rendered_html()

    assert 'class="cb-module-stack"' in html
    assert 'id="module_plot"' in html
    assert 'id="module_violin_plot"' in html


def test_top_markers_has_logfc_and_embedded_enrichment_controls() -> None:
    html = _rendered_html()

    marker_top_n = html.index('id="marker_top_n"')
    assert 'value="25"' in html[marker_top_n : marker_top_n + 220]
    assert 'id="marker_min_logfc"' in html
    assert 'id="run_marker_enrichment"' in html
    assert 'id="marker_enrichment_plot"' in html
    assert 'id="marker_enrichment_table_container"' in html
    assert 'class="cb-top-markers-stack"' in html
    assert "#marker_heatmap_container { display:block" in html
    assert html.index('id="marker_heatmap_container"') < html.index('id="run_marker_enrichment"')


def test_ora_tab_runs_all_cluster_markers_and_pathway_heatmap() -> None:
    html = _rendered_html()

    assert 'id="ora_library"' in html
    assert "GO Biological Process 2025" in html
    assert 'id="run_ora"' in html
    assert 'id="ora_heatmap_container"' in html
    assert 'id="ora_marker_table"' in html
    assert 'id="ora_table"' in html
    assert "Ranks positive markers one cluster versus all remaining cells" in html


def test_refmats_uses_scrollable_stack_and_includes_correlation_table() -> None:
    html = _rendered_html()

    assert 'class="cb-refmats-stack"' in html
    assert 'id="reference_correlation_container"' in html
    assert 'id="reference_correlation_table"' in html
    assert 'id="discard_reference_predictions"' not in html
    assert html.index('id="apply_reference_predictions"') < html.index(
        'id="reference_prediction_table"'
    )
    assert html.index('id="reference_correlation_container"') < html.index(
        'id="reference_correlation_table"'
    )
