"""Shiny application entry point for the AnnData-native MVP."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd
from shiny import App, Inputs, Outputs, Session, reactive, render, req, ui
from shinywidgets import output_widget, render_plotly

from clustbuster import __version__
from clustbuster.config import AppConfig
from clustbuster.core.enrichment import AllClusterOraResult, EnrichmentError, EnrichmentResult
from clustbuster.core.expression import (
    DotPlotResult,
    ExpressionResult,
    aggregate_dotplot,
    extract_expression,
    parse_gene_list,
)
from clustbuster.core.markers import MarkerResult, rank_all_markers, rank_markers
from clustbuster.core.modules import MODULE_PRESETS, ModuleScoreResult, calculate_module_score
from clustbuster.core.workspace import workspace_from_import
from clustbuster.integrations.enrichr import DEFAULT_LIBRARY, EnrichrClient
from clustbuster.integrations.pyclustifyr import (
    SUPPORTED_METHODS,
    PyClustifyrAdapter,
    ReferenceAnnotationParameters,
    ReferenceAnnotationResult,
)
from clustbuster.models import (
    CellTypeSummary,
    ExpressionSource,
    LoadedReference,
    MarkerSet,
    ReferenceFilters,
    Workspace,
)
from clustbuster.plotting.dotplot import dotplot_figure
from clustbuster.plotting.embedding import embedding_figure
from clustbuster.plotting.enrichment import enrichment_figure
from clustbuster.plotting.feature import feature_gene_figure
from clustbuster.plotting.heatmap import marker_heatmap_figure, marker_heatmap_height
from clustbuster.plotting.module import module_score_figure, module_score_violin_figure
from clustbuster.plotting.ora import ora_heatmap_figure, ora_heatmap_height
from clustbuster.plotting.reference import (
    reference_correlation_figure,
    reference_correlation_height,
)
from clustbuster.resources.markers.csv import CsvMarkerProvider
from clustbuster.resources.references.local import LocalReferenceProvider
from clustbuster.services.exports import WorkspaceExportService
from clustbuster.services.imports import ImportService, SessionFiles
from clustbuster.services.workspaces import configure_workspace

config = AppConfig.from_env()
logging.basicConfig(level=config.log_level)
logger = logging.getLogger("clustbuster")
_DEFAULT_GENE_PANEL = (
    "PTPRC, CD3E, CD4, CD8A, MS4A1, CD14, NKG7, EPCAM, PECAM1, COL1A1, ACTA2, MKI67"
)
_ORA_LIBRARIES = {
    DEFAULT_LIBRARY: "GO Biological Process 2025",
    "GO_Molecular_Function_2025": "GO Molecular Function 2025",
    "GO_Cellular_Component_2025": "GO Cellular Component 2025",
    "Reactome_Pathways_2024": "Reactome Pathways 2024",
    "KEGG_2021_Human": "KEGG Human 2021",
    "MSigDB_Hallmark_2020": "MSigDB Hallmark 2020",
}
_upload_accept = [".h5ad", "application/x-hdf5"]
if config.enable_seurat_import:
    _upload_accept.extend([".rds", ".h5seurat"])
elif config.enable_sce_import:
    _upload_accept.append(".rds")


def _styles() -> ui.Tag:
    return ui.tags.style(
        """
        :root { --cb-navy: #19324a; --cb-teal: #2d8c88; --cb-bg: #f4f7f8; }
        body { background: var(--cb-bg); color: var(--cb-navy); }
        .cb-workspace-brand { display:block; padding:.15rem 0 1rem; margin-bottom:1rem;
                              border-bottom:1px solid #dbe4e8; text-align:center; }
        .cb-workspace-brand > div { padding-top:.7rem; }
        .cb-logo { display:block; width:100%; height:auto; aspect-ratio:1; object-fit:cover;
                   border-radius:18px; }
        .cb-title { margin:0; font-size:1.55rem; font-weight:700; }
        .cb-subtitle { margin:0; color:#607180; font-size:.88rem; }
        .cb-version { display:block; margin-top:.2rem; color:#607180; font-size:.78rem;
                      font-weight:600; letter-spacing:.03em; }
        .cb-empty { min-height:420px; display:flex; align-items:center; justify-content:center;
                    flex-direction:column; color:#6a7c89; text-align:center; padding:3rem; }
        .cb-summary { display:grid; grid-template-columns:repeat(2, minmax(0,1fr)); gap:.6rem; }
        .cb-stat { background:#fff; border:1px solid #dbe4e8; border-radius:.6rem; padding:.65rem; }
        .cb-stat strong { display:block; font-size:1.15rem; color:var(--cb-navy); }
        .cb-annotation-sidebar { border-left-color:#cbd8dd; }
        .cb-annotation-sidebar .sidebar-content { min-width:310px; }
        .cb-annotation-table { width:100%; border-collapse:separate; border-spacing:0 .35rem; }
        .cb-annotation-table th { font-size:.78rem; color:#607180; padding:0 .3rem; }
        .cb-annotation-table td { padding:0 .2rem; vertical-align:middle; }
        .cb-annotation-table td:first-child { font-weight:600; width:18%; }
        .cb-annotation-table .form-group { margin:0; }
        .cb-annotation-table input { font-size:.84rem; padding:.3rem .45rem; min-width:0; }
        .cb-annotation-actions { display:flex; justify-content:flex-end; gap:.4rem;
                                 flex:0 0 auto; }
        .cb-icon-button { width:2.25rem; height:2.25rem; padding:.25rem;
                          display:inline-flex; align-items:center; justify-content:center;
                          font-size:1rem; line-height:1; }
        .cb-cluster-summary { background:#fff; border:1px solid #dbe4e8;
                              border-radius:.5rem; padding:.75rem; }
        .cb-cluster-summary p { margin-bottom:0; }
        .cb-feature-stack { display:flex; flex-direction:column; gap:1rem; }
        .cb-module-stack { display:flex; flex-direction:column; gap:1rem; width:100%;
                           padding-bottom:1rem; }
        .cb-module-stack > * { flex:0 0 auto !important; margin-bottom:0 !important; }
        .cb-top-markers-stack { display:flex; flex-direction:column; gap:1rem;
                                width:100%; padding-bottom:1rem; }
        .cb-top-markers-stack > * { flex:0 0 auto !important; margin-bottom:0 !important; }
        #marker_heatmap_container { display:block; width:100%; flex:0 0 auto !important; }
        .cb-marker-heatmap-frame { display:block; width:100%; flex:none !important; }
        #marker_enrichment_plot { display:block; width:100%; flex:0 0 560px !important;
                                  min-height:560px; }
        .cb-ora-stack { display:flex; flex-direction:column; gap:1rem; width:100%;
                        padding-bottom:1rem; }
        .cb-ora-stack > * { flex:0 0 auto !important; margin-bottom:0 !important; }
        #ora_heatmap_container { display:block; width:100%; flex:0 0 auto !important; }
        .cb-ora-heatmap-frame { display:block; width:100%; flex:none !important; }
        .cb-refmats-stack { display:flex; flex-direction:column; gap:1rem; width:100%;
                            padding-bottom:1rem; }
        .cb-refmats-stack > * { flex:0 0 auto !important; margin-bottom:0 !important; }
        #reference_correlation_container { display:block; width:100%;
                                           flex:0 0 auto !important; }
        .cb-reference-heatmap-frame { display:block; width:100%; flex:none !important; }
        #dataset_progress.shiny-file-input-progress { height:1.5rem; min-height:1.5rem;
                                                       margin-top:.65rem; margin-bottom:1.25rem;
                                                       border-radius:.5rem; overflow:hidden; }
        #dataset_progress .progress-bar { min-height:1.5rem; line-height:1.5rem; }
        .cb-overview-stack { display:flex; flex-direction:column; gap:0; }
        .cb-overview-plot { width:100%; height:1160px; min-height:1160px; margin:0; }
        #embedding_plot { width:100% !important; height:1160px !important;
                          min-height:1160px !important; margin:0 !important; }
        .cb-import-report { margin-top:0 !important; }
        .btn-primary { background-color:var(--cb-teal); border-color:var(--cb-teal); }
        """
    )


def _plot_refresh_script() -> ui.Tag:
    return ui.tags.script(
        """
        (function () {
          var refreshTimer = null;
          function refreshPlotlyOutputs() {
            if (document.visibilityState !== "visible") return;
            window.clearTimeout(refreshTimer);
            refreshTimer = window.setTimeout(function () {
              if (window.Shiny && typeof window.Shiny.setInputValue === "function") {
                window.Shiny.setInputValue("plot_refresh", Date.now(), {priority: "event"});
              }
              if (window.Plotly) {
                document.querySelectorAll(".js-plotly-plot").forEach(function (plot) {
                  window.Plotly.Plots.resize(plot);
                });
              }
            }, 100);
          }
          document.addEventListener("visibilitychange", refreshPlotlyOutputs);
          window.addEventListener("pageshow", refreshPlotlyOutputs);
          document.addEventListener("shiny:connected", refreshPlotlyOutputs);
          if (window.jQuery) {
            window.jQuery(document).on("shiny:connected", refreshPlotlyOutputs);
          }
        })();
        """
    )


app_ui = ui.page_fillable(
    _styles(),
    _plot_refresh_script(),
    ui.layout_sidebar(
        ui.sidebar(
            ui.div(
                ui.tags.img(
                    src="/assets/clustbuster-logo.png",
                    class_="cb-logo",
                    alt="ClustBuster logo",
                ),
                ui.div(
                    ui.h1("ClustBuster", class_="cb-title"),
                    ui.p("Assisted single-cell cluster annotation", class_="cb-subtitle"),
                    ui.tags.span(f"Version {__version__}", class_="cb-version"),
                ),
                class_="cb-workspace-brand",
            ),
            ui.input_file(
                "dataset",
                "Upload single-cell object",
                accept=_upload_accept,
                button_label="Choose dataset",
                placeholder="No dataset selected",
            ),
            ui.output_ui("import_panel"),
            ui.input_select(
                "color_by",
                "Color embedding by",
                {"cluster": "Source cluster", "annotation": "Current annotation"},
                selected="cluster",
            ),
            ui.output_ui("initialize_workspace_control"),
            width=330,
            open="desktop",
        ),
        ui.navset_card_tab(
            ui.nav_panel(
                "Overview",
                ui.output_ui("overview_header"),
                ui.div(
                    ui.div(
                        output_widget("embedding_plot", width="100%", height="1160px"),
                        class_="cb-overview-plot",
                    ),
                    ui.card(
                        ui.card_header("Import report"),
                        ui.output_ui("report_summary"),
                        ui.output_data_frame("report_table"),
                        class_="cb-import-report",
                        fill=False,
                    ),
                    class_="cb-overview-stack",
                ),
            ),
            ui.nav_panel(
                "Feature plot",
                ui.card(
                    ui.input_text_area(
                        "feature_genes",
                        "Genes",
                        value=_DEFAULT_GENE_PANEL,
                        placeholder="Comma, space, or newline separated",
                        rows=2,
                    ),
                    ui.input_checkbox(
                        "feature_show_annotations",
                        "Show cluster annotations on plots",
                        value=False,
                    ),
                    ui.help_text("Plots update automatically when the gene list changes."),
                    fill=False,
                ),
                ui.output_ui("feature_feedback"),
                ui.output_ui("feature_plot_container"),
            ),
            ui.nav_panel(
                "Dot plot",
                ui.card(
                    ui.input_text_area(
                        "dot_genes",
                        "Gene panel",
                        value=_DEFAULT_GENE_PANEL,
                        placeholder="Comma, space, or newline separated",
                        rows=2,
                    ),
                    ui.help_text("The dot plot updates automatically when the gene list changes."),
                    fill=False,
                ),
                ui.output_ui("dotplot_feedback"),
                output_widget("dot_plot", height="1100px"),
            ),
            ui.nav_panel(
                "Module scores",
                ui.div(
                    ui.card(
                        ui.layout_columns(
                            ui.input_select(
                                "module_preset",
                                "Gene-set preset",
                                {
                                    "custom": "Custom gene set",
                                    **{preset.key: preset.label for preset in MODULE_PRESETS},
                                },
                                selected="t_cell",
                            ),
                            ui.input_text("module_name", "Score label", value="T cell module"),
                            ui.input_action_button(
                                "run_module", "Calculate score", class_="btn-primary"
                            ),
                            col_widths=(4, 5, 3),
                        ),
                        ui.input_text_area(
                            "module_genes",
                            "Genes",
                            value="CD3D, IL7R",
                            placeholder="Comma, space, or newline separated",
                            rows=3,
                        ),
                        ui.help_text(
                            "Scores are the mean of per-gene standardized expression "
                            "for the selected expression source."
                        ),
                        fill=False,
                    ),
                    ui.output_ui("module_feedback"),
                    output_widget("module_plot", height="580px"),
                    output_widget("module_violin_plot", height="620px"),
                    ui.card(
                        ui.card_header("Cluster summary"),
                        ui.output_data_frame("module_summary"),
                        fill=False,
                    ),
                    class_="cb-module-stack",
                ),
            ),
            ui.nav_panel(
                "Top Markers",
                ui.div(
                    ui.card(
                        ui.layout_columns(
                            ui.input_select("marker_cluster", "Selected cluster", {}),
                            ui.input_numeric("marker_top_n", "Top genes", value=25, min=1, max=50),
                            ui.input_numeric(
                                "marker_min_fraction",
                                "Minimum expressing fraction",
                                value=0.1,
                                min=0,
                                max=1,
                                step=0.05,
                            ),
                            ui.input_numeric(
                                "marker_min_logfc",
                                "Minimum logFC",
                                value=0.25,
                                min=0,
                                step=0.05,
                            ),
                            ui.input_action_button(
                                "run_markers", "Rank markers", class_="btn-primary"
                            ),
                            col_widths=(3, 2, 3, 2, 2),
                        ),
                        ui.help_text(
                            "Ranks positive markers for the selected cluster versus all "
                            "remaining cells using Welch's t-test with "
                            "Benjamini-Hochberg correction. "
                            "logFC is the selected-minus-rest mean on the active expression scale."
                        ),
                        fill=False,
                    ),
                    ui.output_ui("marker_feedback"),
                    ui.card(
                        ui.card_header("Ranked markers"),
                        ui.output_data_frame("marker_table"),
                        fill=False,
                    ),
                    ui.output_ui("marker_heatmap_container"),
                    ui.card(
                        ui.card_header("Marker enrichment"),
                        ui.layout_columns(
                            ui.input_select(
                                "marker_enrichment_library",
                                "Gene-set library",
                                {DEFAULT_LIBRARY: "GO Biological Process 2025"},
                            ),
                            ui.input_numeric(
                                "marker_enrichment_top_n",
                                "Terms to display",
                                value=12,
                                min=3,
                                max=30,
                            ),
                            ui.input_action_button(
                                "run_marker_enrichment",
                                "Run enrichment",
                                class_="btn-primary",
                            ),
                            col_widths=(5, 3, 4),
                        ),
                        ui.help_text(
                            "Uses the ranked marker genes above. Running enrichment sends only "
                            "those gene symbols to the Ma'ayan Lab Enrichr service; expression "
                            "values and cell metadata are not transmitted."
                        ),
                        fill=False,
                    ),
                    ui.output_ui("marker_enrichment_feedback"),
                    output_widget("marker_enrichment_plot", height="560px"),
                    ui.card(
                        ui.card_header("Enriched terms"),
                        ui.output_ui("marker_enrichment_table_container"),
                        fill=False,
                    ),
                    class_="cb-top-markers-stack",
                ),
            ),
            ui.nav_panel(
                "MarkerCodex",
                ui.output_ui("provider_status"),
                ui.card(
                    ui.card_header("Marker catalog"),
                    ui.layout_columns(
                        ui.input_text(
                            "marker_query", "Cell-type search", placeholder="e.g. T cell"
                        ),
                        ui.input_select("marker_species", "Species", {"human": "Human"}),
                        ui.input_select("marker_tissue", "Tissue", {"blood": "Blood"}),
                        ui.input_action_button(
                            "search_markers", "Search catalog", class_="btn-primary"
                        ),
                        col_widths=(5, 2, 2, 3),
                    ),
                    ui.output_data_frame("marker_search_table"),
                    ui.layout_columns(
                        ui.input_select("catalog_cell_type", "Cell type", {}),
                        ui.input_action_button("load_marker_set", "Load marker set"),
                        ui.input_action_button("use_markers_feature", "Use in feature plot"),
                        ui.input_action_button("use_markers_dot", "Use in dot plot"),
                        col_widths=(4, 3, 3, 2),
                    ),
                    ui.output_ui("marker_set_feedback"),
                    ui.output_data_frame("marker_set_table"),
                    fill=False,
                ),
                ui.card(
                    ui.card_header("Local reference matrices"),
                    ui.output_ui("reference_controls"),
                    ui.output_ui("reference_feedback"),
                    ui.output_data_frame("reference_table"),
                    fill=False,
                ),
            ),
            ui.nav_panel(
                "Refmats",
                ui.div(
                    ui.card(
                        ui.card_header("Reference-based cluster annotation"),
                        ui.output_ui("reference_annotation_controls"),
                        ui.help_text(
                            "Scoring runs locally with pyclustifyr. Results remain a preview "
                            "until you explicitly apply them to the current annotation state."
                        ),
                        fill=False,
                    ),
                    ui.output_ui("reference_annotation_feedback"),
                    ui.layout_columns(
                        ui.input_action_button(
                            "apply_reference_predictions",
                            "Apply predictions",
                            class_="btn-primary",
                        ),
                        col_widths=(3,),
                    ),
                    ui.card(
                        ui.card_header("Prediction preview"),
                        ui.output_data_frame("reference_prediction_table"),
                        fill=False,
                    ),
                    ui.output_ui("reference_correlation_container"),
                    ui.card(
                        ui.card_header("Correlation matrix"),
                        ui.output_data_frame("reference_correlation_table"),
                        fill=False,
                    ),
                    class_="cb-refmats-stack",
                ),
            ),
            ui.nav_panel(
                "ORA",
                ui.div(
                    ui.card(
                        ui.layout_columns(
                            ui.input_select(
                                "ora_library",
                                "Pathway library",
                                _ORA_LIBRARIES,
                                selected=DEFAULT_LIBRARY,
                            ),
                            ui.input_numeric(
                                "ora_marker_top_n",
                                "Markers per cluster",
                                value=25,
                                min=2,
                                max=100,
                            ),
                            ui.input_numeric(
                                "ora_min_fraction",
                                "Minimum expressing fraction",
                                value=0.1,
                                min=0,
                                max=1,
                                step=0.05,
                            ),
                            ui.input_numeric(
                                "ora_min_logfc",
                                "Minimum logFC",
                                value=0.25,
                                min=0,
                                step=0.05,
                            ),
                            ui.input_numeric(
                                "ora_pathway_top_n",
                                "Pathways to display",
                                value=30,
                                min=5,
                                max=100,
                            ),
                            ui.input_action_button(
                                "run_ora", "Run all-cluster ORA", class_="btn-primary"
                            ),
                            col_widths=(3, 2, 2, 2, 2, 3),
                        ),
                        ui.help_text(
                            "Ranks positive markers one cluster versus all remaining cells, "
                            "then runs over-representation analysis separately for every source "
                            "cluster. Only marker gene symbols are sent to Enrichr."
                        ),
                        fill=False,
                    ),
                    ui.output_ui("ora_feedback"),
                    ui.output_ui("ora_heatmap_container"),
                    ui.card(
                        ui.card_header("All-cluster markers"),
                        ui.output_data_frame("ora_marker_table"),
                        fill=False,
                    ),
                    ui.card(
                        ui.card_header("Cluster pathway results"),
                        ui.output_data_frame("ora_table"),
                        fill=False,
                    ),
                    class_="cb-ora-stack",
                ),
            ),
            ui.nav_panel(
                "Export",
                ui.output_ui("export_status"),
                ui.card(
                    ui.card_header("Download annotations"),
                    ui.p("The ZIP contains separate cluster-level and cell-level CSV files."),
                    ui.download_button(
                        "download_annotations",
                        "Download annotation ZIP",
                        class_="btn-primary",
                    ),
                    fill=False,
                ),
                ui.card(
                    ui.card_header("Download annotated object"),
                    ui.p(
                        "Creates a copy with namespaced annotation columns and provenance, "
                        "then reopens it for validation before download."
                    ),
                    ui.download_button(
                        "download_h5ad", "Download annotated H5AD", class_="btn-primary"
                    ),
                    fill=False,
                ),
            ),
            sidebar=ui.sidebar(
                ui.div(
                    ui.download_button(
                        "download_annotation_csv",
                        "💾",
                        title="Download annotations as CSV",
                        aria_label="Download annotations as CSV",
                        class_="btn btn-outline-primary btn-sm cb-icon-button",
                    ),
                    ui.input_action_button(
                        "request_reset_annotations",
                        "🗑️",
                        title="Reset all annotations",
                        aria_label="Reset all annotations",
                        class_="btn btn-outline-danger btn-sm cb-icon-button",
                    ),
                    class_="cb-annotation-actions",
                ),
                ui.help_text(
                    "Edit labels or notes directly in the table. Changes save automatically "
                    "when you leave a field."
                ),
                ui.output_ui("annotation_table"),
                ui.output_ui("annotation_sidebar_status"),
                title="Annotations",
                position="right",
                open="always",
                width=390,
                class_="cb-annotation-sidebar",
                fillable=True,
            ),
            full_screen=True,
        ),
        fillable=True,
    ),
    title="ClustBuster",
)


def server(input: Inputs, output: Outputs, session: Session) -> None:
    session_files = SessionFiles.create(config.workspace_root)
    import_service = ImportService(
        config.max_upload_mb,
        enable_seurat_import=config.enable_seurat_import,
        enable_sce_import=config.enable_sce_import,
    )
    export_service = WorkspaceExportService()
    workspace = reactive.Value[Workspace | None](None)
    configured = reactive.Value(False)
    error_message = reactive.Value[str | None](None)
    revision = reactive.Value(0)
    feature_result = reactive.Value[ExpressionResult | None](None)
    feature_error = reactive.Value[str | None](None)
    dotplot_result = reactive.Value[DotPlotResult | None](None)
    dotplot_error = reactive.Value[str | None](None)
    module_score_result = reactive.Value[ModuleScoreResult | None](None)
    module_error = reactive.Value[str | None](None)
    marker_result = reactive.Value[MarkerResult | None](None)
    marker_error = reactive.Value[str | None](None)
    enrichment_result = reactive.Value[EnrichmentResult | None](None)
    enrichment_error = reactive.Value[str | None](None)
    ora_result = reactive.Value[AllClusterOraResult | None](None)
    ora_error = reactive.Value[str | None](None)
    enrichment_client = EnrichrClient(timeout_seconds=config.enrichment_timeout_seconds)
    marker_provider = CsvMarkerProvider(config.marker_catalog_path)
    reference_provider = LocalReferenceProvider(config.reference_root)
    reference_adapter = PyClustifyrAdapter()
    marker_search_results = reactive.Value[list[CellTypeSummary]]([])
    loaded_marker_set = reactive.Value[MarkerSet | None](None)
    loaded_reference = reactive.Value[LoadedReference | None](None)
    marker_resource_error = reactive.Value[str | None](None)
    reference_resource_error = reactive.Value[str | None](None)
    reference_annotation_result = reactive.Value[ReferenceAnnotationResult | None](None)
    reference_annotation_error = reactive.Value[str | None](None)
    reference_annotation_applied = reactive.Value(False)

    session.on_ended(session_files.cleanup)

    def configure_current_workspace(
        current: Workspace,
        *,
        cluster_column: str,
        embedding_key: str,
        expression_source: ExpressionSource,
    ) -> None:
        configure_workspace(
            current,
            cluster_column=cluster_column,
            embedding_key=embedding_key,
            expression_source=expression_source,
        )
        cluster_choices = {
            record.cluster_id.serialized: record.cluster_id.display
            for record in current.annotations.records()
        }
        ui.update_select("marker_cluster", choices=cluster_choices, session=session)
        configured.set(True)
        feature_result.set(None)
        dotplot_result.set(None)
        module_score_result.set(None)
        marker_result.set(None)
        enrichment_result.set(None)
        ora_result.set(None)
        ora_error.set(None)
        reference_annotation_result.set(None)
        reference_annotation_error.set(None)
        reference_annotation_applied.set(False)
        revision.set(revision.get() + 1)

    @reactive.effect
    @reactive.event(input.dataset)
    def import_dataset() -> None:
        upload_value = input.dataset()
        if not upload_value:
            return
        error_message.set(None)
        configured.set(False)
        with ui.Progress(min=0, max=1, session=session) as progress:
            progress.set(0.15, message="Validating object", detail="Reading object metadata")
            try:
                uploaded = cast(dict[str, Any], upload_value[0])
                result = import_service.import_upload(uploaded, session_files)
                current = workspace_from_import(result)
                report = current.import_report
                if report.candidate_cluster_columns and report.embeddings:
                    configure_current_workspace(
                        current,
                        cluster_column=report.candidate_cluster_columns[0],
                        embedding_key=report.embeddings[0],
                        expression_source=report.expression_sources[0],
                    )
                workspace.set(current)
                progress.set(1, message="Import complete")
                ui.notification_show(
                    f"Loaded {result.report.cell_count:,} cells and "
                    f"{result.report.feature_count:,} features",
                    type="message",
                    session=session,
                )
            except Exception as exc:
                logger.exception("Object import failed", extra={"error_type": type(exc).__name__})
                workspace.set(None)
                error_message.set(str(exc))
                ui.notification_show(str(exc), type="error", duration=10, session=session)

    @output
    @render.ui
    def import_panel() -> ui.TagChild:
        error = error_message.get()
        current = workspace.get()
        if error:
            return ui.div(ui.strong("Import failed"), ui.p(error), class_="alert alert-danger")
        if current is None:
            return ui.p("Upload an H5AD or supported Seurat file to begin.", class_="text-muted")
        report = current.import_report
        cluster_choices = {name: name for name in report.candidate_cluster_columns}
        embedding_choices = {name: name for name in report.embeddings}
        expression_choices = {source.label: source.label for source in report.expression_sources}
        if not cluster_choices or not embedding_choices:
            return ui.div(
                ui.strong("Configuration unavailable"),
                ui.p("The object needs a compatible cluster column and 2D embedding."),
                class_="alert alert-warning",
            )
        return ui.div(
            ui.div(
                ui.div(ui.strong(f"{report.cell_count:,}"), "Cells", class_="cb-stat"),
                ui.div(ui.strong(f"{report.feature_count:,}"), "Features", class_="cb-stat"),
                class_="cb-summary mb-3",
            ),
            ui.input_select(
                "cluster_column",
                "Cluster column",
                cluster_choices,
                selected=report.candidate_cluster_columns[0],
            ),
            ui.input_select(
                "embedding_key",
                "Embedding",
                embedding_choices,
                selected=report.embeddings[0],
            ),
            ui.input_select(
                "expression_source",
                "Expression source",
                expression_choices,
                selected=report.expression_sources[0].label,
            ),
        )

    @output
    @render.ui
    def initialize_workspace_control() -> ui.TagChild:
        current = workspace.get()
        if current is None:
            return ui.div()
        report = current.import_report
        if not report.candidate_cluster_columns or not report.embeddings:
            return ui.div()
        return ui.input_action_button(
            "configure", "Reinitialize Workspace", class_="btn-primary w-100 mt-3"
        )

    @reactive.effect
    @reactive.event(input.configure)
    def apply_configuration() -> None:
        current = workspace.get()
        req(current is not None)
        assert current is not None
        try:
            configure_current_workspace(
                current,
                cluster_column=str(input.cluster_column()),
                embedding_key=str(input.embedding_key()),
                expression_source=ExpressionSource.from_label(str(input.expression_source())),
            )
            ui.notification_show("Workspace reinitialized", type="message", session=session)
        except Exception as exc:
            ui.notification_show(str(exc), type="error", duration=8, session=session)

    @output
    @render.ui
    def overview_header() -> ui.TagChild:
        current = workspace.get()
        if current is None:
            return ui.div(
                ui.h3("Start with a single-cell workspace"),
                ui.p(
                    "Upload a dataset, review its discovered fields, and configure the workspace."
                ),
                class_="cb-empty",
            )
        if not configured.get():
            return ui.div(
                ui.h3("Review the detected workspace settings"),
                ui.p("Choose the cluster column, embedding, and expression source in the sidebar."),
                class_="cb-empty",
            )
        return ui.div()

    @output
    @render.ui
    def annotation_sidebar_status() -> ui.TagChild:
        current = workspace.get()
        if current is None:
            return ui.div(
                ui.strong("No workspace"),
                ui.p("Upload and configure a supported object to edit annotations."),
                class_="alert alert-light",
            )
        if not configured.get():
            return ui.div(
                ui.strong("Workspace not configured"),
                ui.p("Choose the source cluster column in the left sidebar."),
                class_="alert alert-warning",
            )
        revision.get()
        return ui.div(
            ui.strong(f"{len(current.annotations):,} source clusters"),
            ui.p("Table edits save automatically and apply to every cell in the cluster."),
            class_="cb-cluster-summary",
        )

    @output
    @render.ui
    def annotation_table() -> ui.TagChild:
        current = workspace.get()
        if current is None or not configured.get():
            return ui.p("No annotation rows are available.", class_="text-muted")
        rows = [
            ui.tags.tr(
                ui.tags.td(record.cluster_id.display),
                ui.tags.td(
                    ui.input_text(
                        f"annotation_cell_{index}",
                        "",
                        value=record.annotation,
                        update_on="blur",
                    )
                ),
                ui.tags.td(
                    ui.input_text(
                        f"notes_cell_{index}",
                        "",
                        value=record.notes or "",
                        update_on="blur",
                    )
                ),
            )
            for index, record in enumerate(current.annotations.records())
        ]
        return ui.div(
            ui.tags.table(
                ui.tags.thead(
                    ui.tags.tr(
                        ui.tags.th("Cluster"),
                        ui.tags.th("Annotation"),
                        ui.tags.th("Notes"),
                    )
                ),
                ui.tags.tbody(*rows),
                class_="cb-annotation-table",
            ),
            class_="overflow-auto",
        )

    def _refresh_annotation_table(current: Workspace) -> None:
        for index, record in enumerate(current.annotations.records()):
            ui.update_text(f"annotation_cell_{index}", value=record.annotation, session=session)
            ui.update_text(f"notes_cell_{index}", value=record.notes or "", session=session)

    @reactive.effect
    def autosave_annotation_table() -> None:
        current = workspace.get()
        req(current is not None and configured.get())
        assert current is not None
        edits: list[tuple[str, str, str]] = []
        for index, record in enumerate(current.annotations.records()):
            annotation = input[f"annotation_cell_{index}"]()
            notes = input[f"notes_cell_{index}"]()
            req(annotation is not None and notes is not None)
            edits.append((record.cluster_id.serialized, str(annotation), str(notes)))
        try:
            updated = current.annotations.apply_table_edits(edits)
            if updated:
                with reactive.isolate():
                    revision.set(revision.get() + 1)
        except (KeyError, ValueError) as exc:
            _refresh_annotation_table(current)
            ui.notification_show(str(exc), type="error", duration=8, session=session)

    @reactive.effect
    @reactive.event(input.request_reset_annotations)
    def request_reset_annotations() -> None:
        current = workspace.get()
        req(current is not None and configured.get())
        ui.modal_show(
            ui.modal(
                ui.p(
                    "This will replace every current annotation with its original source "
                    "cluster label and clear all notes."
                ),
                title="Reset all annotations?",
                footer=ui.TagList(
                    ui.modal_button("Cancel"),
                    ui.input_action_button(
                        "confirm_reset_annotations",
                        "Reset annotations",
                        class_="btn-danger",
                    ),
                ),
                easy_close=True,
            ),
            session=session,
        )

    @reactive.effect
    @reactive.event(input.confirm_reset_annotations)
    def reset_all_annotations() -> None:
        current = workspace.get()
        req(current is not None and configured.get())
        assert current is not None
        current.annotations.reset()
        revision.set(revision.get() + 1)
        _refresh_annotation_table(current)
        ui.modal_remove(session=session)
        ui.notification_show("All annotations reset", type="message", session=session)

    @output
    @render_plotly
    def embedding_plot() -> Any:
        input.plot_refresh()
        current = workspace.get()
        req(current is not None and configured.get())
        assert current is not None
        revision.get()
        color_by = str(input.color_by()) if input.color_by() else "cluster"
        return embedding_figure(current, color_by=color_by)

    @output
    @render.ui
    def report_summary() -> ui.TagChild:
        current = workspace.get()
        if current is None:
            return ui.div("No import report is available yet.", class_="cb-empty")
        report = current.import_report
        report_messages = [
            *report.warnings,
            *(f"Unsupported: {item}" for item in report.unsupported_components),
        ]
        warning_text = (
            ui.tags.ul(*(ui.tags.li(item) for item in report_messages))
            if report_messages
            else ui.p("No import warnings.", class_="text-success")
        )
        return ui.div(
            ui.h5(report.source_filename),
            ui.p(
                f"{report.cell_count:,} cells x {report.feature_count:,} features · "
                f"{'sparse' if report.sparse else 'dense'} expression matrix"
            ),
            warning_text,
        )

    @output
    @render.data_frame
    def report_table() -> pd.DataFrame:
        current = workspace.get()
        req(current is not None)
        assert current is not None
        report = current.import_report
        rows = [
            ("Layers", ", ".join(report.layers) or "None"),
            ("Embeddings", ", ".join(report.embeddings) or "None"),
            ("Candidate clusters", ", ".join(report.candidate_cluster_columns) or "None"),
            ("Expression sources", ", ".join(x.label for x in report.expression_sources)),
        ]
        return pd.DataFrame(rows, columns=["Workspace component", "Detected values"])

    registered_feature_plots: set[int] = set()

    def feature_annotation_overlays(
        current: Workspace, coordinates: np.ndarray
    ) -> tuple[tuple[float, float, str], ...]:
        assert current.cluster_column is not None
        groups: dict[str, list[int]] = {}
        labels: dict[str, str] = {}
        for index, raw_cluster in enumerate(current.adata.obs[current.cluster_column].tolist()):
            record = current.annotations.get(raw_cluster)
            cluster_id = record.cluster_id.serialized
            groups.setdefault(cluster_id, []).append(index)
            labels[cluster_id] = record.annotation
        return tuple(
            (
                float(np.median(coordinates[indices, 0])),
                float(np.median(coordinates[indices, 1])),
                labels[cluster_id],
            )
            for cluster_id, indices in groups.items()
            if labels[cluster_id]
        )

    def register_feature_plot(index: int) -> None:
        if index in registered_feature_plots:
            return

        @output(id=f"feature_gene_plot_{index}")
        @render.plot(alt="Feature expression plot")
        def feature_gene_plot() -> Any:
            current = workspace.get()
            result = feature_result.get()
            req(
                current is not None
                and configured.get()
                and result is not None
                and index < len(result.values.columns)
            )
            assert current is not None and result is not None
            assert current.embedding_key is not None
            gene = str(result.values.columns[index])
            coordinates = np.asarray(current.adata.obsm[current.embedding_key])
            annotations: tuple[tuple[float, float, str], ...] = ()
            if bool(input.feature_show_annotations()):
                revision.get()
                annotations = feature_annotation_overlays(current, coordinates)
            return feature_gene_figure(
                coordinates,
                gene,
                result.values[gene].to_numpy(),
                annotations,
            )

        registered_feature_plots.add(index)

    @reactive.effect
    def update_feature_plots() -> None:
        current = workspace.get()
        req(current is not None and configured.get())
        assert current is not None
        feature_error.set(None)
        try:
            genes = parse_gene_list(str(input.feature_genes()), limit=None)
            result = extract_expression(current.adata, current.expression_source, genes)
            for index in range(len(result.values.columns)):
                register_feature_plot(index)
            feature_result.set(result)
        except Exception as exc:
            feature_result.set(None)
            feature_error.set(str(exc))

    @output
    @render.ui
    def feature_feedback() -> ui.TagChild:
        error = feature_error.get()
        result = feature_result.get()
        if error:
            return ui.div(error, class_="alert alert-danger")
        if result is None:
            return ui.p("Configure a workspace and enter one or more genes.")
        missing = ", ".join((*result.report.missing, *result.report.ambiguous))
        if missing:
            return ui.div(f"Not plotted: {missing}", class_="alert alert-warning")
        return ui.div(
            f"Matched {len(result.report.matched)} gene(s).", class_="alert alert-success"
        )

    @output
    @render.ui
    def feature_plot_container() -> ui.TagChild:
        result = feature_result.get()
        if result is None:
            return ui.div()
        return ui.div(
            *(
                ui.card(
                    ui.card_header(str(gene)),
                    ui.output_plot(
                        f"feature_gene_plot_{index}", width="100%", height="520px"
                    ),
                    fill=False,
                )
                for index, gene in enumerate(result.values.columns)
            ),
            class_="cb-feature-stack",
        )

    @reactive.effect
    def update_dot_plot() -> None:
        current = workspace.get()
        req(current is not None and configured.get() and current.cluster_column is not None)
        assert current is not None and current.cluster_column is not None
        dotplot_error.set(None)
        try:
            genes = parse_gene_list(str(input.dot_genes()), limit=24)
            result = aggregate_dotplot(
                current.adata,
                current.expression_source,
                current.cluster_column,
                genes,
            )
            dotplot_result.set(result)
        except Exception as exc:
            dotplot_result.set(None)
            dotplot_error.set(str(exc))

    @output
    @render.ui
    def dotplot_feedback() -> ui.TagChild:
        error = dotplot_error.get()
        result = dotplot_result.get()
        if error:
            return ui.div(error, class_="alert alert-danger")
        if result is None:
            return ui.p("Configure a workspace and enter a gene panel.")
        missing = ", ".join((*result.report.missing, *result.report.ambiguous))
        if missing:
            return ui.div(f"Not plotted: {missing}", class_="alert alert-warning")
        return ui.div(
            f"Aggregated {len(result.report.matched)} gene(s) across "
            f"{result.values['cluster_id'].nunique()} clusters.",
            class_="alert alert-success",
        )

    @output
    @render_plotly
    def dot_plot() -> Any:
        input.plot_refresh()
        result = dotplot_result.get()
        req(result is not None)
        assert result is not None
        return dotplot_figure(result)

    @output
    @render.ui
    def provider_status() -> ui.TagChild:
        marker_status = marker_provider.status()
        reference_status = reference_provider.status()
        status_class = (
            "alert alert-success"
            if marker_status.available and reference_status.available
            else "alert alert-warning"
        )
        return ui.div(
            ui.strong("Resource status"),
            ui.p(f"Markers: {marker_status.message} · References: {reference_status.message}"),
            class_=status_class,
        )

    @reactive.effect
    @reactive.event(input.search_markers)
    def search_marker_catalog() -> None:
        marker_resource_error.set(None)
        loaded_marker_set.set(None)
        try:
            results = marker_provider.search_cell_types(
                str(input.marker_query()),
                species=str(input.marker_species()),
                tissue=str(input.marker_tissue()),
            )
            marker_search_results.set(results)
            ui.update_select(
                "catalog_cell_type",
                choices={item.cell_type: item.cell_type for item in results},
                session=session,
            )
            if not results:
                marker_resource_error.set("No matching cell types were found")
        except Exception as exc:
            marker_search_results.set([])
            marker_resource_error.set(str(exc))

    @output
    @render.data_frame
    def marker_search_table() -> pd.DataFrame:
        results = marker_search_results.get()
        return pd.DataFrame(
            [
                {
                    "Cell type": item.cell_type,
                    "Species": item.species,
                    "Tissue": item.tissue,
                    "Markers": item.marker_count,
                }
                for item in results
            ]
        )

    @reactive.effect
    @reactive.event(input.load_marker_set)
    def load_catalog_marker_set() -> None:
        marker_resource_error.set(None)
        try:
            result = marker_provider.get_markers(
                str(input.catalog_cell_type()),
                species=str(input.marker_species()),
                tissue=str(input.marker_tissue()),
            )
            loaded_marker_set.set(result)
        except Exception as exc:
            loaded_marker_set.set(None)
            marker_resource_error.set(str(exc))

    @output
    @render.ui
    def marker_set_feedback() -> ui.TagChild:
        error = marker_resource_error.get()
        marker_set = loaded_marker_set.get()
        if error:
            return ui.div(error, class_="alert alert-warning")
        if marker_set is None:
            return ui.p("Search for a cell type and load its marker set.")
        positive = sum(record.direction == "positive" for record in marker_set.records)
        negative = len(marker_set.records) - positive
        return ui.div(
            f"Loaded {positive} positive and {negative} negative marker(s) for "
            f"{marker_set.cell_type}.",
            class_="alert alert-success",
        )

    @output
    @render.data_frame
    def marker_set_table() -> pd.DataFrame:
        marker_set = loaded_marker_set.get()
        if marker_set is None:
            return pd.DataFrame()
        return pd.DataFrame(
            [
                {
                    "Gene": record.gene,
                    "Direction": record.direction,
                    "Confidence": record.confidence,
                    "Evidence": record.evidence,
                }
                for record in marker_set.records
            ]
        )

    def _positive_catalog_genes() -> str:
        marker_set = loaded_marker_set.get()
        req(marker_set is not None)
        assert marker_set is not None
        genes = [record.gene for record in marker_set.records if record.direction == "positive"]
        return ", ".join(genes)

    @reactive.effect
    @reactive.event(input.use_markers_feature)
    def use_catalog_markers_in_feature_plot() -> None:
        ui.update_text_area("feature_genes", value=_positive_catalog_genes(), session=session)
        ui.notification_show("Feature-plot genes updated", type="message", session=session)

    @reactive.effect
    @reactive.event(input.use_markers_dot)
    def use_catalog_markers_in_dot_plot() -> None:
        ui.update_text_area("dot_genes", value=_positive_catalog_genes(), session=session)
        ui.notification_show("Dot-plot genes updated", type="message", session=session)

    @output
    @render.ui
    def reference_controls() -> ui.TagChild:
        try:
            references = reference_provider.list_references(ReferenceFilters())
        except Exception as exc:
            return ui.div(str(exc), class_="alert alert-warning")
        return ui.layout_columns(
            ui.input_select(
                "reference_id",
                "Reference",
                {item.reference_id: item.name for item in references},
            ),
            ui.input_action_button("validate_reference", "Validate reference"),
            col_widths=(8, 4),
        )

    @reactive.effect
    @reactive.event(input.validate_reference)
    def validate_local_reference() -> None:
        reference_resource_error.set(None)
        try:
            loaded_reference.set(reference_provider.load_reference(str(input.reference_id())))
        except Exception as exc:
            loaded_reference.set(None)
            reference_resource_error.set(str(exc))

    @output
    @render.ui
    def reference_feedback() -> ui.TagChild:
        error = reference_resource_error.get()
        reference = loaded_reference.get()
        if error:
            return ui.div(error, class_="alert alert-warning")
        if reference is None:
            return ui.p("Choose a reference to validate its schema and checksum.")
        return ui.div(
            f"Validated {reference.summary.name}: {reference.metadata['gene_count']} genes "
            f"across {len(reference.metadata['cell_types'])} cell types.",
            class_="alert alert-success",
        )

    @output
    @render.data_frame
    def reference_table() -> pd.DataFrame:
        references = reference_provider.list_references(ReferenceFilters())
        return pd.DataFrame(
            [
                {
                    "Reference": item.name,
                    "Species": item.species,
                    "Tissue": item.tissue,
                    "Disease/context": item.disease,
                    "Assay": item.assay,
                    "Version": item.resource_version,
                }
                for item in references
            ]
        )

    @output
    @render.ui
    def reference_annotation_controls() -> ui.TagChild:
        try:
            references = reference_provider.list_references(ReferenceFilters())
        except Exception as exc:
            return ui.div(str(exc), class_="alert alert-warning")
        return ui.layout_columns(
            ui.input_select(
                "annotation_reference_id",
                "Reference",
                {item.reference_id: item.name for item in references},
            ),
            ui.input_select(
                "reference_method",
                "Similarity",
                {method: method.title() for method in SUPPORTED_METHODS},
                selected="spearman",
            ),
            ui.input_numeric(
                "reference_min_overlap", "Minimum shared genes", value=5, min=2, step=1
            ),
            ui.input_numeric(
                "reference_threshold",
                "Minimum score",
                value=0,
                min=-1,
                max=1,
                step=0.05,
            ),
            ui.input_checkbox(
                "reference_query_is_log",
                "Query is log-normalized",
                value=True,
            ),
            ui.input_action_button(
                "run_reference_annotation", "Run annotation", class_="btn-primary"
            ),
            col_widths=(3, 2, 2, 2, 3, 3),
        )

    @reactive.effect
    @reactive.event(input.run_reference_annotation)
    def run_reference_annotation() -> None:
        current = workspace.get()
        req(current is not None and configured.get())
        assert current is not None
        reference_annotation_error.set(None)
        reference_annotation_applied.set(False)
        with ui.Progress(min=0, max=1, session=session) as progress:
            progress.set(0.15, message="Validating reference and gene overlap")
            try:
                reference = reference_provider.load_reference(str(input.annotation_reference_id()))
                parameters = ReferenceAnnotationParameters(
                    compute_method=str(input.reference_method()),
                    minimum_gene_overlap=int(input.reference_min_overlap()),
                    threshold=float(input.reference_threshold()),
                    query_is_log_normalized=bool(input.reference_query_is_log()),
                )
                progress.set(0.45, message="Scoring source clusters")
                result = reference_adapter.annotate(current, reference, parameters)
                reference_annotation_result.set(result)
                progress.set(1, message="Prediction preview ready")
                ui.notification_show(
                    f"Previewed predictions for {len(result.predictions)} clusters",
                    type="message",
                    session=session,
                )
            except Exception as exc:
                reference_annotation_result.set(None)
                reference_annotation_error.set(str(exc))
                logger.exception(
                    "Reference annotation failed", extra={"error_type": type(exc).__name__}
                )
                ui.notification_show(str(exc), type="error", duration=10, session=session)

    @output
    @render.ui
    def reference_annotation_feedback() -> ui.TagChild:
        error = reference_annotation_error.get()
        result = reference_annotation_result.get()
        if error:
            return ui.div(error, class_="alert alert-danger")
        if result is None:
            return ui.p("Configure a workspace, select a reference, and run annotation.")
        status = (
            "Predictions applied to current annotations."
            if reference_annotation_applied.get()
            else "Preview only — review the table and correlations before applying."
        )
        warnings = (
            ui.tags.ul(*(ui.tags.li(item) for item in result.warnings))
            if result.warnings
            else ui.p("No scoring warnings.")
        )
        return ui.div(
            ui.strong(
                f"{len(result.matched_genes)} shared genes · pyclustifyr "
                f"{result.package_version} · {result.parameters.compute_method}"
            ),
            ui.p(status),
            warnings,
            class_=(
                "alert alert-success" if reference_annotation_applied.get() else "alert alert-info"
            ),
        )

    @output
    @render.data_frame
    def reference_prediction_table() -> pd.DataFrame:
        result = reference_annotation_result.get()
        if result is None:
            return pd.DataFrame()
        return pd.DataFrame(
            [
                {
                    "Cluster": prediction.cluster_display,
                    "Predicted annotation": prediction.annotation,
                    "Best score": round(prediction.confidence, 4),
                    "Margin": round(prediction.margin, 4),
                }
                for prediction in result.predictions
            ]
        )

    @output
    @render.ui
    def reference_correlation_container() -> ui.TagChild:
        result = reference_annotation_result.get()
        if result is None:
            return ui.div()
        height = reference_correlation_height(len(result.correlations.index))
        return ui.div(
            output_widget("reference_correlation_plot", width="100%", height=f"{height}px"),
            class_="cb-reference-heatmap-frame",
            style=f"height:{height}px; min-height:{height}px;",
        )

    @output
    @render_plotly
    def reference_correlation_plot() -> Any:
        input.plot_refresh()
        result = reference_annotation_result.get()
        req(result is not None)
        assert result is not None
        return reference_correlation_figure(result)

    @output
    @render.data_frame
    def reference_correlation_table() -> pd.DataFrame:
        result = reference_annotation_result.get()
        req(result is not None)
        assert result is not None
        display_by_id = {
            prediction.cluster_id: prediction.cluster_display for prediction in result.predictions
        }
        table = result.correlations.copy().round(4)
        table.insert(
            0,
            "Source cluster",
            [display_by_id.get(str(cluster_id), str(cluster_id)) for cluster_id in table.index],
        )
        return table.reset_index(drop=True)

    @reactive.effect
    @reactive.event(input.apply_reference_predictions)
    def apply_reference_predictions() -> None:
        current = workspace.get()
        result = reference_annotation_result.get()
        req(current is not None and configured.get() and result is not None)
        assert current is not None and result is not None
        version_label = result.reference.summary.resource_version or "unknown"
        reference_label = f"{result.reference.summary.reference_id}@{version_label}"
        try:
            current.annotations.apply_previewed(
                result.predictions,
                source="pyclustifyr",
                reference_id=reference_label,
            )
            reference_annotation_applied.set(True)
            revision.set(revision.get() + 1)
            _refresh_annotation_table(current)
            ui.notification_show(
                f"Applied {len(result.predictions)} reference predictions",
                type="message",
                session=session,
            )
        except Exception as exc:
            ui.notification_show(str(exc), type="error", duration=8, session=session)

    @reactive.effect
    @reactive.event(input.module_preset, ignore_init=True)
    def apply_module_preset() -> None:
        selected = str(input.module_preset())
        if selected == "custom":
            return
        preset = next((item for item in MODULE_PRESETS if item.key == selected), None)
        if preset is None:
            return
        ui.update_text_area("module_genes", value=", ".join(preset.genes), session=session)
        ui.update_text("module_name", value=f"{preset.label} module", session=session)

    @reactive.effect
    @reactive.event(input.run_module)
    def run_module_score() -> None:
        current = workspace.get()
        req(current is not None and configured.get() and current.cluster_column is not None)
        assert current is not None and current.cluster_column is not None
        module_error.set(None)
        try:
            genes = parse_gene_list(str(input.module_genes()), limit=100)
            result = calculate_module_score(
                current.adata,
                current.expression_source,
                current.cluster_column,
                genes,
                name=str(input.module_name()),
            )
            module_score_result.set(result)
        except Exception as exc:
            module_score_result.set(None)
            module_error.set(str(exc))
            ui.notification_show(str(exc), type="error", duration=8, session=session)

    @output
    @render.ui
    def module_feedback() -> ui.TagChild:
        error = module_error.get()
        result = module_score_result.get()
        if error:
            return ui.div(error, class_="alert alert-danger")
        if result is None:
            return ui.p("Configure a workspace, choose a gene set, and calculate its score.")
        unmatched = ", ".join((*result.report.missing, *result.report.ambiguous))
        message = f"Scored {len(result.report.matched)} matched gene(s)."
        if unmatched:
            return ui.div(f"{message} Not included: {unmatched}", class_="alert alert-warning")
        return ui.div(message, class_="alert alert-success")

    @output
    @render_plotly
    def module_plot() -> Any:
        input.plot_refresh()
        current = workspace.get()
        result = module_score_result.get()
        req(current is not None and configured.get() and result is not None)
        assert current is not None and result is not None and current.embedding_key is not None
        coordinates = np.asarray(current.adata.obsm[current.embedding_key])
        return module_score_figure(
            coordinates, current.adata.obs_names.astype(str).tolist(), result
        )

    @output
    @render_plotly
    def module_violin_plot() -> Any:
        input.plot_refresh()
        result = module_score_result.get()
        req(result is not None)
        assert result is not None
        return module_score_violin_figure(result)

    @output
    @render.data_frame
    def module_summary() -> pd.DataFrame:
        result = module_score_result.get()
        req(result is not None)
        assert result is not None
        summary = result.cluster_summary[["cluster", "cells", "mean_score", "median_score"]]
        return summary.rename(
            columns={
                "cluster": "Cluster",
                "cells": "Cells",
                "mean_score": "Mean score",
                "median_score": "Median score",
            }
        )

    @reactive.effect
    @reactive.event(input.run_markers)
    def run_marker_ranking() -> None:
        current = workspace.get()
        req(current is not None and configured.get() and current.cluster_column is not None)
        assert current is not None and current.cluster_column is not None
        marker_error.set(None)
        try:
            result = rank_markers(
                current.adata,
                current.expression_source,
                current.cluster_column,
                str(input.marker_cluster()),
                top_n=int(input.marker_top_n()),
                min_fraction=float(input.marker_min_fraction()),
                min_log_fold_change=float(input.marker_min_logfc()),
            )
            marker_result.set(result)
            enrichment_result.set(None)
        except Exception as exc:
            marker_result.set(None)
            marker_error.set(str(exc))
            ui.notification_show(str(exc), type="error", duration=8, session=session)

    @output
    @render.ui
    def marker_feedback() -> ui.TagChild:
        error = marker_error.get()
        result = marker_result.get()
        if error:
            return ui.div(error, class_="alert alert-danger")
        if result is None:
            return ui.p("Configure a workspace, select a cluster, and rank its markers.")
        return ui.div(
            f"Ranked {len(result.values)} positive marker gene(s) for cluster "
            f"{result.selected_cluster}.",
            class_="alert alert-success",
        )

    @output
    @render.data_frame
    def marker_table() -> pd.DataFrame:
        result = marker_result.get()
        req(result is not None)
        assert result is not None
        table = result.values[
            [
                "rank",
                "gene",
                "score",
                "p_adjusted",
                "log_fold_change",
                "fraction_selected",
                "fraction_rest",
            ]
        ].copy()
        return table.rename(
            columns={
                "rank": "Rank",
                "gene": "Gene",
                "score": "Welch score",
                "p_adjusted": "Adjusted p-value",
                "log_fold_change": "logFC",
                "fraction_selected": "Fraction selected",
                "fraction_rest": "Fraction rest",
            }
        )

    @output
    @render.ui
    def marker_heatmap_container() -> ui.TagChild:
        result = marker_result.get()
        if result is None:
            return ui.div()
        height = marker_heatmap_height(len(result.heatmap.index))
        return ui.div(
            output_widget("marker_heatmap", width="100%", height=f"{height}px"),
            class_="cb-marker-heatmap-frame",
            style=f"height:{height}px; min-height:{height}px;",
        )

    @output
    @render_plotly
    def marker_heatmap() -> Any:
        input.plot_refresh()
        result = marker_result.get()
        req(result is not None)
        assert result is not None
        return marker_heatmap_figure(result)

    @reactive.effect
    @reactive.event(input.run_marker_enrichment)
    def run_go_enrichment() -> None:
        markers = marker_result.get()
        req(markers is not None)
        assert markers is not None
        enrichment_error.set(None)
        genes = tuple(markers.values["gene"].astype(str))
        with ui.Progress(min=0, max=1, session=session) as progress:
            progress.set(0.15, message="Submitting marker genes to Enrichr")
            try:
                result = enrichment_client.enrich(
                    genes,
                    description=f"ClustBuster cluster {markers.selected_cluster} markers",
                )
                enrichment_result.set(result)
                progress.set(1, message="GO enrichment complete")
            except Exception as exc:
                enrichment_result.set(None)
                enrichment_error.set(str(exc))
                ui.notification_show(str(exc), type="error", duration=10, session=session)

    def _enrichment_feedback_content() -> ui.TagChild:
        error = enrichment_error.get()
        result = enrichment_result.get()
        if error:
            return ui.div(error, class_="alert alert-danger")
        if result is None:
            return ui.p("Rank markers first, then run GO Biological Process enrichment.")
        return ui.div(
            f"Returned {len(result.values)} GO terms from {result.source} using "
            f"{len(result.genes)} marker gene(s).",
            class_="alert alert-success",
        )

    @output
    @render.ui
    def marker_enrichment_feedback() -> ui.TagChild:
        return _enrichment_feedback_content()

    @output
    @render_plotly
    def marker_enrichment_plot() -> Any:
        input.plot_refresh()
        result = enrichment_result.get()
        req(result is not None)
        assert result is not None
        return enrichment_figure(result, top_n=int(input.marker_enrichment_top_n()))

    def _enrichment_table_frame() -> pd.DataFrame:
        result = enrichment_result.get()
        req(result is not None)
        assert result is not None
        table = result.values[
            [
                "rank",
                "term",
                "adjusted_p_value",
                "odds_ratio",
                "combined_score",
                "overlap_genes",
            ]
        ].copy()
        return table.rename(
            columns={
                "rank": "Rank",
                "term": "GO Biological Process",
                "adjusted_p_value": "Adjusted p-value",
                "odds_ratio": "Odds ratio",
                "combined_score": "Combined score",
                "overlap_genes": "Overlapping genes",
            }
        )

    @output
    @render.data_frame
    def marker_enrichment_table() -> pd.DataFrame:
        return _enrichment_table_frame()

    @output
    @render.ui
    def marker_enrichment_table_container() -> ui.TagChild:
        result = enrichment_result.get()
        if result is None:
            return ui.p("Run enrichment to display the enriched-term table here.")
        return ui.TagList(
            ui.p(f"Showing all {len(result.values)} returned enriched terms."),
            ui.output_data_frame("marker_enrichment_table"),
        )

    @reactive.effect
    @reactive.event(input.run_ora)
    def run_all_cluster_ora() -> None:
        current = workspace.get()
        req(current is not None and configured.get() and current.cluster_column is not None)
        assert current is not None and current.cluster_column is not None
        ora_error.set(None)
        ora_result.set(None)
        with ui.Progress(min=0, max=1, session=session) as progress:
            try:
                progress.set(0.05, message="Ranking markers for every source cluster")
                markers = rank_all_markers(
                    current.adata,
                    current.expression_source,
                    current.cluster_column,
                    top_n_per_cluster=int(input.ora_marker_top_n()),
                    min_fraction=float(input.ora_min_fraction()),
                    min_log_fold_change=float(input.ora_min_logfc()),
                )
                library = str(input.ora_library())
                client = EnrichrClient(
                    library=library,
                    timeout_seconds=config.enrichment_timeout_seconds,
                )
                grouped = list(markers.values.groupby("cluster_id", sort=False))
                frames: list[pd.DataFrame] = []
                failures = list(markers.failures)
                for index, (cluster_id, cluster_markers) in enumerate(grouped, start=1):
                    cluster_display = str(cluster_markers.iloc[0]["cluster"])
                    progress.set(
                        0.1 + 0.85 * (index - 1) / len(grouped),
                        message=f"Running ORA for cluster {cluster_display}",
                        detail=f"{index} of {len(grouped)} clusters",
                    )
                    genes = tuple(cluster_markers["gene"].astype(str))
                    try:
                        enrichment = client.enrich(
                            genes,
                            description=(
                                f"ClustBuster cluster {cluster_display} all-cluster markers"
                            ),
                        )
                    except Exception as exc:
                        failures.append(f"Cluster {cluster_display}: {exc}")
                        continue
                    values = enrichment.values.copy()
                    values.insert(0, "cluster_id", str(cluster_id))
                    values.insert(1, "cluster", cluster_display)
                    frames.append(values)
                if not frames:
                    raise EnrichmentError(
                        "ORA did not return pathway results for any source cluster"
                    )
                combined = pd.concat(frames, ignore_index=True)
                ora_result.set(
                    AllClusterOraResult(
                        library=library,
                        source="Enrichr",
                        markers=markers.values,
                        values=combined,
                        failures=tuple(failures),
                    )
                )
                progress.set(1, message="All-cluster ORA complete")
            except Exception as exc:
                ora_error.set(str(exc))
                logger.exception("All-cluster ORA failed")
                ui.notification_show(str(exc), type="error", duration=10, session=session)

    @output
    @render.ui
    def ora_feedback() -> ui.TagChild:
        error = ora_error.get()
        result = ora_result.get()
        if error:
            return ui.div(error, class_="alert alert-danger")
        if result is None:
            return ui.p("Configure a workspace, choose a pathway library, and run all-cluster ORA.")
        cluster_count = result.values["cluster_id"].nunique()
        warning_list = (
            ui.tags.ul(*(ui.tags.li(item) for item in result.failures))
            if result.failures
            else ui.p("All source clusters completed successfully.")
        )
        return ui.div(
            ui.strong(
                f"{cluster_count} clusters · {len(result.markers)} markers · "
                f"{result.values['term'].nunique()} pathways"
            ),
            warning_list,
            class_="alert alert-warning" if result.failures else "alert alert-success",
        )

    @output
    @render.ui
    def ora_heatmap_container() -> ui.TagChild:
        result = ora_result.get()
        if result is None:
            return ui.div()
        pathway_count = min(int(input.ora_pathway_top_n()), result.values["term"].nunique())
        height = ora_heatmap_height(pathway_count)
        return ui.div(
            output_widget("ora_heatmap", width="100%", height=f"{height}px"),
            class_="cb-ora-heatmap-frame",
            style=f"height:{height}px; min-height:{height}px;",
        )

    @output
    @render_plotly
    def ora_heatmap() -> Any:
        input.plot_refresh()
        result = ora_result.get()
        req(result is not None)
        assert result is not None
        return ora_heatmap_figure(result, top_n_pathways=int(input.ora_pathway_top_n()))

    @output
    @render.data_frame
    def ora_marker_table() -> pd.DataFrame:
        result = ora_result.get()
        req(result is not None)
        assert result is not None
        return result.markers[
            [
                "cluster",
                "rank",
                "gene",
                "log_fold_change",
                "p_adjusted",
                "fraction_selected",
                "fraction_rest",
            ]
        ].rename(
            columns={
                "cluster": "Source cluster",
                "rank": "Rank",
                "gene": "Gene",
                "log_fold_change": "logFC",
                "p_adjusted": "Adjusted p-value",
                "fraction_selected": "Fraction selected",
                "fraction_rest": "Fraction rest",
            }
        )

    @output
    @render.data_frame
    def ora_table() -> pd.DataFrame:
        result = ora_result.get()
        req(result is not None)
        assert result is not None
        return result.values[
            [
                "cluster",
                "rank",
                "term",
                "adjusted_p_value",
                "odds_ratio",
                "combined_score",
                "overlap_genes",
            ]
        ].rename(
            columns={
                "cluster": "Source cluster",
                "rank": "Rank",
                "term": "Pathway",
                "adjusted_p_value": "Adjusted p-value",
                "odds_ratio": "Odds ratio",
                "combined_score": "Combined score",
                "overlap_genes": "Overlapping genes",
            }
        )

    @output
    @render.ui
    def export_status() -> ui.TagChild:
        current = workspace.get()
        if current is None or not configured.get():
            return ui.div(
                ui.h3("Configure a workspace before exporting"),
                ui.p("Exports are generated only from explicit download actions."),
                class_="cb-empty",
            )
        return ui.div(
            ui.strong("Ready to export"),
            ui.p(
                f"{len(current.annotations):,} cluster annotations and "
                f"{current.adata.n_obs:,} cell annotations will be included."
            ),
            class_="alert alert-success",
        )

    @output
    @render.download_button(
        filename="clustbuster-cluster-annotations.csv",
        media_type="text/csv",
    )
    def download_annotation_csv() -> str:
        current = workspace.get()
        req(current is not None and configured.get())
        assert current is not None
        result = export_service.export_cluster_annotations_csv(current, session_files.exports)
        return str(result.artifacts[0].path)

    @output
    @render.download_button(
        filename="clustbuster-annotations.zip",
        media_type="application/zip",
    )
    def download_annotations() -> str:
        current = workspace.get()
        req(current is not None and configured.get())
        assert current is not None
        with ui.Progress(min=0, max=1, session=session) as progress:
            progress.set(0.2, message="Preparing annotation tables")
            result = export_service.export_annotation_zip(current, session_files.exports)
            progress.set(1, message="Annotation export ready")
        return str(result.artifacts[0].path)

    @output
    @render.download_button(
        filename="clustbuster-annotated.h5ad",
        media_type="application/x-hdf5",
    )
    def download_h5ad() -> str:
        current = workspace.get()
        req(current is not None and configured.get())
        assert current is not None
        with ui.Progress(min=0, max=1, session=session) as progress:
            progress.set(0.1, message="Creating annotated H5AD copy")
            result = export_service.export_h5ad(current, session_files.exports)
            progress.set(1, message="Annotated H5AD validated")
        return str(result.artifacts[0].path)


app = App(app_ui, server, static_assets={"/assets": Path(__file__).parent / "www"})


def main() -> None:
    import subprocess
    import sys

    raise SystemExit(subprocess.call(["shiny", "run", "clustbuster.app:app", *sys.argv[1:]]))
