"""Shiny application entry point for the AnnData-native MVP."""

from __future__ import annotations

import asyncio
import logging
from contextlib import suppress
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd
from shiny import App, Inputs, Outputs, Session, reactive, render, req, ui
from shinywidgets import output_widget, render_plotly

from clustbuster import __version__
from clustbuster.config import AppConfig
from clustbuster.core.enrichment import (
    AllClusterOraResult,
    EnrichmentError,
    EnrichmentProvider,
    EnrichmentResult,
)
from clustbuster.core.expression import (
    DotPlotResult,
    ExpressionResult,
    aggregate_dotplot,
    expression_matrix,
    extract_expression,
    parse_gene_list,
)
from clustbuster.core.markers import AllMarkerResult, MarkerResult, rank_all_markers, rank_markers
from clustbuster.core.modules import ModuleScoreResult, calculate_module_score
from clustbuster.core.workspace import workspace_from_import
from clustbuster.integrations.enrichr import DEFAULT_LIBRARY, EnrichrClient
from clustbuster.integrations.gene_species import MOUSE_LIBRARIES
from clustbuster.integrations.offline_enrichment import GeneSetCache, OfflineEnrichmentClient
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
    MarkerRecord,
    MarkerSet,
    ReferenceFilters,
    Workspace,
)
from clustbuster.plotting.dotplot import dotplot_figure, dotplot_height, dotplot_width
from clustbuster.plotting.embedding import embedding_figure
from clustbuster.plotting.enrichment import enrichment_figure
from clustbuster.plotting.feature import feature_gene_figure
from clustbuster.plotting.heatmap import (
    all_marker_heatmap_figure,
    all_marker_heatmap_height,
    all_marker_heatmap_width,
    marker_heatmap_figure,
    marker_heatmap_height,
)
from clustbuster.plotting.labels import annotation_labels
from clustbuster.plotting.module import (
    module_score_static_figure,
    module_score_violin_figure,
)
from clustbuster.plotting.ora import ora_heatmap_figure, ora_heatmap_height
from clustbuster.plotting.reference import (
    reference_correlation_figure,
    reference_correlation_height,
)
from clustbuster.resources.providers import (
    marker_provider_from_config,
    reference_provider_from_config,
)
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
        .cb-square-module { position:relative; width:100%; aspect-ratio:1 / 1;
                            flex:0 0 auto !important; }
        .cb-square-module > .shiny-ipywidget-output { position:absolute; inset:0; }
        .cb-catalog-plot-stack { display:flex; flex-direction:column; gap:1rem; width:100%; }
        .cb-catalog-plot-stack > * { flex:0 0 auto !important; }
        .cb-catalog-module-frame { position:relative; width:100%; max-width:1000px;
                                   aspect-ratio:1 / 1; }
        .cb-catalog-module-frame > .shiny-plot-output { position:absolute; inset:0; }
        .cb-catalog-module-frame img { width:100% !important; height:100% !important;
                                       object-fit:contain; }
        .cb-feature-stack { display:flex; flex-direction:column; gap:1rem; }
        .cb-module-stack { display:flex; flex-direction:column; gap:1rem; width:100%;
                           padding-bottom:1rem; }
        .cb-module-stack > * { flex:0 0 auto !important; margin-bottom:0 !important; }
        .cb-top-markers-stack { display:flex; flex-direction:column; gap:1rem;
                                width:100%; padding-bottom:1rem; }
        .cb-top-markers-stack > * { flex:0 0 auto !important; margin-bottom:0 !important; }
        .cb-all-markers-stack { display:flex; flex-direction:column; gap:1rem;
                                width:100%; padding-bottom:1rem; }
        .cb-all-markers-stack > * { flex:0 0 auto !important; margin-bottom:0 !important; }
        #dot_plot_container { display:block; width:100%; flex:none !important; }
        #all_marker_heatmap_container { display:block; width:100%;
                                        flex:0 0 auto !important; }
        .cb-all-marker-heatmap-frame { display:block; width:100%; flex:none !important; }
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


def _default_enrichment_mode(cache_root: Path) -> str:
    """Prefer local analysis when the default library has been downloaded."""
    for root in (cache_root, cache_root / "mouse"):
        try:
            GeneSetCache(root).load(DEFAULT_LIBRARY)
        except EnrichmentError:
            continue
        return "offline"
    return "online"


def _square_module_output(output_id: str) -> ui.Tag:
    return ui.div(
        ui.output_plot(output_id, width="100%", height="100%", fill=False),
        class_="cb-catalog-module-frame",
    )


def _scrolling_plot_widget(output_id: str, height: int, width: int) -> ui.Tag:
    return ui.div(
        ui.div(
            output_widget(output_id, width="100%", height=f"{height}px", fill=False),
            style=f"height:{height}px; min-height:{height}px; min-width:{width}px;",
        ),
        style="width:100%; overflow-x:auto; flex:none;",
    )


def _catalog_plot_output(index: int, figure: Any, kind: str) -> ui.Tag:
    if hasattr(figure, "to_plotly_json"):
        # Plot builders grow with cluster count; never crop them to a fixed viewport.
        height = int(figure.layout.height or 700)
        if kind == "dot":
            return _scrolling_plot_widget(
                f"catalog_widget_{index}", height, int(figure.layout.width or 800),
            )
        return ui.div(
            output_widget(f"catalog_widget_{index}", width="100%", height=f"{height}px"),
            style=f"height:{height}px; min-height:{height}px;",
        )
    if kind == "module":
        return ui.div(
            ui.output_plot("catalog_module_plot", width="100%", height="100%", fill=False),
            class_="cb-catalog-module-frame",
        )
    return ui.output_plot(f"catalog_feature_{index}", width="100%", height="520px")


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

            ui.output_ui("initialize_workspace_control"),
            ui.input_switch("show_umap_annotations", "Show UMAP annotations", value=True),
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
                ui.output_ui("dot_plot_container"),
            ),
            ui.nav_panel(
                "Module scores",
                ui.div(
                    ui.card(
                        ui.input_action_button(
                            "run_module", "Calculate score", class_="btn-primary"
                        ),
                        ui.input_text_area(
                            "module_genes",
                            "Genes",
                            value="CD8A, CD8B, CD3D, CD3E, CD3G, TRAC",
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
                    _square_module_output("module_plot"),
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
                                _ORA_LIBRARIES,
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
                            "Uses the ranked marker genes above and the method chosen in Settings. "
                            "Online mode sends gene symbols to Enrichr. Offline mode runs locally "
                            "with a downloaded library."
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
                "All Markers",
                ui.div(
                    ui.card(
                        ui.layout_columns(
                            ui.input_numeric(
                                "all_marker_top_n",
                                "Genes per cluster",
                                value=10,
                                min=1,
                                max=100,
                            ),
                            ui.input_numeric(
                                "all_marker_min_fraction",
                                "Minimum expressing fraction",
                                value=0.1,
                                min=0,
                                max=1,
                                step=0.05,
                            ),
                            ui.input_numeric(
                                "all_marker_min_logfc",
                                "Minimum logFC",
                                value=0.25,
                                min=0,
                                step=0.05,
                            ),
                            ui.input_action_button(
                                "run_all_markers",
                                "Find all markers",
                                class_="btn-primary",
                            ),
                            col_widths=(3, 3, 3, 3),
                        ),
                        ui.help_text(
                            "Ranks every source cluster one versus the remaining cells and "
                            "keeps the requested number of positive marker genes per cluster."
                        ),
                        fill=False,
                    ),
                    ui.output_ui("all_marker_feedback"),
                    ui.output_ui("all_marker_heatmap_container"),
                    ui.card(
                        ui.card_header("All-cluster marker table"),
                        ui.output_data_frame("all_marker_table"),
                        fill=False,
                    ),
                    class_="cb-all-markers-stack",
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
                            "cluster using the method chosen in Settings."
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
                "MarkerCodex",
                ui.output_ui("provider_status"),
                ui.card(
                    ui.card_header("Marker catalog"),
                    ui.output_ui("marker_search_controls"),
                    ui.output_data_frame("marker_search_table"),
                    ui.layout_columns(
                        ui.input_select("catalog_cell_type", "Cell type", {}),
                        ui.input_action_button("load_marker_set", "Load marker set"),
                        col_widths=(8, 4),
                    ),
                    ui.output_ui("marker_set_feedback"),
                    ui.output_ui("marker_set_table"),
                    ui.layout_columns(
                        ui.input_action_button("use_markers_feature", "Plot features"),
                        ui.input_action_button("use_markers_dot", "Plot dot plot"),
                        ui.input_action_button("use_markers_module", "Plot module scores"),
                        ui.input_action_button("use_markers_heatmap", "Plot heatmap"),
                        col_widths=(3, 3, 3, 3),
                    ),
                    ui.output_ui("catalog_plot_feedback"),
                    ui.output_ui("catalog_plot_container"),
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
                ui.card(
                    ui.card_header("Export reference matrix"),
                    ui.p(
                        "Average the active expression source by the current annotations "
                        "or by any cell metadata column. The resulting CSV can be used as "
                        "a Refmat."
                    ),
                    ui.output_ui("reference_matrix_export_controls"),
                    fill=False,
                ),
            ),
            ui.nav_panel(
                "Settings",
                ui.card(
                    ui.card_header("Enrichment analysis"),
                    ui.input_select(
                        "enrichment_mode", "Enrichment method",
                        {"online": "Online — Enrichr", "offline": "Offline — local ORA"},
                        selected=_default_enrichment_mode(config.gene_set_cache_root),
                    ),
                    ui.help_text(
                        "Applies to Top Markers enrichment and all-cluster ORA. "
                        "Offline mode uses a downloaded library and the genes in the active "
                        "expression source as background. It makes no enrichment network requests. "
                        "P-values can differ from online Enrichr. "
                        "Its combined score is unavailable offline. "
                        "Changing mode clears existing enrichment results."
                    ),
                    ui.input_select("offline_species", "Dataset species (offline)",
                                    {
                                        "auto": "Detect automatically",
                                        "human": "Human", "mouse": "Mouse",
                                    }),
                    ui.output_ui("offline_species_status"),
                    ui.input_select("download_species", "Library species to download",
                                    {"human": "Human", "mouse": "Mouse"}),
                    ui.input_select("offline_library", "Library to download", _ORA_LIBRARIES),
                    ui.input_action_button("download_gene_library", "Download for offline use"),
                    ui.help_text(
                        "Download each library once while connected, then use it offline. "
                        "Mouse uses separate MSigDB mouse GO/Reactome libraries. "
                        "Gene symbols are matched exactly; no homolog conversion is used."
                    ),
                    ui.output_ui("offline_library_status"),
                    fill=False,
                ),
                ui.card(
                    ui.card_header("Annotation display"),
            ui.input_select(
                "color_by",
                "Annotations to show",
                {"cluster": "Source cluster", "annotation": "Current annotation"},
                selected="annotation",
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
    configuration_revision = reactive.Value(0)
    feature_result = reactive.Value[ExpressionResult | None](None)
    feature_error = reactive.Value[str | None](None)
    dotplot_result = reactive.Value[DotPlotResult | None](None)
    dotplot_error = reactive.Value[str | None](None)
    module_score_result = reactive.Value[ModuleScoreResult | None](None)
    module_error = reactive.Value[str | None](None)
    marker_result = reactive.Value[MarkerResult | None](None)
    marker_error = reactive.Value[str | None](None)
    all_marker_result = reactive.Value[AllMarkerResult | None](None)
    all_marker_error = reactive.Value[str | None](None)
    enrichment_result = reactive.Value[EnrichmentResult | None](None)
    enrichment_error = reactive.Value[str | None](None)
    ora_result = reactive.Value[AllClusterOraResult | None](None)
    ora_error = reactive.Value[str | None](None)
    gene_set_cache = GeneSetCache(config.gene_set_cache_root)

    enrichment_mode_initialized = False

    @reactive.effect
    @reactive.event(input.enrichment_mode)
    def initialize_enrichment_mode() -> None:
        nonlocal enrichment_mode_initialized
        if enrichment_mode_initialized:
            return
        enrichment_mode_initialized = True
        ui.update_select(
            "enrichment_mode", selected=_default_enrichment_mode(config.gene_set_cache_root),
            session=session,
        )

    gene_set_cache_revision = reactive.Value(0)

    @reactive.extended_task
    async def download_offline_library(library: str, species: str) -> tuple[str, int]:
        count = await asyncio.to_thread(gene_set_cache.download, library, species=species)
        return library, count

    session.on_ended(download_offline_library.cancel)

    @reactive.effect
    @reactive.event(input.download_species)
    def update_download_libraries() -> None:
        choices = {
            key: (label if input.download_species() == "human" else
                  f"{MOUSE_LIBRARIES[key]} (MSigDB 2025.1.Mm)")
            for key, label in _ORA_LIBRARIES.items()
            if input.download_species() == "human" or key in MOUSE_LIBRARIES
        }
        ui.update_select("offline_library", choices=choices, session=session)

    @reactive.effect
    @reactive.event(input.download_gene_library)
    def request_offline_library_download() -> None:
        if download_offline_library.status() != "running":
            download_offline_library(str(input.offline_library()), str(input.download_species()))

    @reactive.effect
    @reactive.event(download_offline_library.status)
    def finish_offline_library_download() -> None:
        if download_offline_library.status() != "success":
            return
        try:
            library, count = download_offline_library.result()
        except Exception:
            return  # The status output displays the download error.
        gene_set_cache_revision.set(gene_set_cache_revision.get() + 1)
        ui.notification_show(f"{library}: {count:,} entries ready offline", session=session)

    @output
    @render.ui
    def offline_library_status() -> ui.TagChild:
        gene_set_cache_revision.get()
        status = download_offline_library.status()
        message: ui.TagChild = ui.div()
        if status == "running":
            message = ui.p("Downloading gene-set library…")
        elif status == "error":
            try:
                download_offline_library.result()
            except Exception as exc:
                message = ui.div(str(exc), class_="alert alert-warning")
        return ui.div(
            message,
            ui.tags.ul(*(ui.tags.li(
                f"{label if input.download_species() == 'human' else MOUSE_LIBRARIES[library]}: "
                + ("Available offline" if GeneSetCache(
                    gene_set_cache.root / ("mouse" if input.download_species() == "mouse" else "")
                ).path(library).is_file()
                               else "Not downloaded")
            ) for library, label in _ORA_LIBRARIES.items()
                if input.download_species() == "human" or library in MOUSE_LIBRARIES)),
        )

    def enrichment_provider(library: str) -> EnrichmentProvider:
        if str(input.enrichment_mode()) == "offline":
            current = workspace.get()
            if current is None:
                raise EnrichmentError("Configure a workspace before running offline enrichment")
            _, names = expression_matrix(current.adata, current.expression_source)
            return OfflineEnrichmentClient(
                gene_set_cache, library, tuple(names.astype(str)), str(input.offline_species()),
            )
        return EnrichrClient(library=library, timeout_seconds=config.enrichment_timeout_seconds)

    @output
    @render.ui
    def offline_species_status() -> ui.TagChild:
        gene_set_cache_revision.get()
        if str(input.enrichment_mode()) != "offline" or workspace.get() is None:
            return ui.p("Load a dataset and select Offline to check species and gene matching.")
        try:
            client = enrichment_provider(str(input.offline_library()))
            assert isinstance(client, OfflineEnrichmentClient)
            matched = len(set().union(*client.terms.values()))
            return ui.p(
                f"Species: {client.species}. {len(client.background):,} unique background genes; "
                f"{matched:,} match the selected library."
            )
        except EnrichmentError as exc:
            return ui.p(str(exc), class_="text-warning")

    @reactive.effect
    @reactive.event(input.enrichment_mode, input.offline_species, gene_set_cache_revision)
    def reset_enrichment_mode() -> None:
        enrichment_result.set(None)
        enrichment_error.set(None)
        ora_result.set(None)
        ora_error.set(None)
    marker_provider = marker_provider_from_config(config)
    reference_provider = reference_provider_from_config(config)
    reference_adapter = PyClustifyrAdapter()
    marker_search_results = reactive.Value[list[CellTypeSummary]]([])
    loaded_marker_set = reactive.Value[MarkerSet | None](None)
    loaded_catalog_selection = reactive.Value[CellTypeSummary | None](None)
    catalog_gene_generation = reactive.Value(0)
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
        feature_error.set(None)
        dotplot_result.set(None)
        dotplot_error.set(None)
        module_score_result.set(None)
        module_error.set(None)
        marker_result.set(None)
        marker_error.set(None)
        all_marker_result.set(None)
        all_marker_error.set(None)
        enrichment_result.set(None)
        enrichment_error.set(None)
        ora_result.set(None)
        ora_error.set(None)
        reference_annotation_result.set(None)
        reference_annotation_error.set(None)
        reference_annotation_applied.set(False)
        configuration_revision.set(configuration_revision.get() + 1)
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
            _refresh_annotation_table(current)
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
        configuration_revision.get()
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
        color_by = str(input.color_by()) if input.color_by() else "annotation"
        return embedding_figure(
            current, color_by=color_by, show_annotations=bool(input.show_umap_annotations()),
        )

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

    def plot_labels() -> dict[str, str]:
        revision.get()
        current = workspace.get()
        return annotation_labels(current, str(input.color_by())) if current is not None else {}

    def plot_display_labels() -> dict[str, str]:
        labels = plot_labels()
        current = workspace.get()
        if current is not None:
            labels.update({record.cluster_id.display: labels[record.cluster_id.serialized]
                           for record in current.annotations.records()})
        return labels

    def feature_annotation_overlays(
        current: Workspace, coordinates: np.ndarray
    ) -> tuple[tuple[float, float, str], ...]:
        if not bool(input.show_umap_annotations()):
            return ()
        assert current.cluster_column is not None
        groups: dict[str, list[int]] = {}
        labels = plot_labels()
        for index, raw_cluster in enumerate(current.adata.obs[current.cluster_column].tolist()):
            record = current.annotations.get(raw_cluster)
            cluster_id = record.cluster_id.serialized
            groups.setdefault(cluster_id, []).append(index)
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
        configuration_revision.get()
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
                    ui.output_plot(f"feature_gene_plot_{index}", width="100%", height="520px"),
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
        configuration_revision.get()
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
    @render.ui
    def dot_plot_container() -> ui.TagChild:
        result = dotplot_result.get()
        if result is None:
            return ui.div()
        return _scrolling_plot_widget(
            "dot_plot", dotplot_height(result.values["cluster_id"].nunique()),
            dotplot_width(result.values["gene"].nunique()),
        )

    @output
    @render_plotly
    def dot_plot() -> Any:
        input.plot_refresh()
        result = dotplot_result.get()
        req(result is not None)
        assert result is not None
        return dotplot_figure(result, plot_labels())

    @output
    @render.ui
    def provider_status() -> ui.TagChild:
        marker_status = marker_provider.status()
        status_class = (
            "alert alert-success"
            if marker_status.available
            else "alert alert-warning"
        )
        return ui.div(
            ui.strong("Resource status"),
            ui.p(
                f"Markers: {marker_status.message}"
                + (f" ({marker_status.version})" if marker_status.version else "")
            ),
            class_=status_class,
        )

    @output
    @render.ui
    def marker_search_controls() -> ui.TagChild:
        try:
            facets = marker_provider.list_facets()
        except Exception as exc:
            return ui.div(str(exc), class_="alert alert-warning")
        species_choices = {"": "All"} | {value: value for value in facets.species}
        tissue_choices = {"": "All"} | {value: value for value in facets.tissues}
        return ui.layout_columns(
            ui.input_text("marker_query", "Search", placeholder="Search all fields"),
            ui.input_select("marker_species", "Species", species_choices),
            ui.input_select("marker_tissue", "Tissue", tissue_choices),
            ui.input_select("marker_submitter", "Submitter", {"": "All"}),
            ui.input_select("marker_collection", "Collection", {"": "All"}),
            ui.input_action_button("search_markers", "Search catalog", class_="btn-primary"),
            col_widths=(4, 2, 2, 2, 2, 12),
        )

    catalog_search_generation = reactive.Value(0)
    catalog_search_pending = reactive.Value(False)

    @reactive.extended_task
    async def query_marker_catalog(
        generation: int, query: str, species: str | None, tissue: str | None,
    ) -> tuple[int, list[CellTypeSummary], str | None]:
        try:
            results = await asyncio.to_thread(
                marker_provider.search_cell_types, query, species=species, tissue=tissue,
            )
            return generation, results, None
        except Exception as exc:
            return generation, [], str(exc)

    session.on_ended(query_marker_catalog.cancel)

    @reactive.effect
    @reactive.event(input.search_markers)
    async def search_marker_catalog() -> None:
        query_marker_catalog.cancel()
        generation = catalog_search_generation.get() + 1
        catalog_search_generation.set(generation)
        catalog_search_pending.set(True)
        marker_resource_error.set(None)
        loaded_marker_set.set(None)
        loaded_catalog_selection.set(None)
        reset_catalog_plots()
        # Reset the browser's selection before replacing the table's data.
        with suppress(Exception):
            await marker_search_table.update_cell_selection(None)
        marker_search_results.set([])
        ui.update_select("catalog_cell_type", choices={}, selected="", session=session)
        query_marker_catalog(
            generation, str(input.marker_query()),
            str(input.marker_species()) or None, str(input.marker_tissue()) or None,
        )

    @reactive.effect
    @reactive.event(query_marker_catalog.result)
    def receive_marker_catalog_search() -> None:
        generation, results, error = query_marker_catalog.result()
        if generation != catalog_search_generation.get():
            return
        marker_search_results.set(results)
        catalog_search_pending.set(False)
        marker_resource_error.set(
            error or (None if results else "No matching cell types were found")
        )
        ui.update_select(
            "catalog_cell_type",
            choices={
                str(index): " · ".join(filter(None, (item.cell_type, item.species, item.tissue)))
                for index, item in enumerate(results)
            },
            selected="0" if results else "",
            session=session,
        )

    @output
    @render.data_frame
    def marker_search_table() -> render.DataGrid[pd.DataFrame]:
        results = marker_search_results.get()
        table = pd.DataFrame(
            [
                {
                    "Cell type": item.cell_type,
                    "Species": item.species,
                    "Tissue": item.tissue,
                    "Markers": item.marker_count,
                }
                for item in results
            ],
            columns=["Cell type", "Species", "Tissue", "Markers"],
        )

        return render.DataGrid(table, selection_mode="row")

    def ensure_catalog_marker_set(selected_index: int | None = None) -> MarkerSet:
        marker_resource_error.set(None)
        try:
            if catalog_search_pending.get():
                raise ValueError("Wait for the catalog search to finish before selecting markers.")
            results = marker_search_results.get()
            if selected_index is None:
                selected_index = int(str(input.catalog_cell_type()))
            if selected_index < 0 or selected_index >= len(results):
                raise ValueError("Choose a marker search result first.")
            selected = results[selected_index]
            existing = loaded_marker_set.get()
            if existing is not None and loaded_catalog_selection.get() == selected:
                return existing
            reset_catalog_plots()
            result = marker_provider.get_markers(
                selected.cell_type, species=selected.species, tissue=selected.tissue,
            )
            catalog_gene_generation.set(catalog_gene_generation.get() + 1)
            loaded_marker_set.set(result)
            loaded_catalog_selection.set(selected)
            return result
        except Exception as exc:
            loaded_marker_set.set(None)
            loaded_catalog_selection.set(None)
            reset_catalog_plots()
            marker_resource_error.set(str(exc) or "Choose a marker search result first.")
            raise

    @reactive.effect
    @reactive.event(input.load_marker_set)
    def load_catalog_marker_set() -> None:
        # The shared loader displays failures beside the marker list.
        with suppress(Exception):
            ensure_catalog_marker_set()

    @reactive.effect
    @reactive.event(input.marker_search_table_cell_selection)
    def load_selected_catalog_row() -> None:
        if catalog_search_pending.get():
            return
        selection = marker_search_table.cell_selection()
        rows = selection["rows"]
        if not rows:
            return
        selected_index = rows[0]
        ui.update_select("catalog_cell_type", selected=str(selected_index), session=session)
        # The shared loader displays failures beside the marker list.
        with suppress(Exception):
            ensure_catalog_marker_set(selected_index)

    @output
    @render.ui
    def marker_set_feedback() -> ui.TagChild:
        error = marker_resource_error.get()
        marker_set = loaded_marker_set.get()
        if error:
            return ui.div(error, class_="alert alert-warning")
        if catalog_search_pending.get():
            return ui.p("Searching marker catalog…")
        if marker_set is None:
            return ui.p("Search for a cell type and load its marker set.")
        positive = sum(record.direction == "positive" for record in marker_set.records)
        negative = len(marker_set.records) - positive
        return ui.div(
            f"Loaded {positive} positive and {negative} negative marker(s) for "
            f"{marker_set.cell_type}.",
            class_="alert alert-success",
        )

    def catalog_gene_input_id(index: int) -> str:
        return f"catalog_gene_{catalog_gene_generation.get()}_{index}"

    def selected_catalog_genes(marker_set: MarkerSet) -> tuple[str, ...]:
        genes = tuple(dict.fromkeys(record.gene for record in marker_set.records))
        return tuple(
            gene for index, gene in enumerate(genes)
            if not input[catalog_gene_input_id(index)].is_set()
            or bool(input[catalog_gene_input_id(index)]())
        )

    @output
    @render.ui
    def marker_set_table() -> ui.TagChild:
        marker_set = loaded_marker_set.get()
        if marker_set is None:
            return ui.div()
        rows = []
        genes = tuple(dict.fromkeys(record.gene for record in marker_set.records))
        def values(records: list[MarkerRecord], field: str) -> str:
            return "; ".join(dict.fromkeys(
                str(value) for record in records
                if (value := getattr(record, field)) is not None and str(value) != ""
            ))

        for index, gene in enumerate(genes):
            records = [record for record in marker_set.records if record.gene == gene]
            confidence = "; ".join(dict.fromkeys(
                str(record.confidence_label or record.confidence)
                for record in records
                if record.confidence_label or record.confidence is not None
            ))
            rows.append(ui.tags.tr(
                ui.tags.td(ui.input_checkbox(catalog_gene_input_id(index), gene, value=True)),
                *(ui.tags.td(value) for value in (
                    values(records, "direction"), confidence, values(records, "verified"),
                    values(records, "evidence"), values(records, "citation"),
                )),
            ))
        return ui.div(
            ui.help_text("Click a gene checkbox to include or exclude it from the next plot. "
                         "All genes start enabled."),
            ui.tags.table(
                ui.tags.thead(ui.tags.tr(*(ui.tags.th(label) for label in (
                    "Plot / Gene", "Direction", "Confidence", "Verified", "Evidence", "Citation",
                )))),
                ui.tags.tbody(*rows),
                class_="table table-striped table-hover",
            ),
            style="max-height:420px; overflow:auto;",
        )

    catalog_plots = reactive.Value[list[Any]]([])
    catalog_plotted_genes = reactive.Value[tuple[str, ...]](())
    catalog_module_result = reactive.Value[ModuleScoreResult | None](None)
    catalog_plot_error = reactive.Value[str | None](None)
    catalog_plot_missing = reactive.Value[tuple[str, ...]](())
    registered_catalog_plots: set[int] = set()

    catalog_plot_kind = reactive.Value("")

    def reset_catalog_plots() -> None:
        catalog_plots.set([])
        catalog_module_result.set(None)
        catalog_plot_error.set(None)
        catalog_plot_missing.set(())
        catalog_plot_kind.set("")

    @reactive.effect
    def clear_catalog_plots() -> None:
        workspace.get()
        configuration_revision.get()
        reset_catalog_plots()

    def register_catalog_plot(index: int) -> None:
        if index in registered_catalog_plots:
            return

        @output(id=f"catalog_feature_{index}")
        @render.plot(alt="Marker catalog feature expression")
        def catalog_feature() -> Any:
            plots = catalog_plots.get()
            req(index < len(plots) and not hasattr(plots[index], "to_plotly_json"))
            return plots[index]

        registered_catalog_plots.add(index)

    def plot_catalog(kind: str, plotted_genes: tuple[str, ...] | None = None) -> None:
        catalog_plots.set([])
        catalog_module_result.set(None)
        catalog_plot_error.set(None)
        catalog_plot_missing.set(())
        try:
            current = workspace.get()
            marker_set = ensure_catalog_marker_set()
            if current is None or not configured.get():
                raise ValueError("Import and configure a workspace before plotting.")
            assert current.cluster_column is not None and current.embedding_key is not None
            genes = (selected_catalog_genes(marker_set)
                     if plotted_genes is None else plotted_genes)
            if not genes:
                raise ValueError("Select at least one gene in the marker list before plotting.")
            coordinates = np.asarray(current.adata.obsm[current.embedding_key])
            plots: list[Any]
            if kind == "feature":
                expression = extract_expression(current.adata, current.expression_source, genes)
                report = expression.report
                plots = [
                    feature_gene_figure(
                        coordinates, str(gene), expression.values[gene].to_numpy(),
                        feature_annotation_overlays(current, coordinates),
                    )
                    for gene in expression.values.columns
                ]
                for index in range(len(plots)):
                    register_catalog_plot(index)
            elif kind == "module":
                module = calculate_module_score(
                    current.adata, current.expression_source, current.cluster_column,
                    genes, name=marker_set.cell_type,
                )
                catalog_module_result.set(module)
                report = module.report
                plots = [
                    module_score_static_figure(
                        coordinates, module, feature_annotation_overlays(current, coordinates),
                    ),
                    module_score_violin_figure(module, plot_labels()),
                ]
            else:
                dots = aggregate_dotplot(
                    current.adata, current.expression_source, current.cluster_column, genes,
                )
                report = dots.report
                if kind == "dot":
                    plots = [dotplot_figure(dots, plot_labels())]
                else:
                    means = dots.values.pivot(
                        index="cluster_id", columns="gene", values="mean_expression"
                    )
                    labels = dots.values.drop_duplicates("cluster_id").set_index("cluster_id")
                    means.index = labels.loc[means.index, "cluster_display"].astype(str)
                    figure = marker_heatmap_figure(
                        MarkerResult(marker_set.cell_type, "catalog", pd.DataFrame(), means),
                        plot_display_labels(),
                    )
                    figure.update_layout(title=f"{marker_set.cell_type} marker expression")
                    plots = [figure]
            catalog_plot_missing.set((*report.missing, *report.ambiguous))
            catalog_plotted_genes.set(genes)
            catalog_plot_kind.set(kind)
            catalog_plots.set(plots)
        except Exception as exc:
            catalog_plot_error.set(str(exc))

    @reactive.effect
    @reactive.event(input.use_markers_feature)
    def plot_catalog_features() -> None:
        plot_catalog("feature")

    @reactive.effect
    @reactive.event(input.use_markers_dot)
    def plot_catalog_dot() -> None:
        plot_catalog("dot")

    @reactive.effect
    @reactive.event(input.use_markers_module)
    def plot_catalog_module() -> None:
        plot_catalog("module")

    @reactive.effect
    @reactive.event(input.color_by, input.show_umap_annotations, revision, ignore_init=True)
    def refresh_catalog_annotation_labels() -> None:
        kind = catalog_plot_kind.get()
        if kind and catalog_plots.get():
            plot_catalog(kind, catalog_plotted_genes.get())

    @reactive.effect
    @reactive.event(input.use_markers_heatmap)
    def plot_catalog_heatmap() -> None:
        plot_catalog("heatmap")

    @output
    @render.ui
    def catalog_plot_feedback() -> ui.TagChild:
        error = catalog_plot_error.get()
        if error:
            return ui.div(error, class_="alert alert-danger")
        missing = catalog_plot_missing.get()
        if missing:
            return ui.div("Not plotted: " + ", ".join(missing), class_="alert alert-warning")
        return ui.div()

    @output
    @render.ui
    def catalog_plot_container() -> ui.TagChild:
        return ui.div(
            *(
                _catalog_plot_output(index, figure, catalog_plot_kind.get())
                for index, figure in enumerate(catalog_plots.get())
            ),
            *(
                [ui.card(ui.card_header("Cluster summary"),
                         ui.output_data_frame("catalog_module_summary"), fill=False)]
                if catalog_plot_kind.get() == "module" and catalog_plots.get() else []
            ),
            class_="cb-catalog-plot-stack",
        )

    @output
    @render.data_frame
    def catalog_module_summary() -> pd.DataFrame:
        result = catalog_module_result.get()
        req(result is not None and catalog_plot_kind.get() == "module")
        assert result is not None
        summary = result.cluster_summary.copy()
        summary["cluster"] = summary["cluster_id"].map(plot_labels()).fillna(summary["cluster"])
        return summary[["cluster", "cells", "mean_score", "median_score"]].rename(
            columns={"cluster": "Cluster", "cells": "Cells", "mean_score": "Mean score",
                     "median_score": "Median score"}
        )

    @output
    @render.plot(width=1000, height=1000, alt="Marker catalog module-score UMAP")
    def catalog_module_plot() -> Any:
        plots = catalog_plots.get()
        req(catalog_plot_kind.get() == "module" and plots)
        return plots[0]

    @output
    @render_plotly
    def catalog_widget_0() -> Any:
        plots = catalog_plots.get()
        if not plots or not hasattr(plots[0], "to_plotly_json"):
            return None
        return plots[0]

    @output
    @render_plotly
    def catalog_widget_1() -> Any:
        plots = catalog_plots.get()
        if len(plots) < 2 or not hasattr(plots[1], "to_plotly_json"):
            return None
        return plots[1]

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
        return reference_correlation_figure(result, plot_labels())

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
                name="Gene-set module score",
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
    @render.plot(width=1000, height=1000, alt="Module-score UMAP")
    def module_plot() -> Any:
        current = workspace.get()
        result = module_score_result.get()
        req(current is not None and configured.get() and result is not None)
        assert current is not None and result is not None and current.embedding_key is not None
        coordinates = np.asarray(current.adata.obsm[current.embedding_key])
        return module_score_static_figure(
            coordinates, result, feature_annotation_overlays(current, coordinates),
        )

    @output
    @render_plotly
    def module_violin_plot() -> Any:
        input.plot_refresh()
        result = module_score_result.get()
        req(result is not None)
        assert result is not None
        return module_score_violin_figure(result, plot_labels())

    @output
    @render.data_frame
    def module_summary() -> pd.DataFrame:
        result = module_score_result.get()
        req(result is not None)
        assert result is not None
        summary = result.cluster_summary.copy()
        summary["cluster"] = summary["cluster_id"].map(plot_labels()).fillna(summary["cluster"])
        summary = summary[["cluster", "cells", "mean_score", "median_score"]]
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
        return marker_heatmap_figure(result, plot_display_labels())

    @reactive.effect
    @reactive.event(input.run_all_markers)
    def run_all_marker_ranking() -> None:
        current = workspace.get()
        req(current is not None and configured.get() and current.cluster_column is not None)
        assert current is not None and current.cluster_column is not None
        all_marker_error.set(None)
        all_marker_result.set(None)
        with ui.Progress(min=0, max=1, session=session) as progress:
            try:
                progress.set(0.1, message="Ranking markers for every source cluster")
                result = rank_all_markers(
                    current.adata,
                    current.expression_source,
                    current.cluster_column,
                    top_n_per_cluster=int(input.all_marker_top_n()),
                    min_fraction=float(input.all_marker_min_fraction()),
                    min_log_fold_change=float(input.all_marker_min_logfc()),
                )
                all_marker_result.set(result)
                progress.set(1, message="All-cluster marker ranking complete")
            except Exception as exc:
                all_marker_error.set(str(exc))
                logger.exception("All-cluster marker ranking failed")
                ui.notification_show(str(exc), type="error", duration=10, session=session)

    @output
    @render.ui
    def all_marker_feedback() -> ui.TagChild:
        error = all_marker_error.get()
        result = all_marker_result.get()
        if error:
            return ui.div(error, class_="alert alert-danger")
        if result is None:
            return ui.p("Configure a workspace, then find markers for every source cluster.")
        cluster_count = result.values["cluster_id"].nunique()
        warnings = (
            ui.tags.ul(*(ui.tags.li(item) for item in result.failures))
            if result.failures
            else ui.p("All source clusters completed successfully.")
        )
        return ui.div(
            ui.strong(f"{cluster_count} clusters · {len(result.values)} marker records"),
            warnings,
            class_="alert alert-warning" if result.failures else "alert alert-success",
        )

    @output
    @render.ui
    def all_marker_heatmap_container() -> ui.TagChild:
        result = all_marker_result.get()
        if result is None:
            return ui.div()
        gene_count = result.values["gene"].nunique()
        height = all_marker_heatmap_height(gene_count)
        width = all_marker_heatmap_width(result.values["cluster_id"].nunique())
        return ui.div(
            ui.div(
                ui.output_plot("all_marker_heatmap", width="100%", height=f"{height}px",
                               fill=False),
                class_="cb-all-marker-heatmap-frame",
                style=f"height:{height}px; min-height:{height}px; min-width:{width}px;",
            ),
            style="width:100%; overflow-x:auto; flex:none;",
        )

    @output
    @render.plot(alt="All-cluster marker expression heatmap")
    def all_marker_heatmap() -> Any:
        current = workspace.get()
        result = all_marker_result.get()
        req(
            current is not None
            and configured.get()
            and current.cluster_column is not None
            and result is not None
        )
        assert current is not None and current.cluster_column is not None and result is not None
        return all_marker_heatmap_figure(
            current.adata,
            current.expression_source,
            current.cluster_column,
            result,
            plot_labels(),
        )

    @output
    @render.data_frame
    def all_marker_table() -> pd.DataFrame:
        result = all_marker_result.get()
        req(result is not None)
        assert result is not None
        return result.values[
            [
                "cluster",
                "rank",
                "gene",
                "score",
                "p_adjusted",
                "log_fold_change",
                "fraction_selected",
                "fraction_rest",
            ]
        ].rename(
            columns={
                "cluster": "Source cluster",
                "rank": "Rank",
                "gene": "Gene",
                "score": "Welch score",
                "p_adjusted": "Adjusted p-value",
                "log_fold_change": "logFC",
                "fraction_selected": "Fraction selected",
                "fraction_rest": "Fraction rest",
            }
        )

    @reactive.effect
    @reactive.event(input.run_marker_enrichment)
    def run_go_enrichment() -> None:
        markers = marker_result.get()
        req(markers is not None)
        assert markers is not None
        enrichment_error.set(None)
        genes = tuple(markers.values["gene"].astype(str))
        with ui.Progress(min=0, max=1, session=session) as progress:
            progress.set(0.15, message="Running marker enrichment")
            try:
                result = enrichment_provider(str(input.marker_enrichment_library())).enrich(
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
        if table["combined_score"].isna().all():
            table = table.drop(columns="combined_score")
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
                client = enrichment_provider(library)
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
                        source=enrichment.source,
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
                f"{result.source} · {cluster_count} clusters · {len(result.markers)} markers · "
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
        return ora_heatmap_figure(
            result, top_n_pathways=int(input.ora_pathway_top_n()), labels=plot_labels(),
        )

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
        table = result.values[
            [
                "cluster",
                "rank",
                "term",
                "adjusted_p_value",
                "odds_ratio",
                "combined_score",
                "overlap_genes",
            ]
        ].copy()
        if table["combined_score"].isna().all():
            table = table.drop(columns="combined_score")
        return table.rename(
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
    @render.ui
    def reference_matrix_export_controls() -> ui.TagChild:
        current = workspace.get()
        if current is None or not configured.get():
            return ui.p("Configure a workspace to choose a Refmat grouping.", class_="text-muted")
        choices = {"annotations": "Current annotations"}
        choices.update(
            {f"metadata:{column}": f"Metadata: {column}" for column in current.adata.obs.columns}
        )
        return ui.div(
            ui.input_select(
                "reference_matrix_group",
                "Group cells by",
                choices,
                selected="annotations",
            ),
            ui.download_button(
                "download_reference_matrix",
                "Download Refmat CSV",
                class_="btn-primary",
            ),
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
        filename="clustbuster-reference-matrix.csv",
        media_type="text/csv",
    )
    def download_reference_matrix() -> str:
        current = workspace.get()
        req(current is not None and configured.get())
        assert current is not None
        grouping = str(input.reference_matrix_group())
        metadata_column = (
            grouping.removeprefix("metadata:") if grouping.startswith("metadata:") else None
        )
        with ui.Progress(min=0, max=1, session=session) as progress:
            progress.set(0.2, message="Averaging expression by group")
            result = export_service.export_reference_matrix_csv(
                current,
                session_files.exports,
                metadata_column=metadata_column,
            )
            progress.set(1, message="Refmat export ready")
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
