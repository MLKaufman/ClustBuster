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
from clustbuster.core.enrichment import EnrichmentResult
from clustbuster.core.expression import (
    DotPlotResult,
    ExpressionResult,
    aggregate_dotplot,
    extract_expression,
    parse_gene_list,
)
from clustbuster.core.markers import MarkerResult, rank_markers
from clustbuster.core.modules import MODULE_PRESETS, ModuleScoreResult, calculate_module_score
from clustbuster.core.workspace import workspace_from_import
from clustbuster.integrations.enrichr import DEFAULT_LIBRARY, EnrichrClient
from clustbuster.models import ExpressionSource, Workspace
from clustbuster.plotting.dotplot import dotplot_figure
from clustbuster.plotting.embedding import embedding_figure
from clustbuster.plotting.enrichment import enrichment_figure
from clustbuster.plotting.feature import feature_figure
from clustbuster.plotting.heatmap import marker_heatmap_figure
from clustbuster.plotting.module import module_score_figure
from clustbuster.services.exports import WorkspaceExportService
from clustbuster.services.imports import ImportService, SessionFiles
from clustbuster.services.workspaces import configure_workspace

config = AppConfig.from_env()
logging.basicConfig(level=config.log_level)
logger = logging.getLogger("clustbuster")


def _styles() -> ui.Tag:
    return ui.tags.style(
        """
        :root { --cb-navy: #19324a; --cb-teal: #2d8c88; --cb-bg: #f4f7f8; }
        body { background: var(--cb-bg); color: var(--cb-navy); }
        .cb-header { display:flex; align-items:center; gap:.8rem; padding:.65rem 1rem; }
        .cb-logo { width:62px; height:62px; border-radius:14px; }
        .cb-title { margin:0; font-size:1.55rem; font-weight:700; }
        .cb-subtitle { margin:0; color:#607180; font-size:.88rem; }
        .cb-empty { min-height:420px; display:flex; align-items:center; justify-content:center;
                    flex-direction:column; color:#6a7c89; text-align:center; padding:3rem; }
        .cb-summary { display:grid; grid-template-columns:repeat(2, minmax(0,1fr)); gap:.6rem; }
        .cb-stat { background:#fff; border:1px solid #dbe4e8; border-radius:.6rem; padding:.65rem; }
        .cb-stat strong { display:block; font-size:1.15rem; color:var(--cb-navy); }
        .btn-primary { background-color:var(--cb-teal); border-color:var(--cb-teal); }
        """
    )


app_ui = ui.page_fillable(
    _styles(),
    ui.tags.header(
        ui.tags.img(
            src="/assets/clustbuster-logo.png", class_="cb-logo", alt="ClustBuster logo"
        ),
        ui.tags.div(
            ui.h1("ClustBuster", class_="cb-title"),
            ui.p("Assisted single-cell cluster annotation", class_="cb-subtitle"),
        ),
        class_="cb-header",
    ),
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_file(
                "dataset",
                "Upload AnnData",
                accept=[".h5ad", "application/x-hdf5"],
                button_label="Choose H5AD",
                placeholder="No dataset selected",
            ),
            ui.output_ui("import_panel"),
            title="Workspace",
            width=330,
            open="desktop",
        ),
        ui.navset_card_tab(
            ui.nav_panel(
                "Overview",
                ui.output_ui("overview_header"),
                ui.card(
                    ui.card_header("Annotation controls"),
                    ui.layout_columns(
                        ui.input_select(
                            "color_by",
                            "Color embedding by",
                            {"cluster": "Source cluster", "annotation": "Current annotation"},
                            selected="cluster",
                        ),
                        ui.input_select("annotation_cluster", "Cluster", {}),
                        ui.input_text(
                            "annotation_label",
                            "Annotation",
                            placeholder="e.g. CD4 T cell",
                        ),
                        ui.help_text("Labels save automatically as you edit."),
                        col_widths=(4, 3, 3, 2),
                    ),
                    fill=False,
                ),
                output_widget("embedding_plot", height="580px"),
            ),
            ui.nav_panel(
                "Import report",
                ui.output_ui("report_summary"),
                ui.output_data_frame("report_table"),
            ),
            ui.nav_panel(
                "Feature plot",
                ui.card(
                    ui.layout_columns(
                        ui.input_text_area(
                            "feature_genes",
                            "Genes",
                            value="CD3D, LYZ",
                            placeholder="Comma, space, or newline separated",
                            rows=2,
                        ),
                        ui.input_action_button(
                            "run_feature", "Run feature plot", class_="btn-primary"
                        ),
                        col_widths=(9, 3),
                    ),
                    fill=False,
                ),
                ui.output_ui("feature_feedback"),
                output_widget("feature_plot", height="580px"),
            ),
            ui.nav_panel(
                "Dot plot",
                ui.card(
                    ui.layout_columns(
                        ui.input_text_area(
                            "dot_genes",
                            "Gene panel",
                            value="CD3D, LYZ, MS4A1, NKG7",
                            placeholder="Comma, space, or newline separated",
                            rows=2,
                        ),
                        ui.input_action_button(
                            "run_dotplot", "Run dot plot", class_="btn-primary"
                        ),
                        col_widths=(9, 3),
                    ),
                    fill=False,
                ),
                ui.output_ui("dotplot_feedback"),
                output_widget("dot_plot", height="580px"),
            ),
            ui.nav_panel(
                "Module scores",
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
                        ui.input_text(
                            "module_name", "Score label", value="T cell module"
                        ),
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
                ui.card(
                    ui.card_header("Cluster summary"),
                    ui.output_data_frame("module_summary"),
                ),
            ),
            ui.nav_panel(
                "Markers & heatmap",
                ui.card(
                    ui.layout_columns(
                        ui.input_select("marker_cluster", "Selected cluster", {}),
                        ui.input_numeric(
                            "marker_top_n", "Top genes", value=12, min=1, max=50
                        ),
                        ui.input_numeric(
                            "marker_min_fraction",
                            "Minimum expressing fraction",
                            value=0.1,
                            min=0,
                            max=1,
                            step=0.05,
                        ),
                        ui.input_action_button(
                            "run_markers", "Rank markers", class_="btn-primary"
                        ),
                        col_widths=(3, 2, 4, 3),
                    ),
                    ui.help_text(
                        "Ranks positive markers for the selected cluster versus all "
                        "remaining cells using Welch's t-test with Benjamini-Hochberg correction."
                    ),
                    fill=False,
                ),
                ui.output_ui("marker_feedback"),
                ui.card(
                    ui.card_header("Ranked markers"),
                    ui.output_data_frame("marker_table"),
                ),
                output_widget("marker_heatmap", height="520px"),
            ),
            ui.nav_panel(
                "Enrichment",
                ui.card(
                    ui.layout_columns(
                        ui.input_select(
                            "enrichment_library",
                            "Gene-set library",
                            {DEFAULT_LIBRARY: "GO Biological Process 2025"},
                        ),
                        ui.input_numeric(
                            "enrichment_top_n", "Terms to display", value=12, min=3, max=30
                        ),
                        ui.input_action_button(
                            "run_enrichment", "Run enrichment", class_="btn-primary"
                        ),
                        col_widths=(5, 3, 4),
                    ),
                    ui.help_text(
                        "Uses the current ranked marker genes. Running enrichment sends only "
                        "those gene symbols to the Ma'ayan Lab Enrichr service; expression "
                        "values and cell metadata are not transmitted."
                    ),
                    fill=False,
                ),
                ui.output_ui("enrichment_feedback"),
                output_widget("enrichment_plot", height="560px"),
                ui.card(
                    ui.card_header("Enriched terms"),
                    ui.output_data_frame("enrichment_table"),
                ),
            ),
            ui.nav_panel(
                "Export",
                ui.output_ui("export_status"),
                ui.card(
                    ui.card_header("Download annotations"),
                    ui.p(
                        "The ZIP contains separate cluster-level and cell-level CSV files."
                    ),
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
            title=ui.tags.span(f"AnnData MVP · v{__version__}"),
            full_screen=True,
        ),
        fillable=True,
    ),
    title="ClustBuster",
)


def server(input: Inputs, output: Outputs, session: Session) -> None:
    session_files = SessionFiles.create(config.workspace_root)
    import_service = ImportService(config.max_upload_mb)
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
    enrichment_client = EnrichrClient(
        timeout_seconds=config.enrichment_timeout_seconds
    )

    session.on_ended(session_files.cleanup)

    @reactive.effect
    @reactive.event(input.dataset)
    def import_dataset() -> None:
        upload_value = input.dataset()
        if not upload_value:
            return
        error_message.set(None)
        configured.set(False)
        with ui.Progress(min=0, max=1, session=session) as progress:
            progress.set(0.15, message="Validating H5AD", detail="Reading object metadata")
            try:
                uploaded = cast(dict[str, Any], upload_value[0])
                result = import_service.import_upload(uploaded, session_files)
                workspace.set(workspace_from_import(result))
                progress.set(1, message="Import complete")
                ui.notification_show(
                    f"Loaded {result.report.cell_count:,} cells and "
                    f"{result.report.feature_count:,} features",
                    type="message",
                    session=session,
                )
            except Exception as exc:
                logger.exception("H5AD import failed", extra={"error_type": type(exc).__name__})
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
            return ui.p("Upload an H5AD file to begin.", class_="text-muted")
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
            ui.input_action_button(
                "configure", "Configure workspace", class_="btn-primary w-100"
            ),
        )

    @reactive.effect
    @reactive.event(input.configure)
    def apply_configuration() -> None:
        current = workspace.get()
        req(current is not None)
        assert current is not None
        try:
            configure_workspace(
                current,
                cluster_column=str(input.cluster_column()),
                embedding_key=str(input.embedding_key()),
                expression_source=ExpressionSource.from_label(str(input.expression_source())),
            )
            cluster_choices = {
                record.cluster_id.serialized: record.cluster_id.display
                for record in current.annotations.records()
            }
            ui.update_select(
                "annotation_cluster", choices=cluster_choices, session=session
            )
            ui.update_select("marker_cluster", choices=cluster_choices, session=session)
            configured.set(True)
            feature_result.set(None)
            dotplot_result.set(None)
            module_score_result.set(None)
            marker_result.set(None)
            enrichment_result.set(None)
            revision.set(revision.get() + 1)
            ui.notification_show("Workspace configured", type="message", session=session)
        except Exception as exc:
            ui.notification_show(str(exc), type="error", duration=8, session=session)

    @output
    @render.ui
    def overview_header() -> ui.TagChild:
        current = workspace.get()
        if current is None:
            return ui.div(
                ui.h3("Start with an H5AD workspace"),
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

    @reactive.effect
    @reactive.event(input.annotation_label, ignore_init=True)
    def assign_annotation() -> None:
        annotation_label = str(input.annotation_label()).strip()
        req(annotation_label)
        with reactive.isolate():
            current = workspace.get()
            req(current is not None and configured.get())
            assert current is not None
            try:
                current.annotations.assign_serialized(
                    str(input.annotation_cluster()),
                    annotation_label,
                    source="manual",
                )
                revision.set(revision.get() + 1)
                ui.notification_show("Annotation updated", type="message", session=session)
            except Exception as exc:
                ui.notification_show(str(exc), type="error", duration=8, session=session)

    @output
    @render_plotly
    def embedding_plot() -> Any:
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
        warning_text = (
            ui.tags.ul(*(ui.tags.li(item) for item in report.warnings))
            if report.warnings
            else ui.p("No import warnings.", class_="text-success")
        )
        return ui.card(
            ui.card_header(report.source_filename),
            ui.p(
                f"{report.cell_count:,} cells x {report.feature_count:,} features · "
                f"{'sparse' if report.sparse else 'dense'} expression matrix"
            ),
            warning_text,
            fill=False,
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

    @reactive.effect
    @reactive.event(input.run_feature)
    def run_feature_plot() -> None:
        current = workspace.get()
        req(current is not None and configured.get())
        assert current is not None
        feature_error.set(None)
        try:
            genes = parse_gene_list(str(input.feature_genes()), limit=6)
            result = extract_expression(current.adata, current.expression_source, genes)
            feature_result.set(result)
        except Exception as exc:
            feature_result.set(None)
            feature_error.set(str(exc))
            ui.notification_show(str(exc), type="error", duration=8, session=session)

    @output
    @render.ui
    def feature_feedback() -> ui.TagChild:
        error = feature_error.get()
        result = feature_result.get()
        if error:
            return ui.div(error, class_="alert alert-danger")
        if result is None:
            return ui.p("Configure a workspace, enter genes, and run the plot.")
        missing = ", ".join((*result.report.missing, *result.report.ambiguous))
        if missing:
            return ui.div(f"Not plotted: {missing}", class_="alert alert-warning")
        return ui.div(
            f"Matched {len(result.report.matched)} gene(s).", class_="alert alert-success"
        )

    @output
    @render_plotly
    def feature_plot() -> Any:
        current = workspace.get()
        result = feature_result.get()
        req(current is not None and configured.get() and result is not None)
        assert current is not None and result is not None and current.embedding_key is not None
        coordinates = np.asarray(current.adata.obsm[current.embedding_key])
        return feature_figure(coordinates, current.adata.obs_names.astype(str).tolist(), result)

    @reactive.effect
    @reactive.event(input.run_dotplot)
    def run_dot_plot() -> None:
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
            ui.notification_show(str(exc), type="error", duration=8, session=session)

    @output
    @render.ui
    def dotplot_feedback() -> ui.TagChild:
        error = dotplot_error.get()
        result = dotplot_result.get()
        if error:
            return ui.div(error, class_="alert alert-danger")
        if result is None:
            return ui.p("Configure a workspace, enter a gene panel, and run the plot.")
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
        result = dotplot_result.get()
        req(result is not None)
        assert result is not None
        return dotplot_figure(result)

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
        current = workspace.get()
        result = module_score_result.get()
        req(current is not None and configured.get() and result is not None)
        assert current is not None and result is not None and current.embedding_key is not None
        coordinates = np.asarray(current.adata.obsm[current.embedding_key])
        return module_score_figure(
            coordinates, current.adata.obs_names.astype(str).tolist(), result
        )

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
                "mean_difference",
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
                "mean_difference": "Mean difference",
                "fraction_selected": "Fraction selected",
                "fraction_rest": "Fraction rest",
            }
        )

    @output
    @render_plotly
    def marker_heatmap() -> Any:
        result = marker_result.get()
        req(result is not None)
        assert result is not None
        return marker_heatmap_figure(result)

    @reactive.effect
    @reactive.event(input.run_enrichment)
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

    @output
    @render.ui
    def enrichment_feedback() -> ui.TagChild:
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
    @render_plotly
    def enrichment_plot() -> Any:
        result = enrichment_result.get()
        req(result is not None)
        assert result is not None
        return enrichment_figure(result, top_n=int(input.enrichment_top_n()))

    @output
    @render.data_frame
    def enrichment_table() -> pd.DataFrame:
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
