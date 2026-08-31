"""Minimal Shiny shell used while the AnnData MVP is assembled."""

from __future__ import annotations

from pathlib import Path

from shiny import App, ui

from clustbuster import __version__

app_ui = ui.page_fillable(
    ui.tags.style(
        """
        :root { --cb-navy: #19324a; --cb-teal: #2d8c88; --cb-bg: #f4f7f8; }
        body { background: var(--cb-bg); color: var(--cb-navy); }
        .cb-card { max-width: 760px; margin: 8vh auto; padding: 2rem; text-align: center; }
        .cb-mark { width: 108px; height: 108px; border-radius: 24px; }
        """
    ),
    ui.card(
        ui.tags.div(
            ui.tags.img(
                src="/assets/clustbuster-logo.png",
                class_="cb-mark",
                alt="ClustBuster logo",
            ),
            ui.h1("ClustBuster"),
            ui.p("Assisted annotation for single-cell RNA-sequencing clusters."),
            ui.p(
                "The application foundation is running. H5AD upload and annotation workflows "
                "are the next implementation slice."
            ),
            ui.tags.small(f"Version {__version__}"),
            class_="cb-card",
        )
    ),
    title="ClustBuster",
)


def server(input: object, output: object, session: object) -> None:
    """No reactive behavior is exposed in the Phase 0 shell."""


app = App(app_ui, server, static_assets={"/assets": Path(__file__).parent / "www"})


def main() -> None:
    import subprocess
    import sys

    raise SystemExit(subprocess.call(["shiny", "run", "clustbuster.app:app", *sys.argv[1:]]))
