# ClustBuster

ClustBuster is a Python-first web application for assisted annotation of
single-cell RNA-sequencing clusters. The current implementation is an early
Phase 0 foundation; the supported MVP input will be H5AD.

## Development

Requirements: Python 3.11–3.14 and [uv](https://docs.astral.sh/uv/).

```bash
uv sync --extra dev
uv run pytest
uv run ruff check .
uv run shiny run --reload src/clustbuster/app.py
```

The application is intentionally session-scoped. Uploaded data will live in an
ephemeral workspace and the source upload will never be modified. Download all
exports before the session or container is removed.

## Current scope

- Typed, format-independent workspace and expression-source models
- Independent cluster annotation state with type-safe cluster identifiers
- Import and resource-provider protocols
- Configuration from a single validated startup object
- Minimal Shiny application shell and non-root Docker image

H5AD import, analysis views, and export workflows are being built next. Seurat
and SingleCellExperiment RDS round-tripping is not currently supported.

