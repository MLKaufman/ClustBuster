# ClustBuster

![ClustBuster logo](img/clustbuster-logo.png)

ClustBuster is a Python-first web application for assisted annotation of
single-cell RNA-sequencing clusters. The current Python-native application supports
H5AD plus a tested Seurat v5 RDS/H5Seurat compatibility slice.

## Development

Requirements: Python 3.11–3.14 and [uv](https://docs.astral.sh/uv/).

```bash
uv sync --extra dev
uv run pytest
uv run ruff check .
uv run shiny run --reload src/clustbuster/app.py
```

Container verification is available with:

```bash
docker build -t clustbuster:test .
scripts/container_smoke.sh clustbuster:test
```

The application is intentionally session-scoped. Uploaded data will live in an
ephemeral workspace and the source upload will never be modified. Download all
exports before the session or container is removed.

## Current scope

- Typed, format-independent workspace and expression-source models
- Independent cluster annotation state with type-safe cluster identifiers
- Session-isolated H5AD and Seurat upload, structural validation, and import reports
- Cluster-column, embedding, and expression-source configuration
- Interactive Plotly/WebGL embedding colored by source cluster or annotation
- Sparse-safe multi-gene feature plots and cluster-by-gene dot plots
- Sparse-safe standardized module scores with reusable immune-cell presets,
  embedding visualization, and per-cluster summaries
- Sparse-safe selected-cluster marker ranking with Welch statistics,
  multiple-testing correction, ranked tables, and cluster heatmaps
- GO Biological Process enrichment through a replaceable Enrichr adapter,
  with ranked terms and an interactive significance summary
- Validated local CSV marker catalogs and CSV/TSV/Parquet reference matrices,
  including provider health, filtering, checksums, and bundled demonstration data
- Reference-based cluster annotation through a pinned `pyclustifyr` adapter, with
  shared-gene validation, correlation previews, thresholds, and explicit apply/discard
- Persistent editable annotation/notes sidebar with selected/all reset and bounded undo/redo
- Cluster- and cell-level annotation CSVs packaged as a ZIP
- Annotated H5AD export with provenance and reopen validation
- Import and resource-provider protocols
- Configuration from a single validated startup object
- Shiny application and non-root Docker image
- Automated independent-session and hardened-container smoke tests

Manual fixtures are documented in `testdata/`. Imported Seurat objects export as CSV
and annotated H5AD; native Seurat or SingleCellExperiment RDS write-back is not
supported.

## Seurat compatibility

Set `CLUSTBUSTER_ENABLE_SEURAT_IMPORT=1` to accept `.rds` and `.h5seurat` uploads
(enabled by `compose.dev.yml`). The tested RDS contract is an in-memory Seurat v5
object with one active assay, cell metadata, compatible counts/scale layers, and
standard dimensional reductions. Compatible secondary assays with identical cell and
feature identifiers are mapped to namespaced layers; incompatible assays are reported.
The adapter pins `readseurat==0.1.0` and repairs its known Assay5 feature-coordinate
conversion defect. On-disk layers, custom S4 extensions, and native RDS export remain
unsupported.

GO enrichment is the only current workflow that sends data outside the local
container. An explicit Run action submits only the displayed marker gene symbols
and a non-identifying cluster description to the
[Ma'ayan Lab Enrichr API](https://maayanlab.cloud/Enrichr/); expression values,
cell identifiers, and observation metadata remain local.

## Local biological resources

The bundled demonstration marker catalog is in `resources/marker_sets/`, and
reference matrices plus JSON sidecars are in `resources/reference_matrices/`.
Override them with `CLUSTBUSTER_MARKER_CATALOG_PATH` and
`CLUSTBUSTER_REFERENCE_ROOT`. Provider files are read-only, validated before use,
and unavailable resources degrade only the Resources panel.

The bundled records are synthetic workflow fixtures, not an authoritative marker
database or a substitute for biological review.

## Reference annotation

After configuring a workspace, open **Reference annotation**, choose a validated
local reference, and run scoring. ClustBuster shows the full cluster-by-cell-type
similarity matrix and best-call preview. Nothing changes in annotation state until
**Apply predictions** is selected; applied records retain the score,
`source=pyclustifyr`, and the reference identifier/version. The bundled PBMC
reference is synthetic and intended only to test the workflow with
`testdata/clustbuster-demo.h5ad`.

The adapter is pinned to pyclustifyr revision
`db8761a87072b814f95ce7e0767540b4d10689bd`; API findings and current limitations
are recorded in `docs/research/pyclustifyr.md`.
