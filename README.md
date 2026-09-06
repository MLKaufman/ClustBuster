# ClustBuster

![ClustBuster logo](img/clustbuster-logo.png)

ClustBuster is a Python-first web application for assisted annotation of
single-cell RNA-sequencing clusters. The current Python-native application supports
H5AD, tested Seurat v4/v5 RDS and H5Seurat objects, and tested in-memory
SingleCellExperiment RDS objects.

## Development

Requirements: Python 3.11–3.14 and [uv](https://docs.astral.sh/uv/).

```bash
uv sync --extra dev
uv run pytest
uv run ruff check .
CLUSTBUSTER_ENABLE_SEURAT_IMPORT=1 \
CLUSTBUSTER_ENABLE_SCE_IMPORT=1 \
uv run shiny run --reload src/clustbuster/app.py
```

Container verification is available with:

```bash
docker build -t clustbuster:test .
scripts/container_smoke.sh clustbuster:test
```

Production deployment templates for Traefik and ShinyProxy are under `deployment/`.
See [`docs/OPERATIONS.md`](docs/OPERATIONS.md) for TLS, authentication, resource
limits, session lifecycle, upgrades, troubleshooting, and cleanup.

The application is intentionally session-scoped. Uploaded data will live in an
ephemeral workspace and the source upload will never be modified. Download all
exports before the session or container is removed.

## Current scope

- Typed, format-independent workspace and expression-source models
- Independent cluster annotation state with type-safe cluster identifiers
- Session-isolated H5AD, Seurat, and SingleCellExperiment upload, structural validation,
  and import reports
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

## SingleCellExperiment compatibility

Set `CLUSTBUSTER_ENABLE_SCE_IMPORT=1` to accept `.rds` files containing tested,
in-memory SingleCellExperiment objects (enabled by `compose.dev.yml`). All aligned
assays are preserved as AnnData layers, `logcounts` is preferred as `X`, `colData`
and `rowData` map to observation and feature metadata, and `reducedDims` map to
AnnData embeddings. Delayed/on-disk assays, alternative experiments, pair data,
and custom S4 extensions are outside the current compatibility contract and are
rejected or reported explicitly. Exports remain CSV and annotated H5AD.

## Seurat compatibility

Set `CLUSTBUSTER_ENABLE_SEURAT_IMPORT=1` to accept `.rds` and `.h5seurat` uploads
(enabled by `compose.dev.yml`). The tested RDS contract covers in-memory Seurat v4
and v5 objects with one active assay, cell metadata, compatible expression layers,
and standard dimensional reductions. Compatible secondary assays with identical cell and
feature identifiers are mapped to namespaced layers; incompatible assays are reported.
The adapter pins `readseurat==0.1.0` and repairs its known Assay5 feature-coordinate
conversion defect. On-disk layers, custom S4 extensions, and native RDS export remain
unsupported.

The repository includes genuine Seurat v4 and v5 RDS fixtures plus a deterministic
H5Seurat compatibility fixture under `testdata/`.

`testdata/so-large.rds` is an optional 1.9 GB real-world stress fixture. The
development Compose configuration sizes Shiny's temporary upload area and the
session workspace for this file, but Docker Desktop must also have sufficient memory
and disk available. Its measured import peak was approximately 18.3 GiB, so 24 GiB
per session is the production default. It is intentionally excluded from the routine
automated suite.

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

### Sovereign Atlas / MarkerCodex

ClustBuster can instead consume a Sovereign Atlas checkout or immutable release
bundle directly. The two DuckDB catalogs are always opened read-only, and reference
matrix paths are resolved beneath a separately configured read-only files directory.
For a sibling checkout, start ClustBuster with:

```bash
CLUSTBUSTER_MARKER_DB_PATH=../sovereign-atlas/data/markercodex.duckdb \
CLUSTBUSTER_REFERENCE_CATALOG_PATH=../sovereign-atlas/data/reference_matrices.duckdb \
CLUSTBUSTER_REFERENCE_FILES_ROOT=../sovereign-atlas/data/reference_matrices/files \
CLUSTBUSTER_ATLAS_VERSION=preview-local \
uv run shiny run clustbuster.app:app
```

Alternatively, `./run_local.sh` automatically detects `../sovereign-atlas`, configures
all three data paths, and labels the resource with that checkout's current commit.
Set `CLUSTBUSTER_ATLAS_REPO=/another/path/sovereign-atlas` before running the script
when the checkout is elsewhere.

Setting a database path selects its DuckDB provider; leaving it unset preserves the
bundled CSV/local provider. `CLUSTBUSTER_ATLAS_VERSION` should be a release tag or
source commit. If it is omitted, the marker provider reports a short checksum of the
database instead.

The MarkerCodex provider targets the `marker_atlas` consumer view and recognizes a
small set of documented prerelease column aliases. It fails only the Resources panel
with a compatibility message when required consumer fields disappear. This keeps
schema adaptation at the provider boundary while Sovereign Atlas is still evolving.

DuckDB and the operating-system page cache handle marker queries. ClustBuster caches
only resolved schemas and successful reference checksums for a provider's lifetime;
it does not copy either catalog, cache query-expression data between sessions, or
download atlas data at runtime. Production deployments should mount a pinned atlas
bundle read-only or bake that exact bundle into the session image.

#### Updating Sovereign Atlas data

Atlas updates are intentionally not applied automatically. ClustBuster treats one
Sovereign Atlas revision as an immutable snapshot so a running analysis cannot mix
database metadata from one revision with reference-matrix bytes from another.

For a local sibling checkout:

1. Stop the running ClustBuster process.
2. Update and validate Sovereign Atlas according to that repository's instructions:

   ```bash
   git -C ../sovereign-atlas pull --ff-only
   git -C ../sovereign-atlas rev-parse --short HEAD
   ```

3. Restart ClustBuster with the reported commit as the resource version:

   ```bash
   CLUSTBUSTER_MARKER_DB_PATH=../sovereign-atlas/data/markercodex.duckdb \
   CLUSTBUSTER_REFERENCE_CATALOG_PATH=../sovereign-atlas/data/reference_matrices.duckdb \
   CLUSTBUSTER_REFERENCE_FILES_ROOT=../sovereign-atlas/data/reference_matrices/files \
   CLUSTBUSTER_ATLAS_VERSION=NEW_COMMIT \
   uv run shiny run clustbuster.app:app
   ```

   Replace `NEW_COMMIT` with the revision printed in step 2.

No manual cache deletion is required. Restarting ClustBuster discards its in-process
schema and checksum caches, while the operating system manages its own file cache.
Do not run the Sovereign Atlas curator or replace its databases and matrix files while
ClustBuster is using that checkout.

For production, publish or stage each atlas revision in a new versioned directory,
for example:

```text
/opt/sovereign-atlas/releases/2026.09.01/
/opt/sovereign-atlas/releases/2026.09.15/
```

Each directory should contain the two DuckDB databases and the complete
`reference_matrices/files/` directory from the same revision. Validate the new
bundle, mount it read-only into new ClustBuster session containers, set
`CLUSTBUSTER_ATLAS_VERSION` to its release tag or commit, and restart the application.
Keep the preceding directory unchanged until the new deployment has been verified;
rollback then consists of restoring the previous image or read-only mount.

Compatible schema changes require no ClustBuster changes. If Sovereign Atlas removes
or renames a required consumer field, the affected MarkerCodex or Refmats panel shows
a compatibility warning instead of preventing the rest of ClustBuster from starting.
Update the corresponding DuckDB provider's centralized column aliases or contract,
run the provider tests, and redeploy before adopting that atlas revision.

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
