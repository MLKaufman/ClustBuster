# ClustBuster Implementation Plan

## 1. Purpose

ClustBuster is a Python-first web application for assisted annotation of single-cell RNA-sequencing clusters. It will reproduce and improve the core workflow of the `mlkaufman/SOWhat` R Shiny application while accepting Scanpy/AnnData, Seurat, and SingleCellExperiment inputs through one format-independent workspace.

This document is the implementation contract for the project. Codex should use it to scaffold the repository, make architectural decisions, implement features in the stated phases, and evaluate phase acceptance criteria. Where a library API or file-format capability is uncertain, implement a narrow adapter and a tested research spike before coupling the rest of the application to it.

## 2. Product goals

ClustBuster must:

- Let a scientist load an `.h5ad`, Seurat `.rds`/`.h5seurat`, or SingleCellExperiment `.rds` object.
- Normalize supported inputs to an `anndata.AnnData` workspace at the import boundary.
- Let the user select a clustering column, embedding, and expression source rather than assuming format-specific names.
- Support interactive cluster annotation using dimensional-reduction, feature, dot, marker, heatmap, module-score, enrichment, and reference-correlation views.
- Integrate `pyclustifyr` for reference-based cluster annotation.
- Preserve the source clustering and maintain editable annotation state separately.
- Always export annotations as CSV and export an annotated H5AD object.
- Eventually support safe native Seurat and SingleCellExperiment write-back, without making native RDS export a prerequisite for the MVP.
- Run in Docker and support isolated concurrent user sessions on a local compute network.
- Expose stable provider interfaces for a future marker database and reference-matrix database.

### Non-goals for the MVP

- Reproducing every possible Seurat or Bioconductor object variant.
- Writing fully valid Seurat or SingleCellExperiment RDS objects from Python.
- Running R in the main application image.
- Recomputing an entire single-cell analysis pipeline from raw reads.
- Browser-side matrix analysis or WebAssembly optimization.
- Shared, real-time editing of one workspace by multiple users.
- Permanent server-side storage of user uploads after a session ends.

## 3. Guiding decisions

1. **Python is the application language.** Do not add an R runtime to the primary application container.
2. **AnnData is the canonical internal model.** Format-specific behavior belongs in `io/`; all analysis, plotting, UI, and annotation code operates on AnnData or derived typed models.
3. **H5AD is the native interchange and object-export format.** CSV is the universal annotation output.
4. **Scanpy is the analysis backend, not the main browser plotting layer.** Use Plotly/WebGL for interactive plots.
5. **Shiny for Python is the initial web framework.** Preserve clean boundaries so UI modules can be replaced later without rewriting core analysis.
6. **Source data is immutable by default.** Annotation changes affect session state and explicit exports, never the uploaded source file.
7. **ClustBuster owns analysis and interaction; external databases own biological knowledge.** Marker and reference resources are accessed through provider protocols.
8. **Optimize only after profiling.** Preserve sparse matrices, cache expensive derived results, and consider Rust/PyO3 only for demonstrated server-side bottlenecks.

## 4. System architecture

```text
Browser
  |
  v
Traefik                  reverse proxy, TLS, routing, optional auth entry point
  |
  v
ShinyProxy               authentication, per-user container lifecycle, quotas
  |
  +---- ClustBuster container: session A ----+
  +---- ClustBuster container: session B ----+-- isolated temporary workspaces
  +---- ClustBuster container: session C ----+

Inside each ClustBuster container

Shiny for Python UI
  |
  +--> reactive session/application services
          |
          +--> Workspace + AnnotationStore
          +--> Scanpy/NumPy/SciPy analysis services
          +--> pyclustifyr adapter
          +--> Plotly figure builders
          +--> MarkerProvider
          +--> ReferenceProvider
          +--> import/export adapters
                  |
                  v
             canonical AnnData
```

Traefik routes traffic but does not create user containers. ShinyProxy owns application-container creation and teardown. Each user session receives a separate ClustBuster process, memory space, and temporary workspace.

### Application boundaries

- `ui/` owns layout, input bindings, reactive effects, notifications, and download handlers.
- `services/` coordinates use cases and is the only layer the UI calls for substantive work.
- `core/` owns format-independent models and analysis logic.
- `plotting/` turns validated, compact data into Plotly figures; it does not mutate AnnData.
- `io/` detects, imports, validates, and exports files.
- `resources/` provides marker and reference-matrix access through protocols.
- `integrations/` wraps external packages such as `pyclustifyr`.

No UI module should call `readseurat`, `rdata`, DuckDB, or database-specific code directly.

## 5. Canonical AnnData model

All successful imports produce an `AnnData` object with these conventions:

| AnnData location | ClustBuster meaning |
|---|---|
| `adata.X` | Default expression matrix when explicitly selected |
| `adata.layers[...]` | Counts, normalized, scaled, log-count, or other named matrices |
| `adata.raw` | Optional immutable raw expression representation |
| `adata.obs` | Cell metadata, cluster labels, and exported annotation columns |
| `adata.var` | Feature metadata; `var_names` is the canonical feature identifier |
| `adata.obsm[...]` | Embeddings such as `X_umap`, `X_tsne`, `X_pca`, and `X_harmony` |
| `adata.uns["clustbuster"]` | Namespaced export metadata only |

Do not silently densify sparse matrices. Core functions must accept SciPy sparse and dense inputs and return small derived tables when practical.

### Import normalization rules

- Cell identifiers and feature identifiers must be non-null and unique. Use AnnData's uniqueness helpers only after recording that names changed.
- Preserve original metadata column names.
- Preserve all safely mapped expression layers and embeddings.
- Record import warnings and a mapping report; never silently discard an assay, layer, embedding, or identifier collision.
- Never guess that counts are normalized or that normalized data are counts. Store what the source provides and expose the choice to the user.
- Do not run normalization, scaling, PCA, neighbor graph construction, or clustering automatically on upload.
- Validate that every matrix has dimensions compatible with `obs` and `var`.

### Workspace model

Implement a typed session model similar to:

```python
@dataclass
class Workspace:
    adata: AnnData
    source_format: ObjectFormat
    source_filename: str
    cluster_column: str | None
    embedding_key: str | None
    expression_source: ExpressionSource
    annotations: AnnotationStore
    import_report: ImportReport
```

`ExpressionSource` must distinguish `X`, `raw`, and a named `layer`; do not encode all three as ambiguous strings. `Workspace` is session-scoped and must not be stored in module globals.

## 6. Input strategy

Implement a common importer contract:

```python
class ObjectImporter(Protocol):
    def probe(self, path: Path) -> ProbeResult: ...
    def load(self, path: Path, options: ImportOptions) -> ImportResult: ...
```

`ImportResult` contains the AnnData object plus source-format metadata, warnings, mapping decisions, and unsupported source components. File-extension checks are advisory; importers must report useful parse and validation errors.

Uploads first enter a random, session-specific directory with sanitized server-generated filenames. Enforce configurable upload-size, free-disk, and extraction limits. Never interpolate an uploaded filename into a shell command.

### H5AD

- Use `anndata.read_h5ad()`.
- Treat H5AD as the reference path and implement it first.
- Initially load normal MVP-sized files into memory.
- Add a configurable threshold for backed/lazy strategies after measuring actual workflows; do not claim full backed-mode compatibility until every enabled analysis path is tested.
- Validate annotations, layers, `raw`, embeddings, and sparse matrices with fixture files.

### Seurat

- Prefer a pure-Python Seurat adapter built around the maintained `readseurat` API for supported `.rds` and `.h5seurat` inputs.
- Keep `readseurat` behind `io/seurat.py`; no other package may depend on its data structures.
- Map cells to `obs`, features to `var`, selected/default assay data to `X`, other usable assay representations to named layers where dimensions permit, reductions to `obsm`, and cell metadata to `obs`.
- Preserve assay and layer provenance in the import report.
- Create versioned fixtures covering representative Seurat v4 and v5 objects, sparse data, multiple assays, metadata, and reductions.
- If a source structure is unsupported, fail with an actionable compatibility message; do not return a partially corrupted workspace.

The Seurat importer is a compatibility boundary, not a promise that every serialized Seurat object can be decoded without R. Before declaring support, pin and verify the chosen package/API against repository fixtures.

### SingleCellExperiment

Treat direct SCE `.rds` support as a dedicated compatibility milestone:

1. Evaluate `rdata` parsing plus current BiocPy `SingleCellExperiment` interoperability.
2. Prefer conversion through a maintained BiocPy object if it preserves assays, `colData`, `rowData`, and `reducedDims` correctly.
3. Add narrow custom converters only for missing, documented R S4 classes.
4. Map assays to `X`/layers, `colData` to `obs`, `rowData` to `var`, and `reducedDims` to `obsm`.
5. Test dense and sparse assays, multiple assays, missing reduced dimensions, categorical metadata, and common delayed representations.

SCE support is complete only for the explicitly tested compatibility matrix. Unsupported delayed/on-disk or custom S4 content must be reported clearly. A separate, documented offline R-to-H5AD conversion path may be offered as a fallback, but the main ClustBuster image remains Python-only.

### Object detection and selection

After import:

- Suggest cluster columns from categorical, integer-like, and low-cardinality `obs` columns; let the user choose any compatible column.
- Discover embeddings from compatible two-or-more-column arrays in `obsm`; prefer `X_umap`, then `X_tsne`, but require user confirmation.
- Enumerate `X`, `raw` when available, and every compatible named layer as expression choices.
- Show an import summary containing cells, genes, sparsity, layers, candidate cluster columns, embeddings, and warnings before analysis.

## 7. Analysis backend

Use Scanpy, AnnData, NumPy, SciPy, and pandas for server-side processing. Core functions should be deterministic, independently testable, and accept explicit cluster, layer, and parameter inputs.

### Required analysis services

- Expression extraction for one or more genes without densifying the whole matrix.
- Cluster-wise fraction-expressing and mean-expression aggregation for dot plots.
- Module scoring through `scanpy.tl.score_genes()` or an equivalently validated method.
- Differential marker ranking through `scanpy.tl.rank_genes_groups()` with explicit method and comparison parameters.
- Cluster-by-gene matrices for heatmaps.
- GO Biological Process enrichment through a replaceable enrichment adapter, initially using a Python-native library such as GSEApy or GOATOOLS.
- Gene identifier matching with explicit case/alias behavior and a report of missing or duplicate genes.

Never mutate shared state merely to calculate a view. When a Scanpy operation writes into `uns` or `obs`, run it on an appropriate copy or under a unique namespaced result key, then return a typed result.

## 8. Shiny for Python and Plotly frontend

Build the UI with modular Shiny for Python components and Plotly figures. Use `shinywidgets` or the current supported Shiny/Plotly bridge selected during scaffolding.

### Initial application flow

1. Upload and validate an object.
2. Review the import report.
3. Select cluster column, embedding, and expression source.
4. Explore clusters and evidence.
5. Edit annotations and optional metadata.
6. Export CSV and annotated H5AD.

### Primary views

- **Overview/annotation:** WebGL embedding colored by source cluster or current annotation, cluster selection, editable annotation table, and assignment controls.
- **Feature plot:** one or more genes, expression color scale, missing-gene feedback, and optional selected-cluster emphasis.
- **Dot plot:** configurable cluster grouping and gene panel, with dot size for fraction expressing and color for aggregated expression.
- **Module scores:** gene-list input/presets, score calculation, and embedding/table visualization.
- **Reference annotation:** reference selection, `pyclustifyr` parameters, correlation results, and explicit apply/reject controls.
- **Markers/heatmap:** differential-marker parameters, ranked results, export, and heatmap.
- **Enrichment:** selected marker set, organism/database settings, results table, and Plotly summary.
- **Export:** validation summary and downloads.

### Interaction rules

- Use `scattergl` for ordinary large embeddings and preserve stable cell identifiers in Plotly `customdata`.
- Cluster assignment acts on cluster identifiers, not the set of points currently rendered after downsampling.
- Long-running actions show progress and do not trigger on every reactive input change; require an explicit Run action.
- Cache derived data by workspace identity, expression source, cluster column, gene set, and analysis parameters.
- Display actionable errors in the relevant panel while logging full exception details server-side.
- Disable analysis controls until import and workspace validation succeed.

## 9. `pyclustifyr` integration

Wrap the external project in `integrations/pyclustifyr.py`. The wrapper must:

- Accept the current AnnData workspace, selected cluster column, expression source, a loaded reference matrix, and explicit parameters.
- Validate gene overlap and report overlap counts before running.
- Avoid relying on undocumented mutation of `adata.obs`.
- Return a typed result containing per-cluster predictions, correlations/scores, reference metadata, warnings, and package version.
- Let the user preview results before applying them to annotation state.
- Record `source="pyclustifyr"` and the reference identifier/version when results are applied.

Pin `pyclustifyr` to a known-compatible release or Git revision until its public API is stable. Add a small contract test that runs against a fixture AnnData/reference pair.

## 10. External biological-resource providers

The first release uses local files, but application code must target provider protocols so the future databases can be added without changing analysis or UI behavior.

### Marker provider

```python
class MarkerProvider(Protocol):
    def search_cell_types(
        self, query: str, *, species: str | None = None,
        tissue: str | None = None, limit: int = 50
    ) -> list[CellTypeSummary]: ...

    def get_markers(
        self, cell_type: str, *, species: str | None = None,
        tissue: str | None = None
    ) -> MarkerSet: ...
```

Canonical marker records should support:

- cell type
- gene identifier/symbol
- species
- tissue
- positive or negative marker direction
- evidence/source
- citation
- confidence
- resource version

Implement `CsvMarkerProvider` first. Add `DuckDbMarkerProvider` later using parameterized queries and read-only connections. Provider errors and empty results must be distinguishable.

### Reference-matrix provider

```python
class ReferenceProvider(Protocol):
    def list_references(self, filters: ReferenceFilters) -> list[ReferenceSummary]: ...
    def load_reference(self, reference_id: str) -> LoadedReference: ...
```

Reference metadata should support:

- stable ID and name
- species
- tissue
- disease/context
- assay/platform
- source dataset and citation
- cell-type coverage
- matrix location, checksum, schema version, and resource version

Implement `LocalReferenceProvider` first for validated CSV/TSV/Parquet matrices plus sidecar metadata. A future database-backed provider may resolve metadata and return local or streamed artifacts, but `pyclustifyr` must receive the same `LoadedReference` model either way.

### Provider configuration

Select implementations using configuration/dependency injection. Do not scatter environment-variable lookups through core modules. Providers must expose a health/status result so unavailable optional resources degrade the relevant panel rather than breaking the app.

## 11. Annotation state

Never overwrite the chosen source cluster column. Maintain an independent table keyed by the exact serialized cluster identifier:

| Field | Required | Meaning |
|---|---:|---|
| `cluster_id` | yes | Source cluster value |
| `annotation` | yes | Current human-readable cell-type label |
| `notes` | no | User notes |
| `source` | yes | `manual`, `pyclustifyr`, `imported`, etc. |
| `confidence` | no | Controlled category or validated numeric score |
| `reference_id` | no | Reference used for a prediction |
| `updated_at` | yes | UTC ISO-8601 timestamp |

Initialize one row per observed cluster. Preserve the distinction between numeric `1`, string `"1"`, and missing values during validation and serialization.

`AnnotationStore` should provide typed operations for assigning one cluster, assigning several clusters, applying a previewed prediction set, resetting selected/all annotations, and materializing cell-level annotations. Keep a bounded undo/redo history for the session if feasible after the basic editor is stable.

Only an explicit export materializes the current labels into:

```text
adata.obs["clustbuster_cluster"]
adata.obs["clustbuster_annotation"]
```

Store provenance under `adata.uns["clustbuster"]`, including ClustBuster version, selected source cluster, expression source, timestamp, and annotation-table schema version. Do not overwrite an existing ClustBuster export column without an explicit replace/version choice.

## 12. Export behavior

### Annotation CSV (required for every supported input)

Export both:

1. A cluster-level annotation table containing all annotation-state fields.
2. A cell-level table containing `cell_id`, `cluster_id`, and `annotation` (plus selected optional provenance fields).

If the UI exposes one primary CSV download, package both CSVs in a ZIP or label the two downloads unambiguously. CSV output must round-trip commas, quotes, Unicode, missing values, and nonnumeric cluster IDs.

### Annotated H5AD (required)

- Create an export copy or temporary on-disk artifact; do not mutate the source upload.
- Add the two namespaced `obs` columns and namespaced `uns` provenance described above.
- Preserve existing matrices, layers, embeddings, and metadata.
- Re-open and validate the written H5AD before exposing it for download.
- Generate collision-safe output names and checksums for logged diagnostics.

### Native R output (post-MVP)

Python's ability to serialize an RDS payload does not by itself guarantee a valid Seurat or SCE S4 object. Native `.rds` write-back is a separate milestone requiring versioned round-trip tests in R. Until those tests pass:

- H5AD input -> CSV + annotated H5AD
- Seurat input -> CSV + annotated H5AD
- SCE input -> CSV + annotated H5AD

Do not label a generic Python-written RDS file as a native Seurat or SCE export.

## 13. Session storage and security

Use a session directory such as:

```text
/workspace/<opaque-session-id>/
  uploads/
  cache/
  state/
  exports/
```

Requirements:

- Generate the session ID server-side; do not use usernames or uploaded filenames as paths.
- Restrict file operations to the session root and reject traversal.
- Apply upload size, memory, CPU, disk, and session-duration limits through configuration and deployment.
- Default to ephemeral storage and remove the per-session workspace when the user container ends.
- Do not log expression values, annotation contents, access tokens, or identifying cell/sample metadata.
- Document the privacy model and the fact that exports must be downloaded before session teardown.
- Use non-root application containers, a read-only root filesystem where practical, a writable workspace mount, and pinned base-image/dependency versions.

## 14. Deployment

### Container image

Build one production ClustBuster image that:

- installs the locked Python environment (prefer `uv` with `pyproject.toml` and `uv.lock`);
- runs the Shiny for Python app on a configurable host/port;
- includes a health endpoint or health command;
- runs as a non-root user;
- writes only to configured temporary/cache/workspace locations;
- contains no R runtime unless a future, separately approved compatibility image is introduced.

### ShinyProxy

Configure ShinyProxy to launch one ClustBuster container per authenticated user/session. Define:

- image and network;
- CPU and memory limits;
- maximum upload/workspace disk policy;
- idle and absolute timeouts;
- maximum concurrent sessions;
- cleanup behavior;
- allowed mounts and environment variables;
- authentication/authorization appropriate to the home network.

Start with a measured resource profile, for example 4 CPU and 16 GB RAM per large-data session, but keep it configurable rather than hard-coded.

### Traefik

Traefik terminates TLS and routes requests to ShinyProxy. Configure forwarded headers and WebSocket/long-lived connection behavior required by Shiny. Do not give Traefik responsibility for per-user container lifecycle. Do not expose the Docker socket broadly; use the least-privileged supported integration or a protected socket proxy where applicable.

### Deployment deliverables

- `Dockerfile`
- `.dockerignore`
- local development Compose file for ClustBuster
- deployment Compose/configuration for Traefik and ShinyProxy
- example environment file with no secrets
- health checks, resource limits, persistent configuration, and ephemeral workspace volumes
- operator documentation for build, upgrade, backup of configuration, troubleshooting, and cleanup

## 15. Proposed repository structure

```text
ClustBuster/
├── implementation.md
├── README.md
├── CONTRIBUTING.md
├── LICENSE
├── pyproject.toml
├── uv.lock
├── Dockerfile
├── .dockerignore
├── compose.dev.yml
├── deploy/
│   ├── compose.yml
│   ├── shinyproxy/
│   │   └── application.yml
│   └── traefik/
│       └── dynamic.yml
├── src/clustbuster/
│   ├── __init__.py
│   ├── app.py
│   ├── config.py
│   ├── models.py
│   ├── io/
│   │   ├── detect.py
│   │   ├── h5ad.py
│   │   ├── seurat.py
│   │   ├── sce.py
│   │   ├── validation.py
│   │   └── export.py
│   ├── core/
│   │   ├── workspace.py
│   │   ├── annotations.py
│   │   ├── expression.py
│   │   ├── markers.py
│   │   ├── modules.py
│   │   ├── enrichment.py
│   │   └── cache.py
│   ├── services/
│   │   ├── imports.py
│   │   ├── analysis.py
│   │   ├── annotation.py
│   │   └── exports.py
│   ├── integrations/
│   │   └── pyclustifyr.py
│   ├── resources/
│   │   ├── markers/
│   │   │   ├── base.py
│   │   │   ├── csv.py
│   │   │   └── duckdb.py
│   │   └── references/
│   │       ├── base.py
│   │       ├── local.py
│   │       └── database.py
│   ├── plotting/
│   │   ├── embedding.py
│   │   ├── feature.py
│   │   ├── dotplot.py
│   │   ├── heatmap.py
│   │   └── enrichment.py
│   └── ui/
│       ├── layout.py
│       ├── overview.py
│       ├── annotation.py
│       ├── feature.py
│       ├── dotplot.py
│       ├── modules.py
│       ├── references.py
│       ├── markers.py
│       ├── enrichment.py
│       └── export.py
├── resources/
│   ├── marker_sets/
│   └── reference_matrices/
├── tests/
│   ├── unit/
│   ├── integration/
│   ├── ui/
│   ├── deployment/
│   └── data/
│       ├── h5ad/
│       ├── seurat/
│       ├── sce/
│       └── references/
└── scripts/
    └── create_test_fixtures/
```

Modules may be consolidated while small, but dependency direction must remain UI -> services -> core/adapters, never the reverse.

## 16. Configuration and observability

Define one typed configuration object populated at startup. Include upload limits, workspace root, allowed input types, cache limits, provider implementation/locations, session constraints, log level, and feature flags for non-MVP importers.

Use structured logs with request/session correlation IDs that are random and non-identifying. Record operation names, durations, shapes, sparsity, package versions, warnings, and exception types. Do not record matrix values or sensitive metadata. Add timing around import, aggregation, marker ranking, `pyclustifyr`, and export so performance work is evidence-based.

## 17. Testing strategy

### Unit tests

- Annotation assignment, reset, provenance, materialization, collisions, and serialization.
- Expression-source selection and sparse-safe extraction/aggregation.
- Candidate cluster/embedding detection.
- Gene matching and missing-gene reporting.
- Provider filtering, schema validation, and unavailable-provider behavior.
- Figure builders with small deterministic inputs.

### Import compatibility tests

Maintain small, redistributable fixtures with documented generator versions and expected mappings:

- dense and sparse H5AD;
- H5AD with `raw`, multiple layers, categorical metadata, and multiple embeddings;
- representative Seurat v4/v5 inputs with multiple assays/reductions;
- representative SCE inputs with multiple assays/reduced dimensions;
- malformed, truncated, oversized, and unsupported files.

Assert shapes, identifiers, selected values, sparse/dense representation, metadata categories, embeddings, layers, and warning reports. Golden expectations should focus on semantic mapping rather than unstable serialized bytes.

### Export and round-trip tests

- Write annotated H5AD, reopen it, and compare source content plus added namespaced fields.
- Verify both cluster- and cell-level CSVs with special characters and missing values.
- Verify the uploaded file remains byte-for-byte unchanged.
- Native RDS export, when attempted later, must be opened and semantically validated by the target R package in a separate compatibility test environment.

### Integration tests

- Upload -> configure workspace -> annotate -> CSV/H5AD export.
- Marker provider -> gene panel -> dot/feature plot.
- Reference provider -> `pyclustifyr` -> preview -> apply.
- Marker ranking -> heatmap -> enrichment.
- Concurrent sessions cannot see or modify each other's files or state.

### UI and end-to-end tests

Automate the critical browser workflow against a small fixture: upload, choose cluster/embedding/layer, select a cluster, change a label, render feature/dot plots, and download/validate exports. Test progress, validation, missing genes, invalid uploads, and analysis errors.

### Deployment and nonfunctional tests

- Container health, non-root execution, clean shutdown, and workspace cleanup.
- ShinyProxy lifecycle and Traefik routing/TLS in a staging network.
- Two or more simultaneous sessions with enforced quotas and isolation.
- Benchmarks at representative cell/gene sizes; capture import time, peak memory, plot latency, analysis latency, export time, and browser payload size.
- Dependency, image, and secret scanning in CI.

CI should run formatting/linting, static typing, unit/integration tests, image build, and a minimal container smoke test. Keep slow format-compatibility and large performance tests in scheduled or explicitly triggered jobs.

## 18. Performance strategy

1. Establish representative benchmark datasets and budgets before optimizing.
2. Preserve SciPy sparse matrices and operate on selected genes/clusters only.
3. Precompute compact embedding data and cluster membership indexes after workspace configuration.
4. Cache dot-plot aggregates, module scores, marker results, and reference correlations using complete parameter-aware keys.
5. Use Plotly `scattergl` below a measured threshold; downsample only visualization points above it while preserving annotation semantics.
6. Evaluate Datashader/rasterization for hundreds of thousands of cells when Plotly payload or interaction latency becomes unacceptable.
7. Limit concurrent expensive jobs per session and enforce container memory/CPU quotas.
8. Use backed H5AD or Zarr only after auditing which operations remain safe and performant in lazy mode.
9. Profile CPU, memory, serialization, and network payloads in production-like containers.
10. Introduce a Rust extension through PyO3 only for a stable, benchmarked hotspot such as sparse aggregation, gene matching, parsing, or serialization. Do not use WASM merely to move server-side matrix work into the browser.

Suggested initial performance targets for the agreed MVP benchmark fixture should be set after the fixture is selected. Record them as testable budgets in `tests/performance/README.md`; avoid claiming arbitrary dataset-scale guarantees without hardware and fixture definitions.

## 19. Phased implementation plan

### Phase 0: repository and technical spikes

Deliver:

- Python package, locked dependencies, quality tooling, CI, configuration, and basic container.
- Small canonical AnnData fixtures.
- Verified Shiny for Python + Plotly integration.
- Research spikes documenting the exact current APIs and limitations of `readseurat`, `rdata`/BiocPy, and `pyclustifyr`.
- Import/export and provider protocols with typed models.

Exit gate: a containerized placeholder app runs; CI passes; the compatibility decisions and pinned versions are recorded.

### Phase 1: AnnData-native MVP

Deliver:

- H5AD upload, validation, and import report.
- Cluster-column, embedding, and expression-source discovery/selection.
- Interactive Plotly embedding, feature plot, and dot plot.
- Editable cluster-level annotation table and application to the view.
- Cluster-level and cell-level CSV exports.
- Annotated H5AD export with provenance and re-open validation.
- Session isolation, progress/error states, and core automated tests.

Exit gate: all MVP acceptance criteria in section 20 pass for the reference H5AD fixture.

### Phase 2: SOWhat functional parity

Deliver:

- Module scores and reusable gene lists.
- Differential markers and heatmap.
- GO:BP enrichment.
- `pyclustifyr` integration with local reference matrices, preview, and apply workflow.
- Local CSV marker provider and local reference provider.
- Reset workflow and, if feasible, undo/redo.

Exit gate: every named analysis can be run, reviewed, and exported without mutating the original clustering.

### Phase 3: Seurat import

Deliver:

- Pure-Python Seurat importer for the documented compatibility matrix.
- Assay/layer, cell/feature metadata, and reduction mapping.
- Seurat v4/v5 regression fixtures and actionable unsupported-input errors.
- CSV and annotated H5AD export from imported Seurat objects.

Exit gate: all semantic fixture assertions and the full annotation/export workflow pass for every supported Seurat fixture.

### Phase 4: SingleCellExperiment import

Deliver:

- Tested `rdata`/BiocPy or narrowly customized SCE adapter.
- Assay, `colData`, `rowData`, and `reducedDims` mapping.
- Explicit compatibility matrix and fallback conversion documentation.
- CSV and annotated H5AD export from imported SCE objects.

Exit gate: the full workflow passes for every documented SCE fixture and unsupported structures fail safely.

### Phase 5: production multi-user deployment

Deliver:

- Hardened image and Compose/configuration for ShinyProxy + Traefik.
- Authentication, TLS, quotas, timeouts, health checks, cleanup, and operator guide.
- Concurrency/isolation tests and measured resource profile.

Exit gate: simultaneous staged users receive isolated containers and workspaces; teardown releases memory and deletes ephemeral data; routing and downloads work through Traefik.

### Phase 6: database providers and optional native R export

Deliver independently:

- Read-only DuckDB marker provider.
- Reference-database provider.
- Resource versioning/checksums and UI status.
- Native Seurat/SCE write-back only if versioned semantic round-trip tests prove it safe.

Exit gate: provider contract suites pass for local and database implementations; any native R output opens and validates in its target R ecosystem.

### Phase 7: measured scaling work

Deliver only as justified by benchmarks:

- Smarter caching and job controls.
- Backed/Zarr workflows.
- Datashader or alternative large-point rendering.
- PyO3/Rust acceleration for confirmed hotspots.

## 20. MVP acceptance criteria

The MVP is complete only when all of the following are demonstrated in the production container using a documented reference H5AD fixture:

1. A user can upload a valid H5AD file and receives a summary of cell/gene counts, candidate cluster columns, embeddings, expression sources, and import warnings.
2. Invalid or unsupported files fail without crashing the session and show an actionable message.
3. The user can choose any compatible cluster column, embedding, and `X`/`raw`/layer expression source.
4. The overview renders an interactive WebGL embedding colored by cluster and by current annotation.
5. The user can select a cluster and assign/edit its annotation; all cells in that cluster update consistently in the view.
6. The source clustering column remains unchanged throughout the session and in the uploaded file.
7. A requested existing gene renders in a feature plot; a missing gene produces clear feedback rather than an empty or misleading plot.
8. A multi-gene dot plot correctly represents per-cluster fraction expressing and aggregate expression for the selected expression source.
9. Annotation state retains cluster ID, annotation, source, notes/confidence when supplied, and timestamp.
10. Cluster-level and cell-level CSV outputs contain correct identifiers and annotations and handle nonnumeric cluster IDs and quoted/Unicode labels.
11. Annotated H5AD export adds `clustbuster_cluster`, `clustbuster_annotation`, and namespaced provenance, preserves the original object content, and successfully reopens in AnnData.
12. Export does not alter the uploaded source file.
13. Two concurrent sessions have independent reactive state and cannot access each other's uploads, caches, annotations, or exports.
14. The application runs as a non-root process in Docker, passes its health check, and cleans its ephemeral session workspace when the container is removed.
15. Automated unit, integration, end-to-end critical-path, and container smoke tests pass in CI.
16. The README documents local startup, supported inputs for the MVP (H5AD only), known limitations, privacy/session cleanup, and export behavior without implying native RDS round-trip support.

Seurat import, SCE import, `pyclustifyr`, enrichment, database-backed providers, and ShinyProxy/Traefik production orchestration are planned capabilities but are not allowed to weaken or ambiguously redefine this MVP gate.

## 21. Definition of done for every feature

A feature is done when:

- its behavior and failure modes are represented by typed interfaces;
- core logic is independent of Shiny reactive code where practical;
- tests cover the happy path, missing/invalid inputs, and sparse data where relevant;
- the UI provides progress, validation, and actionable errors;
- no input file is mutated;
- relevant operations are timed without logging sensitive data;
- documentation and configuration examples are updated;
- the production image builds and the existing end-to-end workflow still passes.

## 22. First coding tasks

Codex should begin in this order:

1. Scaffold the package, development tooling, CI, and non-root Docker image.
2. Define `Workspace`, `ExpressionSource`, `AnnotationStore`, import/export result models, and provider protocols.
3. Create tiny dense and sparse H5AD fixtures and expected semantic mappings.
4. Implement H5AD import validation and the import-summary service.
5. Implement annotation-state operations and CSV/H5AD exporters with round-trip tests.
6. Build the minimal Shiny flow: upload, configure, embedding, annotation editor, and downloads.
7. Add sparse-safe feature and dot-plot services plus Plotly figure builders.
8. Add concurrency/session-isolation and container smoke tests.
9. Only after the MVP gate passes, proceed to SOWhat parity and format-compatibility phases.

When implementation details conflict with this plan, preserve the guiding decisions and document the proposed deviation in an architecture decision record before changing a public interface or scope boundary.
