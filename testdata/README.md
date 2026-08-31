# Manual test data

This directory contains small, synthetic objects for exercising ClustBuster without
using biological or identifying data.

## Available now

- `clustbuster-demo.h5ad` — 90 cells, 10 genes, five clusters, sparse `X`, `raw`,
  `counts` and `log1p` layers, UMAP/PCA embeddings, and categorical cell metadata.
- `so.rds` — genuine Seurat v5 object produced with SeuratObject 5.4.0: 40 cells,
  eight genes, sparse counts, four `seurat_clusters`, cell metadata, a UMAP
  reduction, a compatible secondary assay, and an intentionally incompatible ADT
  assay for mapping/reporting tests.
- `so-v4.rds` — genuine legacy Seurat object produced with SeuratObject 4.1.4:
  24 cells, six genes, sparse counts, three clusters, UMAP, a compatible ALT
  assay, and an intentionally incompatible ADT assay.
- `so.h5seurat` — deterministic H5Seurat compatibility fixture produced by the
  pinned `readseurat` writer: 12 cells, five genes, sparse expression and counts,
  three `seurat_clusters`, and a UMAP reduction.
- `sce.rds` — genuine SingleCellExperiment 1.34.0 object produced with
  Bioconductor 3.23: 40 cells, eight genes, sparse counts, dense log-counts,
  four `sce_clusters` groups, a UMAP reduced dimension, and an explicitly reported
  ADT alternative experiment.

Regenerate the fixtures with:

```bash
uv run python scripts/create_test_fixtures/generate_h5ad.py
uv run python scripts/create_test_fixtures/generate_h5seurat.py
Rscript scripts/create_test_fixtures/generate_rds.R seurat
Rscript scripts/create_test_fixtures/generate_rds.R sce
```

The v4 fixture requires an isolated SeuratObject 4.1.4 library so it cannot be
accidentally regenerated with the v5 object model:

```bash
R_LIBS=/path/to/seurat-v4-library Rscript scripts/create_test_fixtures/generate_seurat_v4.R
```

The R generator accepts `seurat`, `sce`, or `all` and refuses to create substitute
files when the corresponding real R package is absent.
