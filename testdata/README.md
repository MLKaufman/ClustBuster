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

Regenerate the fixtures with:

```bash
uv run python scripts/create_test_fixtures/generate_h5ad.py
Rscript scripts/create_test_fixtures/generate_rds.R seurat
```

## Planned compatibility fixtures

- `sce.rds` — genuine SingleCellExperiment object; added when the Phase 4 SCE adapter
  is built.

The R generator accepts `seurat`, `sce`, or `all` and refuses to create substitute
files when the corresponding real R package is absent.
