# Manual test data

This directory contains small, synthetic objects for exercising ClustBuster without
using biological or identifying data.

## Available now

- `clustbuster-demo.h5ad` — 90 cells, 10 genes, five clusters, sparse `X`, `raw`,
  `counts` and `log1p` layers, UMAP/PCA embeddings, and categorical cell metadata.

Regenerate it with:

```bash
uv run python scripts/create_test_fixtures/generate_h5ad.py
```

## Planned compatibility fixtures

- `so.rds` — genuine Seurat object; added when the Phase 3 Seurat adapter is built.
- `sce.rds` — genuine SingleCellExperiment object; added when the Phase 4 SCE adapter
  is built.

The reproducible R generator already lives at
`scripts/create_test_fixtures/generate_rds.R`. It deliberately refuses to create
substitute files when the real Seurat and SingleCellExperiment packages are absent.
This machine currently has R but not those packages, so no misleading placeholder RDS
files are committed.

