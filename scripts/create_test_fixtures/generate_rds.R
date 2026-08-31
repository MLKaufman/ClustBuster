# Generate genuine Seurat and SingleCellExperiment compatibility fixtures.
# Run this only in an R environment with both target packages installed.

required <- c("Seurat", "SingleCellExperiment", "S4Vectors")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install required R packages before generating RDS fixtures: ", paste(missing, collapse = ", "))
}

set.seed(20260830)
genes <- c("CD3D", "IL7R", "LYZ", "S100A8", "MS4A1", "CD79A", "NKG7", "GNLY")
cells <- sprintf("demo-cell-%03d", seq_len(40))
counts <- matrix(rpois(length(genes) * length(cells), lambda = 1), nrow = length(genes),
                 dimnames = list(genes, cells))
clusters <- factor(rep(0:3, each = 10))
embedding <- cbind(UMAP_1 = rep(c(-2, 2, 0, 0), each = 10) + rnorm(40, sd = 0.3),
                   UMAP_2 = rep(c(0, 0, -2, 2), each = 10) + rnorm(40, sd = 0.3))
rownames(embedding) <- cells

seurat <- Seurat::CreateSeuratObject(counts = counts, project = "ClustBusterFixture")
seurat$clustbuster_cluster <- clusters
seurat[["umap"]] <- Seurat::CreateDimReducObject(
  embeddings = embedding, key = "UMAP_", assay = Seurat::DefaultAssay(seurat)
)

sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts, logcounts = log1p(counts)),
  colData = S4Vectors::DataFrame(clustbuster_cluster = clusters)
)
SingleCellExperiment::reducedDim(sce, "UMAP") <- embedding

root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "..", ".."), mustWork = TRUE)
dir.create(file.path(root, "testdata"), showWarnings = FALSE, recursive = TRUE)
saveRDS(seurat, file.path(root, "testdata", "so.rds"), version = 3)
saveRDS(sce, file.path(root, "testdata", "sce.rds"), version = 3)
message("Wrote testdata/so.rds and testdata/sce.rds")

