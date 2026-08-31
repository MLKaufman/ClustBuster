# Generate a genuine SeuratObject 4.1.4 regression fixture.
# Usage: R_LIBS=/path/to/seurat-v4-library Rscript generate_seurat_v4.R

if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop("Install SeuratObject 4.1.4 into an isolated R library first")
}
installed_version <- as.character(utils::packageVersion("SeuratObject"))
if (installed_version != "4.1.4") {
  stop("Expected SeuratObject 4.1.4, found ", installed_version)
}

set.seed(20260831)
genes <- c("CD3D", "IL7R", "LYZ", "S100A8", "MS4A1", "NKG7")
cells <- sprintf("v4-cell-%03d", seq_len(24))
counts <- matrix(
  rpois(length(genes) * length(cells), lambda = 1),
  nrow = length(genes),
  dimnames = list(genes, cells)
)
clusters <- factor(rep(0:2, each = 8))
embedding <- cbind(
  UMAP_1 = rep(c(-2, 2, 0), each = 8) + rnorm(24, sd = 0.25),
  UMAP_2 = rep(c(0, 0, 2), each = 8) + rnorm(24, sd = 0.25)
)
rownames(embedding) <- cells

seurat <- SeuratObject::CreateSeuratObject(
  counts = counts,
  project = "ClustBusterV4Fixture"
)
seurat$seurat_clusters <- clusters
seurat[["ALT"]] <- SeuratObject::CreateAssayObject(counts = counts * 2)
protein_counts <- matrix(
  rpois(3 * length(cells), lambda = 2),
  nrow = 3,
  dimnames = list(c("CD3", "CD14", "CD19"), cells)
)
seurat[["ADT"]] <- SeuratObject::CreateAssayObject(counts = protein_counts)
seurat[["umap"]] <- SeuratObject::CreateDimReducObject(
  embeddings = embedding,
  key = "UMAP_",
  assay = SeuratObject::DefaultAssay(seurat)
)

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(file_arg)) stop("Run this generator with Rscript")
script_path <- sub("^--file=", "", file_arg[[1]])
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
output <- file.path(root, "testdata", "so-v4.rds")
saveRDS(seurat, output, version = 3)
message("Wrote testdata/so-v4.rds with SeuratObject ", installed_version)
