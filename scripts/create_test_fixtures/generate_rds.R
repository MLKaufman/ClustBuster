# Generate genuine Seurat and/or SingleCellExperiment compatibility fixtures.
# Usage: Rscript generate_rds.R [all|seurat|sce]

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args)) tolower(args[[1]]) else "all"
if (!target %in% c("all", "seurat", "sce")) {
  stop("Target must be one of: all, seurat, sce")
}

required <- character()
if (target %in% c("all", "seurat")) required <- c(required, "SeuratObject")
if (target %in% c("all", "sce")) {
  required <- c(required, "SingleCellExperiment", "S4Vectors")
}
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

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(file_arg)) stop("Run this generator with Rscript")
script_path <- sub("^--file=", "", file_arg[[1]])
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
dir.create(file.path(root, "testdata"), showWarnings = FALSE, recursive = TRUE)

if (target %in% c("all", "seurat")) {
  seurat <- SeuratObject::CreateSeuratObject(
    counts = counts, project = "ClustBusterFixture"
  )
  seurat$seurat_clusters <- clusters
  seurat[["ALT"]] <- SeuratObject::CreateAssay5Object(counts = counts * 2)
  protein_counts <- matrix(
    rpois(3 * length(cells), lambda = 2),
    nrow = 3,
    dimnames = list(c("CD3", "CD14", "CD19"), cells)
  )
  seurat[["ADT"]] <- SeuratObject::CreateAssay5Object(counts = protein_counts)
  seurat[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(seurat)
  )
  saveRDS(seurat, file.path(root, "testdata", "so.rds"), version = 3)
  message("Wrote testdata/so.rds")
}

if (target %in% c("all", "sce")) {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts, logcounts = log1p(counts)),
    colData = S4Vectors::DataFrame(clustbuster_cluster = clusters)
  )
  SingleCellExperiment::reducedDim(sce, "UMAP") <- embedding
  saveRDS(sce, file.path(root, "testdata", "sce.rds"), version = 3)
  message("Wrote testdata/sce.rds")
}
