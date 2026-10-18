# =============================================================================
# tf.R | DoRothEA + VIPER transcription-factor activity analysis (human)
# Input : data/af1_Figure7_FINAL_reannotated_noNA.rds
#         Seurat metadata must contain 'Celltype_revised'.
# Output: results/dorothea/
# Run   : Rscript scripts/tf.R  (from the repository root)
# =============================================================================

set.seed(1234)
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dorothea)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
})

# 1. Configuration ------------------------------------------------------------
input_file <- file.path("data", "af1_Figure7_FINAL_reannotated_noNA.rds")
out_dir <- file.path("results", "dorothea")
celltype_column <- "Celltype_revised"
input_assay <- "RNA"
regulon_confidence <- c("A", "B", "C")
n_top_tfs <- 20
max_dense_elements <- 5e7    # Safety limit; as.matrix() needs ~8 bytes/value
save_full_seurat <- TRUE      # Large RDS output

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Input and Seurat v5 layer handling ---------------------------------------
if (!file.exists(input_file)) stop("Missing Seurat object: ", input_file)
sce <- readRDS(input_file)
if (!celltype_column %in% names(sce@meta.data)) {
  stop("Missing required metadata: ", celltype_column)
}
if (!input_assay %in% Assays(sce)) stop("Missing assay: ", input_assay)
DefaultAssay(sce) <- input_assay
layer_names <- Layers(sce[[input_assay]])
if (length(grep("^data", layer_names)) > 1L ||
    length(grep("^counts", layer_names)) > 1L) {
  sce <- JoinLayers(sce, assay = input_assay)
}
if (!"data" %in% Layers(sce[[input_assay]])) {
  sce <- NormalizeData(sce, assay = input_assay, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)
}
expr <- LayerData(sce, assay = input_assay, layer = "data")
if (prod(dim(expr)) > max_dense_elements) {
  stop("The normalized matrix is too large for safe dense VIPER conversion. ",
       "Subset cell types/cells or raise max_dense_elements if RAM allows.")
}
expr <- as.matrix(expr)

# 3. Regulons and inferred TF activity ---------------------------------------
# Use data() to load the dataset; get(data(...)) is not a valid loading pattern.
data("dorothea_hs", package = "dorothea")
regulon <- dorothea_hs %>%
  filter(confidence %in% regulon_confidence)
if (!nrow(regulon)) stop("The selected DoRothEA regulon is empty.")

scores <- dorothea::run_viper(
  input = expr,
  regulons = regulon,
  options = list(method = "scale", minsize = 4,
                 eset.filter = FALSE, cores = 1, verbose = FALSE)
)
scores <- as.matrix(scores)

# run_viper(matrix) returns TFs x cells. Verify and align explicitly.
cells_in_rows <- sum(rownames(scores) %in% colnames(sce))
cells_in_cols <- sum(colnames(scores) %in% colnames(sce))
if (cells_in_rows > cells_in_cols) scores <- t(scores)
shared_cells <- intersect(colnames(sce), colnames(scores))
if (!length(shared_cells)) stop("VIPER score matrix has no matching cell IDs.")
sce <- subset(sce, cells = shared_cells)
scores <- scores[, colnames(sce), drop = FALSE]
stopifnot(identical(colnames(scores), colnames(sce)))
if (nrow(scores) < 3L) stop("Not enough inferred TFs for PCA/heatmap analysis.")

# 4. Seurat TF assay, PCA, clustering and t-SNE -------------------------------
sce[["dorothea"]] <- CreateAssayObject(data = scores)
DefaultAssay(sce) <- "dorothea"
sce <- ScaleData(sce, assay = "dorothea", features = rownames(scores),
                 verbose = FALSE)
sce <- RunPCA(sce, assay = "dorothea", features = rownames(scores),
              npcs = min(30L, nrow(scores) - 1L), verbose = FALSE)
n_pcs <- min(10L, ncol(Embeddings(sce, reduction = "pca")))
if (n_pcs < 2L) stop("Too few PCs for TF activity clustering.")
sce <- FindNeighbors(sce, reduction = "pca", dims = seq_len(n_pcs),
                     verbose = FALSE)
sce <- FindClusters(sce, resolution = 0.5, verbose = FALSE)
sce <- RunTSNE(sce, reduction = "pca", dims = seq_len(n_pcs),
               check_duplicates = FALSE, verbose = FALSE)
Idents(sce) <- celltype_column
p <- DimPlot(sce, reduction = "tsne", label = TRUE, pt.size = 0.5) +
  NoLegend()
ggsave(file.path(out_dir, "tf_activity_tsne.pdf"), p, width = 9, height = 7)

# 5. Compute TF activity summaries without a giant long-format data frame ----
scaled <- as.matrix(LayerData(sce, assay = "dorothea", layer = "scale.data"))
labels <- as.character(sce@meta.data[colnames(scaled), celltype_column])
if (anyNA(labels) || any(!nzchar(labels))) {
  stop("Cell-type annotations contain missing or blank values.")
}
celltypes <- unique(labels)
if (length(celltypes) < 2L) {
  stop("At least two cell types are needed to rank variable TF activities.")
}

mean_mat <- vapply(celltypes, function(g) {
  rowMeans(scaled[, labels == g, drop = FALSE], na.rm = TRUE)
}, numeric(nrow(scaled)))
sd_mat <- vapply(celltypes, function(g) {
  apply(scaled[, labels == g, drop = FALSE], 1, sd, na.rm = TRUE)
}, numeric(nrow(scaled)))
rownames(mean_mat) <- rownames(scaled)
rownames(sd_mat) <- rownames(scaled)
colnames(mean_mat) <- celltypes
colnames(sd_mat) <- celltypes

summary_df <- data.frame(
  tf = rep(rownames(mean_mat), times = ncol(mean_mat)),
  cell_type = rep(colnames(mean_mat), each = nrow(mean_mat)),
  avg = as.vector(mean_mat),
  std = as.vector(sd_mat),
  stringsAsFactors = FALSE
)
tf_variance <- apply(mean_mat, 1, var, na.rm = TRUE)
tf_variance <- sort(tf_variance[is.finite(tf_variance)], decreasing = TRUE)
if (!length(tf_variance)) stop("No finite TF variance across cell types.")
top_tfs <- head(names(tf_variance), n_top_tfs)
top_table <- data.frame(tf = top_tfs,
                        variance = as.numeric(tf_variance[top_tfs]))

# 6. Save tables, heatmap and Seurat object ----------------------------------
write.csv(t(mean_mat), file.path(out_dir, "tf.csv"), row.names = TRUE)
write.csv(top_table, file.path(out_dir, "top20_variable_TFs.csv"),
          row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "all_TF_activity_by_celltype.csv"),
          row.names = FALSE)

heatmap_mat <- mean_mat[top_tfs, , drop = FALSE]
heatmap_mat <- heatmap_mat[apply(heatmap_mat, 1, function(x) all(is.finite(x))),
                           , drop = FALSE]
if (!nrow(heatmap_mat)) stop("No finite TF activity rows to plot.")
max_abs <- suppressWarnings(max(abs(heatmap_mat), na.rm = TRUE))
if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
palette <- grDevices::colorRampPalette(c("Darkblue", "white", "red"))(100)
pheatmap::pheatmap(
  heatmap_mat,
  color = palette,
  breaks = seq(-max_abs, max_abs, length.out = 101),
  main = "DoRothEA (ABC) TF Activity",
  fontsize_row = 8, fontsize_col = 10, angle_col = 45,
  cluster_rows = nrow(heatmap_mat) >= 2L,
  cluster_cols = ncol(heatmap_mat) >= 2L,
  border_color = NA,
  filename = file.path(out_dir, "tf_activity_heatmap.pdf"),
  width = 11, height = 8
)
if (save_full_seurat) {
  saveRDS(sce, file.path(out_dir, "pbmc_dorothea_seuratV5.rds"))
}
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("DoRothEA/VIPER analysis complete: ", normalizePath(out_dir))
