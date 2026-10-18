# =============================================================================
# monocle_main.R | Monocle 2 pseudotime analysis (DDRTree)
# Input : data/af_Activated_Macrophages.rds (Seurat; raw RNA counts)
# Output: results/monocle2/
# Run   : Rscript scripts/monocle_main.R  (from the repository root)
# IMPORTANT: This is Monocle 2 ('monocle'), NOT Monocle 3 ('monocle3').
# =============================================================================

set.seed(1234)
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(monocle)
  library(Biobase)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# 1. Configuration ------------------------------------------------------------
input_file <- file.path("data", "af_Activated_Macrophages.rds")
out_dir <- file.path("results", "monocle2")
key_genes <- c("M6PR", "BLOC1S1", "SORT1", "KIF1B", "CTSO")
root_state <- NULL              # Optional: specify a biologically justified state
max_heatmap_genes <- 100
make_full_heatmap <- FALSE      # Expensive on large datasets
num_heatmap_clusters <- 4
num_cores <- 1
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, filename, w = 10, h = 8) {
  ggsave(file.path(out_dir, filename), plot = plot, width = w, height = h)
}
save_pheatmap <- function(ht, filename, w = 8, h = 10) {
  grDevices::pdf(file.path(out_dir, filename), width = w, height = h)
  grid::grid.newpage()
  grid::grid.draw(ht$gtable)
  grDevices::dev.off()
}

# 2. Build a Monocle 2 CellDataSet from Seurat v5 raw counts -------------------
if (!file.exists(input_file)) stop("Missing Seurat object: ", input_file)
sce <- readRDS(input_file)
DefaultAssay(sce) <- "RNA"
if (length(grep("^counts", Layers(sce[["RNA"]]))) > 1L) {
  sce <- JoinLayers(sce, assay = "RNA")
}
if (!"counts" %in% Layers(sce[["RNA"]])) {
  stop("RNA counts layer is required for Monocle 2.")
}
counts <- LayerData(sce, assay = "RNA", layer = "counts")
if (is.null(rownames(counts)) || is.null(colnames(counts))) {
  stop("Raw count matrix has no gene/cell names.")
}
meta <- sce@meta.data[colnames(counts), , drop = FALSE]
stopifnot(identical(rownames(meta), colnames(counts)))
gene_meta <- data.frame(gene_short_name = rownames(counts),
                        row.names = rownames(counts))

cds <- newCellDataSet(counts,
                      phenoData = new("AnnotatedDataFrame", data = meta),
                      featureData = new("AnnotatedDataFrame", data = gene_meta),
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

# 3. Ordering genes, trajectory and pseudotime -------------------------------
disp <- dispersionTable(cds)
ordering_genes <- subset(disp, mean_expression >= 0.1 &
                         dispersion_empirical >= dispersion_fit)$gene_id
ordering_genes <- intersect(ordering_genes, rownames(cds))
if (length(ordering_genes) < 3L) {
  stop("Too few ordering genes; inspect dispersion-based cutoffs.")
}
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- if (is.null(root_state)) orderCells(cds) else {
  orderCells(cds, root_state = root_state)
}

# State/Pseudotime are created by orderCells(). Other fields require metadata.
color_fields <- c("seurat_clusters", "State", "Celltype", "Pseudotime")
for (field in color_fields) {
  if (!field %in% colnames(pData(cds))) {
    message("Missing metadata column; skipped trajectory: ", field)
    next
  }
  p <- plot_cell_trajectory(cds, color_by = field) + theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12))
  save_plot(p, paste0("trajectory_", field, ".pdf"), 11, 8)
}

if ("Celltype" %in% colnames(pData(cds))) {
  plot_df <- pData(cds)
  p_density <- ggplot(plot_df, aes(x = Pseudotime,
                                   colour = Celltype, fill = Celltype)) +
    geom_density(bw = 0.5, linewidth = 1, alpha = 0.25) +
    theme_classic()
  save_plot(p_density, "celltype_pseudotime_density.pdf", 12, 7)
}
write.csv(pData(cds), file.path(out_dir, "pseudotime_metadata.csv"),
          row.names = TRUE)
saveRDS(cds, file.path(out_dir, "monocle2_cds.rds"))

# 4. Key gene expression along pseudotime ------------------------------------
key_genes_present <- intersect(key_genes, rownames(cds))
missing_genes <- setdiff(key_genes, key_genes_present)
if (length(missing_genes)) {
  message("Missing key genes: ", paste(missing_genes, collapse = ", "))
}
if (length(key_genes_present)) {
  cds_key <- cds[key_genes_present, ]
  panels <- list()
  for (field in c("State", "Celltype", "Pseudotime")) {
    if (field %in% colnames(pData(cds))) {
      panels[[field]] <- plot_genes_in_pseudotime(cds_key, color_by = field)
    }
  }
  if (length(panels)) {
    p_gene <- patchwork::wrap_plots(panels, ncol = length(panels))
    save_plot(p_gene, "key_genes_pseudotime.pdf", 16, 8)
  }
}

# 5. Test trajectory-associated genes ----------------------------------------
diff_test <- differentialGeneTest(cds[ordering_genes, ],
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                  cores = num_cores)
diff_test$gene_id <- rownames(diff_test)
diff_test <- diff_test %>% arrange(qval)
write.csv(diff_test, file.path(out_dir, "pseudotime_differential_genes.csv"),
          row.names = FALSE)

ranked_genes <- diff_test$gene_id[!is.na(diff_test$qval)]
ranked_genes <- intersect(ranked_genes, rownames(cds))
top_genes <- head(ranked_genes, max_heatmap_genes)
if (length(top_genes) >= 2L) {
  k_top <- min(num_heatmap_clusters, length(top_genes))
  ht <- plot_pseudotime_heatmap(cds[top_genes, ],
                                num_clusters = k_top,
                                show_rownames = TRUE,
                                return_heatmap = TRUE,
                                cores = num_cores)
  save_pheatmap(ht, "top_pseudotime_genes_heatmap.pdf", 8, 11)
  if (!is.null(ht$tree_row)) {
    clusters <- stats::cutree(ht$tree_row, k = k_top)
    write.csv(data.frame(gene = names(clusters), cluster = clusters),
              file.path(out_dir, "top_pseudotime_gene_clusters.csv"),
              row.names = FALSE)
  }
}

if (make_full_heatmap && length(ranked_genes) >= 2L) {
  ht_all <- plot_pseudotime_heatmap(cds[ranked_genes, ],
                                    num_clusters = min(num_heatmap_clusters,
                                                       length(ranked_genes)),
                                    show_rownames = FALSE,
                                    return_heatmap = TRUE,
                                    cores = num_cores)
  save_pheatmap(ht_all, "all_pseudotime_genes_heatmap.pdf", 9, 12)
}

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Monocle 2 analysis complete: ", normalizePath(out_dir))
