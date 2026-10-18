# =============================================================================
# main.R | Sepsis scRNA-seq preprocessing, clustering and manual annotation
# Input : data/output_folder/<sample>/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}
# Output: results/main/ and data/af.rds
# Run   : Rscript scripts/main.R  (from the repository root)
# =============================================================================

set.seed(1234)
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(harmony)
})

# 1. Configuration ------------------------------------------------------------
data_dir <- file.path("data", "output_folder")
annotation_file <- file.path("config", "cluster_annotations.csv")
sample_metadata_file <- file.path("config", "sample_metadata.csv")
out_dir <- file.path("results", "main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

qc_min_counts <- 1000
qc_max_counts <- 40000
qc_min_features <- 200
qc_max_features <- 5000
qc_max_mito <- 25                         # Percent, not fraction
n_variable_features <- 2000
n_pcs <- 9
selected_resolution <- 1.0
cluster_resolutions <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8,
                         1, 1.2, 1.5, 2, 2.5, 3)
excluded_celltypes <- c("Ribo-high / Low-info Cells", "Doublets")

save_plot <- function(plot, name, width = 8, height = 6) {
  ggsave(file.path(out_dir, name), plot = plot, width = width, height = height)
}

# 2. Read sample-specific Cell Ranger matrices --------------------------------
if (!dir.exists(data_dir)) {
  stop("Input folder missing: ", data_dir)
}
sample_dirs <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)
sample_names <- basename(sample_dirs)
if (!length(sample_dirs)) stop("No sample subdirectories found in ", data_dir)
if (anyDuplicated(sample_names)) stop("Sample directory names must be unique.")

read_sample <- function(path, sample_id) {
  message("Reading sample: ", sample_id)
  x <- Read10X(data.dir = path)
  if (is.list(x)) {
    if (!"Gene Expression" %in% names(x)) {
      stop("No 'Gene Expression' matrix in sample: ", sample_id)
    }
    x <- x[["Gene Expression"]]
  }
  obj <- CreateSeuratObject(counts = x, project = sample_id,
                            min.cells = 3, min.features = 200)
  obj$Sample <- sample_id
  obj
}
objs <- Map(read_sample, sample_dirs, sample_names)
af <- if (length(objs) == 1L) {
  objs[[1L]]
} else {
  merge(objs[[1L]], y = objs[-1L], add.cell.ids = sample_names)
}
rm(objs)
DefaultAssay(af) <- "RNA"

# Merged Seurat v5 objects can contain counts.* layers; join before QC / DEG.
if (length(grep("^counts", Layers(af[["RNA"]]))) > 1L) {
  af <- JoinLayers(af, assay = "RNA")
}

# 3. Quality control ----------------------------------------------------------
af[["percent.MT"]] <- PercentageFeatureSet(af, pattern = "^MT-")
af[["percent.rb"]] <- PercentageFeatureSet(af, pattern = "^RP[SL]")
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.MT", "percent.rb")
save_plot(VlnPlot(af, features = qc_features, ncol = 4),
          "01_qc_before.pdf", 12, 5)

qc_scatter <- FeatureScatter(af, "nCount_RNA", "percent.MT") +
  FeatureScatter(af, "nCount_RNA", "percent.rb") +
  FeatureScatter(af, "nCount_RNA", "nFeature_RNA")
save_plot(qc_scatter, "01_qc_correlations_before.pdf", 12, 5)

keep <- af$nCount_RNA >= qc_min_counts & af$nCount_RNA <= qc_max_counts &
  af$nFeature_RNA >= qc_min_features & af$nFeature_RNA <= qc_max_features &
  af$percent.MT <= qc_max_mito
if (!any(keep)) stop("All cells failed QC; review thresholds and input data.")
af <- af[, keep]
save_plot(VlnPlot(af, features = qc_features, ncol = 4),
          "01_qc_after.pdf", 12, 5)

# 4. Normalize, reduce dimensions, remove sample effects ---------------------
af <- NormalizeData(af, normalization.method = "LogNormalize",
                    scale.factor = 10000)
af <- FindVariableFeatures(af, selection.method = "vst",
                           nfeatures = n_variable_features)
af <- ScaleData(af, verbose = FALSE)
af <- RunPCA(af, features = VariableFeatures(af), verbose = FALSE)

p_var <- VariableFeaturePlot(af)
top_genes <- head(VariableFeatures(af), 10)
save_plot(LabelPoints(plot = p_var, points = top_genes, repel = TRUE),
          "02_variable_features.pdf", 7, 6)
save_plot(DimPlot(af, reduction = "pca", group.by = "Sample"),
          "02_pca_before_harmony.pdf", 8, 6)
n_pcs_initial <- min(20L, ncol(Embeddings(af, reduction = "pca")))
if (n_pcs_initial < 2L) stop("Insufficient principal components for t-SNE.")
af <- RunTSNE(af, reduction = "pca", dims = seq_len(n_pcs_initial))
save_plot(DimPlot(af, reduction = "tsne", group.by = "Sample"),
          "02_tsne_before_harmony.pdf", 8, 6)

if (length(sample_names) > 1L) {
  af <- RunHarmony(af, group.by.vars = "Sample")
  reduction_used <- "harmony"
  save_plot(ElbowPlot(af, ndims = 30, reduction = "harmony"),
            "03_harmony_elbow.pdf", 6, 5)
  save_plot(DimPlot(af, reduction = "harmony", group.by = "Sample"),
            "03_harmony_embedding.pdf", 8, 6)
} else {
  reduction_used <- "pca"
  message("Only one sample: Harmony correction is not needed; using PCA.")
}
save_plot(ElbowPlot(af, ndims = 30, reduction = "pca"),
          "03_pca_elbow.pdf", 6, 5)

# 5. Cluster at several resolutions, then choose the original resolution -----
if (ncol(Embeddings(af, reduction = reduction_used)) < n_pcs) {
  stop("Not enough PCs; reduce n_pcs.")
}
af <- FindNeighbors(af, reduction = reduction_used, dims = seq_len(n_pcs))
for (r in cluster_resolutions) {
  af <- FindClusters(af, graph.name = "RNA_snn", resolution = r,
                     verbose = FALSE)
}
if (requireNamespace("clustree", quietly = TRUE)) {
  save_plot(clustree::clustree(af@meta.data, prefix = "RNA_snn_res."),
            "04_resolution_clustree.pdf", 12, 10)
}
af <- FindClusters(af, graph.name = "RNA_snn",
                   resolution = selected_resolution, verbose = FALSE)
af <- RunUMAP(af, reduction = reduction_used, dims = seq_len(n_pcs))
af <- RunTSNE(af, reduction = reduction_used, dims = seq_len(n_pcs),
              check_duplicates = FALSE)
save_plot(DimPlot(af, reduction = "umap", label = TRUE),
          "04_clusters_umap.pdf")
save_plot(DimPlot(af, reduction = "tsne", label = TRUE),
          "04_clusters_tsne.pdf")
saveRDS(af, file.path(out_dir, "af_before_annotation.rds"))

# 6. Cluster markers and manual annotation -----------------------------------
Idents(af) <- "seurat_clusters"
cluster_markers <- FindAllMarkers(af, test.use = "roc", return.thresh = 0.25,
                                  logfc.threshold = 0.5, min.pct = 0.3,
                                  only.pos = TRUE)
write.csv(cluster_markers, file.path(out_dir, "05_cluster_markers.csv"),
          row.names = FALSE)
if ("power" %in% names(cluster_markers)) {
  top_markers <- cluster_markers %>%
    group_by(cluster) %>%
    slice_max(order_by = power, n = 20, with_ties = FALSE) %>%
    ungroup()
  write.csv(top_markers, file.path(out_dir, "05_top_cluster_markers.csv"),
            row.names = FALSE)
}

if (!file.exists(annotation_file)) {
  stop("Annotation CSV not found: ", annotation_file,
       ". Review your cluster markers and create the mapping file.")
}
annotation <- read.csv(annotation_file, stringsAsFactors = FALSE,
                       fileEncoding = "UTF-8")
if (!all(c("cluster", "cell_type") %in% names(annotation))) {
  stop("Annotation CSV needs 'cluster' and 'cell_type' columns.")
}
if (anyDuplicated(annotation$cluster) || anyNA(annotation$cell_type) ||
    any(!nzchar(annotation$cell_type))) {
  stop("Cluster annotations must be unique and non-empty.")
}
actual_clusters <- sort(unique(as.character(af$seurat_clusters)))
if (!setequal(actual_clusters, as.character(annotation$cluster))) {
  stop("Cluster labels do not match config/cluster_annotations.csv. Found: ",
       paste(actual_clusters, collapse = ", "),
       ". Please review markers and update the annotation mapping.")
}
annotation_map <- setNames(annotation$cell_type, as.character(annotation$cluster))
af$Celltype <- unname(annotation_map[as.character(af$seurat_clusters)])
Idents(af) <- "Celltype"
af1 <- subset(af, idents = excluded_celltypes, invert = TRUE)
if (!ncol(af1)) stop("No cells remain after excluding unwanted types.")

save_plot(DimPlot(af1, reduction = "umap", group.by = "Celltype",
                  label = TRUE, repel = TRUE),
          "06_annotated_umap.pdf", 10, 7)
if ("PHLDA3" %in% rownames(af1)) {
  save_plot(DotPlot(af1, features = "PHLDA3") + RotatedAxis(),
            "06_PHLDA3.pdf", 6, 5)
}

# 7. Optional group-level proportion comparison ------------------------------
# Metadata CSV must map actual sample folder names to HC / Sepsis conditions.
if (file.exists(sample_metadata_file)) {
  sample_meta <- read.csv(sample_metadata_file, stringsAsFactors = FALSE)
  if (!all(c("sample", "condition") %in% names(sample_meta))) {
    stop("sample_metadata.csv needs 'sample' and 'condition' columns.")
  }
  if (anyDuplicated(sample_meta$sample) ||
      anyNA(sample_meta$condition) || any(!nzchar(sample_meta$condition)) ||
      !all(unique(af1$Sample) %in% sample_meta$sample)) {
    stop("sample_metadata.csv must uniquely map every input sample.")
  }
  sample_to_condition <- setNames(sample_meta$condition, sample_meta$sample)
  af1$Condition <- unname(sample_to_condition[af1$Sample])
  cell_ratio <- prop.table(table(af1$Celltype, af1$Condition), margin = 2)
  cell_ratio <- as.data.frame(cell_ratio)
  names(cell_ratio) <- c("Celltype", "Condition", "Ratio")
  save_plot(ggplot(cell_ratio, aes(x = Condition, y = Ratio, fill = Celltype)) +
              geom_col(width = 0.75) + coord_flip() + theme_classic(),
            "07_celltype_proportions.pdf", 9, 6)
  write.csv(cell_ratio, file.path(out_dir, "07_celltype_proportions.csv"),
            row.names = FALSE)
} else {
  message("No config/sample_metadata.csv: skipping condition-level proportions.")
}

# 8. Cell-type markers and expression visualizations -------------------------
Idents(af1) <- "Celltype"
cell_markers <- FindAllMarkers(af1, test.use = "roc", return.thresh = 0.25,
                               logfc.threshold = 0.5, min.pct = 0.3,
                               only.pos = TRUE)
write.csv(cell_markers, file.path(out_dir, "08_celltype_markers.csv"),
          row.names = FALSE)

marker_list <- list(
  "B Cells" = c("CD79A", "IGHM"),
  "B Cells (Naïve)" = c("CD79A", "IGHM"),
  "Plasma Cells" = c("JCHAIN", "MZB1"),
  "T Cells (General)" = c("CD3D", "TRAC"),
  "Memory CD4+ T" = c("IL7R", "CCR7"),
  "Effector CD8+ T" = c("CD3E", "CCL5"),
  "NK Cells" = c("NKG7", "GNLY"),
  "CD14+ Monocytes" = c("LYZ", "FCN1"),
  "Inflammatory Monocytes" = c("S100A8", "S100A12"),
  "FCGR3A+ Monocytes" = c("FCGR3A", "LST1"),
  "Activated Myeloid Cells" = c("CYBB", "TYROBP"),
  "Myeloid" = c("LYZ", "TYROBP"),
  "cDC-like APC" = c("HLA-DRA", "CD74"),
  "Neutrophils" = c("CSF3R", "CD177"),
  "Platelets" = c("PPBP", "PF4"),
  "Erythroid Precursors" = c("HBA1", "HBB"),
  "Proliferating Cells" = c("MKI67", "STMN1")
)
features_use <- intersect(unique(unlist(marker_list)), rownames(af1))
if (length(features_use)) {
  p_dot <- DotPlot(af1, features = features_use, group.by = "Celltype") +
    RotatedAxis() + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_dot, "08_marker_dotplot.pdf", 13, 7)
}

# ROC marker outputs differ from Wilcoxon outputs. Label the actual effect metric.
marker_pairs <- data.frame(
  cluster = rep(names(marker_list), lengths(marker_list)),
  gene = unlist(marker_list, use.names = FALSE), stringsAsFactors = FALSE
)
if (nrow(cell_markers) > 0 && all(c("cluster", "gene", "pct.1") %in% names(cell_markers))) {
  effect_col <- intersect(c("avg_log2FC", "avg_diff", "myAUC", "power"),
                          names(cell_markers))
  if (length(effect_col)) {
    bubble_data <- dplyr::inner_join(cell_markers, marker_pairs,
                                    by = c("cluster", "gene"))
    if (nrow(bubble_data)) {
      p_bubble <- ggplot(bubble_data,
                         aes(x = gene, y = cluster, size = .data[["pct.1"]],
                             color = .data[[effect_col[1]]])) +
        geom_point(alpha = 0.85) + theme_bw() +
        labs(x = "Marker gene", y = "Cell type", size = "pct.1",
             color = effect_col[1]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      save_plot(p_bubble, "08_marker_bubble.pdf", 11, 7)
    }
  }
}

saveRDS(af1, file.path("data", "af.rds"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Finished. Outputs: ", normalizePath(out_dir), " and data/af.rds")
