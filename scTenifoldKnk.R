# =============================================================================
# scTenifoldKnk.R | Independent in-silico hub-gene knockouts
# Input : data/af_Activated_Macrophages.rds (Seurat; Celltype metadata)
# Output: results/scTenifoldKnk/
# Run   : Rscript scripts/scTenifoldKnk.R  (from the repository root)
# Note  : Original script actually enables 3 hub genes, not 5.
# =============================================================================

set.seed(1234)
suppressPackageStartupMessages({
  library(scTenifoldKnk)
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# 1. Configuration ------------------------------------------------------------
input_file <- file.path("data", "af_Activated_Macrophages.rds")
out_dir <- file.path("results", "scTenifoldKnk")
target_celltype <- "FCGR3A+ macrophages"
hub_genes <- c("KIF1B", "CTSO", "BLOC1S1")
# To reproduce the originally described 5-gene design, validate and add:
# hub_genes <- c("SORT1", "M6PR", "KIF1B", "CTSO", "BLOC1S1")
n_hvgs <- 2000
z_cutoff <- 2
adjusted_p_cutoff <- 0.01
base_seed <- 1234
n_networks <- 10
max_cells_per_network <- 500
save_combined_objects <- TRUE       # Can consume considerable RAM

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Input cell subset --------------------------------------------------------
if (!file.exists(input_file)) stop("Missing Seurat object: ", input_file)
sc <- readRDS(input_file)
if (!"Celltype" %in% names(sc@meta.data)) {
  stop("Expected Seurat metadata column 'Celltype'.")
}
Idents(sc) <- "Celltype"
if (!target_celltype %in% levels(Idents(sc))) {
  stop("Target cell type not found: '", target_celltype,
       "'. Available: ", paste(levels(Idents(sc)), collapse = ", "))
}
sc <- subset(sc, idents = target_celltype)
DefaultAssay(sc) <- "RNA"
if (length(grep("^counts", Layers(sc[["RNA"]]))) > 1L) {
  sc <- JoinLayers(sc, assay = "RNA")
}
if (!"counts" %in% Layers(sc[["RNA"]])) stop("Raw RNA counts layer missing.")
counts_full <- LayerData(sc, assay = "RNA", layer = "counts")

# 3. HVGs plus explicitly included knockout targets --------------------------
found_targets <- intersect(hub_genes, rownames(counts_full))
missing_targets <- setdiff(hub_genes, found_targets)
if (length(missing_targets)) {
  warning("Gene(s) not found in expression data: ",
          paste(missing_targets, collapse = ", "))
}
if (!length(found_targets)) stop("None of the selected hub genes is present.")
if (!"data" %in% Layers(sc[["RNA"]])) {
  sc <- NormalizeData(sc, assay = "RNA", verbose = FALSE)
}
sc <- FindVariableFeatures(sc, assay = "RNA", selection.method = "vst",
                            nfeatures = n_hvgs, verbose = FALSE)
hvgs <- VariableFeatures(sc)
network_genes <- intersect(unique(c(hvgs, found_targets)),
                           rownames(counts_full))
counts <- counts_full[network_genes, , drop = FALSE]
stopifnot(all(found_targets %in% rownames(counts)))
write.csv(data.frame(gene = rownames(counts),
                     is_hvg = rownames(counts) %in% hvgs,
                     is_hub_gene = rownames(counts) %in% found_targets),
          file.path(out_dir, "network_features.csv"), row.names = FALSE)
message("Network input: ", nrow(counts), " genes x ", ncol(counts), " cells.")

# 4. Current and legacy scTenifoldKnk API compatibility -----------------------
# Recent versions use qc_maxMTratio/qc_minLibSize; older versions use
# qc_mtThreshold/qc_minLSize. Only pass options accepted by the installed API.
formal_args <- names(formals(scTenifoldKnk::scTenifoldKnk))
num_cores <- parallel::detectCores(logical = TRUE)
if (is.na(num_cores) || num_cores < 1L) num_cores <- 1L
num_cores <- max(1L, min(8L, num_cores - 2L))

ko_args <- list(
  countMatrix = counts,
  qc = TRUE,
  nc_lambda = 0,
  nc_nNet = n_networks,
  nc_nCells = min(max_cells_per_network, ncol(counts)),
  nc_nComp = 3,
  nc_scaleScores = TRUE,
  nc_symmetric = FALSE,
  nc_q = 0.9,
  td_K = 3,
  td_maxIter = 1000,
  td_maxError = 1e-5,
  td_nDecimal = 3,
  ma_nDim = 2
)
if ("qc_maxMTratio" %in% formal_args) {
  ko_args$qc_maxMTratio <- 0.1
} else if ("qc_mtThreshold" %in% formal_args) {
  ko_args$qc_mtThreshold <- 0.1
}
if ("qc_minLibSize" %in% formal_args) {
  ko_args$qc_minLibSize <- 1000
} else if ("qc_minLSize" %in% formal_args) {
  ko_args$qc_minLSize <- 1000
}
if ("nCores" %in% formal_args) ko_args$nCores <- num_cores
unknown <- setdiff(names(ko_args), formal_args)
if (length(unknown)) stop("Unsupported package arguments: ",
                          paste(unknown, collapse = ", "))

# 5. Plot functions -----------------------------------------------------------
plot_z_vs_q <- function(dat, gene) {
  dat <- dat %>%
    mutate(log_p = -log10(pmax(p.adj, 1e-300)),
           significant = !is.na(Z) & !is.na(p.adj) &
             abs(Z) >= z_cutoff & p.adj < adjusted_p_cutoff)
  labels <- dat %>% filter(significant) %>%
    arrange(p.adj) %>% slice_head(n = 50)
  ggplot(dat, aes(x = Z, y = log_p, colour = significant)) +
    geom_point(size = 1.6, alpha = 0.75, na.rm = TRUE) +
    geom_vline(xintercept = c(-z_cutoff, z_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(adjusted_p_cutoff), linetype = "dashed") +
    ggrepel::geom_text_repel(data = labels, aes(label = gene),
                             size = 3, seed = 123, max.overlaps = Inf) +
    scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "#D73027")) +
    labs(title = paste(gene, "virtual knockout"), x = "Z-score",
         y = expression(-log[10](adjusted~P)), colour = "Significant") +
    theme_classic()
}

plot_top_fc <- function(dat, gene) {
  ggplot(dat, aes(x = reorder(gene, FC), y = FC, fill = FC)) +
    geom_col(width = 0.7) + coord_flip() +
    scale_fill_gradient(low = "#9ecae1", high = "#08519c") +
    labs(title = paste(gene, "virtual knockout"),
         x = "Gene", y = "Fold change (FC)") +
    theme_minimal() + guides(fill = "none")
}

# 6. Independent knockout runs ------------------------------------------------
results <- list()
tables <- list()
summaries <- list()
for (target in found_targets) {
  message("Running virtual knockout: ", target)
  set.seed(base_seed)  # Reproducible initialization for each independent run
  args <- c(ko_args, list(gKO = target))
  run <- tryCatch(do.call(scTenifoldKnk::scTenifoldKnk, args),
                  error = function(e) {
                    message("KO failed [", target, "]: ", conditionMessage(e))
                    NULL
                  })
  if (is.null(run)) {
    summaries[[target]] <- data.frame(KO_gene = target, Status = "FAILED",
                                     N_tested = NA_integer_,
                                     N_significant = NA_integer_)
    next
  }

  df <- as.data.frame(run$diffRegulation)
  if (!"gene" %in% names(df) && !is.null(rownames(df))) {
    df$gene <- rownames(df)
  }
  needed <- c("gene", "FC", "Z", "p.adj")
  if (!all(needed %in% names(df))) {
    message("Unexpected differential-regulation table for ", target,
            "; missing: ", paste(setdiff(needed, names(df)), collapse = ", "))
    summaries[[target]] <- data.frame(KO_gene = target,
                                     Status = "INVALID_OUTPUT",
                                     N_tested = NA_integer_,
                                     N_significant = NA_integer_)
    next
  }
  df <- df %>% filter(gene != target) %>% mutate(KO_gene = target)
  sig <- df %>% filter(!is.na(Z), !is.na(p.adj),
                       abs(Z) >= z_cutoff, p.adj < adjusted_p_cutoff) %>%
    arrange(p.adj, desc(abs(Z)))
  top <- df %>% filter(!is.na(FC)) %>%
    arrange(desc(FC)) %>% slice_head(n = 20)

  saveRDS(run, file.path(out_dir, paste0(target, "_KO_result.rds")))
  write.csv(df, file.path(out_dir, paste0(target, "_all_genes.csv")),
            row.names = FALSE)
  write.csv(sig, file.path(out_dir, paste0(target, "_significant.csv")),
            row.names = FALSE)
  write.csv(top, file.path(out_dir, paste0(target, "_top20_FC.csv")),
            row.names = FALSE)
  p_z <- plot_z_vs_q(df, target)
  ggsave(file.path(out_dir, paste0(target, "_Z_vs_adjustedP.pdf")),
         p_z, width = 8, height = 6)
  ggsave(file.path(out_dir, paste0(target, "_Z_vs_adjustedP.png")),
         p_z, width = 8, height = 6, dpi = 600)
  if (nrow(top)) {
    p_fc <- plot_top_fc(top, target)
    ggsave(file.path(out_dir, paste0(target, "_top20_FC.pdf")),
           p_fc, width = 7, height = 6)
    ggsave(file.path(out_dir, paste0(target, "_top20_FC.png")),
           p_fc, width = 7, height = 6, dpi = 600)
  }
  summaries[[target]] <- data.frame(KO_gene = target, Status = "SUCCESS",
                                   N_tested = nrow(df),
                                   N_significant = nrow(sig))
  tables[[target]] <- df
  if (save_combined_objects) results[[target]] <- run
  rm(run, df, sig, top, p_z)
  gc(verbose = FALSE)
}

# 7. Summary and reproducibility ----------------------------------------------
summary_df <- dplyr::bind_rows(summaries)
summary_df$KO_gene <- factor(summary_df$KO_gene, levels = hub_genes)
summary_df <- dplyr::arrange(summary_df, KO_gene)
write.csv(summary_df, file.path(out_dir, "HubGene_KO_summary.csv"),
          row.names = FALSE)
if (save_combined_objects) {
  saveRDS(results, file.path(out_dir, "All_HubGene_KO_results.rds"))
}
if (length(tables)) {
  all_genes <- dplyr::bind_rows(tables)
  all_sig <- all_genes %>%
    filter(!is.na(Z), !is.na(p.adj), abs(Z) >= z_cutoff,
           p.adj < adjusted_p_cutoff)
  write.csv(all_genes, file.path(out_dir, "All_HubGene_KO_all_genes.csv"),
            row.names = FALSE)
  write.csv(all_sig, file.path(out_dir, "All_HubGene_KO_significant.csv"),
            row.names = FALSE)
}
plot_data <- summary_df %>% filter(Status == "SUCCESS")
if (nrow(plot_data)) {
  p <- ggplot(plot_data, aes(x = KO_gene, y = N_significant)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = N_significant), vjust = -0.4) +
    labs(title = paste("Virtual knockouts in", target_celltype),
         x = NULL, y = "Significantly perturbed genes") +
    theme_classic()
  ggsave(file.path(out_dir, "HubGene_KO_comparison.pdf"),
         p, width = 8, height = 5)
  ggsave(file.path(out_dir, "HubGene_KO_comparison.png"),
         p, width = 8, height = 5, dpi = 600)
}
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
print(summary_df)
message("Knockout analysis complete: ", normalizePath(out_dir))
