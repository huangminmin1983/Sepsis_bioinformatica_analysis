# =============================================================================
# cell.R | CellChat ligand-receptor communication analysis (human)
# Input : data/af_Activated_Macrophages.rds (Seurat; Celltype metadata)
# Output: results/cellchat/
# Run   : Rscript scripts/cell.R  (from the repository root)
# Note  : Written for the CellChat v1/v2 workflow; not Spatial CellChat v3.
# =============================================================================

set.seed(1234)
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(CellChat)
  library(ggplot2)
  library(dplyr)
  library(igraph)
  library(ggraph)
})

# 1. Configuration ------------------------------------------------------------
input_file <- file.path("data", "af_Activated_Macrophages.rds")
out_dir <- file.path("results", "cellchat")
celltype_col <- "Celltype"
min_cells_per_group <- 3
signaling_database <- "Secreted Signaling"   # Human-specific source database
pathways_to_plot <- c("GALECTIN", "RESISTIN", "ANNEXIN", "VISFATIN",
                      "VEGF", "CALCR")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, filename, w = 8, h = 6) {
  ggsave(file.path(out_dir, filename), plot = p, width = w, height = h)
}

# 2. Input and expression layers ----------------------------------------------
if (!file.exists(input_file)) stop("Missing Seurat object: ", input_file)
sce <- readRDS(input_file)
if (!celltype_col %in% colnames(sce@meta.data)) {
  stop("Missing cell-type metadata column: ", celltype_col)
}
DefaultAssay(sce) <- "RNA"
if (length(grep("^data", Layers(sce[["RNA"]]))) > 1L) {
  sce <- JoinLayers(sce, assay = "RNA")
}
if (!"data" %in% Layers(sce[["RNA"]])) {
  sce <- NormalizeData(sce, assay = "RNA", verbose = FALSE)
}

identity <- data.frame(celltype = as.character(sce@meta.data[[celltype_col]]),
                       row.names = rownames(sce@meta.data))
if (anyNA(identity$celltype) || any(!nzchar(identity$celltype))) {
  stop("Celltype contains missing or blank entries.")
}
expr <- LayerData(sce, assay = "RNA", layer = "data")
if (!identical(colnames(expr), rownames(identity))) {
  identity <- identity[colnames(expr), , drop = FALSE]
}
stopifnot(identical(colnames(expr), rownames(identity)))

# 3. Inference (CellChat v1/v2 API) -------------------------------------------
cellchat <- createCellChat(object = expr, meta = identity, group.by = "celltype")
cellchat@DB <- subsetDB(CellChatDB.human, search = signaling_database)
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = min_cells_per_group)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

write.csv(subsetCommunication(cellchat),
          file.path(out_dir, "ligand_receptor_interactions.csv"),
          row.names = FALSE)
write.csv(subsetCommunication(cellchat, slot.name = "netP"),
          file.path(out_dir, "pathway_interactions.csv"),
          row.names = FALSE)
saveRDS(cellchat, file.path(out_dir, "cellchat.rds"))

# 4. Network circle visualization --------------------------------------------
# Node sizes are matched by *cell-type name* (not table order).
plot_network <- function(network, metadata, title, edge_legend) {
  nodes <- data.frame(name = rownames(network), stringsAsFactors = FALSE)
  sizes <- table(metadata$celltype)
  nodes$size <- as.numeric(sizes[nodes$name])
  nodes$size[is.na(nodes$size)] <- 0

  edges <- as.data.frame(as.table(network), stringsAsFactors = FALSE)
  names(edges) <- c("Source", "Target", "Weight")
  edges <- edges[!is.na(edges$Weight) & edges$Weight > 0, , drop = FALSE]
  if (!nrow(edges)) {
    return(ggplot() + theme_void() +
             labs(title = title) +
             annotate("text", x = 0, y = 0, label = "No inferred interactions"))
  }
  g <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  ggraph(g, layout = "linear", circular = TRUE) +
    geom_edge_arc(aes(edge_width = Weight, edge_colour = node1.name),
                  arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed"),
                  alpha = 0.8, show.legend = TRUE) +
    geom_node_point(aes(size = size, colour = name)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    scale_edge_width(range = c(0.2, 2), name = edge_legend) +
    scale_edge_colour_discrete(name = "Source cell") +
    scale_colour_discrete(name = "Cell type") +
    scale_size_continuous(range = c(3, 10), name = "Cell count") +
    coord_fixed() + theme_void() +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")
}

save_plot(plot_network(cellchat@net$count, identity,
                       "Number of Interactions", "Count"),
          "network_interaction_count.pdf", 10, 10)
save_plot(plot_network(cellchat@net$weight, identity,
                       "Interaction Strength", "Weight"),
          "network_interaction_weight.pdf", 10, 10)

# 5. Bubble plots -------------------------------------------------------------
# Original analysis selected group indices (1, 2, 3, 4). Review the printed
# order below before biological interpretation; indices depend on this dataset.
group_names <- levels(cellchat@idents)
message("CellChat groups: ", paste(seq_along(group_names), group_names,
                                    sep = "=", collapse = "; "))
if (length(group_names) >= 4L) {
  p_source <- netVisual_bubble(cellchat, sources.use = 2,
                               targets.use = c(1, 3, 4), remove.isolate = FALSE)
  p_target <- netVisual_bubble(cellchat, sources.use = c(1, 3, 4),
                               targets.use = 2, remove.isolate = FALSE)
  save_plot(p_source, "bubble_source.pdf", 7, 5)
  save_plot(p_target, "bubble_target.pdf", 7, 5)
} else {
  message("Fewer than four groups: skipping index-based bubble comparisons.")
}

# 6. Signaling pathway expression and centrality ------------------------------
available_pathways <- cellchat@netP$pathways
for (pathway in pathways_to_plot) {
  if (!pathway %in% available_pathways) {
    message("Pathway not inferred; skipped: ", pathway)
    next
  }
  p <- plotGeneExpression(cellchat, signaling = pathway,
                          enriched.only = FALSE)
  save_plot(p, paste0("expression_", pathway, ".pdf"), 8, 6)
}

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
save_plot(netAnalysis_signalingRole_scatter(cellchat),
          "centrality_sender_receiver.pdf", 7, 6)

# Heatmaps are ComplexHeatmap objects, not ggplot objects: draw to PDF directly.
for (direction in c("outgoing", "incoming")) {
  ht <- netAnalysis_signalingRole_heatmap(cellchat, pattern = direction)
  grDevices::pdf(file.path(out_dir, paste0("centrality_", direction, ".pdf")),
                 width = 10, height = 7)
  ComplexHeatmap::draw(ht)
  grDevices::dev.off()
}
if ("MK" %in% available_pathways) {
  grDevices::pdf(file.path(out_dir, "MK_signaling_network.pdf"),
                 width = 9, height = 5)
  netAnalysis_signalingRole_network(cellchat, signaling = "MK",
                                    width = 8, height = 2.5, font.size = 10)
  grDevices::dev.off()
}

saveRDS(cellchat, file.path(out_dir, "cellchat_with_centrality.rds"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("CellChat analysis complete: ", normalizePath(out_dir))
