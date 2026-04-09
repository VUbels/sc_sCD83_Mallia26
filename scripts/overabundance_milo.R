#!/usr/bin/env Rscript

# =============================================================================
# Milo DA Pipeline with Consistent Seurat-Based Embeddings
# =============================================================================
# Per mapping_cell_type:
#   - Subset Seurat object
#   - cluster_subcluster() -> PCA + UMAP in Seurat
#   - Plot split UMAP and single UMAP from the Seurat reduction
#   - Convert to SCE, inject the Seurat PCA + UMAP into reducedDim
#   - Build Milo on the SAME PCA -> DA testing
#   - Plot Milo DA neighbourhood graph in the SAME UMAP space
#   - Full diagnostic suite: beeswarm, volcano, p-value hist, nhood sizes
#   - Extract DA cells, save summary stats
#
# Key design decision: Milo never recomputes PCA/UMAP via scater/scran.
# The Seurat subclustering reduction is THE reduction for everything.
# This guarantees that the split UMAP, single UMAP, and neighbourhood
# graph all share identical coordinates.
#
# Assumes Seurat v5 (layers, not assays).
# =============================================================================

library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scuttle)
library(dplyr)
library(ggplot2)
library(patchwork)
library(BiocParallel)
library(igraph)

# Source helpers
source("./scripts/helper_functions.R")
source("./scripts/milo_plot_helpers.R")


###################################################
# PARAMETERS
###################################################

main_folder <- "./"
input_file  <- paste0(main_folder, "fine_annotated_obj.rds")
output_dir  <- file.path(main_folder, "milo_results/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Subclustering parameters per cell type.
# These match the clustering script and control Seurat PCA dims.
# The same dims are passed to Milo's buildGraph(d = ...).
subcluster_params <- list(
  Keratinocytes = list(dims = 1:30),
  Fibroblasts   = list(dims = 1:20, n_genes = 2000, features = 2000, conserve.memory = FALSE),
  Immune        = list(dims = 1:20, n_genes = 2000, features = 2000, conserve.memory = FALSE),
  Endothelial   = list(dims = 1:20, n_genes = 2000, features = 2000, conserve.memory = FALSE),
  Remaining     = list(dims = 1:20, n_genes = 2000, features = 2000, conserve.memory = FALSE)
)


###################################################
# DATA LOADING
###################################################

obj <- readRDS(input_file)

# Force RNA so SCE conversion pulls raw counts, not SCT residuals
DefaultAssay(obj) <- "RNA"

# Derive mapping_cell_type from broad_cluster if absent
if (!"mapping_cell_type" %in% colnames(obj@meta.data)) {
  obj$mapping_cell_type <- sub("\\.\\d+$", "", obj$broad_cluster)
}

# Normalize the RNA layer
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Resolve cell types to loop over
cell_types <- unique(obj@meta.data$mapping_cell_type)
cell_types <- cell_types[!is.na(cell_types)]

cat("Cell types to analyze:\n")
print(cell_types)


###################################################
# MAIN LOOP
###################################################

bpparam <- SerialParam()
register(bpparam)

for (ct in cell_types) {
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("  Processing:", ct, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  cell_type_dir <- file.path(output_dir, ct)
  dir.create(cell_type_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Subset
  obj_sub <- subset(obj, subset = mapping_cell_type == ct)
  cat("Cells:", ncol(obj_sub), "\n")
  
  # Subcluster: Seurat PCA + UMAP.
  # After this, obj_sub has "pca" and "umap" reductions.
  # These are THE reductions used everywhere downstream.
  params <- subcluster_params[[ct]]
  if (!is.null(params)) {
    obj_sub <- do.call(
      cluster_subcluster,
      c(list(obj_sub, output_dir = cell_type_dir), params)
    )
  } else {
    cat("  No preset params for", ct, "-- using defaults (dims 1:20)\n")
    obj_sub <- cluster_subcluster(obj_sub, output_dir = cell_type_dir, dims = 1:20)
    params <- list(dims = 1:20)
  }
  
  # Split UMAP from Seurat
  p_split <- plot_split_dimred(
    obj_sub,
    split_by  = "treatment",
    color_by  = "fine_clust",
    reduction = "umap",
    pt_size   = 0.3,
    title     = paste0(ct, " | ", ncol(obj_sub), " cells")
  )
  ggsave(file.path(cell_type_dir, "UMAP_split_treatment.png"),
         p_split, width = 10, height = 5, dpi = 300)
  
  # Single UMAP colored by fine_clust
  p_single <- plot_single_dimred(
    obj_sub,
    color_by  = "fine_clust",
    reduction = "umap",
    pt_size   = 0.3,
    title     = paste0(ct, " | ", ncol(obj_sub), " cells")
  )
  ggsave(file.path(cell_type_dir, "UMAP_fine_clust.png"),
         p_single, width = 7, height = 6, dpi = 300)
  
  # Single UMAP colored by treatment
  p_treat <- plot_single_dimred(
    obj_sub,
    color_by  = "treatment",
    reduction = "umap",
    pt_size   = 0.3,
    title     = paste0(ct, " | Treatment")
  )
  ggsave(file.path(cell_type_dir, "UMAP_treatment.png"),
         p_treat, width = 7, height = 6, dpi = 300)
  
  # Run Milo pipeline (using Seurat PCA/UMAP).
  # d is set to the number of PCA dims from subclustering.
  # k and prop are adaptive based on dataset size.
  milo_d   <- length(params$dims)
  adaptive <- get_milo_params(ncol(obj_sub),
                              length(unique(obj_sub$orig.ident)))
  
  milo_res <- run_milo_pipeline(
    obj_sub,
    k             = adaptive$k,
    d             = milo_d,
    prop          = adaptive$prop,
    sample_col    = "orig.ident",
    treatment_col = "treatment",
    pca_name      = "pca",
    umap_name     = "umap"
  )
  
  # Milo DA neighbourhood graph.
  # layout="UMAP" uses the injected Seurat UMAP -> same coordinates.
  p_milo <- plot_milo_da(
    milo_res$milo,
    milo_res$da_results,
    alpha  = 0.05,
    layout = "UMAP",
    title  = paste0(ct, " | Milo DA")
  )
  ggsave(file.path(cell_type_dir, "nhood_graph_DA.png"),
         p_milo, width = 8, height = 8, dpi = 300)
  
  # Nhood graph with p-value filtering
  png(file.path(cell_type_dir, "nhood_graph_pval.png"),
      width = 8, height = 8, units = "in", res = 300)
  print(plotNhoodGraphDA_pval(milo_res$milo, milo_res$da_results,
                              alpha = 0.05, use_pvalue = TRUE,
                              layout = "UMAP"))
  dev.off()
  
  # Compose: [Split UMAP] | [Milo DA graph]
  p_row <- compose_row(p_split, p_milo, widths = c(3, 2),
                       row_label = ct)
  ggsave(file.path(cell_type_dir, "combined_umap_milo.png"),
         p_row, width = 16, height = 6, dpi = 300)
  
  # Neighbourhood size histogram
  p_nhood <- plot_nhood_size_hist(
    milo_res$milo,
    title = paste0(ct, " - Neighbourhood Sizes")
  )
  ggsave(file.path(cell_type_dir, "neighborhood_size_hist.png"),
         p_nhood, width = 6, height = 5, dpi = 300)
  
  # P-value histogram
  p_pval <- plot_pvalue_hist(
    milo_res$da_results,
    title = paste0(ct, " - P-value Distribution")
  )
  ggsave(file.path(cell_type_dir, "pvalue_histogram.png"),
         p_pval, width = 6, height = 5, dpi = 300)
  
  # Volcano plot
  p_volc <- plot_volcano(
    milo_res$da_results,
    pval_thresh = 0.05,
    title = paste0(ct, " - Volcano Plot")
  )
  ggsave(file.path(cell_type_dir, "volcano_plot.png"),
         p_volc, width = 7, height = 6, dpi = 300)
  
  # Beeswarm by fine_clust
  if ("fine_clust" %in% colnames(milo_res$da_results)) {
    p_bee <- plot_da_beeswarm(
      milo_res$da_results,
      group_by          = "fine_clust",
      title             = ct,
      use_pvalue_as_fdr = TRUE
    )
    ggsave(file.path(cell_type_dir, "DA_beeswarm.png"),
           p_bee, width = 10, height = 6, dpi = 300)
  }
  
  # Extract DA cells
  cells_up <- extract_DA_cells(milo_res$milo, milo_res$da_results,
                               alpha = 0.05, direction = "up",
                               use_pvalue = TRUE)
  cells_down <- extract_DA_cells(milo_res$milo, milo_res$da_results,
                                 alpha = 0.05, direction = "down",
                                 use_pvalue = TRUE)
  
  cat("  DA cells up:", length(cells_up),
      "| DA cells down:", length(cells_down), "\n")
  
  write.csv(data.frame(barcode = cells_up),
            file.path(cell_type_dir, "DA_cells_upregulated.csv"),
            row.names = FALSE)
  write.csv(data.frame(barcode = cells_down),
            file.path(cell_type_dir, "DA_cells_downregulated.csv"),
            row.names = FALSE)
  
  # Summary stats
  summary_stats <- summarise_da_by_cluster(milo_res$da_results,
                                           cluster_col = "fine_clust")
  if (!is.null(summary_stats)) {
    write.csv(summary_stats,
              file.path(cell_type_dir, "DA_summary_by_fineclust.csv"),
              row.names = FALSE)
  }
  
  # Save raw results + Milo object
  write.csv(milo_res$da_results,
            file.path(cell_type_dir, "milo_da_results_full.csv"),
            row.names = FALSE)
  saveRDS(milo_res$milo,
          file.path(cell_type_dir, "milo_object.rds"))
  
  # Cleanup
  rm(obj_sub, milo_res)
  gc()
  
  cat("  Completed:", ct, "\n")
}

cat("\nDone. All outputs in:", output_dir, "\n")