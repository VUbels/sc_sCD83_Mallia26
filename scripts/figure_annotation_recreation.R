#!/usr/bin/env Rscript

# =============================================================================
# Fig S3 Recreation: Split UMAPs + MiloR DA + Marker Dot Plots
# =============================================================================
# Per mapping_cell_type:
#   1. Subset Seurat object
#   2. cluster_subcluster() → Seurat PCA/UMAP on the subset
#   3. Convert to SCE, inject the Seurat PCA + UMAP into reducedDim
#   4. Build Milo on the SAME reduction → DA testing
#   5. Split UMAP + Milo DA graph + marker dot plot all use the same space
#
# This avoids recomputing PCA/UMAP through scater/scran. MiloR only needs
# a reducedDim slot to build its KNN graph — it doesn't care whether PCA
# came from Seurat or scater.
# =============================================================================

library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scuttle)        # only for logNormCounts (MiloR needs logcounts assay)
library(dplyr)
library(ggplot2)
library(patchwork)
library(BiocParallel)
library(igraph)

source("./scripts/helper_functions.R")

# =============================================================================
# 1. PLOT FUNCTIONS
# =============================================================================

#' Run MiloR DA using Seurat reductions
#'
#' Takes a Seurat subset that already has PCA/UMAP (from cluster_subcluster),
#' converts to SCE, injects those reductions, and runs the Milo pipeline.
#'
#' @param obj_sub        Seurat object with PCA + UMAP reductions
#' @param k              KNN k (default 50)
#' @param d              PCA dims to use (default 50, capped to available)
#' @param prop           makeNhoods proportion (default 0.1)
#' @param sample_col     Sample column (default "orig.ident")
#' @param treatment_col  Treatment column (default "treatment")
#' @param pca_name       Name of PCA reduction in Seurat obj (default "pca")
#' @param umap_name      Name of UMAP reduction in Seurat obj (default "umap")
#' @return list: milo (Milo object), da_results (data.frame)
run_milo_pipeline <- function(obj_sub,
                              k             = 50,
                              d             = 50,
                              prop          = 0.1,
                              sample_col    = "orig.ident",
                              treatment_col = "treatment",
                              pca_name      = "pca",
                              umap_name     = "umap") {
  
  # --- Convert to SCE ---
  DefaultAssay(obj_sub) <- "RNA"
  obj.sce <- as.SingleCellExperiment(obj_sub)
  
  # MiloR needs a logcounts assay. The Seurat→SCE conversion should
  # produce one from the RNA "data" layer, but if it's missing, add it.
  if (!"logcounts" %in% assayNames(obj.sce)) {
    cat("  logcounts missing — running logNormCounts...\n")
    obj.sce <- scuttle::logNormCounts(obj.sce)
  }
  
  # --- Inject Seurat reductions ---
  seurat_pca  <- Embeddings(obj_sub, reduction = pca_name)
  seurat_umap <- Embeddings(obj_sub, reduction = umap_name)
  
  # Cap d to available PCA components
  if (ncol(seurat_pca) < d) {
    d <- ncol(seurat_pca)
    cat("  PCA has", ncol(seurat_pca), "components — using d =", d, "\n")
  }
  
  reducedDim(obj.sce, "PCA")  <- seurat_pca[, 1:d]
  reducedDim(obj.sce, "UMAP") <- seurat_umap
  
  # --- Build Milo ---
  cat("  Building Milo object...\n")
  obj.milo <- Milo(obj.sce)
  obj.milo <- buildGraph(obj.milo, k = k, d = d, reduced.dim = "PCA")
  obj.milo <- makeNhoods(obj.milo, prop = prop, k = k, d = d,
                         refined = TRUE, reduced_dims = "PCA")
  
  cat("  Neighborhoods:", ncol(nhoods(obj.milo)),
      " | Mean size:", round(mean(colSums(nhoods(obj.milo))), 1), "\n")
  
  # --- Count + distances ---
  obj.milo <- countCells(obj.milo,
                         meta.data = data.frame(colData(obj.milo)),
                         samples   = sample_col)
  obj.milo <- calcNhoodDistance(obj.milo, d = d, reduced.dim = "PCA")
  
  # --- Design matrix ---
  sample_meta <- colData(obj.milo) %>%
    as.data.frame() %>%
    dplyr::select(all_of(c(sample_col, treatment_col))) %>%
    distinct()
  
  design_matrix <- data.frame(
    Sample    = sample_meta[[sample_col]],
    Treatment = sample_meta[[treatment_col]]
  )
  rownames(design_matrix) <- design_matrix$Sample
  design_matrix <- design_matrix[colnames(nhoodCounts(obj.milo)), ]
  
  cat("  Design matrix:\n")
  print(design_matrix)
  
  # --- DA testing ---
  cat("  Testing for DA...\n")
  da_results <- testNhoods(obj.milo,
                           design    = ~ Treatment,
                           design.df = design_matrix)
  
  cat("  Sig neighborhoods (P < 0.05):",
      sum(da_results$PValue < 0.05, na.rm = TRUE),
      "| Up:", sum(da_results$PValue < 0.05 & da_results$logFC > 0, na.rm = TRUE),
      "| Down:", sum(da_results$PValue < 0.05 & da_results$logFC < 0, na.rm = TRUE), "\n")
  
  # --- Nhood graph + annotate ---
  obj.milo <- buildNhoodGraph(obj.milo)
  
  if ("fine_clust" %in% colnames(colData(obj.milo))) {
    da_results <- annotateNhoods(obj.milo, da_results,
                                 coldata_col = "fine_clust")
  }
  
  return(list(milo = obj.milo, da_results = da_results))
}


#' Split UMAP from Seurat object, colored by fine_clust
#'
#' Uses the same UMAP reduction that was injected into the Milo object.
#'
#' @param obj_sub   Seurat object (already has UMAP from cluster_subcluster)
#' @param split_by  Facet column (default "treatment")
#' @param color_by  Color column (default "fine_clust")
#' @param reduction Reduction name (default "umap")
#' @param pt_size   Point size (default 0.4)
#' @param title     Plot title
#' @return ggplot
plot_split_umap <- function(obj_sub,
                            split_by  = "treatment",
                            color_by  = "fine_clust",
                            reduction = "umap",
                            pt_size   = 0.4,
                            title     = NULL) {
  
  emb <- Embeddings(obj_sub, reduction = reduction)
  df  <- data.frame(
    UMAP_1    = emb[, 1],
    UMAP_2    = emb[, 2],
    split_var = obj_sub@meta.data[[split_by]],
    color_var = obj_sub@meta.data[[color_by]],
    stringsAsFactors = FALSE
  )
  
  set.seed(42)
  df <- df[sample(nrow(df)), ]
  
  if (is.null(title)) title <- paste0(ncol(obj_sub), " cells")
  
  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = color_var)) +
    geom_point(size = pt_size, stroke = 0) +
    facet_wrap(~ split_var, ncol = 2) +
    labs(title = title, color = color_by) +
    theme_void(base_size = 11) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
      strip.text      = element_text(face = "bold", size = 12),
      legend.position = "right",
      legend.title    = element_text(size = 9),
      legend.text     = element_text(size = 8),
      legend.key.size = unit(0.35, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))
}


#' MiloR DA neighborhood graph
#'
#' layout="UMAP" pulls from the same reducedDim injected from Seurat.
plot_milo_da <- function(milo_obj, da_results, alpha = 0.05, title = NULL) {
  
  p <- plotNhoodGraphDA(milo_obj, da_results, layout = "UMAP", alpha = alpha)
  if (!is.null(title)) p <- p + ggtitle(title)
  
  p + theme_void(base_size = 11) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    )
}


#' Marker gene dot plot
plot_marker_dotplot <- function(obj_sub,
                                markers,
                                cluster_order = NULL,
                                col_low       = "steelblue",
                                col_mid       = "white",
                                col_high      = "darkred",
                                dot_range     = c(0.5, 6)) {
  
  Idents(obj_sub) <- "fine_clust"
  
  if (!is.null(cluster_order)) {
    present <- intersect(cluster_order, unique(obj_sub$fine_clust))
    obj_sub@meta.data$fine_clust <- factor(obj_sub@meta.data$fine_clust,
                                           levels = present)
    Idents(obj_sub) <- "fine_clust"
  }
  
  all_genes <- unlist(markers, use.names = FALSE)
  missing   <- setdiff(all_genes, rownames(obj_sub))
  if (length(missing) > 0) {
    warning("Genes not found (skipped): ", paste(missing, collapse = ", "))
    all_genes <- intersect(all_genes, rownames(obj_sub))
    markers   <- lapply(markers, function(g) g[g %in% all_genes])
    markers   <- markers[sapply(markers, length) > 0]
    all_genes <- unlist(markers, use.names = FALSE)
  }
  
  p <- DotPlot(obj_sub, features = all_genes) +
    scale_color_gradient2(low = col_low, mid = col_mid, high = col_high,
                          midpoint = 0, name = "Scaled\nExpression") +
    scale_size_continuous(range = dot_range, name = "% Expressed",
                          breaks = c(25, 50, 75, 100)) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                  size = 9, face = "bold"),
      axis.text.y  = element_text(face = "bold", size = 10),
      axis.title   = element_blank(),
      panel.grid.major = element_line(linewidth = 0.15, color = "grey85"),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      legend.key.height = unit(0.3, "cm"),
      legend.key.width  = unit(0.3, "cm"),
      legend.title = element_text(size = 8.5),
      legend.text  = element_text(size = 7.5)
    ) +
    coord_flip()
  
  gene_group_sizes <- sapply(markers, length)
  v_breaks <- cumsum(gene_group_sizes)
  v_breaks <- v_breaks[-length(v_breaks)] + 0.5
  if (length(v_breaks) > 0) {
    p <- p + geom_vline(xintercept = v_breaks, linewidth = 0.5, color = "grey30")
  }
  
  return(p)
}


#' Compose row: [split UMAP] | [Milo DA] | [dot plot]
compose_row <- function(p_umap, p_milo, p_dot,
                        widths = c(3, 2, 3), row_label = NULL) {
  row <- p_umap + p_milo + p_dot + plot_layout(ncol = 3, widths = widths)
  if (!is.null(row_label)) {
    row <- row + plot_annotation(
      title = row_label,
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))
    )
  }
  return(row)
}


# =============================================================================
# 2. MARKER GENES AND CLUSTER ORDERS
# =============================================================================
# From *_cluster_identification in the clustering script.
# Corrections: CPE3→CPA3, NFKB→NFKBIA, CALD→CALD1
# =============================================================================

marker_genes <- list(
  
  Keratinocytes = list(
    ORS.1          = c("CXCL8", "MMP3", "CXCL1"),
    ORS.2          = c("CXCL14", "S100A2"),
    ORS.Suprabasal = c("KRT6A", "KRT16", "KRT17"),
    ORS.3          = c("KRT5", "KRT14", "GAS5"),
    Bulge          = c("EGFR", "CDH13", "RUNX1"),
    Eccrine.cells  = c("KRT19", "AQP5", "MUCL1"),
    Ishtmus        = c("KLK10", "MXD1", "FOXO3", "PRDM1"),
    Matrix         = c("KRT35", "KRT85")
  ),
  
  Fibroblasts = list(
    FBs.Activated          = c("CEMIP", "MT1X", "IGFBP2"),
    FBs.Interstitial       = c("COL1A1", "COL6A3", "GAS1"),
    FBs.Infl.Myofibroblast = c("MMP3", "MMP1", "CXCL5"),
    FBs.Pericytes          = c("PI15", "ESAM", "ABCC9", "RGS5"),
    FBs.DP.like            = c("VCAM1", "NGFR", "OXTR"),
    FBs.Perifolicular      = c("NTN1", "PLXDC2", "LAMA2"),
    FBs.Cycling            = c("DDX21", "WDR43", "ODC1"),
    FBs.NF.kB              = c("CXCL2", "CXCL3", "NFKBIA"),
    FBs.Dermal.Sheath      = c("COL11A1", "TAGLN", "PPP1R14A")
  ),
  
  Immune = list(
    T.Naive.ANK3    = c("PTPRC", "CD96", "IL7R", "ANK3"),
    T.Naive         = c("CD3D", "CD3E", "IL7R"),
    M2.Macrophages  = c("MRC1", "CCL3", "CD83"),
    MAIT.cells      = c("KLRB1", "BLK"),
    Mast.cells      = c("CPA3", "KIT"),
    T.reg           = c("FOXP3", "IKZF2"),
    T.cyto          = c("CD8A", "NKG7"),
    Dendritic.cells = c("HLA-DRA", "LAMP3"),
    Plasma.cells    = c("IGKC", "JCHAIN", "MZB1"),
    B.cells         = c("CD19", "MS4A1", "PAX5")
  ),
  
  Endothelial = list(
    Endo.General    = c("PLVAP", "PECAM1", "VWF"),
    Endo.Vascular   = c("VEGFA", "PLVAP"),
    Endo.Angiogenic = c("BICC1"),
    Endo.Lymphatic  = c("CCL21")
  ),
  
  Remaining = list(
    Melanocytes   = c("MLANA", "SOX5", "DCT"),
    Schwann.cells = c("NGFR", "PMP2", "SOX10")
  )
)


# Actual fine_clust strings (hyphens where applicable)
cluster_order_map <- list(
  
  Keratinocytes = c(
    "Matrix", "Ishtmus", "Eccrine.cells", "Bulge",
    "ORS.3", "ORS.Suprabasal", "ORS.2", "ORS.1",
    "FBs.Reticular"
  ),
  
  Fibroblasts = c(
    "FBs.Dermal.Sheath", "FBs.NF-kB", "FBs.Cycling",
    "FBs.Perifolicular", "FBs.DP-like", "FBs.Pericytes",
    "FBs.Infl.Myofibroblast", "FBs.Interstitial", "FBs.Activated"
  ),
  
  Immune = c(
    "B.cells", "Plasma.cells", "Dendritic.cells",
    "T.cyto", "T.reg", "Mast.cells", "MAIT.cells",
    "M2.Macrophages", "T.Naive", "T.Naive-ANK3",
    "FBs.Dermal.Sheath"
  ),
  
  Endothelial = c(
    "Endo.Lymphatic", "Endo.Angiogenic",
    "Endo.Vascular", "Endo.General"
  ),
  
  Remaining = c(
    "Schwann.cells", "Melanocytes"
  )
)


# Subclustering parameters per cell type (matching the clustering script)
subcluster_params <- list(
  Keratinocytes = list(dims = 1:30),
  Fibroblasts   = list(dims = 1:20, n_genes = 2000, features = 2000, conserve.memory = FALSE),
  Immune        = list(dims = 1:20, n_genes = 2000, features = 2000, conserve.memory = FALSE),
  Endothelial   = list(dims = 1:20, n_genes = 2000, features = 2000, conserve.memory = FALSE),
  Remaining     = list(dims = 1:20, n_genes = 2000, features = 2000, conserve.memory = FALSE)
)


# =============================================================================
# 3. MASTER PIPELINE
# =============================================================================

main_folder <- "./"
input_file  <- paste0(main_folder, "fine_annotated_obj.rds")
output_dir  <- file.path(main_folder, "fig_s3_panels/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(input_file)
DefaultAssay(obj) <- "RNA"

if (!"mapping_cell_type" %in% colnames(obj@meta.data)) {
  obj$mapping_cell_type <- sub("\\.\\d+$", "", obj$broad_cluster)
}

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

# --- Verify fine_clust labels ---
cat("=== fine_clust per mapping_cell_type ===\n")
for (ct in names(marker_genes)) {
  present  <- unique(obj$fine_clust[obj$mapping_cell_type == ct])
  expected <- cluster_order_map[[ct]]
  cat("\n", ct, ":\n  Present:", paste(sort(present), collapse = ", "), "\n")
  miss  <- setdiff(expected, present)
  extra <- setdiff(present, expected)
  if (length(miss) > 0)  cat("  MISSING:", paste(miss, collapse = ", "), "\n")
  if (length(extra) > 0) cat("  EXTRA:", paste(extra, collapse = ", "), "\n")
}


# =============================================================================
# 4. LOOP
# =============================================================================

cell_types <- names(marker_genes)
all_rows   <- list()

for (ct in cell_types) {
  
  cat("\n===== Processing:", ct, "=====\n")
  
  # --- Subset ---
  obj_sub <- subset(obj, subset = mapping_cell_type == ct)
  cat("Cells:", ncol(obj_sub), "\n")
  
  # --- Subcluster: Seurat PCA + UMAP ---
  # Uses the same function and parameters as the clustering script.
  # After this, obj_sub has "pca" and "umap" reductions.
  params <- subcluster_params[[ct]]
  obj_sub <- do.call(cluster_subcluster, c(list(obj_sub, output_dir = file.path(output_dir, ct)), params))
  
  # --- A: Split UMAP (from Seurat reduction) ---
  p_umap <- plot_split_umap(
    obj_sub,
    split_by  = "treatment",
    color_by  = "fine_clust",
    reduction = "umap",
    pt_size   = 0.3,
    title     = paste0(ncol(obj_sub), " cells")
  )
  
  # --- B: MiloR DA (same Seurat PCA/UMAP injected into Milo) ---
  # The d parameter is capped to the number of PCA dims from cluster_subcluster
  milo_d <- length(params$dims)  # e.g., 30 for keratinocytes, 20 for others
  
  milo_res <- run_milo_pipeline(
    obj_sub,
    k         = 50,
    d         = milo_d,
    prop      = 0.1,
    sample_col    = "orig.ident",
    treatment_col = "treatment",
    pca_name  = "pca",
    umap_name = "umap"
  )
  
  p_milo <- plot_milo_da(
    milo_res$milo,
    milo_res$da_results,
    alpha = 0.05
  )
  
  # --- C: Marker dot plot ---
  p_dot <- plot_marker_dotplot(
    obj_sub,
    markers       = marker_genes[[ct]],
    cluster_order = cluster_order_map[[ct]]
  )
  
  # --- Compose + save ---
  row <- compose_row(p_umap, p_milo, p_dot, row_label = ct)
  all_rows[[ct]] <- row
  
  ggsave(file.path(output_dir, paste0(ct, "_umap.png")),    p_umap, width = 10, height = 5,  dpi = 300)
  ggsave(file.path(output_dir, paste0(ct, "_milo.png")),    p_milo, width = 7,  height = 6,  dpi = 300)
  ggsave(file.path(output_dir, paste0(ct, "_dotplot.png")), p_dot,  width = 10, height = 6,  dpi = 300)
  
  # DA results + diagnostics
  cell_type_dir <- file.path(output_dir, ct)
  dir.create(cell_type_dir, showWarnings = FALSE, recursive = TRUE)
  
  write.csv(milo_res$da_results,
            file.path(cell_type_dir, "milo_da_results.csv"), row.names = FALSE)
  
  # Beeswarm
  if ("fine_clust" %in% colnames(milo_res$da_results)) {
    da_plot <- milo_res$da_results
    da_plot$SpatialFDR <- da_plot$PValue
    p_bee <- plotDAbeeswarm(da_plot, group.by = "fine_clust") +
      ggtitle(ct) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
    ggsave(file.path(cell_type_dir, "DA_beeswarm.png"), p_bee, width = 10, height = 6, dpi = 300)
  }
  
  # DA cells
  cells_up   <- extract_DA_cells(milo_res$milo, milo_res$da_results, alpha = 0.05, direction = "up", use_pvalue = TRUE)
  cells_down <- extract_DA_cells(milo_res$milo, milo_res$da_results, alpha = 0.05, direction = "down", use_pvalue = TRUE)
  write.csv(data.frame(barcode = cells_up),   file.path(cell_type_dir, "DA_cells_up.csv"),   row.names = FALSE)
  write.csv(data.frame(barcode = cells_down), file.path(cell_type_dir, "DA_cells_down.csv"), row.names = FALSE)
  
  # Summary stats
  if ("fine_clust" %in% colnames(milo_res$da_results)) {
    summary_stats <- milo_res$da_results %>%
      group_by(fine_clust) %>%
      summarise(
        n_nhoods     = n(),
        n_sig_p0.05  = sum(PValue < 0.05, na.rm = TRUE),
        n_up         = sum(PValue < 0.05 & logFC > 0, na.rm = TRUE),
        n_down       = sum(PValue < 0.05 & logFC < 0, na.rm = TRUE),
        mean_logFC   = mean(logFC, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(n_sig_p0.05))
    write.csv(summary_stats, file.path(cell_type_dir, "DA_summary.csv"), row.names = FALSE)
  }
  
  saveRDS(milo_res$milo, file.path(cell_type_dir, "milo_object.rds"))
  
  rm(obj_sub, milo_res)
  gc()
}

# =============================================================================
# 5. FULL FIGURE
# =============================================================================

full_fig <- wrap_plots(all_rows, ncol = 1)

ggsave(file.path(output_dir, "Fig_S3_full.png"), full_fig,
       width = 22, height = 8 * length(cell_types), dpi = 300, limitsize = FALSE)
ggsave(file.path(output_dir, "Fig_S3_full.pdf"), full_fig,
       width = 22, height = 8 * length(cell_types), limitsize = FALSE)

cat("\nDone. Output:", output_dir, "\n")