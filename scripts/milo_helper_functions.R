#!/usr/bin/env Rscript

# =============================================================================
# HELPER FUNCTIONS: Milo DA Pipeline + Consistent Plotting
# =============================================================================
# Standalone helper functions for running the Milo DA pipeline using
# Seurat reductions and producing consistent plots where the split
# UMAP/TSNE, single UMAP/TSNE, and Milo neighbourhood graph all
# share identical coordinates.
#
# Assumes Seurat v5 (layers, not assays).
# External dependency: helper_functions.R for cluster_subcluster()
# and theme_BOR() if used.
# =============================================================================


###################################################
# ADAPTIVE MILO PARAMETER SELECTION
###################################################
get_milo_params <- function(n_cells, n_samples) {
  # Select k, d, prop for Milo based on dataset size.
  # Returns a named list with k, d, prop.
  
  if (n_cells < 5000) {
    list(k = 15, d = 20, prop = 0.15)
  } else if (n_cells < 10000) {
    list(k = 20, d = 25, prop = 0.1)
  } else {
    list(k = 30, d = 30, prop = 0.1)
  }
}


###################################################
# RUN MILO PIPELINE USING SEURAT REDUCTIONS
###################################################
run_milo_pipeline <- function(obj_sub,
                              k             = 50,
                              d             = 50,
                              prop          = 0.1,
                              sample_col    = "orig.ident",
                              treatment_col = "treatment",
                              pca_name      = "pca",
                              umap_name     = "umap") {
  # Converts Seurat -> SCE, injects the Seurat PCA + UMAP into
  # reducedDim slots so Milo builds its KNN graph on the SAME
  # reduction used for the cell-level UMAP. Returns a list with
  # milo (Milo object) and da_results (data.frame).
  
  # Convert to SCE from RNA layer
  DefaultAssay(obj_sub) <- "RNA"
  obj.sce <- as.SingleCellExperiment(obj_sub)
  
  # MiloR requires a logcounts assay; add if missing
  if (!"logcounts" %in% assayNames(obj.sce)) {
    cat("  logcounts missing -- running logNormCounts...\n")
    obj.sce <- scuttle::logNormCounts(obj.sce)
  }
  
  # Inject Seurat reductions into SCE
  seurat_pca  <- Embeddings(obj_sub, reduction = pca_name)
  seurat_umap <- Embeddings(obj_sub, reduction = umap_name)
  
  # Cap d to available PCA components
  if (ncol(seurat_pca) < d) {
    d <- ncol(seurat_pca)
    cat("  PCA has", ncol(seurat_pca), "components -- capping d =", d, "\n")
  }
  
  reducedDim(obj.sce, "PCA")  <- seurat_pca[, 1:d]
  reducedDim(obj.sce, "UMAP") <- seurat_umap
  
  # Build Milo
  cat("  Building Milo object...\n")
  obj.milo <- Milo(obj.sce)
  obj.milo <- buildGraph(obj.milo, k = k, d = d, reduced.dim = "PCA")
  obj.milo <- makeNhoods(obj.milo, prop = prop, k = k, d = d,
                         refined = TRUE, reduced_dims = "PCA")
  
  cat("  Neighborhoods:", ncol(nhoods(obj.milo)),
      " | Mean size:", round(mean(colSums(nhoods(obj.milo))), 1), "\n")
  
  # Count cells + distances
  obj.milo <- countCells(obj.milo,
                         meta.data = data.frame(colData(obj.milo)),
                         samples   = sample_col)
  obj.milo <- calcNhoodDistance(obj.milo, d = d, reduced.dim = "PCA")
  
  # Design matrix
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
  
  # DA testing
  cat("  Testing for DA...\n")
  da_results <- testNhoods(obj.milo,
                           design    = ~ Treatment,
                           design.df = design_matrix)
  
  cat("  Sig neighborhoods (P < 0.05):",
      sum(da_results$PValue < 0.05, na.rm = TRUE),
      "| Up:", sum(da_results$PValue < 0.05 & da_results$logFC > 0, na.rm = TRUE),
      "| Down:", sum(da_results$PValue < 0.05 & da_results$logFC < 0, na.rm = TRUE), "\n")
  cat("  Sig neighborhoods (SpatialFDR < 0.2):",
      sum(da_results$SpatialFDR < 0.2, na.rm = TRUE), "\n")
  
  # Build nhood graph + annotate
  obj.milo <- buildNhoodGraph(obj.milo)
  
  if ("fine_clust" %in% colnames(colData(obj.milo))) {
    da_results <- annotateNhoods(obj.milo, da_results,
                                 coldata_col = "fine_clust")
  }
  
  return(list(milo = obj.milo, da_results = da_results))
}


###################################################
# CUSTOM PLOTTING FUNCTION - P-VALUES INSTEAD OF FDR
###################################################
plotNhoodGraphDA_pval <- function(x, milo_res, alpha = 0.1, res_column = "logFC", 
                                  use_pvalue = TRUE, layout = "UMAP", ...) {
  
  # Check if neighborhood graph exists
  if (is.null(nhoodGraph(x)) || length(igraph::E(nhoodGraph(x))) == 0) {
    stop("Not a valid Milo object - neighbourhood graph is missing. Please run buildNhoodGraph() first.")
  }
  
  # Check if layout is valid
  if (is.character(layout)) {
    if (!layout %in% names(reducedDims(x))) {
      stop(layout, " is not in reducedDim(x) - choose a different layout")
    }
  }
  
  # Add milo results to colData
  signif_res <- milo_res
  
  # Use PValue instead of SpatialFDR if specified
  if (use_pvalue) {
    signif_res$test_stat <- signif_res$PValue
  } else {
    signif_res$test_stat <- signif_res$SpatialFDR
  }
  
  # Handle NAs
  signif_res$test_stat[is.na(signif_res$test_stat)] <- 1
  
  # Set logFC to 0 for non-significant neighborhoods
  signif_res[signif_res$test_stat > alpha, res_column] <- 0
  
  # Add results to colData
  colData(x)[res_column] <- NA
  
  # Handle nhood subsetting
  if (any(names(list(...)) %in% c("subset.nhoods"))) {
    subset.nhoods <- list(...)$subset.nhoods
    sub.indices <- nhoodIndex(x)[subset.nhoods]
    colData(x)[unlist(sub.indices[signif_res$Nhood]), res_column] <- signif_res[, res_column]
  } else {
    colData(x)[unlist(nhoodIndex(x)[signif_res$Nhood]), res_column] <- signif_res[, res_column]
  }
  
  # Check for res_column in graph vertex attributes
  g_atts <- names(igraph::vertex_attr(nhoodGraph(x)))
  if (isFALSE(res_column %in% g_atts)) {
    message("Adding nhood effect sizes to neighbourhood graph attributes")
    
    if (any(names(list(...)) %in% c("subset.nhoods"))) {
      nh.v <- igraph::V(nhoodGraph(x))
      drop.v <- setdiff(nh.v, sub.indices)
      nhgraph <- nhoodGraph(x)
      nhgraph <- igraph::subgraph(nhgraph, sub.indices)
      nhgraph <- igraph::set_vertex_attr(nhgraph,
                                         name = res_column, value = signif_res[, res_column])
      nhoodGraph(x) <- nhgraph
    } else {
      nhoodGraph(x) <- igraph::set_vertex_attr(nhoodGraph(x), 
                                               name = res_column, 
                                               value = signif_res[, res_column])
    }
  }
  
  # Plot logFC - pass layout explicitly
  plotNhoodGraph(x, colour_by = res_column, layout = layout, 
                 size_range = c(1, 5), is.da = TRUE, ...) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
}


###################################################
# CELL EXTRACTION FUNCTION
###################################################
extract_DA_cells <- function(milo_obj, da_results, alpha = 0.05, 
                             direction = "both", use_pvalue = TRUE) {
  # direction: "up" (logFC > 0), "down" (logFC < 0), or "both"
  
  if (use_pvalue) {
    sig_col <- "PValue"
  } else {
    sig_col <- "SpatialFDR"
  }
  
  # Filter significant neighborhoods
  if (direction == "up") {
    sig_nhoods <- da_results$Nhood[da_results[[sig_col]] < alpha & da_results$logFC > 0]
  } else if (direction == "down") {
    sig_nhoods <- da_results$Nhood[da_results[[sig_col]] < alpha & da_results$logFC < 0]
  } else {
    sig_nhoods <- da_results$Nhood[da_results[[sig_col]] < alpha]
  }
  
  # Extract cells from significant neighborhoods
  cell_barcodes <- c()
  for (i in sig_nhoods) {
    nhood_cells <- colnames(milo_obj)[nhoods(milo_obj)[, i] == 1]
    cell_barcodes <- c(cell_barcodes, nhood_cells)
  }
  
  return(unique(cell_barcodes))
}


###################################################
# SPLIT UMAP/TSNE FROM SEURAT OBJECT
###################################################
plot_split_dimred <- function(obj_sub,
                              split_by  = "treatment",
                              color_by  = "fine_clust",
                              reduction = "umap",
                              pt_size   = 0.4,
                              title     = NULL) {
  # Faceted (split) UMAP or TSNE from a Seurat object.
  # Extracts embeddings directly from the Seurat reduction slot,
  # so coordinates are identical to those injected into Milo.
  
  emb <- Embeddings(obj_sub, reduction = reduction)
  
  # Determine axis labels from reduction type
  red_upper <- toupper(reduction)
  ax1 <- paste0(red_upper, "_1")
  ax2 <- paste0(red_upper, "_2")
  
  df <- data.frame(
    dim_1     = emb[, 1],
    dim_2     = emb[, 2],
    split_var = obj_sub@meta.data[[split_by]],
    color_var = obj_sub@meta.data[[color_by]],
    stringsAsFactors = FALSE
  )
  
  # Shuffle to avoid overplotting bias
  set.seed(42)
  df <- df[sample(nrow(df)), ]
  
  if (is.null(title)) title <- paste0(ncol(obj_sub), " cells")
  
  ggplot(df, aes(x = dim_1, y = dim_2, color = color_var)) +
    geom_point(size = pt_size, stroke = 0) +
    facet_wrap(~ split_var, ncol = 2) +
    labs(title = title, x = ax1, y = ax2, color = color_by) +
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


###################################################
# UNSPLIT (SINGLE) UMAP/TSNE COLORED BY METADATA
###################################################
plot_single_dimred <- function(obj_sub,
                               color_by  = "fine_clust",
                               reduction = "umap",
                               pt_size   = 0.4,
                               title     = NULL) {
  # Single UMAP or TSNE colored by a metadata column.
  
  emb <- Embeddings(obj_sub, reduction = reduction)
  red_upper <- toupper(reduction)
  ax1 <- paste0(red_upper, "_1")
  ax2 <- paste0(red_upper, "_2")
  
  df <- data.frame(
    dim_1     = emb[, 1],
    dim_2     = emb[, 2],
    color_var = obj_sub@meta.data[[color_by]],
    stringsAsFactors = FALSE
  )
  
  set.seed(42)
  df <- df[sample(nrow(df)), ]
  
  if (is.null(title)) title <- paste0(ncol(obj_sub), " cells")
  
  ggplot(df, aes(x = dim_1, y = dim_2, color = color_var)) +
    geom_point(size = pt_size, stroke = 0) +
    labs(title = title, x = ax1, y = ax2, color = color_by) +
    theme_void(base_size = 11) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.position = "right",
      legend.title    = element_text(size = 9),
      legend.text     = element_text(size = 8),
      legend.key.size = unit(0.35, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))
}


###################################################
# MILO DA NEIGHBOURHOOD GRAPH (CONSISTENT EMBEDDING)
###################################################
plot_milo_da <- function(milo_obj, da_results,
                         alpha  = 0.05,
                         layout = "UMAP",
                         title  = NULL) {
  # Plots the MiloR DA neighbourhood graph using layout="UMAP"
  # which pulls from the same Seurat UMAP injected during
  # run_milo_pipeline(). Same coordinate space as the cell UMAPs.
  
  p <- plotNhoodGraphDA(milo_obj, da_results,
                        layout = layout, alpha = alpha)
  
  if (!is.null(title)) p <- p + ggtitle(title)
  
  p + theme_void(base_size = 11) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    )
}


###################################################
# P-VALUE HISTOGRAM
###################################################
plot_pvalue_hist <- function(da_results, title = NULL) {
  
  if (is.null(title)) title <- "P-value Distribution"
  
  ggplot(da_results, aes(x = PValue)) +
    geom_histogram(bins = 50, fill = "#272E6A", color = "white") +
    theme_classic(base_size = 11) +
    labs(title = title, x = "P-value", y = "Count") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
    )
}


###################################################
# VOLCANO PLOT
###################################################
plot_volcano <- function(da_results, pval_thresh = 0.05, title = NULL) {
  
  if (is.null(title)) title <- "Volcano Plot"
  
  ggplot(da_results, aes(x = logFC, y = -log10(PValue))) +
    geom_point(aes(color = PValue < pval_thresh), size = 2) +
    scale_color_manual(
      values = c("grey60", "#D51F26"),
      labels = c("NS", paste0("P < ", pval_thresh)),
      name   = "Significance"
    ) +
    geom_hline(yintercept = -log10(pval_thresh),
               linetype = "dashed", color = "red") +
    theme_classic(base_size = 11) +
    labs(title = title, x = "Log Fold Change", y = "-log10(P-value)") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
    )
}


###################################################
# DA BEESWARM BY FINE_CLUST
###################################################
plot_da_beeswarm <- function(da_results, group_by = "fine_clust",
                             title = NULL, use_pvalue_as_fdr = TRUE) {
  # If use_pvalue_as_fdr = TRUE, copies PValue into SpatialFDR
  # before calling plotDAbeeswarm(). This is a deliberate workaround
  # when SpatialFDR is too conservative for the sample size.
  
  da_plot <- da_results
  
  if (use_pvalue_as_fdr) {
    da_plot$SpatialFDR <- da_plot$PValue
  }
  
  p <- plotDAbeeswarm(da_plot, group.by = group_by) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y  = element_text(size = 10),
      axis.title   = element_text(size = 12),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 13)
    )
  
  if (!is.null(title)) p <- p + ggtitle(title)
  
  return(p)
}


###################################################
# NEIGHBOURHOOD SIZE HISTOGRAM
###################################################
plot_nhood_size_hist <- function(milo_obj, title = NULL) {
  
  p <- plotNhoodSizeHist(milo_obj)
  
  subtitle_text <- paste("Total neighborhoods:", ncol(nhoods(milo_obj)),
                         "| Mean size:", round(mean(colSums(nhoods(milo_obj))), 1))
  
  if (!is.null(title)) {
    p <- p + labs(title = title, subtitle = subtitle_text)
  } else {
    p <- p + labs(subtitle = subtitle_text)
  }
  
  p + theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )
}


###################################################
# COMPOSE ROW FOR PATCHWORK
###################################################
compose_row <- function(..., widths = NULL, row_label = NULL) {
  # Compose multiple ggplots into a single patchwork row.
  
  plots <- list(...)
  n <- length(plots)
  if (is.null(widths)) widths <- rep(1, n)
  
  row <- Reduce(`+`, plots) + plot_layout(ncol = n, widths = widths)
  
  if (!is.null(row_label)) {
    row <- row + plot_annotation(
      title = row_label,
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))
    )
  }
  
  return(row)
}


###################################################
# SUMMARISE DA RESULTS BY CLUSTER
###################################################
summarise_da_by_cluster <- function(da_results, cluster_col = "fine_clust") {
  
  if (!cluster_col %in% colnames(da_results)) {
    warning(cluster_col, " not found in da_results. Returning NULL.")
    return(NULL)
  }
  
  da_results %>%
    group_by(.data[[cluster_col]]) %>%
    summarise(
      n_nhoods     = n(),
      n_sig_p0.05  = sum(PValue < 0.05, na.rm = TRUE),
      n_sig_fdr0.2 = sum(SpatialFDR < 0.2, na.rm = TRUE),
      n_up         = sum(PValue < 0.05 & logFC > 0, na.rm = TRUE),
      n_down       = sum(PValue < 0.05 & logFC < 0, na.rm = TRUE),
      mean_logFC   = mean(logFC, na.rm = TRUE),
      median_logFC = median(logFC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_sig_p0.05))
}