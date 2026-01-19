###################################################
# CUSTOM PLOTTING FUNCTION - P-VALUES INSTEAD OF FDR
###################################################

plotNhoodGraphDA_pval <- function(x, milo_res, alpha = 0.1, res_column = "logFC", 
                                  use_pvalue = TRUE, layout = "UMAP", ...) {
  
  # Check if neighborhood graph exists
  if(is.null(nhoodGraph(x)) || length(igraph::E(nhoodGraph(x))) == 0){
    stop("Not a valid Milo object - neighbourhood graph is missing. Please run buildNhoodGraph() first.")
  }
  
  # Check if layout is valid
  if (is.character(layout)) {
    if (!layout %in% names(reducedDims(x))) {
      stop(layout, " is not in reducedDim(x) - choose a different layout")
    }
  }
  
  ## Add milo results to colData
  signif_res <- milo_res
  
  # Use PValue instead of SpatialFDR if specified
  if(use_pvalue) {
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
  if(any(names(list(...)) %in% c("subset.nhoods"))){
    subset.nhoods <- list(...)$subset.nhoods
    sub.indices <- nhoodIndex(x)[subset.nhoods]
    colData(x)[unlist(sub.indices[signif_res$Nhood]), res_column] <- signif_res[,res_column]
  } else{
    colData(x)[unlist(nhoodIndex(x)[signif_res$Nhood]), res_column] <- signif_res[,res_column]
  }
  
  # Check for res_column in graph vertex attributes
  g_atts <- names(igraph::vertex_attr(nhoodGraph(x)))
  if(isFALSE(res_column %in% g_atts)){
    message("Adding nhood effect sizes to neighbourhood graph attributes")
    
    if(any(names(list(...)) %in% c("subset.nhoods"))){
      nh.v <- igraph::V(nhoodGraph(x))
      drop.v <- setdiff(nh.v, sub.indices)
      nhgraph <- nhoodGraph(x)
      nhgraph <- igraph::subgraph(nhgraph, sub.indices)
      nhgraph <- igraph::set_vertex_attr(nhgraph,
                                         name = res_column, value = signif_res[, res_column])
      nhoodGraph(x) <- nhgraph
    } else{
      nhoodGraph(x) <- igraph::set_vertex_attr(nhoodGraph(x), 
                                               name = res_column, 
                                               value = signif_res[, res_column])
    }
  }
  
  ## Plot logFC - pass layout explicitly
  plotNhoodGraph(x, colour_by = res_column, layout = layout, is.da = TRUE, ...) +
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
  
  if(use_pvalue) {
    sig_col <- "PValue"
  } else {
    sig_col <- "SpatialFDR"
  }
  
  # Filter significant neighborhoods
  if(direction == "up") {
    sig_nhoods <- da_results$Nhood[da_results[[sig_col]] < alpha & da_results$logFC > 0]
  } else if(direction == "down") {
    sig_nhoods <- da_results$Nhood[da_results[[sig_col]] < alpha & da_results$logFC < 0]
  } else {
    sig_nhoods <- da_results$Nhood[da_results[[sig_col]] < alpha]
  }
  
  # Extract cells from significant neighborhoods
  cell_barcodes <- c()
  for(i in sig_nhoods) {
    nhood_cells <- colnames(milo_obj)[nhoods(milo_obj)[, i] == 1]
    cell_barcodes <- c(cell_barcodes, nhood_cells)
  }
  
  return(unique(cell_barcodes))
}


# ============================================================================
# ENHANCED SPATIAL OVERLAP WITH EXPRESSION THRESHOLDS
# ============================================================================

plot_spatial_overlap_thresholded <- function(spatial_obj, gene1, gene2,
                                             sigma = 50,
                                             assay = "SCT",
                                             overlap_method = "geometric_mean",
                                             # NEW THRESHOLD PARAMETERS
                                             min_expr_percentile = 20,  # Only consider top 80% of expression
                                             min_expr_absolute = NULL,   # Optional: absolute cutoff (e.g., 0.5)
                                             min_overlap_percentile = 10, # Only show top 90% of overlap
                                             require_both_genes = TRUE,   # Both genes must pass threshold
                                             output_path = NULL,
                                             pt.size.factor = 1.5,
                                             show_diagnostics = TRUE) {
  "
  Visualize spatial overlap with expression thresholds to remove background.
  
  NEW Parameters:
  ---------------
  min_expr_percentile : numeric (0-100)
      Percentile threshold for individual gene expression.
      e.g., 20 means only spots in top 80% of expression are considered positive.
      Higher = stricter (less background, fewer false positives)
      Recommended: 10-30
      
  min_expr_absolute : numeric or NULL
      Absolute expression threshold (applied AFTER smoothing).
      e.g., 0.5 means smoothed expression must be > 0.5
      NULL = use percentile only
      
  min_overlap_percentile : numeric (0-100)  
      Percentile threshold for final overlap scores.
      Only display overlap values above this percentile.
      Recommended: 5-20
      
  require_both_genes : logical
      If TRUE, both genes must exceed threshold for overlap to be counted.
      If FALSE, overlap is calculated then thresholded separately.
  "
  
  # Get expression data
  assay_data <- GetAssayData(spatial_obj, assay = assay, layer = "data")
  
  # Check genes exist
  if(!gene1 %in% rownames(assay_data)) {
    stop(sprintf("Gene '%s' not found in assay '%s'", gene1, assay))
  }
  if(!gene2 %in% rownames(assay_data)) {
    stop(sprintf("Gene '%s' not found in assay '%s'", gene2, assay))
  }
  
  # Apply Gaussian smoothing
  smoothed <- gaussian_smooth_spatial(spatial_obj, c(gene1, gene2), 
                                      sigma = sigma, assay = assay)
  
  gene1_expr <- as.numeric(smoothed[gene1, ])
  gene2_expr <- as.numeric(smoothed[gene2, ])
  
  # ========================================================================
  # APPLY EXPRESSION THRESHOLDS
  # ========================================================================
  
  # Calculate percentile thresholds for each gene
  gene1_threshold <- quantile(gene1_expr[gene1_expr > 0], 
                              probs = min_expr_percentile/100, 
                              na.rm = TRUE)
  gene2_threshold <- quantile(gene2_expr[gene2_expr > 0], 
                              probs = min_expr_percentile/100, 
                              na.rm = TRUE)
  
  # Override with absolute threshold if provided
  if(!is.null(min_expr_absolute)) {
    gene1_threshold <- max(gene1_threshold, min_expr_absolute)
    gene2_threshold <- max(gene2_threshold, min_expr_absolute)
  }
  
  if(show_diagnostics) {
    message("\n=== EXPRESSION THRESHOLDS ===")
    message(sprintf("%s threshold: %.4f (raw range: %.4f - %.4f)", 
                    gene1, gene1_threshold, min(gene1_expr), max(gene1_expr)))
    message(sprintf("%s threshold: %.4f (raw range: %.4f - %.4f)", 
                    gene2, gene2_threshold, min(gene2_expr), max(gene2_expr)))
    message(sprintf("%s spots above threshold: %d / %d (%.1f%%)", 
                    gene1, sum(gene1_expr > gene1_threshold), 
                    length(gene1_expr),
                    100 * sum(gene1_expr > gene1_threshold) / length(gene1_expr)))
    message(sprintf("%s spots above threshold: %d / %d (%.1f%%)", 
                    gene2, sum(gene2_expr > gene2_threshold), 
                    length(gene2_expr),
                    100 * sum(gene2_expr > gene2_threshold) / length(gene2_expr)))
  }
  
  # Create thresholded expression vectors
  gene1_thresh <- gene1_expr
  gene2_thresh <- gene2_expr
  
  gene1_thresh[gene1_expr < gene1_threshold] <- 0
  gene2_thresh[gene2_expr < gene2_threshold] <- 0
  
  # Calculate overlap using thresholded expression
  if(require_both_genes) {
    # Both genes must be positive for overlap to exist
    overlap_raw <- ifelse(gene1_thresh > 0 & gene2_thresh > 0,
                          switch(overlap_method,
                                 "geometric_mean" = sqrt(gene1_thresh * gene2_thresh),
                                 "minimum" = pmin(gene1_thresh, gene2_thresh),
                                 "product" = gene1_thresh * gene2_thresh,
                                 "average" = (gene1_thresh + gene2_thresh) / 2,
                                 stop("Invalid overlap_method")),
                          0)
  } else {
    # Calculate overlap, then threshold
    overlap_raw <- switch(overlap_method,
                          "geometric_mean" = sqrt(gene1_thresh * gene2_thresh),
                          "minimum" = pmin(gene1_thresh, gene2_thresh),
                          "product" = gene1_thresh * gene2_thresh,
                          "average" = (gene1_thresh + gene2_thresh) / 2,
                          stop("Invalid overlap_method")
    )
  }
  
  # Apply overlap percentile threshold
  nonzero_overlap <- overlap_raw[overlap_raw > 0]
  if(length(nonzero_overlap) > 0) {
    overlap_threshold <- quantile(nonzero_overlap, 
                                  probs = min_overlap_percentile/100, 
                                  na.rm = TRUE)
    overlap_raw[overlap_raw < overlap_threshold] <- 0
  } else {
    overlap_threshold <- 0
  }
  
  # Align with Seurat object
  names(overlap_raw) <- colnames(smoothed)
  overlap_aligned <- overlap_raw[Cells(spatial_obj)]
  overlap_aligned[is.na(overlap_aligned)] <- 0
  
  # Calculate quantiles for color scaling (only from positive values)
  nonzero_vals <- overlap_aligned[overlap_aligned > 0]
  
  if(length(nonzero_vals) == 0) {
    warning("No overlap detected after thresholding. Try lower thresholds.")
    q01 <- 0
    q50 <- 0
    q99 <- max(overlap_aligned)
  } else {
    q01 <- quantile(nonzero_vals, 0.01, na.rm = TRUE)
    q50 <- quantile(nonzero_vals, 0.50, na.rm = TRUE)
    q99 <- quantile(nonzero_vals, 0.99, na.rm = TRUE)
  }
  
  # Diagnostics
  if(show_diagnostics) {
    message("\n=== OVERLAP DIAGNOSTICS ===")
    message(sprintf("Overlap threshold: %.4f", overlap_threshold))
    message(sprintf("Overlap range: %.4f to %.4f", 
                    min(overlap_aligned), max(overlap_aligned)))
    message(sprintf("Spots with overlap > 0: %d / %d (%.1f%%)", 
                    sum(overlap_aligned > 0), 
                    length(overlap_aligned),
                    100 * sum(overlap_aligned > 0) / length(overlap_aligned)))
    message(sprintf("Quantiles - 1%%: %.4f, 50%%: %.4f, 99%%: %.4f", 
                    q01, q50, q99))
    message("============================\n")
  }
  
  # Add to metadata
  spatial_obj$temp_overlap <- overlap_aligned
  
  # Create spatial plot
  # Use discrete colors to emphasize true overlap regions
  p <- SpatialFeaturePlot(
    spatial_obj, 
    features = "temp_overlap",
    image.alpha = 0,
    pt.size.factor = pt.size.factor
  ) +
    scale_fill_gradientn(
      colors = c("grey95", "grey90", "#FFCCCC", "#FF9999", 
                 "#FF5555", "red", "#CC0000", "darkred"),
      values = scales::rescale(c(0, 0.001, q01, q01*1.5, q50*0.8, 
                                 q50*1.2, q99*0.9, q99)),
      limits = c(0, q99),
      na.value = "white",
      name = "Overlap\nScore"
    ) +
    ggtitle(paste0(gene1, " & ", gene2, " spatial overlap"),
            subtitle = sprintf("σ=%dμm | Gene thresh=%d%% | Overlap thresh=%d%%", 
                               sigma, min_expr_percentile, min_overlap_percentile)) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.key.height = unit(1.2, "cm"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 9)
    )
  
  # Clean up
  spatial_obj$temp_overlap <- NULL
  
  # Save if requested
  if(!is.null(output_path)) {
    ggsave(output_path, p, width = 10, height = 10, dpi = 300, bg = "white")
    message(sprintf("Plot saved to: %s", output_path))
  }
  
  return(list(
    plot = p,
    overlap_scores = overlap_aligned,
    gene1_smoothed = gene1_expr,
    gene2_smoothed = gene2_expr,
    gene1_thresholded = gene1_thresh,
    gene2_thresholded = gene2_thresh,
    thresholds = list(
      gene1 = gene1_threshold,
      gene2 = gene2_threshold,
      overlap = overlap_threshold
    ),
    sigma = sigma,
    overlap_method = overlap_method
  ))
}

# ============================================================================
# SIDE-BY-SIDE COMPARISON: RAW vs THRESHOLDED
# ============================================================================

compare_thresholding <- function(spatial_obj, gene1, gene2,
                                 sigma = 50,
                                 threshold_levels = c(0, 10, 20, 30)) {
  "
  Compare overlap with different threshold stringency levels.
  
  Parameters:
  -----------
  threshold_levels : numeric vector
      Percentile thresholds to test (0 = no threshold, 30 = very strict)
  "
  
  library(patchwork)
  
  plots <- list()
  
  for(thresh in threshold_levels) {
    message(sprintf("\n--- Testing threshold = %d%% ---", thresh))
    
    result <- plot_spatial_overlap_thresholded(
      spatial_obj = spatial_obj,
      gene1 = gene1,
      gene2 = gene2,
      sigma = sigma,
      min_expr_percentile = thresh,
      min_overlap_percentile = max(5, thresh/2),  # Scale overlap threshold too
      show_diagnostics = TRUE,
      pt.size.factor = 1.2
    )
    
    label <- ifelse(thresh == 0, "No threshold", sprintf("%d%% threshold", thresh))
    plots[[as.character(thresh)]] <- result$plot + 
      ggtitle(label)
  }
  
  combined <- wrap_plots(plots, ncol = 2)
  
  return(combined)
}

# ============================================================================
# ADAPTIVE THRESHOLDING BASED ON DISTRIBUTION
# ============================================================================

plot_spatial_overlap_adaptive <- function(spatial_obj, gene1, gene2,
                                          sigma = 50,
                                          assay = "SCT",
                                          overlap_method = "geometric_mean",
                                          sensitivity = "medium",  # "low", "medium", "high"
                                          output_path = NULL,
                                          pt.size.factor = 1.5) {
  "
  Automatically determine thresholds based on expression distribution.
  
  Parameters:
  -----------
  sensitivity : character
      - -low: Strict thresholds (30th percentile, fewer false positives)
      - -medium: Moderate thresholds (20th percentile) [DEFAULT]
      - -high: Lenient thresholds (10th percentile, more sensitive)
  "
  
  # Map sensitivity to threshold values
  threshold_map <- list(
    "low" = list(expr = 30, overlap = 15),
    "medium" = list(expr = 20, overlap = 10),
    "high" = list(expr = 10, overlap = 5)
  )
  
  if(!sensitivity %in% names(threshold_map)) {
    stop("sensitivity must be 'low', 'medium', or 'high'")
  }
  
  thresholds <- threshold_map[[sensitivity]]
  
  message(sprintf("Using %s sensitivity (expr: %d%%, overlap: %d%%)", 
                  sensitivity, thresholds$expr, thresholds$overlap))
  
  result <- plot_spatial_overlap_thresholded(
    spatial_obj = spatial_obj,
    gene1 = gene1,
    gene2 = gene2,
    sigma = sigma,
    assay = assay,
    overlap_method = overlap_method,
    min_expr_percentile = thresholds$expr,
    min_overlap_percentile = thresholds$overlap,
    require_both_genes = TRUE,
    output_path = output_path,
    pt.size.factor = pt.size.factor,
    show_diagnostics = TRUE
  )
  
  return(result)
}
