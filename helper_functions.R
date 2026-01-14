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