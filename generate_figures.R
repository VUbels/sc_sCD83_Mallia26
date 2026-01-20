library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# ============================================
# BLACKLIST FUNCTION
# ============================================

get_blacklist_genes <- function(seurat_obj) {
  
  all_genes <- rownames(seurat_obj)
  
  mt_genes <- grep(pattern = "^MT-", x = all_genes, value = TRUE)
  rps_genes <- grep(pattern = "^RPS", x = all_genes, value = TRUE)
  rpl_genes <- grep(pattern = "^RPL", x = all_genes, value = TRUE)
  
  s_genes <- cc.genes.updated.2019$s.genes
  g2m_genes <- cc.genes.updated.2019$g2m.genes
  
  sex_genes <- tryCatch({
    require(TxDb.Hsapiens.UCSC.hg38.knownGene)
    require(org.Hs.eg.db)
    require(GenomicFeatures)
    
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    geneGR <- GenomicFeatures::genes(txdb)
    sexGenesGR <- geneGR[seqnames(geneGR) %in% c("chrY", "chrX")]
    matchedGeneSymbols <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = sexGenesGR$gene_id,
      columns = c("ENTREZID", "SYMBOL"),
      keytype = "ENTREZID"
    )
    matchedGeneSymbols$SYMBOL[!is.na(matchedGeneSymbols$SYMBOL)]
  }, error = function(e) {
    message("Could not retrieve sex chromosome genes: ", e$message)
    character(0)
  })
  
  blacklist <- unique(c(mt_genes, rps_genes, rpl_genes, s_genes, g2m_genes, sex_genes))
  blacklist <- blacklist[blacklist %in% all_genes]
  
  message(paste0("Blacklist contains ", length(blacklist), " genes"))
  return(blacklist)
}

# ============================================
# CELL FILTERING FUNCTION
# ============================================

#' Filter out cells based on cluster/group labels
#' @param seurat_obj Seurat object
#' @param clusters_to_remove Character vector of cluster labels to remove (e.g., c("ORS.1", "ORS.4"))
#' @param cluster_col Metadata column containing cluster labels (e.g., "fine_clust")
#' @param cells_to_remove Optional: specific cell barcodes to remove directly
#' @param verbose Print filtering statistics

filter_cells <- function(
    seurat_obj,
    clusters_to_remove = NULL,
    cluster_col = "fine_clust",
    cells_to_remove = NULL,
    verbose = TRUE
) {
  
  n_before <- ncol(seurat_obj)
  
  if (verbose) {
    message("\n========== CELL FILTERING ==========")
    message(paste0("Cells before filtering: ", n_before))
  }
  
  cells_to_exclude <- character(0)
  
  # Filter by cluster/group labels
  if (!is.null(clusters_to_remove) && length(clusters_to_remove) > 0) {
    
    if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
      stop(paste0("Column '", cluster_col, "' not found in metadata.\nAvailable columns: ",
                  paste(colnames(seurat_obj@meta.data), collapse = ", ")))
    }
    
    if (verbose) {
      message(paste0("\nFiltering column: ", cluster_col))
      message("Distribution before filtering:")
      cluster_counts <- table(seurat_obj@meta.data[[cluster_col]])
      for (cl in names(cluster_counts)) {
        marker <- ifelse(cl %in% clusters_to_remove, " [REMOVING]", "")
        message(paste0("  ", cl, ": ", cluster_counts[cl], " cells", marker))
      }
    }
    
    # Find cells belonging to clusters to remove
    cells_in_clusters <- colnames(seurat_obj)[seurat_obj@meta.data[[cluster_col]] %in% clusters_to_remove]
    cells_to_exclude <- c(cells_to_exclude, cells_in_clusters)
    
    if (verbose) {
      message(paste0("\nClusters to remove: ", paste(clusters_to_remove, collapse = ", ")))
      message(paste0("Cells in these clusters: ", length(cells_in_clusters)))
    }
  }
  
  # Filter by specific cell barcodes
  if (!is.null(cells_to_remove) && length(cells_to_remove) > 0) {
    cells_found <- cells_to_remove[cells_to_remove %in% colnames(seurat_obj)]
    cells_to_exclude <- c(cells_to_exclude, cells_found)
    
    if (verbose) {
      message(paste0("\nSpecific cells to remove: ", length(cells_to_remove)))
      message(paste0("Found in object: ", length(cells_found)))
    }
  }
  
  cells_to_exclude <- unique(cells_to_exclude)
  
  # Perform filtering
  if (length(cells_to_exclude) > 0) {
    cells_to_keep <- setdiff(colnames(seurat_obj), cells_to_exclude)
    seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
  }
  
  n_after <- ncol(seurat_obj)
  
  if (verbose) {
    message(paste0("\n--- Filtering Summary ---"))
    message(paste0("Cells removed: ", n_before - n_after))
    message(paste0("Cells remaining: ", n_after))
    
    if (!is.null(clusters_to_remove) && cluster_col %in% colnames(seurat_obj@meta.data)) {
      message("\nDistribution after filtering:")
      cluster_counts <- table(seurat_obj@meta.data[[cluster_col]])
      for (cl in names(cluster_counts)) {
        message(paste0("  ", cl, ": ", cluster_counts[cl], " cells"))
      }
    }
    message("=========================================\n")
  }
  
  return(seurat_obj)
}


#' Show condition balance across cell types
#' @param seurat_obj Seurat object
#' @param condition_col Condition column (e.g., "treatment")
#' @param celltype_col Cell type column for DE analysis (e.g., "mapping_cell_type")
#' @param subcluster_col Optional: finer cluster column (e.g., "fine_clust")

show_sample_balance <- function(
    seurat_obj,
    condition_col = "treatment",
    celltype_col = "mapping_cell_type",
    subcluster_col = NULL
) {
  
  message("\n========== SAMPLE BALANCE REPORT ==========\n")
  
  # Cells per condition
  message("--- Cells per Condition ---")
  cond_counts <- table(seurat_obj@meta.data[[condition_col]])
  for (cond in names(cond_counts)) {
    message(paste0("  ", cond, ": ", cond_counts[cond], " cells"))
  }
  
  # Cells per condition per cell type
  message("\n--- Cells per Condition per Cell Type ---")
  ct_cond_table <- table(
    seurat_obj@meta.data[[celltype_col]],
    seurat_obj@meta.data[[condition_col]]
  )
  print(ct_cond_table)
  
  # Balance check
  message("\n--- Balance Check (max/min ratio) ---")
  for (ct in rownames(ct_cond_table)) {
    vals <- ct_cond_table[ct, ]
    if (all(vals > 0)) {
      ratio <- max(vals) / min(vals)
      flag <- ifelse(ratio > 5, " [IMBALANCED]", "")
      message(paste0("  ", ct, ": ", round(ratio, 2), "x", flag))
    } else {
      message(paste0("  ", ct, ": Missing condition [SKIP]"))
    }
  }
  
  # Subcluster breakdown if provided
  if (!is.null(subcluster_col) && subcluster_col %in% colnames(seurat_obj@meta.data)) {
    message("\n--- Subclusters per Condition ---")
    sub_cond_table <- table(
      seurat_obj@meta.data[[subcluster_col]],
      seurat_obj@meta.data[[condition_col]]
    )
    print(sub_cond_table)
  }
  
  message("\n============================================\n")
}


# ============================================
# DE ANALYSIS FUNCTION
# ============================================

run_de_analysis <- function(
    seurat_obj,
    condition_col = "treatment",
    celltype_col = "mapping_cell_type",
    ident.1 = "sCD83",
    ident.2 = "PBS",
    min.pct = 0.1,
    test.use = "wilcox",
    blacklist_genes = NULL
) {
  
  if (is.null(blacklist_genes)) blacklist_genes <- character(0)
  
  cell_types <- sort(unique(seurat_obj@meta.data[[celltype_col]]))
  cell_types <- cell_types[!is.na(cell_types)]
  
  message(paste0("Found ", length(cell_types), " cell types"))
  
  # --- Per cell type DE ---
  de_results_by_celltype <- list()
  
  for (ct in cell_types) {
    message(paste0("\nProcessing: ", ct))
    
    seurat_subset <- subset(seurat_obj, subset = !!sym(celltype_col) == ct)
    
    if (ncol(seurat_subset) < 10) {
      message("  Skipping - too few cells")
      next
    }
    
    conditions_present <- unique(seurat_subset@meta.data[[condition_col]])
    if (!(ident.1 %in% conditions_present) || !(ident.2 %in% conditions_present)) {
      message("  Skipping - missing condition")
      next
    }
    
    n1 <- sum(seurat_subset@meta.data[[condition_col]] == ident.1)
    n2 <- sum(seurat_subset@meta.data[[condition_col]] == ident.2)
    
    if (n1 < 3 || n2 < 3) {
      message("  Skipping - insufficient cells per condition")
      next
    }
    
    message(paste0("  ", ident.1, ": ", n1, ", ", ident.2, ": ", n2))
    
    Idents(seurat_subset) <- condition_col
    features_to_test <- setdiff(rownames(seurat_subset), blacklist_genes)
    
    tryCatch({
      de <- FindMarkers(
        seurat_subset,
        ident.1 = ident.1,
        ident.2 = ident.2,
        features = features_to_test,
        logfc.threshold = 0,
        min.pct = min.pct,
        test.use = test.use
      )
      de$gene <- rownames(de)
      de$cell_type <- ct
      de_results_by_celltype[[ct]] <- de
      
    }, error = function(e) {
      message(paste0("  Error: ", e$message))
    })
  }
  
  # --- Global DE (all cells) ---
  message("\n\nRunning global DE (all cells)...")
  Idents(seurat_obj) <- condition_col
  features_to_test <- setdiff(rownames(seurat_obj), blacklist_genes)
  
  n1_global <- sum(seurat_obj@meta.data[[condition_col]] == ident.1)
  n2_global <- sum(seurat_obj@meta.data[[condition_col]] == ident.2)
  message(paste0("Global: ", ident.1, ": ", n1_global, ", ", ident.2, ": ", n2_global))
  
  de_global <- FindMarkers(
    seurat_obj,
    ident.1 = ident.1,
    ident.2 = ident.2,
    features = features_to_test,
    logfc.threshold = 0,
    min.pct = min.pct,
    test.use = test.use
  )
  de_global$gene <- rownames(de_global)
  de_global$cell_type <- "All"
  
  return(list(
    by_celltype = de_results_by_celltype,
    global = de_global,
    combined = dplyr::bind_rows(de_results_by_celltype)
  ))
}

# ============================================
# PLOTTING FUNCTIONS
# ============================================

create_volcano_plot <- function(
    de_df,
    logfc_thresh = 0.5,
    pval_thresh = 0.05,
    title = "Volcano Plot",
    top_n_label = 20,
    cap_logfc = TRUE,
    cap_pval = TRUE,
    logfc_cap_quantile = 0.99,
    pval_cap_quantile = 0.99
) {
  
  de_df <- de_df %>%
    dplyr::mutate(
      neg_log10_pval = -log10(p_val_adj),
      significance = dplyr::case_when(
        p_val_adj < pval_thresh & avg_log2FC > logfc_thresh ~ "Up",
        p_val_adj < pval_thresh & avg_log2FC < -logfc_thresh ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # Store original values for labeling before capping
  de_df$orig_log2FC <- de_df$avg_log2FC
  de_df$orig_neg_log10_pval <- de_df$neg_log10_pval
  
  # --- Handle extreme values ---
  
  # Cap infinite p-values
  if (any(is.infinite(de_df$neg_log10_pval))) {
    max_finite <- max(de_df$neg_log10_pval[is.finite(de_df$neg_log10_pval)])
    de_df$neg_log10_pval[is.infinite(de_df$neg_log10_pval)] <- max_finite * 1.1
  }
  
  # Cap extreme log2FC values (outliers)
  if (cap_logfc) {
    logfc_upper <- quantile(abs(de_df$avg_log2FC), logfc_cap_quantile, na.rm = TRUE)
    logfc_cap <- max(logfc_upper, logfc_thresh * 2)  # At least 2x threshold
    
    de_df <- de_df %>%
      dplyr::mutate(
        is_logfc_capped = abs(avg_log2FC) > logfc_cap,
        avg_log2FC_plot = dplyr::case_when(
          avg_log2FC > logfc_cap ~ logfc_cap,
          avg_log2FC < -logfc_cap ~ -logfc_cap,
          TRUE ~ avg_log2FC
        )
      )
  } else {
    de_df$avg_log2FC_plot <- de_df$avg_log2FC
    de_df$is_logfc_capped <- FALSE
  }
  
  # Cap extreme p-values (outliers)
  if (cap_pval) {
    pval_upper <- quantile(de_df$neg_log10_pval[is.finite(de_df$neg_log10_pval)], 
                           pval_cap_quantile, na.rm = TRUE)
    pval_cap <- max(pval_upper, -log10(pval_thresh) * 2)  # At least 2x threshold line
    
    de_df <- de_df %>%
      dplyr::mutate(
        is_pval_capped = neg_log10_pval > pval_cap,
        neg_log10_pval_plot = dplyr::if_else(neg_log10_pval > pval_cap, pval_cap, neg_log10_pval)
      )
  } else {
    de_df$neg_log10_pval_plot <- de_df$neg_log10_pval
    de_df$is_pval_capped <- FALSE
  }
  
  # Mark capped points (for different shape)
  de_df$is_capped <- de_df$is_logfc_capped | de_df$is_pval_capped
  
  # --- Select top genes to label: top N UP + top N DOWN ---
  
  top_up <- de_df %>%
    dplyr::filter(significance == "Up") %>%
    dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
    dplyr::slice_head(n = top_n_label)
  
  top_down <- de_df %>%
    dplyr::filter(significance == "Down") %>%
    dplyr::arrange(p_val_adj, avg_log2FC) %>%
    dplyr::slice_head(n = top_n_label)
  
  top_genes <- dplyr::bind_rows(top_up, top_down)
  
  # --- Counts ---
  n_up <- sum(de_df$significance == "Up")
  n_down <- sum(de_df$significance == "Down")
  n_capped <- sum(de_df$is_capped)
  
  # --- Build plot ---
  p <- ggplot(de_df, aes(x = avg_log2FC_plot, y = neg_log10_pval_plot)) +
    # Non-capped points (circles)
    geom_point(
      data = dplyr::filter(de_df, !is_capped),
      aes(color = significance),
      alpha = 0.5,
      size = 1.2
    ) +
    # Capped points (triangles)
    geom_point(
      data = dplyr::filter(de_df, is_capped),
      aes(color = significance),
      alpha = 0.7,
      size = 1.5,
      shape = 17  # Triangle
    ) +
    scale_color_manual(
      values = c("Up" = "#D51F26", "Down" = "#272E6A", "NS" = "grey70"),
      labels = c(
        "Up" = paste0("Up (", n_up, ")"),
        "Down" = paste0("Down (", n_down, ")"),
        "NS" = "NS"
      )
    ) +
    geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 2.5,
      max.overlaps = 30,
      segment.color = "grey50",
      segment.alpha = 0.6,
      min.segment.length = 0.2,
      box.padding = 0.3,
      point.padding = 0.2
    ) +
    labs(
      title = title,
      x = ifelse(cap_logfc && any(de_df$is_logfc_capped), 
                 "log2 Fold Change (capped)", "log2 Fold Change"),
      y = ifelse(cap_pval && any(de_df$is_pval_capped),
                 "-log10(adj. p-value) (capped)", "-log10(adj. p-value)"),
      color = "",
      caption = ifelse(n_capped > 0, 
                       paste0("▲ = ", n_capped, " capped outliers"), "")
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      legend.position = "top",
      plot.caption = element_text(hjust = 0, size = 8, color = "grey40")
    )
  
  # Symmetric x-axis
  max_x <- max(abs(de_df$avg_log2FC_plot), na.rm = TRUE) * 1.05
  p <- p + xlim(-max_x, max_x)
  
  return(p)
}


create_density_plot <- function(
    seurat_obj,
    genes,
    condition_col,
    ident.1,
    ident.2,
    title = "Expression Distribution"
) {
  
  expr_data <- GetAssayData(seurat_obj, layer = "data")
  genes <- genes[genes %in% rownames(expr_data)]
  
  if (length(genes) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No genes") + theme_void())
  }
  
  expr_mat <- expr_data[genes, , drop = FALSE]
  mean_expr <- Matrix::colMeans(expr_mat)
  
  plot_df <- data.frame(
    expression = mean_expr,
    condition = seurat_obj@meta.data[[condition_col]]
  ) %>%
    dplyr::filter(condition %in% c(ident.1, ident.2))
  
  ggplot(plot_df, aes(x = expression, fill = condition)) +
    geom_density(alpha = 0.6, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = setNames(c("#D51F26", "#272E6A"), c(ident.1, ident.2))) +
    labs(title = title, x = "Mean Expression (DE genes)", y = "Density", fill = "") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
          legend.position = "top")
}


create_heatmap <- function(
    seurat_obj,
    genes,
    condition_col,
    ident.1,
    ident.2,
    max_cells = 150,
    title = "Heatmap"
) {
  
  expr_data <- GetAssayData(seurat_obj, layer = "data")
  genes <- genes[genes %in% rownames(expr_data)]
  
  if (length(genes) == 0) return(NULL)
  
  cells_1 <- colnames(seurat_obj)[seurat_obj@meta.data[[condition_col]] == ident.1]
  cells_2 <- colnames(seurat_obj)[seurat_obj@meta.data[[condition_col]] == ident.2]
  
  if (length(cells_1) > max_cells) cells_1 <- sample(cells_1, max_cells)
  if (length(cells_2) > max_cells) cells_2 <- sample(cells_2, max_cells)
  
  cells_ordered <- c(cells_1, cells_2)
  expr_mat <- as.matrix(expr_data[genes, cells_ordered])
  
  expr_scaled <- t(scale(t(expr_mat)))
  expr_scaled[is.nan(expr_scaled)] <- 0
  expr_scaled[expr_scaled > 2] <- 2
  expr_scaled[expr_scaled < -2] <- -2
  
  condition_vec <- c(rep(ident.1, length(cells_1)), rep(ident.2, length(cells_2)))
  
  ha <- HeatmapAnnotation(
    Condition = condition_vec,
    col = list(Condition = setNames(c("#D51F26", "#272E6A"), c(ident.1, ident.2))),
    show_annotation_name = FALSE,
    simple_anno_size = unit(3, "mm")
  )
  
  col_fun <- colorRamp2(c(-2, 0, 2), c("#272E6A", "white", "#D51F26"))
  
  Heatmap(
    expr_scaled,
    name = "Z-score",
    col = col_fun,
    top_annotation = ha,
    column_split = factor(condition_vec, levels = c(ident.1, ident.2)),
    column_gap = unit(1, "mm"),
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_column_names = FALSE,
    show_row_names = length(genes) <= 40,
    row_names_gp = gpar(fontsize = 6),
    column_title = title,
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    border = TRUE,
    use_raster = TRUE
  )
}


# ============================================
# MAIN WRAPPER FUNCTION
# ============================================

run_full_de_pipeline <- function(
    seurat_obj,
    condition_col = "treatment",
    celltype_col = "mapping_cell_type",
    ident.1 = "sCD83",
    ident.2 = "PBS",
    blacklist_genes = NULL,
    clusters_to_remove = NULL,
    cluster_col = "fine_clust",
    cells_to_remove = NULL,
    top_n_genes = 50,
    logfc_thresh = 0.5,
    pval_thresh = 0.05,
    output_dir = "DE_results"
) {
  
  # Create output directory
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # ============================================
  # Step 0: Filter cells/clusters
  # ============================================
  
  if (!is.null(clusters_to_remove) || !is.null(cells_to_remove)) {
    seurat_obj <- filter_cells(
      seurat_obj,
      clusters_to_remove = clusters_to_remove,
      cluster_col = cluster_col,
      cells_to_remove = cells_to_remove,
      verbose = TRUE
    )
  }
  
  # Show balance after filtering
  show_sample_balance(
    seurat_obj,
    condition_col = condition_col,
    celltype_col = celltype_col,
    subcluster_col = cluster_col
  )
  
  # Get blacklist if not provided
  if (is.null(blacklist_genes)) {
    blacklist_genes <- get_blacklist_genes(seurat_obj)
  }
  
  # ============================================
  # Step 1: Run DE analysis
  # ============================================
  message("\n========== RUNNING DE ANALYSIS ==========\n")
  
  de_results <- run_de_analysis(
    seurat_obj = seurat_obj,
    condition_col = condition_col,
    celltype_col = celltype_col,
    ident.1 = ident.1,
    ident.2 = ident.2,
    blacklist_genes = blacklist_genes
  )
  
  # Save DE tables
  write.csv(de_results$global, file.path(output_dir, "DE_global_all_cells.csv"), row.names = FALSE)
  write.csv(de_results$combined, file.path(output_dir, "DE_all_celltypes_combined.csv"), row.names = FALSE)
  
  for (ct in names(de_results$by_celltype)) {
    ct_clean <- gsub("[^A-Za-z0-9]", "_", ct)
    write.csv(de_results$by_celltype[[ct]], 
              file.path(output_dir, paste0("DE_", ct_clean, ".csv")), 
              row.names = FALSE)
  }
  
  # Store all plots
  all_plots <- list()
  
  # ============================================
  # Step 2: Global plots (all cells)
  # ============================================
  message("\n========== GENERATING GLOBAL PLOTS ==========\n")
  
  sig_genes_global <- de_results$global %>%
    dplyr::filter(p_val_adj < pval_thresh, abs(avg_log2FC) > logfc_thresh) %>%
    dplyr::arrange(desc(abs(avg_log2FC))) %>%
    dplyr::slice_head(n = top_n_genes) %>%
    dplyr::pull(gene)
  
  message(paste0("Global significant genes: ", length(sig_genes_global)))
  
  volcano_global <- create_volcano_plot(
    de_results$global,
    logfc_thresh = logfc_thresh,
    pval_thresh = pval_thresh,
    title = paste0("Global: ", ident.1, " vs ", ident.2)
  )
  
  density_global <- create_density_plot(
    seurat_obj,
    genes = sig_genes_global,
    condition_col = condition_col,
    ident.1 = ident.1,
    ident.2 = ident.2,
    title = "Global DE Gene Expression"
  )
  
  heatmap_global <- create_heatmap(
    seurat_obj,
    genes = sig_genes_global,
    condition_col = condition_col,
    ident.1 = ident.1,
    ident.2 = ident.2,
    title = "Global"
  )
  
  all_plots[["Global"]] <- list(
    volcano = volcano_global,
    density = density_global,
    heatmap = heatmap_global,
    sig_genes = sig_genes_global
  )
  
  # Save global plots
  ggsave(file.path(output_dir, "Global_volcano.pdf"), volcano_global, width = 7, height = 6)
  ggsave(file.path(output_dir, "Global_density.pdf"), density_global, width = 6, height = 5)
  
  combined_global <- density_global + volcano_global + plot_layout(ncol = 2)
  ggsave(file.path(output_dir, "Global_density_volcano.pdf"), combined_global, width = 12, height = 5)
  ggsave(file.path(output_dir, "Global_density_volcano.png"), combined_global, width = 12, height = 5, dpi = 150)
  
  if (!is.null(heatmap_global)) {
    pdf(file.path(output_dir, "Global_heatmap.pdf"), width = 10, height = 8)
    draw(heatmap_global)
    dev.off()
    
    png(file.path(output_dir, "Global_heatmap.png"), width = 1000, height = 800, res = 100)
    draw(heatmap_global)
    dev.off()
  }
  
  # ============================================
  # Step 3: Per cell type plots
  # ============================================
  message("\n========== GENERATING PER-CELLTYPE PLOTS ==========\n")
  
  for (ct in names(de_results$by_celltype)) {
    message(paste0("\nPlotting: ", ct))
    
    ct_clean <- gsub("[^A-Za-z0-9]", "_", ct)
    de_ct <- de_results$by_celltype[[ct]]
    
    seurat_ct <- subset(seurat_obj, subset = !!sym(celltype_col) == ct)
    
    sig_genes_ct <- de_ct %>%
      dplyr::filter(p_val_adj < pval_thresh, abs(avg_log2FC) > logfc_thresh) %>%
      dplyr::arrange(desc(abs(avg_log2FC))) %>%
      dplyr::slice_head(n = top_n_genes) %>%
      dplyr::pull(gene)
    
    message(paste0("  Significant genes: ", length(sig_genes_ct)))
    
    # FIX 1: Use de_ct (not de_results$global)
    # FIX 2: Name variable volcano_ct (not volcano)
    volcano_ct <- create_volcano_plot(
      de_ct,
      logfc_thresh = logfc_thresh,
      pval_thresh = pval_thresh,
      top_n_label = 20,
      title = paste0(ct, ": ", ident.1, " vs ", ident.2)
    )
    
    density_ct <- create_density_plot(
      seurat_ct,
      genes = sig_genes_ct,
      condition_col = condition_col,
      ident.1 = ident.1,
      ident.2 = ident.2,
      title = paste0(ct, " DE Genes")
    )
    
    heatmap_ct <- NULL
    if (length(sig_genes_ct) > 0) {
      heatmap_ct <- create_heatmap(
        seurat_ct,
        genes = sig_genes_ct,
        condition_col = condition_col,
        ident.1 = ident.1,
        ident.2 = ident.2,
        title = ct
      )
    }
    
    all_plots[[ct]] <- list(
      volcano = volcano_ct,
      density = density_ct,
      heatmap = heatmap_ct,
      sig_genes = sig_genes_ct
    )
    
    # Save
    ggsave(file.path(output_dir, paste0(ct_clean, "_volcano.pdf")), volcano_ct, width = 7, height = 6)
    ggsave(file.path(output_dir, paste0(ct_clean, "_density.pdf")), density_ct, width = 6, height = 5)
    
    combined_ct <- density_ct + volcano_ct + plot_layout(ncol = 2)
    ggsave(file.path(output_dir, paste0(ct_clean, "_density_volcano.pdf")), combined_ct, width = 12, height = 5)
    ggsave(file.path(output_dir, paste0(ct_clean, "_density_volcano.png")), combined_ct, width = 12, height = 5, dpi = 150)
    
    if (!is.null(heatmap_ct)) {
      pdf(file.path(output_dir, paste0(ct_clean, "_heatmap.pdf")), width = 10, height = 8)
      draw(heatmap_ct)
      dev.off()
      
      png(file.path(output_dir, paste0(ct_clean, "_heatmap.png")), width = 1000, height = 800, res = 100)
      draw(heatmap_ct)
      dev.off()
    }
  }
  
  # ============================================
  # Step 4: Summary
  # ============================================
  summary_df <- data.frame(
    cell_type = c("Global", names(de_results$by_celltype))
  )
  
  summary_df$n_sig_up <- sapply(summary_df$cell_type, function(ct) {
    if (ct == "Global") {
      sum(de_results$global$p_val_adj < pval_thresh & de_results$global$avg_log2FC > logfc_thresh)
    } else {
      sum(de_results$by_celltype[[ct]]$p_val_adj < pval_thresh & 
            de_results$by_celltype[[ct]]$avg_log2FC > logfc_thresh)
    }
  })
  
  summary_df$n_sig_down <- sapply(summary_df$cell_type, function(ct) {
    if (ct == "Global") {
      sum(de_results$global$p_val_adj < pval_thresh & de_results$global$avg_log2FC < -logfc_thresh)
    } else {
      sum(de_results$by_celltype[[ct]]$p_val_adj < pval_thresh & 
            de_results$by_celltype[[ct]]$avg_log2FC < -logfc_thresh)
    }
  })
  
  summary_df$n_sig_total <- summary_df$n_sig_up + summary_df$n_sig_down
  
  write.csv(summary_df, file.path(output_dir, "DE_summary.csv"), row.names = FALSE)
  
  message("\n========== DONE ==========")
  message(paste0("Results saved to: ", output_dir))
  print(summary_df)
  
  return(list(
    de_results = de_results,
    plots = all_plots,
    summary = summary_df,
    seurat_filtered = seurat_obj
  ))
}

###########################################################################
# MAIN FUNCTION TO RUN
##########################################################################

obj <- readRDS("/mnt/d/scrna_output/mallia_25/Mallia_25_LogN_Annotated.RDS")

# Clusters to remove from fine_clust column (imbalanced between conditions)
clusters_to_remove <- c("ORS.1", "ORS.4", "ORS.5", "ORS.6", "Immune-associated.KCs", "Immune-associated.FBs")

# Check balance BEFORE running (optional)
show_sample_balance(
  obj,
  condition_col = "treatment",
  celltype_col = "mapping_cell_type",  # Broad cell types (keratinocytes, endothelial, etc.)
  subcluster_col = "fine_clust"        # Fine clusters (ORS.1, ORS.2, etc.)
)

# Run the full pipeline
results <- run_full_de_pipeline(
  seurat_obj = obj,
  condition_col = "treatment",
  celltype_col = "mapping_cell_type",     # DE analysis grouped by this (broad types)
  ident.1 = "sCD83",
  ident.2 = "PBS",
  blacklist_genes = NULL,
  clusters_to_remove = clusters_to_remove, # Remove these fine_clust labels
  cluster_col = "fine_clust",              # Column containing ORS.1, ORS.2, etc.
  cells_to_remove = NULL,
  top_n_genes = 2000,
  logfc_thresh = 0.5,
  pval_thresh = 0.05,
  output_dir = "DE_results"
)

# View results
print(results$plots$Global$volcano)
draw(results$plots$Global$heatmap)

# View specific cell type (use exact name from mapping_cell_type)
print(results$plots[["Keratinocyte"]]$volcano)