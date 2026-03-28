library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(DESeq2)

# ============================================
# EXPANDED BLACKLIST FUNCTION
# ============================================

get_blacklist_genes <- function(seurat_obj) {
  
  all_genes <- rownames(seurat_obj)
  
  # Mitochondrial
  mt_genes <- grep(pattern = "^MT-", x = all_genes, value = TRUE)
  
  # Ribosomal
  rps_genes <- grep(pattern = "^RPS", x = all_genes, value = TRUE)
  rpl_genes <- grep(pattern = "^RPL", x = all_genes, value = TRUE)
  
  # Cell cycle
  s_genes <- cc.genes.updated.2019$s.genes
  g2m_genes <- cc.genes.updated.2019$g2m.genes
  
  # Hemoglobin genes (common contaminant)
  hb_genes <- grep(pattern = "^HB[^(P)]", x = all_genes, value = TRUE)
  
  # MALAT1 and other high-abundance non-coding RNAs
  noncoding_genes <- intersect(
    c("MALAT1", "NEAT1", "XIST", "TSIX", "KCNQ1OT1", "MEG3"),
    all_genes
  )
  
  # Sex chromosome genes (optional, retrieve if annotation packages available)
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
  
  blacklist <- unique(c(
    mt_genes, rps_genes, rpl_genes, 
    s_genes, g2m_genes, 
    hb_genes, noncoding_genes,
    sex_genes
  ))
  blacklist <- blacklist[blacklist %in% all_genes]
  
  message(paste0("Blacklist contains ", length(blacklist), " genes"))
  message(paste0("  MT: ", length(mt_genes), 
                 " | RPS/RPL: ", length(c(rps_genes, rpl_genes)),
                 " | CC: ", length(intersect(c(s_genes, g2m_genes), all_genes)),
                 " | HB: ", length(hb_genes),
                 " | Noncoding: ", length(noncoding_genes),
                 " | Sex: ", length(intersect(sex_genes, all_genes))))
  return(blacklist)
}


# ============================================
# CELL FILTERING FUNCTION (unchanged)
# ============================================

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
    
    cells_in_clusters <- colnames(seurat_obj)[seurat_obj@meta.data[[cluster_col]] %in% clusters_to_remove]
    cells_to_exclude <- c(cells_to_exclude, cells_in_clusters)
  }
  
  if (!is.null(cells_to_remove) && length(cells_to_remove) > 0) {
    cells_found <- cells_to_remove[cells_to_remove %in% colnames(seurat_obj)]
    cells_to_exclude <- c(cells_to_exclude, cells_found)
  }
  
  cells_to_exclude <- unique(cells_to_exclude)
  
  if (length(cells_to_exclude) > 0) {
    cells_to_keep <- setdiff(colnames(seurat_obj), cells_to_exclude)
    seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
  }
  
  n_after <- ncol(seurat_obj)
  
  if (verbose) {
    message(paste0("\n--- Filtering Summary ---"))
    message(paste0("Cells removed: ", n_before - n_after))
    message(paste0("Cells remaining: ", n_after))
    message("=========================================\n")
  }
  
  return(seurat_obj)
}


# ============================================
# SAMPLE BALANCE REPORT
# ============================================

show_sample_balance <- function(
    seurat_obj,
    condition_col = "treatment",
    celltype_col = "mapping_cell_type",
    sample_col = NULL,
    subcluster_col = NULL
) {
  
  message("\n========== SAMPLE BALANCE REPORT ==========\n")
  
  # Cells per condition
  message("--- Cells per Condition ---")
  cond_counts <- table(seurat_obj@meta.data[[condition_col]])
  for (cond in names(cond_counts)) {
    message(paste0("  ", cond, ": ", cond_counts[cond], " cells"))
  }
  
  # Cells per sample (critical for pseudobulk)
  if (!is.null(sample_col) && sample_col %in% colnames(seurat_obj@meta.data)) {
    message("\n--- Cells per Sample ---")
    sample_cond_table <- table(
      seurat_obj@meta.data[[sample_col]],
      seurat_obj@meta.data[[condition_col]]
    )
    print(sample_cond_table)
    
    message(paste0("\nTotal samples: ", nrow(sample_cond_table)))
    message(paste0("Samples per condition:"))
    for (cond in colnames(sample_cond_table)) {
      n_samples <- sum(sample_cond_table[, cond] > 0)
      message(paste0("  ", cond, ": ", n_samples, " samples"))
    }
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
  
  # Per-sample per-celltype breakdown (important for pseudobulk feasibility)
  if (!is.null(sample_col) && sample_col %in% colnames(seurat_obj@meta.data)) {
    message("\n--- Cells per Sample per Cell Type (pseudobulk feasibility) ---")
    sample_ct_table <- table(
      seurat_obj@meta.data[[sample_col]],
      seurat_obj@meta.data[[celltype_col]]
    )
    
    message("Minimum cells per sample-celltype combination:")
    for (ct in colnames(sample_ct_table)) {
      min_cells <- min(sample_ct_table[, ct])
      n_zero <- sum(sample_ct_table[, ct] == 0)
      flag <- ""
      if (min_cells < 10) flag <- " [LOW - consider excluding]"
      if (n_zero > 0) flag <- paste0(flag, " [", n_zero, " samples with 0 cells]")
      message(paste0("  ", ct, ": min=", min_cells, flag))
    }
  }
  
  message("\n============================================\n")
}


# ============================================
# PSEUDOBULK DE ANALYSIS (PRIMARY METHOD)
# ============================================
#' Pseudobulk DE using AggregateExpression + DESeq2
#' Requires: sample_col with biological replicates (>= 2 per condition)
#' 
#' This aggregates raw counts per sample per cell type, then runs DESeq2
#' treating each sample as an independent observation (not each cell).

run_pseudobulk_de <- function(
    seurat_obj,
    condition_col = "treatment",
    celltype_col = "mapping_cell_type",
    sample_col = "orig.ident",
    ident.1 = "sCD83",
    ident.2 = "PBS",
    min_cells_per_sample = 10,
    min_samples_per_condition = 2,
    blacklist_genes = NULL
) {
  
  if (is.null(blacklist_genes)) blacklist_genes <- character(0)
  
  cell_types <- sort(unique(seurat_obj@meta.data[[celltype_col]]))
  cell_types <- cell_types[!is.na(cell_types)]
  
  message(paste0("Found ", length(cell_types), " cell types"))
  message(paste0("Using sample column: ", sample_col))
  message(paste0("Samples: ", paste(unique(seurat_obj@meta.data[[sample_col]]), collapse = ", ")))
  
  de_results_by_celltype <- list()
  
  for (ct in cell_types) {
    message(paste0("\n--- Processing: ", ct, " ---"))
    
    # Subset to this cell type
    seurat_ct <- subset(seurat_obj, subset = !!sym(celltype_col) == ct)
    
    if (ncol(seurat_ct) < 10) {
      message("  Skipping - too few cells total")
      next
    }
    
    # Check sample-level adequacy
    sample_condition_table <- table(
      seurat_ct@meta.data[[sample_col]],
      seurat_ct@meta.data[[condition_col]]
    )
    
    # Filter to samples with enough cells
    adequate_samples <- list()
    for (cond in c(ident.1, ident.2)) {
      if (!cond %in% colnames(sample_condition_table)) {
        message(paste0("  Skipping - condition '", cond, "' missing"))
        next
      }
      cond_samples <- rownames(sample_condition_table)[sample_condition_table[, cond] >= min_cells_per_sample]
      adequate_samples[[cond]] <- cond_samples
      message(paste0("  ", cond, ": ", length(cond_samples), " samples with >= ", 
                     min_cells_per_sample, " cells"))
    }
    
    if (length(adequate_samples) < 2) next
    if (length(adequate_samples[[ident.1]]) < min_samples_per_condition ||
        length(adequate_samples[[ident.2]]) < min_samples_per_condition) {
      message(paste0("  Skipping - need >= ", min_samples_per_condition, 
                     " samples per condition for pseudobulk"))
      next
    }
    
    # Keep only adequate samples
    keep_samples <- c(adequate_samples[[ident.1]], adequate_samples[[ident.2]])
    cells_keep <- colnames(seurat_ct)[seurat_ct@meta.data[[sample_col]] %in% keep_samples]
    seurat_ct <- subset(seurat_ct, cells = cells_keep)
    
    # ---- Pseudobulk aggregation ----
    # AggregateExpression sums raw counts per sample within this cell type
    # group.by must include sample_col and condition_col
    tryCatch({
      
      pseudo <- AggregateExpression(
        seurat_ct,
        assays = "RNA",
        return.seurat = TRUE,
        group.by = c(sample_col, condition_col)
      )
      
      # Fix known Seurat v5 issue: AggregateExpression can return non-integer 
      # counts when SCTransform models are present (GitHub issue #9771).
      # DESeq2 requires integers, so round the counts layer.
      counts_mat <- GetAssayData(pseudo, layer = "counts")
      if (!all(counts_mat@x == floor(counts_mat@x))) {
        message("  Rounding non-integer pseudobulk counts (known Seurat v5 issue)")
        counts_mat@x <- round(counts_mat@x)
        pseudo <- SetAssayData(pseudo, layer = "counts", new.data = counts_mat)
      }
      
      # Remove blacklisted genes from features to test
      features_to_test <- setdiff(rownames(pseudo), blacklist_genes)
      
      # Set idents to condition for FindMarkers
      Idents(pseudo) <- condition_col
      
      de <- FindMarkers(
        pseudo,
        ident.1 = ident.1,
        ident.2 = ident.2,
        test.use = "DESeq2",
        features = features_to_test,
        min.pct = 0,           # Don't pre-filter for pseudobulk - let DESeq2 handle it
        logfc.threshold = 0    # Get all genes, filter post-hoc
      )
      
      de$gene <- rownames(de)
      de$cell_type <- ct
      de$method <- "pseudobulk_DESeq2"
      de$n_samples_ident1 <- length(adequate_samples[[ident.1]])
      de$n_samples_ident2 <- length(adequate_samples[[ident.2]])
      
      # DESeq2 already uses BH correction (not Bonferroni)
      # but FindMarkers may override with Bonferroni. Recalculate BH if needed.
      de$p_val_BH <- p.adjust(de$p_val, method = "BH")
      
      de_results_by_celltype[[ct]] <- de
      
      n_sig <- sum(de$p_val_adj < 0.05 & abs(de$avg_log2FC) > 0.5, na.rm = TRUE)
      message(paste0("  DESeq2 complete: ", nrow(de), " genes tested, ", n_sig, " significant"))
      
    }, error = function(e) {
      message(paste0("  Pseudobulk DESeq2 failed: ", e$message))
      message("  This can happen with too few samples. Consider MAST fallback.")
    })
  }
  
  # ---- Global pseudobulk (all cell types combined) ----
  message("\n--- Running global pseudobulk DE (all cells) ---")
  
  tryCatch({
    pseudo_global <- AggregateExpression(
      seurat_obj,
      assays = "RNA",
      return.seurat = TRUE,
      group.by = c(sample_col, condition_col)
    )
    
    counts_mat <- GetAssayData(pseudo_global, layer = "counts")
    if (!all(counts_mat@x == floor(counts_mat@x))) {
      counts_mat@x <- round(counts_mat@x)
      pseudo_global <- SetAssayData(pseudo_global, layer = "counts", new.data = counts_mat)
    }
    
    features_to_test <- setdiff(rownames(pseudo_global), blacklist_genes)
    Idents(pseudo_global) <- condition_col
    
    de_global <- FindMarkers(
      pseudo_global,
      ident.1 = ident.1,
      ident.2 = ident.2,
      test.use = "DESeq2",
      features = features_to_test,
      min.pct = 0,
      logfc.threshold = 0
    )
    
    de_global$gene <- rownames(de_global)
    de_global$cell_type <- "All"
    de_global$method <- "pseudobulk_DESeq2"
    de_global$p_val_BH <- p.adjust(de_global$p_val, method = "BH")
    
  }, error = function(e) {
    message(paste0("  Global pseudobulk failed: ", e$message))
    de_global <- data.frame()
  })
  
  return(list(
    by_celltype = de_results_by_celltype,
    global = de_global,
    combined = dplyr::bind_rows(de_results_by_celltype)
  ))
}


# ============================================
# MAST DE ANALYSIS (FALLBACK / COMPLEMENT)
# ============================================
#' MAST with latent variables as fallback when pseudobulk isn't possible
#' (e.g., single replicate per condition, or as a complementary method)
#' 
#' MAST models the cellular detection rate (CDR) as a covariate, which
#' partially accounts for technical variation. Including sample_col as
#' a latent variable helps control for within-sample correlation.

run_mast_de <- function(
    seurat_obj,
    condition_col = "treatment",
    celltype_col = "mapping_cell_type",
    sample_col = NULL,
    ident.1 = "sCD83",
    ident.2 = "PBS",
    min.pct = 0.1,
    blacklist_genes = NULL,
    latent_vars = NULL
) {
  
  if (!requireNamespace("MAST", quietly = TRUE)) {
    stop("MAST package required. Install with: BiocManager::install('MAST')")
  }
  
  if (is.null(blacklist_genes)) blacklist_genes <- character(0)
  
  # Build latent variables list
  # CDR (cellular detection rate) is automatically included by MAST via Seurat
  # Additional latent vars control for technical/biological confounders
  if (is.null(latent_vars)) {
    latent_vars <- c()
    # Include nCount_RNA if available (library size)
    if ("nCount_RNA" %in% colnames(seurat_obj@meta.data)) {
      latent_vars <- c(latent_vars, "nCount_RNA")
    }
    # Include percent.mt if available
    if ("percent.mt" %in% colnames(seurat_obj@meta.data)) {
      latent_vars <- c(latent_vars, "percent.mt")
    }
    # Include sample as latent variable to partially account for 
    # within-sample correlation (not as powerful as pseudobulk but better than nothing)
    if (!is.null(sample_col) && sample_col %in% colnames(seurat_obj@meta.data)) {
      latent_vars <- c(latent_vars, sample_col)
    }
  }
  
  message(paste0("MAST latent variables: ", paste(latent_vars, collapse = ", ")))
  
  cell_types <- sort(unique(seurat_obj@meta.data[[celltype_col]]))
  cell_types <- cell_types[!is.na(cell_types)]
  
  de_results_by_celltype <- list()
  
  for (ct in cell_types) {
    message(paste0("\nProcessing (MAST): ", ct))
    
    seurat_ct <- subset(seurat_obj, subset = !!sym(celltype_col) == ct)
    
    if (ncol(seurat_ct) < 10) {
      message("  Skipping - too few cells")
      next
    }
    
    conditions_present <- unique(seurat_ct@meta.data[[condition_col]])
    if (!(ident.1 %in% conditions_present) || !(ident.2 %in% conditions_present)) {
      message("  Skipping - missing condition")
      next
    }
    
    n1 <- sum(seurat_ct@meta.data[[condition_col]] == ident.1)
    n2 <- sum(seurat_ct@meta.data[[condition_col]] == ident.2)
    
    if (n1 < 3 || n2 < 3) {
      message("  Skipping - insufficient cells per condition")
      next
    }
    
    message(paste0("  ", ident.1, ": ", n1, " cells, ", ident.2, ": ", n2, " cells"))
    
    Idents(seurat_ct) <- condition_col
    features_to_test <- setdiff(rownames(seurat_ct), blacklist_genes)
    
    # Validate latent vars exist in this subset
    valid_latent <- latent_vars[latent_vars %in% colnames(seurat_ct@meta.data)]
    
    # Remove latent vars with zero variance (causes MAST to fail)
    for (lv in valid_latent) {
      if (is.numeric(seurat_ct@meta.data[[lv]])) {
        if (var(seurat_ct@meta.data[[lv]], na.rm = TRUE) == 0) {
          valid_latent <- setdiff(valid_latent, lv)
          message(paste0("  Removing latent var '", lv, "' - zero variance"))
        }
      } else {
        if (length(unique(seurat_ct@meta.data[[lv]])) <= 1) {
          valid_latent <- setdiff(valid_latent, lv)
          message(paste0("  Removing latent var '", lv, "' - single level"))
        }
      }
    }
    
    tryCatch({
      de <- FindMarkers(
        seurat_ct,
        ident.1 = ident.1,
        ident.2 = ident.2,
        test.use = "MAST",
        features = features_to_test,
        min.pct = min.pct,
        logfc.threshold = 0,
        latent.vars = if (length(valid_latent) > 0) valid_latent else NULL
      )
      
      de$gene <- rownames(de)
      de$cell_type <- ct
      de$method <- "MAST"
      de$latent_vars_used <- paste(valid_latent, collapse = ",")
      
      # Recalculate BH-adjusted p-values (FindMarkers uses Bonferroni by default)
      de$p_val_BH <- p.adjust(de$p_val, method = "BH")
      
      de_results_by_celltype[[ct]] <- de
      
      n_sig <- sum(de$p_val_BH < 0.05 & abs(de$avg_log2FC) > 0.5, na.rm = TRUE)
      message(paste0("  MAST complete: ", nrow(de), " genes tested, ", n_sig, " significant (BH)"))
      
    }, error = function(e) {
      message(paste0("  MAST error: ", e$message))
    })
  }
  
  # Global MAST
  message("\n--- Running global MAST DE (all cells) ---")
  Idents(seurat_obj) <- condition_col
  features_to_test <- setdiff(rownames(seurat_obj), blacklist_genes)
  
  valid_latent_global <- latent_vars[latent_vars %in% colnames(seurat_obj@meta.data)]
  
  de_global <- tryCatch({
    de <- FindMarkers(
      seurat_obj,
      ident.1 = ident.1,
      ident.2 = ident.2,
      test.use = "MAST",
      features = features_to_test,
      min.pct = min.pct,
      logfc.threshold = 0,
      latent.vars = if (length(valid_latent_global) > 0) valid_latent_global else NULL
    )
    de$gene <- rownames(de)
    de$cell_type <- "All"
    de$method <- "MAST"
    de$p_val_BH <- p.adjust(de$p_val, method = "BH")
    de
  }, error = function(e) {
    message(paste0("  Global MAST error: ", e$message))
    data.frame()
  })
  
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
    pval_col = "p_val_BH",      # Use BH-adjusted by default instead of Bonferroni
    title = "Volcano Plot",
    top_n_label = 20,
    cap_logfc = TRUE,
    cap_pval = TRUE,
    logfc_cap_quantile = 0.99,
    pval_cap_quantile = 0.99
) {
  
  # Use specified p-value column, fall back to p_val_adj if BH not available
  if (!pval_col %in% colnames(de_df)) {
    pval_col <- "p_val_adj"
  }
  
  de_df <- de_df %>%
    dplyr::mutate(
      pval_use = .data[[pval_col]],
      neg_log10_pval = -log10(pval_use),
      significance = dplyr::case_when(
        pval_use < pval_thresh & avg_log2FC > logfc_thresh ~ "Up",
        pval_use < pval_thresh & avg_log2FC < -logfc_thresh ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # Handle infinite p-values
  if (any(is.infinite(de_df$neg_log10_pval))) {
    max_finite <- max(de_df$neg_log10_pval[is.finite(de_df$neg_log10_pval)])
    de_df$neg_log10_pval[is.infinite(de_df$neg_log10_pval)] <- max_finite * 1.1
  }
  
  # Cap extreme log2FC
  if (cap_logfc) {
    logfc_upper <- quantile(abs(de_df$avg_log2FC), logfc_cap_quantile, na.rm = TRUE)
    logfc_cap <- max(logfc_upper, logfc_thresh * 2)
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
  
  # Cap extreme p-values
  if (cap_pval) {
    pval_upper <- quantile(de_df$neg_log10_pval[is.finite(de_df$neg_log10_pval)], 
                           pval_cap_quantile, na.rm = TRUE)
    pval_cap <- max(pval_upper, -log10(pval_thresh) * 2)
    de_df <- de_df %>%
      dplyr::mutate(
        is_pval_capped = neg_log10_pval > pval_cap,
        neg_log10_pval_plot = dplyr::if_else(neg_log10_pval > pval_cap, pval_cap, neg_log10_pval)
      )
  } else {
    de_df$neg_log10_pval_plot <- de_df$neg_log10_pval
    de_df$is_pval_capped <- FALSE
  }
  
  de_df$is_capped <- de_df$is_logfc_capped | de_df$is_pval_capped
  
  # Top genes to label
  top_up <- de_df %>%
    dplyr::filter(significance == "Up") %>%
    dplyr::arrange(pval_use, desc(avg_log2FC)) %>%
    dplyr::slice_head(n = top_n_label)
  
  top_down <- de_df %>%
    dplyr::filter(significance == "Down") %>%
    dplyr::arrange(pval_use, avg_log2FC) %>%
    dplyr::slice_head(n = top_n_label)
  
  top_genes <- dplyr::bind_rows(top_up, top_down)
  
  n_up <- sum(de_df$significance == "Up")
  n_down <- sum(de_df$significance == "Down")
  n_capped <- sum(de_df$is_capped)
  
  # Determine method label for subtitle
  method_label <- if ("method" %in% colnames(de_df)) unique(de_df$method)[1] else "unknown"
  
  p <- ggplot(de_df, aes(x = avg_log2FC_plot, y = neg_log10_pval_plot)) +
    geom_point(
      data = dplyr::filter(de_df, !is_capped),
      aes(color = significance),
      alpha = 0.5, size = 1.2
    ) +
    geom_point(
      data = dplyr::filter(de_df, is_capped),
      aes(color = significance),
      alpha = 0.7, size = 1.5, shape = 17
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
      subtitle = paste0("Method: ", method_label, " | p-adj: ", pval_col),
      x = ifelse(cap_logfc && any(de_df$is_logfc_capped), 
                 "log2 Fold Change (capped)", "log2 Fold Change"),
      y = ifelse(cap_pval && any(de_df$is_pval_capped),
                 "-log10(adj. p-value) (capped)", "-log10(adj. p-value)"),
      color = "",
      caption = ifelse(n_capped > 0, paste0("▲ = ", n_capped, " capped outliers"), "")
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey40"),
      legend.position = "top",
      plot.caption = element_text(hjust = 0, size = 8, color = "grey40")
    )
  
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
# METHOD CONCORDANCE PLOT
# ============================================
#' Compare pseudobulk vs MAST results when both are run

plot_method_concordance <- function(
    de_pseudobulk,
    de_mast,
    cell_type = "All",
    logfc_thresh = 0.5,
    pval_thresh = 0.05
) {
  
  # Get appropriate dataframes
  if (cell_type == "All") {
    pb <- de_pseudobulk$global
    mast <- de_mast$global
  } else {
    pb <- de_pseudobulk$by_celltype[[cell_type]]
    mast <- de_mast$by_celltype[[cell_type]]
  }
  
  if (is.null(pb) || is.null(mast) || nrow(pb) == 0 || nrow(mast) == 0) {
    message("Cannot compare - missing results for one or both methods")
    return(NULL)
  }
  
  # Merge on gene
  merged <- merge(
    pb %>% dplyr::select(gene, avg_log2FC, p_val_adj) %>% dplyr::rename(lfc_pb = avg_log2FC, padj_pb = p_val_adj),
    mast %>% dplyr::select(gene, avg_log2FC, p_val_BH) %>% dplyr::rename(lfc_mast = avg_log2FC, padj_mast = p_val_BH),
    by = "gene"
  )
  
  merged <- merged %>%
    dplyr::mutate(
      sig_pb = padj_pb < pval_thresh & abs(lfc_pb) > logfc_thresh,
      sig_mast = padj_mast < pval_thresh & abs(lfc_mast) > logfc_thresh,
      category = dplyr::case_when(
        sig_pb & sig_mast ~ "Both",
        sig_pb & !sig_mast ~ "Pseudobulk only",
        !sig_pb & sig_mast ~ "MAST only",
        TRUE ~ "Neither"
      )
    )
  
  # Log fold change correlation
  cor_val <- cor(merged$lfc_pb, merged$lfc_mast, use = "complete.obs", method = "spearman")
  
  p <- ggplot(merged, aes(x = lfc_pb, y = lfc_mast, color = category)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c(
      "Both" = "#D51F26",
      "Pseudobulk only" = "#272E6A",
      "MAST only" = "#F47D2B",
      "Neither" = "grey80"
    )) +
    labs(
      title = paste0("Method Concordance: ", cell_type),
      subtitle = paste0("Spearman rho = ", round(cor_val, 3)),
      x = "log2FC (Pseudobulk DESeq2)",
      y = "log2FC (MAST)",
      color = "Significant in"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      legend.position = "top"
    )
  
  # Print concordance summary
  message(paste0("\n--- Concordance: ", cell_type, " ---"))
  message(paste0("  Both significant: ", sum(merged$category == "Both")))
  message(paste0("  Pseudobulk only: ", sum(merged$category == "Pseudobulk only")))
  message(paste0("  MAST only: ", sum(merged$category == "MAST only")))
  message(paste0("  Spearman rho (logFC): ", round(cor_val, 3)))
  
  return(p)
}


# ============================================
# MAIN WRAPPER FUNCTION
# ============================================

run_full_de_pipeline <- function(
    seurat_obj,
    condition_col = "treatment",
    celltype_col = "mapping_cell_type",
    sample_col = "orig.ident",          # CRITICAL: must be biological replicate ID
    ident.1 = "sCD83",
    ident.2 = "PBS",
    blacklist_genes = NULL,
    clusters_to_remove = NULL,
    cluster_col = "fine_clust",
    cells_to_remove = NULL,
    top_n_genes = 50,
    logfc_thresh = 0.5,
    pval_thresh = 0.05,
    pval_col = "p_val_BH",             # BH correction by default
    run_pseudobulk = FALSE,
    run_mast = TRUE,                    # Run both for concordance
    min_cells_per_sample = 10,
    min_samples_per_condition = 2,
    output_dir = "de_results"
) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # ============================================
  # Step 0: Validate sample column
  # ============================================
  
  if (run_pseudobulk) {
    if (!sample_col %in% colnames(seurat_obj@meta.data)) {
      warning(paste0("sample_col '", sample_col, "' not found in metadata.\n",
                     "Available columns: ", paste(colnames(seurat_obj@meta.data), collapse = ", "),
                     "\nPseudobulk requires a column identifying biological replicates.\n",
                     "Falling back to MAST only."))
      run_pseudobulk <- FALSE
    } else {
      # Check if there are enough samples
      sample_per_cond <- table(
        unique(seurat_obj@meta.data[, c(sample_col, condition_col)])[[condition_col]]
      )
      message("\n--- Sample counts per condition ---")
      print(sample_per_cond)
      
      for (cond in c(ident.1, ident.2)) {
        if (!cond %in% names(sample_per_cond) || sample_per_cond[cond] < min_samples_per_condition) {
          warning(paste0("Condition '", cond, "' has fewer than ", min_samples_per_condition,
                         " samples. Pseudobulk may fail for some cell types."))
        }
      }
    }
  }
  
  # ============================================
  # Step 1: Filter cells/clusters
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
  
  show_sample_balance(
    seurat_obj,
    condition_col = condition_col,
    celltype_col = celltype_col,
    sample_col = if (run_pseudobulk) sample_col else NULL,
    subcluster_col = cluster_col
  )
  
  # Get blacklist
  if (is.null(blacklist_genes)) {
    blacklist_genes <- get_blacklist_genes(seurat_obj)
  }
  
  # ============================================
  # Step 2: Run DE analysis
  # ============================================
  
  de_pseudobulk <- NULL
  de_mast <- NULL
  
  # Primary: Pseudobulk DESeq2
  if (run_pseudobulk) {
    message("\n========== RUNNING PSEUDOBULK DESeq2 ==========\n")
    de_pseudobulk <- run_pseudobulk_de(
      seurat_obj = seurat_obj,
      condition_col = condition_col,
      celltype_col = celltype_col,
      sample_col = sample_col,
      ident.1 = ident.1,
      ident.2 = ident.2,
      min_cells_per_sample = min_cells_per_sample,
      min_samples_per_condition = min_samples_per_condition,
      blacklist_genes = blacklist_genes
    )
  }
  
  # Secondary/fallback: MAST with latent variables
  if (run_mast) {
    message("\n========== RUNNING MAST ==========\n")
    de_mast <- run_mast_de(
      seurat_obj = seurat_obj,
      condition_col = condition_col,
      celltype_col = celltype_col,
      sample_col = sample_col,
      ident.1 = ident.1,
      ident.2 = ident.2,
      blacklist_genes = blacklist_genes
    )
  }
  
  # Pick primary results (prefer pseudobulk)
  de_results <- if (!is.null(de_pseudobulk)) de_pseudobulk else de_mast
  
  # ============================================
  # Step 3: Save DE tables
  # ============================================
  
  if (!is.null(de_pseudobulk)) {
    write.csv(de_pseudobulk$global, file.path(output_dir, "DE_pseudobulk_global.csv"), row.names = FALSE)
    write.csv(de_pseudobulk$combined, file.path(output_dir, "DE_pseudobulk_all_celltypes.csv"), row.names = FALSE)
    for (ct in names(de_pseudobulk$by_celltype)) {
      ct_clean <- gsub("[^A-Za-z0-9]", "_", ct)
      write.csv(de_pseudobulk$by_celltype[[ct]], 
                file.path(output_dir, paste0("DE_pseudobulk_", ct_clean, ".csv")), row.names = FALSE)
    }
  }
  
  if (!is.null(de_mast)) {
    write.csv(de_mast$global, file.path(output_dir, "DE_MAST_global.csv"), row.names = FALSE)
    write.csv(de_mast$combined, file.path(output_dir, "DE_MAST_all_celltypes.csv"), row.names = FALSE)
    for (ct in names(de_mast$by_celltype)) {
      ct_clean <- gsub("[^A-Za-z0-9]", "_", ct)
      write.csv(de_mast$by_celltype[[ct]], 
                file.path(output_dir, paste0("DE_MAST_", ct_clean, ".csv")), row.names = FALSE)
    }
  }
  
  # ============================================
  # Step 4: Generate plots (using primary results)
  # ============================================
  
  all_plots <- list()
  
  # --- Global plots ---
  message("\n========== GENERATING PLOTS ==========\n")
  
  sig_genes_global <- de_results$global %>%
    dplyr::filter(
      (if (pval_col %in% colnames(.)) .data[[pval_col]] else p_val_adj) < pval_thresh, 
      abs(avg_log2FC) > logfc_thresh
    ) %>%
    dplyr::arrange(desc(abs(avg_log2FC))) %>%
    dplyr::slice_head(n = top_n_genes) %>%
    dplyr::pull(gene)
  
  volcano_global <- create_volcano_plot(
    de_results$global,
    logfc_thresh = logfc_thresh,
    pval_thresh = pval_thresh,
    pval_col = pval_col,
    title = paste0("Global: ", ident.1, " vs ", ident.2)
  )
  
  density_global <- create_density_plot(
    seurat_obj, genes = sig_genes_global,
    condition_col = condition_col,
    ident.1 = ident.1, ident.2 = ident.2,
    title = "Global DE Gene Expression"
  )
  
  heatmap_global <- create_heatmap(
    seurat_obj, genes = sig_genes_global,
    condition_col = condition_col,
    ident.1 = ident.1, ident.2 = ident.2,
    title = "Global"
  )
  
  all_plots[["Global"]] <- list(
    volcano = volcano_global,
    density = density_global,
    heatmap = heatmap_global,
    sig_genes = sig_genes_global
  )
  
  # Save global
  ggsave(file.path(output_dir, "Global_volcano.pdf"), volcano_global, width = 7, height = 6)
  ggsave(file.path(output_dir, "Global_density.pdf"), density_global, width = 6, height = 5)
  combined_global <- density_global + volcano_global + plot_layout(ncol = 2)
  ggsave(file.path(output_dir, "Global_density_volcano.pdf"), combined_global, width = 12, height = 5)
  ggsave(file.path(output_dir, "Global_density_volcano.png"), combined_global, width = 12, height = 5, dpi = 150)
  
  if (!is.null(heatmap_global)) {
    pdf(file.path(output_dir, "Global_heatmap.pdf"), width = 10, height = 8)
    draw(heatmap_global)
    dev.off()
  }
  
  # --- Per cell type plots ---
  for (ct in names(de_results$by_celltype)) {
    message(paste0("\nPlotting: ", ct))
    ct_clean <- gsub("[^A-Za-z0-9]", "_", ct)
    de_ct <- de_results$by_celltype[[ct]]
    seurat_ct <- subset(seurat_obj, subset = !!sym(celltype_col) == ct)
    
    sig_genes_ct <- de_ct %>%
      dplyr::filter(
        (if (pval_col %in% colnames(.)) .data[[pval_col]] else p_val_adj) < pval_thresh,
        abs(avg_log2FC) > logfc_thresh
      ) %>%
      dplyr::arrange(desc(abs(avg_log2FC))) %>%
      dplyr::slice_head(n = top_n_genes) %>%
      dplyr::pull(gene)
    
    volcano_ct <- create_volcano_plot(
      de_ct, logfc_thresh = logfc_thresh, pval_thresh = pval_thresh,
      pval_col = pval_col, top_n_label = 20,
      title = paste0(ct, ": ", ident.1, " vs ", ident.2)
    )
    
    density_ct <- create_density_plot(
      seurat_ct, genes = sig_genes_ct,
      condition_col = condition_col,
      ident.1 = ident.1, ident.2 = ident.2,
      title = paste0(ct, " DE Genes")
    )
    
    heatmap_ct <- NULL
    if (length(sig_genes_ct) > 0) {
      heatmap_ct <- create_heatmap(
        seurat_ct, genes = sig_genes_ct,
        condition_col = condition_col,
        ident.1 = ident.1, ident.2 = ident.2,
        title = ct
      )
    }
    
    all_plots[[ct]] <- list(
      volcano = volcano_ct, density = density_ct,
      heatmap = heatmap_ct, sig_genes = sig_genes_ct
    )
    
    ggsave(file.path(output_dir, paste0(ct_clean, "_volcano.pdf")), volcano_ct, width = 7, height = 6)
    ggsave(file.path(output_dir, paste0(ct_clean, "_density.pdf")), density_ct, width = 6, height = 5)
    combined_ct <- density_ct + volcano_ct + plot_layout(ncol = 2)
    ggsave(file.path(output_dir, paste0(ct_clean, "_density_volcano.pdf")), combined_ct, width = 12, height = 5)
    ggsave(file.path(output_dir, paste0(ct_clean, "_density_volcano.png")), combined_ct, width = 12, height = 5, dpi = 150)
    
    if (!is.null(heatmap_ct)) {
      pdf(file.path(output_dir, paste0(ct_clean, "_heatmap.pdf")), width = 10, height = 8)
      draw(heatmap_ct)
      dev.off()
    }
  }
  
  # --- Method concordance plots (if both methods run) ---
  concordance_plots <- list()
  if (!is.null(de_pseudobulk) && !is.null(de_mast)) {
    message("\n========== METHOD CONCORDANCE ==========\n")
    
    # Global concordance
    p_conc_global <- plot_method_concordance(
      de_pseudobulk, de_mast, cell_type = "All",
      logfc_thresh = logfc_thresh, pval_thresh = pval_thresh
    )
    if (!is.null(p_conc_global)) {
      concordance_plots[["Global"]] <- p_conc_global
      ggsave(file.path(output_dir, "concordance_Global.pdf"), p_conc_global, width = 7, height = 6)
    }
    
    # Per cell type concordance
    shared_cts <- intersect(names(de_pseudobulk$by_celltype), names(de_mast$by_celltype))
    for (ct in shared_cts) {
      p_conc <- plot_method_concordance(
        de_pseudobulk, de_mast, cell_type = ct,
        logfc_thresh = logfc_thresh, pval_thresh = pval_thresh
      )
      if (!is.null(p_conc)) {
        concordance_plots[[ct]] <- p_conc
        ct_clean <- gsub("[^A-Za-z0-9]", "_", ct)
        ggsave(file.path(output_dir, paste0("concordance_", ct_clean, ".pdf")), p_conc, width = 7, height = 6)
      }
    }
  }
  
  # ============================================
  # Step 5: Summary
  # ============================================
  
  summary_list <- list()
  
  for (method_name in c("pseudobulk", "MAST")) {
    de_res <- if (method_name == "pseudobulk") de_pseudobulk else de_mast
    if (is.null(de_res)) next
    
    pval_col_use <- if (pval_col %in% colnames(de_res$global)) pval_col else "p_val_adj"
    
    summary_df <- data.frame(
      cell_type = c("Global", names(de_res$by_celltype)),
      method = method_name
    )
    
    summary_df$n_sig_up <- sapply(summary_df$cell_type, function(ct) {
      df <- if (ct == "Global") de_res$global else de_res$by_celltype[[ct]]
      if (is.null(df) || nrow(df) == 0) return(0)
      pv <- if (pval_col_use %in% colnames(df)) df[[pval_col_use]] else df$p_val_adj
      sum(pv < pval_thresh & df$avg_log2FC > logfc_thresh, na.rm = TRUE)
    })
    
    summary_df$n_sig_down <- sapply(summary_df$cell_type, function(ct) {
      df <- if (ct == "Global") de_res$global else de_res$by_celltype[[ct]]
      if (is.null(df) || nrow(df) == 0) return(0)
      pv <- if (pval_col_use %in% colnames(df)) df[[pval_col_use]] else df$p_val_adj
      sum(pv < pval_thresh & df$avg_log2FC < -logfc_thresh, na.rm = TRUE)
    })
    
    summary_df$n_sig_total <- summary_df$n_sig_up + summary_df$n_sig_down
    summary_list[[method_name]] <- summary_df
  }
  
  summary_all <- dplyr::bind_rows(summary_list)
  write.csv(summary_all, file.path(output_dir, "DE_summary_all_methods.csv"), row.names = FALSE)
  
  message("\n========== DONE ==========")
  message(paste0("Results saved to: ", output_dir))
  print(summary_all)
  
  return(list(
    de_pseudobulk = de_pseudobulk,
    de_mast = de_mast,
    de_primary = de_results,
    plots = all_plots,
    concordance_plots = concordance_plots,
    summary = summary_all,
    seurat_filtered = seurat_obj
  ))
}


###########################################################################
# USAGE EXAMPLE
###########################################################################

obj <- readRDS("./fine_annotated_obj.rds")
obj$mapping_cell_type <- sub("\\.\\d+$", "", obj$broad_cluster)

# ---------------------------------------------------------------
# CRITICAL: Verify you have a sample/replicate ID column.
# This must identify biological replicates (e.g., individual mice,
# individual patients, individual wells from different donors).
#
# Check what's available:
print(colnames(obj@meta.data))
#
# Common column names: "orig.ident", "sample_id", "donor_id", "mouse_id"
# If orig.ident maps 1:1 to treatment (i.e., 1 sample per condition),
# pseudobulk cannot be used — set run_pseudobulk = FALSE.
#
# Verify with:
# table(obj$orig.ident, obj$treatment)
# ---------------------------------------------------------------
# Check balance BEFORE running
show_sample_balance(
  obj,
  condition_col = "treatment",
  celltype_col = "mapping_cell_type",
  sample_col = "orig.ident",        
  subcluster_col = "fine_clust"
)

# Run the improved pipeline
results <- run_full_de_pipeline(
  seurat_obj = obj,
  condition_col = "treatment",
  celltype_col = "mapping_cell_type",
  sample_col = "orig.ident",            
  ident.1 = "sCD83",
  ident.2 = "PBS",
  blacklist_genes = NULL,               
  clusters_to_remove = NULL,
  cluster_col = "fine_clust",
  cells_to_remove = NULL,
  top_n_genes = 2000,
  logfc_thresh = 0.5,
  pval_thresh = 0.05,
  pval_col = "p_val_BH",               
  run_pseudobulk = FALSE,                
  run_mast = TRUE,                      
  min_cells_per_sample = 10,
  min_samples_per_condition = 2,
  output_dir = "de_results"
)

# View results
print(results$plots$Global$volcano)
draw(results$plots$Global$heatmap)

# View concordance between methods
print(results$concordance_plots$Global)

# Access primary DE results
head(results$de_primary$global)

# View summary across methods
print(results$summary)