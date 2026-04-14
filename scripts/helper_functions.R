library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(grid)
library(ggrepel)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(DESeq2)
library(tidyr)
library(speckle)

###################################################
# SIMPLE FUNCTION TO PLOT QC METRICS FOR ALL DATA
###################################################

show_qc_metrics <- function(input_folder) {
  
  post_aRNA_removal_dirs <- list.dirs(
    file.path(input_folder, "preprocessing"),
    recursive = FALSE
  )
  
  h5_files <- list.files(
    post_aRNA_removal_dirs,
    pattern = "postcellbender_filtered_seurat.h5",
    full.names = TRUE,
    recursive = TRUE
  )
  
  object.list <- list()
  
  for (i in seq_along(h5_files)) {
    
    object <- h5_files[[i]]  
    stage <- dataset_names[[i]]
    
    # Load CellBender-corrected data
    data.arna_corrected <- Read10X_h5(object, use.names = TRUE)  
    obj <- CreateSeuratObject(counts = data.arna_corrected, project = stage)
    obj$orig.ident <- stage
    
    # Calculate QC metrics
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    
    # QC visualization
    gg1 <- print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
    output_folder = paste0(input_folder, "/preprocessing/qc_metrics/")
    dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
    
    file_name = paste0(output_folder, "qc_VlnPlot_", stage, "_.png")
    ggsave(file_name)
    
    plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    print(plot1 + plot2)
    
    file_name = paste0(output_folder, "qc_Scatterplot_", stage, "_.png")
    ggsave(file_name, width = 12)
    
    object.list[[i]] <- obj
    
    rm(data.arna_corrected)
    rm(obj)
  }
  
  return(object.list)
  
}

###################################################
# FILTER DATA AFTER QC METRIC VISUALIZATION
###################################################

filter_by_qc <- function(input_folder, project_names, min_feature = NULL, max_feature = NULL, min_count = NULL, max_count = NULL, max_percent_mt = NULL) {
  
  post_aRNA_removal_dirs <- list.dirs(
    file.path(input_folder, "preprocessing"),
    recursive = FALSE
  )
  
  h5_files <- list.files(
    post_aRNA_removal_dirs,
    pattern = "postcellbender_filtered_seurat\\.h5$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  # DEBUG: Show file order
  cat("Files found (in order):\n")
  print(basename(h5_files))
  cat("\nProject names (in order):\n")
  print(project_names)
  cat("\n")
  
  object.list <- list()
  
  for (i in seq_along(h5_files)) {
    
    cat("Processing file:", basename(h5_files[i]), "with name:", project_names[i], "\n")
    
    counts <- Seurat::Read10X_h5(h5_files[i])
    
    obj <- Seurat::CreateSeuratObject(
      counts = counts,
      project = as.character(project_names[i])
    )
    
    obj$orig.ident <- project_names[i]
    cat("Set orig.ident to:", unique(obj$orig.ident), "\n")
    
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^mt-|^MT-")
    
    obj <- subset(
      obj, 
      subset = nFeature_RNA > min_feature & 
        nFeature_RNA < max_feature & 
        nCount_RNA > min_count &  
        nCount_RNA < max_count & 
        percent.mt < max_percent_mt
    )
    
    # Rename cells BEFORE adding to list, using project name
    obj <- RenameCells(obj, add.cell.id = project_names[i])
    
    object.list[[i]] <- obj
    
    cat("Remaining cells after QC for", project_names[i], "is", ncol(obj), "cells\n\n")
  }
  
  return(object.list)
}

##################################################
# FILTER DOUBLETS
##################################################

filter_doublets <- function(object_list) {
  for (i in seq_along(object_list)) {
    obj <- object_list[[i]]
    obj <- subset(obj, subset = scDblFinder.class == "singlet")
    object_list[[i]] <- obj
  }
  return(object_list)
}

##############################################################
# BLACKLIST GENES FOR VARIABLE FEATURES
##############################################################

get_blacklist_genes <- function(seurat_obj, remove_cc_genes = FALSE) {
  
  allGenes <- rownames(seurat_obj)
  
  # Mitochondrial
  mt.genes <- grep(pattern = "^MT-", x = allGenes, value = TRUE)
  
  # Ribosomal
  RPS.genes <- grep(pattern = "^RPS", x = allGenes, value = TRUE)
  RPL.genes <- grep(pattern = "^RPL", x = allGenes, value = TRUE)
  
  # X/Y chromosome genes - attempt DB lookup, fall back to pattern matching
  sexChr.genes <- tryCatch({
    library(GenomicFeatures)
    library(org.Hs.eg.db)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    geneGR <- GenomicFeatures::genes(txdb)
    sexGenesGR <- geneGR[seqnames(geneGR) %in% c("chrY", "chrX")]
    matchedGeneSymbols <- AnnotationDbi::select(org.Hs.eg.db,
                                                keys = sexGenesGR$gene_id,
                                                columns = c("ENTREZID", "SYMBOL"),
                                                keytype = "ENTREZID")
    message("Sex chromosome genes retrieved from TxDb")
    matchedGeneSymbols$SYMBOL
  }, error = function(e) {
    message("TxDb unavailable, using pattern-based sex chromosome filtering")
    grep(pattern = "^XIST$|^TSIX$|^RPS4Y|^DDX3Y|^USP9Y|^UTY$|^KDM5D|^EIF1AY|^ZFY$|^SRY$|^NLGN4Y",
         x = allGenes, value = TRUE)
  })
  
  
  s.genes <- character(0)
  g2m.genes <- character(0)
  
  # Cell cycle genes
  if (remove_cc_genes == TRUE) {
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  }  
    
  blacklist.genes <- unique(c(
    mt.genes,
    sexChr.genes,
    s.genes,
    g2m.genes,
    RPS.genes,
    RPL.genes
  ))
  
  # Only return genes present in the data
  blacklist.genes <- blacklist.genes[blacklist.genes %in% allGenes]
  
  message(paste("Blacklisted", length(blacklist.genes), "genes",
                "(MT:", length(mt.genes),
                "RPL:", length(RPL.genes),
                "RPS:", length(RPS.genes),
                "Sex:", length(intersect(sexChr.genes, allGenes)),
                "CC:", length(intersect(c(s.genes, g2m.genes), allGenes)), ")"))
  
  return(blacklist.genes)
}

###############################################################
# SIMPLE CLUSTERIZATION FOR SUBCLUSTERS
###############################################################

cluster_subcluster <- function(obj, output_dir, n_genes = 2000, features = 2000, dims = 1:40, conserve.memory = FALSE) {
  
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(clustree)
  library(glmGamPoi)
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
  obj <- PercentageFeatureSet(obj, pattern = "^RP[LS]", col.name = "percent.ribo")
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  obj$CC.Difference <- obj$S.Score - obj$G2M.Score
  
  obj <- SCTransform(
    obj,
    n_genes = 2000,
    vars.to.regress = c("percent.mt", "percent.ribo", "CC.Difference"), 
    variable.features.n = 2000,
    conserve.memory = FALSE,
    method = "glmGamPoi",
    seed.use = 123
    
  )
  
  # Remove blacklist genes from variable features so they don't
  # drive PCA/clustering. Genes remain in the object for downstream use.
  blacklist <- get_blacklist_genes(obj)
  var_feats <- VariableFeatures(obj)
  VariableFeatures(obj) <- setdiff(var_feats, blacklist)
  message(paste("Variable features:", length(var_feats), "->", 
                length(VariableFeatures(obj)), 
                "(removed", length(intersect(var_feats, blacklist)), "blacklisted)"))
  
  obj <- RunPCA(obj, features = VariableFeatures(obj), seed.use = 123)
  
  obj <- RunUMAP(obj, dims = dims)
  obj <- FindNeighbors(obj, dims = dims)
  obj <- FindClusters(obj, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), algorithm = 1)
  
  dir.create(output_dir, recursive = TRUE)
  
  png(
    filename = file.path(output_dir, ("clustree.png")),
    width = 5,
    height = 5,
    units = "in",
    res = 300
  )
  
  p <- clustree(obj, prefix = "SCT_snn_res.")
  print(p)
  
  dev.off()
  
  return(obj)
  
}

subcluster_and_markers <- function(
    obj,
    cluster_name,
    cluster_col = "fine_clust",
    resolution = 0.5,
    dims = 1:20,
    npcs = 30,
    min_pct = 0.1,
    logfc_threshold = 0.5
) {
  
  message("Subsetting cluster: ", cluster_name)
  
  # identify cells
  cells_use <- colnames(obj)[obj@meta.data[[cluster_col]] == cluster_name]
  
  if (length(cells_use) == 0) {
    stop(paste("No cells found for", cluster_name))
  }
  
  # subset object
  sub_obj <- subset(obj, cells = cells_use)
  
  # recluster
  sub_obj <- RunPCA(sub_obj, npcs = npcs, verbose = FALSE)
  sub_obj <- FindNeighbors(sub_obj, dims = dims)
  sub_obj <- FindClusters(sub_obj, resolution = resolution)
  sub_obj <- RunUMAP(sub_obj, dims = dims)
  
  # markers
  sub_markers <- FindAllMarkers(
    sub_obj,
    only.pos = TRUE,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold
  )
  
  return(list(
    sub_obj = sub_obj,
    sub_markers = sub_markers
  ))
}

insert_subclusters <- function(
    obj,
    sub_obj,
    sub_map,
    new_col = "fine_clust"
) {
  
  # get subcluster IDs
  sub_ids <- as.character(Idents(sub_obj))
  
  # convert to biological labels
  new_labels <- sub_map[sub_ids]
  
  if (any(is.na(new_labels))) {
    stop("Some subclusters are missing from sub_map")
  }
  
  # insert back into original object
  obj[[new_col]][colnames(sub_obj), 1] <- new_labels
  
  return(obj)
}
##################################################
# MARKER GENE DENSITY PLOTS (NEBULOSA)
##################################################

plot_marker_genes <- function(obj,
                              genes,
                              cluster_col = "seurat_clusters",
                              reduction = "umap",
                              output_dir = "./marker_genes",
                              pt_size = 0.3,
                              show_labels = TRUE,
                              label_size = 4) {
  
  library(Nebulosa)
  library(ggplot2)
  library(viridis)
  library(shadowtext)
  library(dplyr)
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  available_genes <- rownames(obj)
  inferno_cols <- c("black", viridis::inferno(100))
  
  # Cluster label coordinates
  umap_coords <- as.data.frame(Embeddings(obj, reduction = reduction))
  colnames(umap_coords) <- c("x", "y")
  umap_coords$cluster <- as.factor(obj@meta.data[[cluster_col]])
  
  label_coords <- umap_coords %>%
    group_by(cluster) %>%
    summarize(x = median(x), y = median(y), .groups = "drop")
  
  for (gene in genes) {
    
    if (!gene %in% available_genes) {
      message(sprintf("Gene not found: %s", gene))
      next
    }
    
    tryCatch({
      
      p <- plot_density(obj, gene, reduction = reduction, size = pt_size) +
        scale_color_gradientn(colors = inferno_cols) +
        ggtitle(gene) +
        theme(
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)
        )
      
      if (show_labels) {
        p <- p +
          geom_shadowtext(
            data = label_coords,
            aes(x = x, y = y, label = cluster),
            color = "white",
            bg.color = "black",
            fontface = "bold",
            size = label_size,
            inherit.aes = FALSE
          )
      }
      
      ggsave(
        filename = file.path(output_dir, paste0(gene, "_density.png")),
        plot = p,
        width = 5,
        height = 5,
        dpi = 300
      )
      
      message(sprintf("Saved: %s", gene))
      
    }, error = function(e) {
      message(sprintf("Failed to plot %s: %s", gene, e$message))
    })
  }
}

#############################################
# CELL FILTERING FUNCTION
#############################################

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


#############################################
# SAMPLE BALANCE REPORT
#############################################

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


#############################################
# PSEUDOBULK DE ANALYSIS (PRIMARY METHOD)
#############################################

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
    
    # Pseudobulk aggregation 
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
  
  # Global pseudobulk (all cell types combined)
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


#############################################
# MAST DE ANALYSIS (FALLBACK / COMPLEMENT)
#############################################

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


#############################################
# PLOTTING FUNCTIONS
#############################################

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
    logfc_cap <- max(logfc_upper, logfc_thresh * 2.5)
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
      values = c("Up" = "#00BEC4", "Down" = "#F8766C", "NS" = "grey70"),
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
      size = 4,
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
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      legend.position = "top",
      plot.caption = element_text(hjust = 0, size = 10, color = "grey40")
    )
  
  max_x <- max(abs(de_df$avg_log2FC_plot), na.rm = TRUE) * 1.05
  p <- p + xlim(-max_x, max_x)
  
  return(p)
}

#############################################
# DENSITY PLOT 
#############################################

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

#############################################
# GENERIC HEATMAP FOR EXPRESSION DATA
#############################################

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


#############################################
# METHOD CONCORDANCE PLOT
#############################################

# Compare pseudobulk vs MAST results when both are run

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


#############################################
# MAIN DE WRAPPER FUNCTION
#############################################

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
    pval_col = "p_val_BH",            
    run_pseudobulk = FALSE,
    run_mast = TRUE,                   
    min_cells_per_sample = 10,
    min_samples_per_condition = 2,
    output_dir = "de_results"
) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  #############################################
  # Validate sample column
  #############################################
  
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
  
  ##############################################
  # Filter cells/clusters
  ##############################################
  
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
  
  #############################################
  # Run DE analysis
  #############################################
  
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
  
  #############################################
  # Save DE tables
  #############################################
  
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
  
  #############################################
  # GENERATE DE PLOTS
  #############################################
  
  all_plots <- list()
  
  # Global plots
  
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
  ggsave(file.path(output_dir, "Global_density_volcano.png"), combined_global, width = 12, height = 5, dpi = 330)
  
  if (!is.null(heatmap_global)) {
    pdf(file.path(output_dir, "Global_heatmap.pdf"), width = 10, height = 8)
    draw(heatmap_global)
    dev.off()
  }
  
  # Per cell type plots 
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
  
  # Method concordance plots
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
  
  #############################################
  # Summary
  #############################################
  
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

####################################################################
# MARKER GENE PLOT SIMILAR TO PYTHON SCANPY
####################################################################

plot_gene_bars_seurat <- function(
    seurat_obj,
    gene_groups,
    large_cluster_col,
    fine_cluster_col,
    outdir = "./marker_genes"
) {
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Loop over each large cluster
  for (grp in names(gene_groups)) {
    message("Processing: ", grp)
    
    # Get cells belonging to the large cluster
    cells_use <- rownames(seurat_obj@meta.data)[
      seurat_obj@meta.data[[large_cluster_col]] == grp
    ]
    
    # Skip if no cells found
    if (length(cells_use) == 0) {
      warning("No cells found for: ", grp)
      next
    }
    
    # Subset Seurat object to the relevant cells
    obj_sub <- subset(seurat_obj, cells = cells_use)
    genes_use <- intersect(gene_groups[[grp]], rownames(seurat_obj))
    
    # Get expression data
    expr <- GetAssayData(obj_sub, layer = "data")[genes_use, , drop = FALSE]
    df <- as.data.frame(t(expr))
    df$cell <- rownames(df)
    
    # Merge with metadata
    meta <- obj_sub@meta.data
    meta$cell <- rownames(meta)
    df <- left_join(df, meta, by = "cell")
    
    # Convert to long format for ggplot
    df_long <- df %>%
      pivot_longer(cols = all_of(genes_use), names_to = "gene", values_to = "expression") %>%
      arrange(.data[[fine_cluster_col]])
    
    # Factorize cells for consistent plotting
    df_long$cell <- factor(df_long$cell, levels = unique(df_long$cell))
    df_long$cell_index <- as.numeric(df_long$cell)
    
    # Calculate mean expression per gene per fine cluster
    df_long <- df_long %>%
      group_by(gene, .data[[fine_cluster_col]]) %>%
      mutate(mean_expr = mean(expression)) %>%
      ungroup()
    
    # Determine cluster boundaries for top bar and bottom labels
    cluster_bounds <- df_long %>%
      distinct(cell_index, fine_cluster = .data[[fine_cluster_col]]) %>%
      group_by(fine_cluster) %>%
      summarise(
        xmin = min(cell_index) - 0.5,
        xmax = max(cell_index) + 0.5,
        xmid = mean(c(min(cell_index), max(cell_index))),
        .groups = "drop"
      )
    
    # Generate color palette for fine clusters
    palette_colors <- hue_pal()(nrow(cluster_bounds))
    names(palette_colors) <- cluster_bounds$fine_cluster
    
    # Top colored bar with reduced height
    cluster_bar_top <- ggplot(cluster_bounds) +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 0.025, fill = fine_cluster), color = NA) +
      scale_x_continuous(expand = c(0,0), limits = c(0, max(cluster_bounds$xmax))) +
      scale_y_continuous(expand = c(0,0), limits = c(0,0.15)) +
      theme_void() +
      theme(legend.position = "none", plot.margin = margin(-10,2,0,2))
    
    # Determine cell boundaries for vertical lines
    df_boundaries <- df_long %>%
      distinct(cell_index, .data[[fine_cluster_col]]) %>%
      group_by(.data[[fine_cluster_col]]) %>%
      summarise(max_x = max(cell_index), .groups = "drop")
    
    # Main faceted gene plot with vertical cluster lines and mean expression
    p <- ggplot(df_long, aes(x = cell_index, y = expression)) +
      geom_col(aes(fill = .data[[fine_cluster_col]]), width = 1) +
      geom_step(aes(y = mean_expr), color = "black", linewidth = 0.5) +
      geom_vline(data = df_boundaries, aes(xintercept = max_x), color = "black", linewidth = 0.3) +
      facet_grid(gene ~ ., scales = "free_y", switch = "y") +
      scale_fill_manual(values = palette_colors) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(breaks = function(x) round(c(min(x), max(x)))) +
      labs(x = NULL, y = NULL) +
      theme_bw() +
      theme(
        axis.text.x = element_blank(), # remove x-axis numbers
        axis.ticks.x = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.right = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 10, face = "bold"), # horizontal labels
        panel.spacing = unit(0,"lines"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,5,0,5)
      )
    
    # Bottom cluster labels closer to gene plot
    cluster_bar_bottom <- ggplot(cluster_bounds, aes(x = xmid, y = 0, label = fine_cluster)) +
      geom_text(angle = 45, size = 3, fontface = "bold", vjust = 0.2, hjust = 1) +
      scale_x_continuous(expand = c(0,0), limits = c(0, max(cluster_bounds$xmax))) +
      scale_y_continuous(expand = c(0,0), limits = c(-0.15,0.05)) +
      theme_void() +
      coord_cartesian(clip = "off")
    
    # Combine top bar, main plot, and bottom labels
    final_plot <- cluster_bar_top / p / cluster_bar_bottom + plot_layout(heights = c(0.5, 4, 0.5))
    
    # Create subdirectory for the group
    subdir <- file.path(outdir, tolower(grp))
    if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)
    
    # Save plot as PDF
    fname <- paste0(tolower(grp), "_markers.pdf")
    ggsave(filename = file.path(subdir, fname),
           plot = final_plot,
           width = 5,
           height = max(6, length(genes_use) * 0.7))
  }
}



###############################################################
# INTEGRATE SCRNA DATA USING SCANORAMA (PYTHON SCRIPT)
###############################################################

scrna_integrate <- function(object.list, output_folder = "./", dataset_names, 
                            python_script_path = "./integrate_scanorama.py",
                            python_path = NULL) { 
  
  cat("\n##################################\n")
  cat("Normalizing and saving filtered objects...\n")
  cat("##################################\n")
  
  sc <- reticulate::import("scanpy")
  temp_dir <- file.path(output_folder, "preprocessing", "filtered_for_integration")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (i in seq_along(object.list)) {
    obj <- object.list[[i]]
    cat("Processing", dataset_names[i], "...\n")
    
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- JoinLayers(obj)
    data_matrix <- LayerData(obj, assay = "RNA", layer = "data")
    
    gene_names <- rownames(data_matrix)
    
    h5_path <- file.path(temp_dir, paste0(dataset_names[i], "_qc_filtered.h5"))
    DropletUtils::write10xCounts(
      path = h5_path, 
      x = data_matrix, 
      gene.id = gene_names,
      gene.symbol = gene_names,
      version = "3", 
      overwrite = TRUE
    )
    
    cat("  Saved:", ncol(obj), "cells x", nrow(obj), "genes\n")
  }
  
  cat("\n##################################\n")
  cat("Running Scanorama via Python...\n")
  cat("##################################\n")
  
  if (is.null(python_path)) {
    python_cmd <- "python"
  } else {
    python_cmd <- python_path
  }
  
  cmd <- sprintf(
    "%s %s %s %s %s",
    shQuote(python_cmd),
    shQuote(python_script_path),
    shQuote(temp_dir),
    shQuote(file.path(output_folder, "preprocessing")),
    paste(shQuote(dataset_names), collapse = " ")
  )
  
  cat("Running:\n", cmd, "\n\n")
  exit_code <- system(cmd)
  
  if (exit_code != 0) {
    stop("Python script failed with exit code ", exit_code)
  }
  
  cat("\n##################################\n")
  cat("Loading integrated data into R...\n")
  cat("##################################\n")
  
  h5ad_path <- file.path(output_folder, "preprocessing", "integrated_scanorama.h5ad")
  adata <- sc$read_h5ad(h5ad_path)
  
  # Convert expression matrix (stays sparse)
  counts <- reticulate::py_to_r(adata$X)
  if (inherits(counts, "scipy.sparse.base.spmatrix") || 
      grepl("sparse", class(counts)[1], ignore.case = TRUE)) {
    counts <- as(counts, "CsparseMatrix")
  }
  counts <- Matrix::t(counts)
  
  gene_names <- adata$var_names$to_list()
  cell_names <- adata$obs_names$to_list()
  
  rownames(counts) <- gene_names
  colnames(counts) <- cell_names
  
  cat("First 5 genes:", paste(head(gene_names, 5), collapse = ", "), "\n")
  
  # Get metadata with leiden clusters
  metadata <- reticulate::py_to_r(adata$obs)
  rownames(metadata) <- cell_names
  
  if ("leiden" %in% colnames(metadata)) {
    metadata$seurat_clusters <- factor(metadata$leiden)
  }
  
  # Create Seurat object
  integrated_seurat <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata
  )
  
  # Transfer UMAP
  umap_coords <- reticulate::py_to_r(adata$obsm$get("X_umap"))
  colnames(umap_coords) <- paste0("UMAP_", 1:2)
  rownames(umap_coords) <- colnames(integrated_seurat)
  integrated_seurat[["umap"]] <- CreateDimReducObject(
    embeddings = umap_coords,
    key = "UMAP_",
    assay = "RNA"
  )
  
  # Transfer scanorama embedding (use as "pca" equivalent for downstream)
  scanorama_coords <- reticulate::py_to_r(adata$obsm$get("X_scanorama"))
  colnames(scanorama_coords) <- paste0("PC_", 1:ncol(scanorama_coords))
  rownames(scanorama_coords) <- colnames(integrated_seurat)
  integrated_seurat[["pca"]] <- CreateDimReducObject(
    embeddings = scanorama_coords,
    key = "PC_",
    assay = "RNA"
  )
  
  # Transfer neighbor graph
  connectivities <- reticulate::py_to_r(adata$obsp$get("connectivities"))
  if (!is.null(connectivities)) {
    connectivities <- as(connectivities, "CsparseMatrix")
    rownames(connectivities) <- colnames(integrated_seurat)
    colnames(connectivities) <- colnames(integrated_seurat)
    integrated_seurat@graphs$RNA_snn <- as.Graph(connectivities)
  }
  
  Idents(integrated_seurat) <- "seurat_clusters"
  
  # Normalize for downstream (FeaturePlot etc.)
  integrated_seurat <- NormalizeData(integrated_seurat, verbose = FALSE)
  
  cat("\n##################################\n")
  cat("Integration complete!\n")
  cat(sprintf("Cells: %d | Genes: %d | Clusters: %d\n", 
              ncol(integrated_seurat), 
              nrow(integrated_seurat),
              length(unique(integrated_seurat$seurat_clusters))))
  cat("##################################\n\n")
  
  return(integrated_seurat)
}

###################################################
# GROUPED DOTPLOT FOR INITIAL IDENTIFICATION
###################################################

grouped_dotplot <- function(obj,
                            group_def,
                            gene_groups,
                            ident_col   = NULL,
                            cluster_map = NULL,
                            bar_colors  = NULL,
                            col_low     = "steelblue",
                            col_mid     = "white",
                            col_high    = "darkred",
                            dot_range   = c(0.5, 6),
                            base_size   = 11,
                            bar_height  = 1.2) {
  
  suppressPackageStartupMessages({
    require(Seurat)
    require(ggplot2)
    require(dplyr)
  })
  
  if (!is.null(cluster_map) && is.null(ident_col)) {
    mapping <- unlist(cluster_map)
    obj$grouped_id <- mapping[as.character(Idents(obj))]
    ident_col <- "grouped_id"
  }
  if (is.null(ident_col)) {
    stop("Supply either `ident_col` or `cluster_map`.")
  }
  if (!ident_col %in% colnames(obj@meta.data)) {
    stop(paste0("Column '", ident_col, "' not found in obj@meta.data. ",
                "Either create it first or use `cluster_map` instead."))
  }
  
  cluster_order <- unlist(group_def, use.names = FALSE)
  obj@meta.data[[ident_col]] <- factor(obj@meta.data[[ident_col]],
                                       levels = cluster_order)
  Idents(obj) <- ident_col
  
  all_genes <- unlist(gene_groups, use.names = FALSE)
  n_genes   <- length(all_genes)
  
  p_tmp  <- DotPlot(obj, features = all_genes) + RotatedAxis()
  dot_df <- p_tmp$data
  dot_df <- dot_df %>%
    dplyr::rename(gene = features.plot, cluster = id)
  
  gene_to_group <- stack(gene_groups) %>%
    dplyr::rename(gene = values, gene_group = ind)
  cluster_to_group <- stack(group_def) %>%
    dplyr::rename(cluster = values, cluster_group = ind)
  
  dot_df <- dot_df %>%
    dplyr::left_join(gene_to_group,    by = "gene") %>%
    dplyr::left_join(cluster_to_group, by = "cluster")
  
  dot_df$gene    <- factor(dot_df$gene,    levels = rev(all_genes))
  dot_df$cluster <- factor(dot_df$cluster, levels = cluster_order)
  dot_df$gene_group    <- factor(dot_df$gene_group,
                                 levels = names(gene_groups))
  dot_df$cluster_group <- factor(dot_df$cluster_group,
                                 levels = names(group_def))
  
  if (is.null(bar_colors)) {
    palette <- c("#E8A0BF","#A0C4E8","#A8E6CF","#FFD3B6","#D5AAFF",
                 "#FFE0AC","#B5EAD7","#C7CEEA","#FFDAC1","#E2F0CB")
    bar_colors <- setNames(palette[seq_along(group_def)],
                           names(group_def))
  }
  
  x_bar <- data.frame(
    group = names(group_def),
    xmin  = sapply(group_def, function(x)
      which(cluster_order == x[1]) - 0.5),
    xmax  = sapply(group_def, function(x)
      which(cluster_order == x[length(x)]) + 0.5),
    fill  = bar_colors[names(group_def)],
    stringsAsFactors = FALSE
  )
  
  v_breaks <- cumsum(sapply(group_def, length))
  v_breaks <- v_breaks[-length(v_breaks)] + 0.5
  
  gene_group_sizes <- sapply(gene_groups, length)
  h_breaks <- cumsum(rev(gene_group_sizes))
  h_breaks <- h_breaks[-length(h_breaks)] + 0.5
  
  bar_y_lo <- n_genes + 0.6
  bar_y_hi <- bar_y_lo + bar_height
  
  p <- ggplot(dot_df, aes(x = cluster, y = gene)) +
    geom_point(aes(size = pct.exp, fill = avg.exp.scaled),
               shape = 21, color = "grey30", stroke = 0.3) +
    scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                         midpoint = 0, name = "Scaled\nExpression") +
    scale_size_continuous(range = dot_range, name = "% Expressed",
                          breaks = c(10, 25, 50, 75, 100)) +
    {if (length(v_breaks) > 0)
      geom_vline(xintercept = v_breaks,
                 linewidth = 0.6, color = "grey20")} +
    {if (length(h_breaks) > 0)
      geom_hline(yintercept = h_breaks,
                 linewidth = 0.6, color = "grey20")} +
    annotate("rect",
             xmin = x_bar$xmin, xmax = x_bar$xmax,
             ymin = bar_y_lo,   ymax = bar_y_hi,
             fill = x_bar$fill, color = "grey30", linewidth = 0.3) +
    annotate("text",
             x = (x_bar$xmin + x_bar$xmax) / 2,
             y = (bar_y_lo + bar_y_hi) / 2,
             label = x_bar$group,
             fontface = "bold", size = 5) +
    coord_cartesian(ylim = c(0.5, bar_y_hi + 0.15), clip = "off") +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 12,
                                  face = "bold"),
      axis.text.y  = element_text(face = "bold", size = 10),
      axis.title   = element_blank(),
      panel.grid.major = element_line(linewidth = 0.15,
                                      color = "grey85"),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      legend.key.height = unit(0.35, "cm"),
      legend.key.width  = unit(0.35, "cm"),
      legend.title = element_text(size = 8.5),
      legend.text  = element_text(size = 7.5),
      plot.margin  = margin(t = 25, r = 10, b = 5, l = 5)
    ) +
    guides(
      fill = guide_colorbar(order = 1, barheight = unit(3, "cm")),
      size = guide_legend(order = 2,
                          override.aes = list(fill = "grey70"))
    )
  
  return(p)
}

##################################################
# MERGE AND PREPROCESS FOR INTEGRATION
##################################################

merge_and_preprocess <- function(object_list, n_variable_features = 3000, n_pcs = 40) {
  
  cat("Merging", length(object_list), "objects...\n")
  merged <- merge(
    x = object_list[[1]],
    y = object_list[2:length(object_list)],
    add.cell.ids = sapply(object_list, function(obj) unique(obj$orig.ident))
  )
  cat("Merged object:", ncol(merged), "cells,", nrow(merged), "genes\n")
  cat("Layers present:", paste(Layers(merged), collapse = ", "), "\n")
  
  # Normalize on split layers 
  merged <- NormalizeData(merged, verbose = FALSE)
  merged <- FindVariableFeatures(merged, nfeatures = n_variable_features, verbose = FALSE)
  merged <- ScaleData(merged, verbose = FALSE)
  merged <- RunPCA(merged, npcs = n_pcs, verbose = FALSE)
  
  # Pre-integration UMAP and clusters from raw PCA
  merged <- RunUMAP(merged, dims = 1:n_pcs, reduction = "pca",
                    reduction.name = "umap.unintegrated", verbose = FALSE)
  merged <- FindNeighbors(merged, dims = 1:n_pcs, reduction = "pca", verbose = FALSE)
  merged <- FindClusters(merged, resolution = 0.5, cluster.name = "unintegrated_clusters", verbose = FALSE)
  
  cat("Preprocessing complete.\n")
  return(merged)
}

##################################################
# HARMONY INTEGRATION VIA SEURAT V5
##################################################

run_harmony_integration <- function(merged, n_pcs = 30, theta = 2, resolution = 0.5) {
  
  cat("Running Harmony via IntegrateLayers...\n")
  merged <- IntegrateLayers(
    object = merged,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    theta = theta,
    verbose = FALSE
  )
  
  # Join layers after integration for downstream use
  merged[["RNA"]] <- JoinLayers(merged[["RNA"]])
  
  # Post-integration neighbors, clusters, UMAP from harmony reduction
  merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
  merged <- FindClusters(merged, resolution = resolution, cluster.name = "harmony_clusters", verbose = FALSE)
  merged <- RunUMAP(merged, reduction = "harmony", dims = 1:n_pcs,
                    reduction.name = "umap.harmony", verbose = FALSE)
  
  cat("Harmony integration complete.\n")
  return(merged)
}

#################################################################
# CALCULTING AVERAGE EXPRESSION AND PERCENT EXPRESSION PER GROUP
#################################################################

#' For each gene and each group (cluster), computes:
#' - avgExpr: mean expression across cells in that group (optionally normalized)
#' - pctExpr: percentage of cells with expression > 0
#'
#' @param mat Expression/score matrix (genes x cells)
#' @param groups Vector of group assignments per cell
#' @param feature_normalize Scale each gene to 0-1 range across groups
#' @param min_pct Minimum percent expressing to include
#' 
avg_prct_expression <- function(mat, groups, feature_normalize = TRUE, min_pct = 0) {
  
  # Ensure groups align with matrix columns
  groups <- as.character(groups)
  unique_groups <- unique(groups) %>% sort()
  
  results <- lapply(unique_groups, function(grp) {
    
    # Cells in this group
    cells_in_grp <- which(groups == grp)
    
    if (length(cells_in_grp) == 0) return(NULL)
    
    # Subset matrix to this group
    sub_mat <- mat[, cells_in_grp, drop = FALSE]
    
    # Average expression per gene
    avg_expr <- rowMeans(sub_mat, na.rm = TRUE)
    
    # Percent of cells expressing (value > 0)
    pct_expr <- rowSums(sub_mat > 0, na.rm = TRUE) / ncol(sub_mat) * 100
    
    data.frame(
      feature = rownames(mat),
      grp = grp,
      avgExpr = avg_expr,
      pctExpr = pct_expr,
      stringsAsFactors = FALSE
    )
  })
  
  df <- do.call(rbind, results)
  
  # Feature normalization: scale each gene's avgExpr to 0-1 range
  # This allows comparison across genes with different expression magnitudes
  if (feature_normalize) {
    df <- df %>%
      group_by(feature) %>%
      mutate(avgExpr = (avgExpr - min(avgExpr)) / (max(avgExpr) - min(avgExpr) + 1e-9)) %>%
      ungroup()
  }
  
  # Filter by minimum percent expressing
  if (min_pct > 0) {
    df <- df %>% filter(pctExpr >= min_pct)
  }
  
  return(as.data.frame(df))
}

#################################################################
# DOTPLOT TO VISUALIZE GENE MARKER EXPRESSION
#################################################################

plot_dotplot <- function(df, xcol, ycol, color_col, size_col, 
                         xorder = NULL, yorder = NULL, cmap = NULL, 
                         color_label = NULL, size_label = NULL, 
                         aspectRatio = NULL, sizeLims = NULL, colorLims = NULL) {
  
  # Set order of axes (important for interpretability)
  if (is.null(xorder)) xorder <- unique(df[, xcol]) %>% sort()
  if (is.null(yorder)) yorder <- unique(df[, ycol]) %>% sort()
  if (is.null(aspectRatio)) aspectRatio <- length(yorder) / length(xorder)
  
  df[, xcol] <- factor(df[, xcol], levels = xorder)
  df[, ycol] <- factor(df[, ycol], levels = yorder)
  df <- df[order(df[, xcol], df[, ycol]), ]
  
  # Build plot
  p <- ggplot(df, aes(
    x = .data[[xcol]], 
    y = .data[[ycol]], 
    color = .data[[color_col]], 
    size = ifelse(.data[[size_col]] > 0, .data[[size_col]], NA)
  )) +
    geom_point() +
    xlab(xcol) +
    ylab(ycol) +
    theme_BOR(border = TRUE) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.25, 0.5, 0.25, 1), "cm"),
      aspect.ratio = aspectRatio,
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    guides(
      colour = guide_colorbar(title = color_label),
      size = guide_legend(title = size_label)
    )
  
  # Color scale
  if (!is.null(cmap)) {
    if (!is.null(colorLims)) {
      p <- p + scale_color_gradientn(colors = cmap, limits = colorLims, oob = scales::squish)
    } else {
      p <- p + scale_color_gradientn(colors = cmap)
    }
  } else {
    p <- p + scale_color_viridis_c(option = "magma")
  }
  
  # Size scale
  if (!is.null(sizeLims)) {
    p <- p + scale_size_continuous(limits = sizeLims, range = c(0.5, 5))
  } else {
    p <- p + scale_size_continuous(range = c(0.5, 5))
  }
  
  return(p)
}

##################################################
# PRE VS POST INTEGRATION UMAP COMPARISON
##################################################

plot_integration_comparison <- function(merged, group_by = "orig.ident", output_path = NULL) {
  
  p1 <- DimPlot(merged, reduction = "umap.unintegrated", group.by = group_by, pt.size = 0.1) +
    ggtitle("Pre-integration — by sample") +
    theme(legend.position = "bottom")
  
  p2 <- DimPlot(merged, reduction = "umap.unintegrated", group.by = "unintegrated_clusters", pt.size = 0.1) +
    ggtitle("Pre-integration — clusters") +
    theme(legend.position="none")
  
  p3 <- DimPlot(merged, reduction = "umap.harmony", group.by = group_by, pt.size = 0.1) +
    ggtitle("Post-integration (Harmony) — by sample") +
    theme(legend.position = "bottom")
  
  p4 <- DimPlot(merged, reduction = "umap.harmony", group.by = "harmony_clusters", pt.size = 0.1) +
    ggtitle("Post-integration (Harmony) — clusters") + 
    theme(legend.position="none")
  
  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(theme = theme(plot.title = element_text(size = 16, face = "bold")))
  
  if (!is.null(output_path)) {
    ggsave(filename = output_path, plot = combined, width = 14, height = 12, dpi = 300)
    cat("Saved plot to:", output_path, "\n")
  }
  
  print(combined)
  return(combined)
}

###################################################
# BINOMIAL ENRICHMENT TEST PER CLUSTER
###################################################

#' Tests whether each cluster is statistically enriched in one of two
#' conditions using a binomial distribution, adapted from Joost et al.,
#' Cell Stem Cell 2020 (STAR Methods).
#'
#' Under the null, cells from each condition distribute across clusters
#' proportionally to their global representation.
#'
#' Two modes for computing expected fraction (p):
#'
#'   use_global_fraction = TRUE (default):
#'     p = total_condA_cells / total_cells
#'     Appropriate when samples were loaded at roughly equal density
#'     and no tissue-level cellularity correction is needed.
#'
#'   use_global_fraction = FALSE:
#'     Uses the Joost et al. %_a formula which cross-references
#'     sequenced counts with external cell isolation counts (cell_counts).
#'     Required when conditions differ in tissue cellularity and capture
#'     does not reflect that difference.
#'
#' Outlier clusters that are known to be condition-specific and large
#' enough to skew the global baseline can be excluded from p_expected
#' computation via exclude_clusters. These clusters are still tested
#' and reported in the output, but do not influence the baseline.
#'
#' @param obj                 Seurat v5 object
#' @param treatment_col       Column in meta.data with two condition labels
#' @param cluster_col         Column in meta.data with cluster identities
#' @param sample_col          Column in meta.data with sample/replicate IDs
#' @param condition_a         Label for condition A. If NULL, first sorted
#'                            unique value is used.
#' @param condition_b         Label for condition B. If NULL, second sorted
#'                            unique value is used.
#' @param exclude_clusters    Character vector of cluster names to exclude
#'                            from p_expected computation. These clusters
#'                            are still tested and included in results.
#' @param use_global_fraction Logical (default TRUE). If TRUE, p is computed
#'                            as the simple global fraction of condition_a
#'                            cells. If FALSE, the Joost et al. cross-sample
#'                            formula is used (supply cell_counts for full
#'                            correction).
#' @param cell_counts         Named numeric vector of total isolated cells per
#'                            sample. Only used when use_global_fraction = FALSE.
#'                            Names must match sample_col values.
#' @param alpha               Significance threshold (default 0.001)
#'

test_cluster_enrichment <- function(obj,
                                    treatment_col       = "treatment",
                                    cluster_col         = "fine_clust",
                                    sample_col          = "orig.ident",
                                    condition_a         = NULL,
                                    condition_b         = NULL,
                                    exclude_clusters    = NULL,
                                    use_global_fraction = TRUE,
                                    cell_counts         = NULL,
                                    alpha               = 0.001) {
  
  suppressPackageStartupMessages({
    require(Seurat)
    require(dplyr)
  })
  
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  md <- obj@meta.data
  
  for (col in c(treatment_col, cluster_col, sample_col)) {
    if (!col %in% colnames(md)) {
      stop(paste0("Column '", col, "' not found in meta.data."))
    }
  }
  
  # Resolve condition labels
  conditions <- if (is.factor(md[[treatment_col]])) {
    levels(md[[treatment_col]])
  } else {
    sort(unique(md[[treatment_col]]))
  }
  if (length(conditions) != 2) {
    stop("Expected exactly 2 conditions in '", treatment_col,
         "', found: ", paste(conditions, collapse = ", "))
  }
  if (is.null(condition_a)) condition_a <- conditions[1]
  if (is.null(condition_b)) condition_b <- conditions[2]
  
  # Baseline metadata: exclude outlier clusters from p computation
  if (!is.null(exclude_clusters)) {
    baseline_idx <- !md[[cluster_col]] %in% exclude_clusters
    md_baseline  <- md[baseline_idx, ]
    n_excluded   <- sum(!baseline_idx)
    cat("Excluding", length(exclude_clusters), "cluster(s) from baseline:",
        paste(exclude_clusters, collapse = ", "), "\n")
    cat("  ->", n_excluded, "cells excluded from p_expected computation\n")
  } else {
    md_baseline <- md
  }
  
  # Per-sample sequenced cell counts (from baseline cells only)
  sample_tab <- md_baseline %>%
    group_by(across(all_of(c(sample_col, treatment_col)))) %>%
    summarise(n_seq = n(), .groups = "drop")
  
  samples_a <- sample_tab %>% filter(.data[[treatment_col]] == condition_a)
  samples_b <- sample_tab %>% filter(.data[[treatment_col]] == condition_b)
  
  S_a <- samples_a$n_seq
  S_b <- samples_b$n_seq
  
  n_total_baseline <- nrow(md_baseline)
  n_a_baseline     <- sum(md_baseline[[treatment_col]] == condition_a)
  
  # Compute p_expected
  if (use_global_fraction) {
    p_expected <- n_a_baseline / n_total_baseline
    cat("Mode: global fraction",
        if (!is.null(exclude_clusters)) "(baseline-corrected)" else "",
        "\n")
  } else {
    if (!is.null(cell_counts)) {
      C_a <- cell_counts[samples_a[[sample_col]]]
      C_b <- cell_counts[samples_b[[sample_col]]]
      if (any(is.na(C_a)) || any(is.na(C_b))) {
        stop("cell_counts names must match all sample IDs in '",
             sample_col, "'.")
      }
    } else {
      C_a <- S_a
      C_b <- S_b
    }
    p_expected <- (sum(S_a) * sum(C_b)) /
      (sum(S_a) * sum(C_b) + sum(S_b) * sum(C_a))
    cat("Mode: Joost et al. cross-sample correction",
        if (!is.null(exclude_clusters)) "(baseline-corrected)" else "",
        "\n")
  }
  
  cat("Condition A:", condition_a, "| Condition B:", condition_b, "\n")
  cat("Baseline cells —",
      condition_a, ":", n_a_baseline, "|",
      condition_b, ":", n_total_baseline - n_a_baseline,
      "(total:", n_total_baseline, ")\n")
  cat("Per-sample baseline cells —",
      condition_a, ":", paste(S_a, collapse = ", "), "|",
      condition_b, ":", paste(S_b, collapse = ", "), "\n")
  cat("p_expected (fraction", condition_a, "):", round(p_expected, 4), "\n")
  
  # Per-cluster binomial test (ALL clusters, including excluded ones)
  clusters <- sort(unique(md[[cluster_col]]))
  
  results <- do.call(rbind, lapply(clusters, function(cl) {
    
    idx <- md[[cluster_col]] == cl
    n   <- sum(idx)
    n_a <- sum(md[[treatment_col]][idx] == condition_a)
    n_b <- n - n_a
    
    # P(X >= n_a): if small, cluster enriched in condition_a
    pval_a <- pbinom(n_a - 1, size = n, prob = p_expected,
                     lower.tail = FALSE)
    # P(X <= n_a): if small, cluster enriched in condition_b
    pval_b <- pbinom(n_a, size = n, prob = p_expected)
    
    enrichment <- "none"
    if (pval_a < alpha) enrichment <- condition_a
    if (pval_b < alpha) enrichment <- condition_b
    if (pval_a < alpha && pval_b < alpha) {
      enrichment <- ifelse(pval_a < pval_b, condition_a, condition_b)
    }
    
    excluded <- cl %in% exclude_clusters
    
    data.frame(
      cluster             = cl,
      n_total             = n,
      n_condA             = n_a,
      n_condB             = n_b,
      frac_condA          = n_a / n,
      p_expected          = p_expected,
      pval_enriched_condA = pval_a,
      pval_enriched_condB = pval_b,
      enrichment          = enrichment,
      excluded_from_baseline = excluded,
      stringsAsFactors    = FALSE
    )
  }))
  
  rownames(results) <- NULL
  attr(results, "condition_a") <- condition_a
  attr(results, "condition_b") <- condition_b
  attr(results, "alpha")       <- alpha
  
  n_a_enr <- sum(results$enrichment == condition_a)
  n_b_enr <- sum(results$enrichment == condition_b)
  n_none  <- sum(results$enrichment == "none")
  cat(n_a_enr, "cluster(s) enriched in", condition_a, "|",
      n_b_enr, "cluster(s) enriched in", condition_b, "|",
      n_none, "not enriched (alpha =", alpha, ")\n")
  
  if (!is.null(exclude_clusters)) {
    excl_res <- results %>% filter(excluded_from_baseline)
    cat("Excluded cluster results:\n")
    for (i in seq_len(nrow(excl_res))) {
      cat("  ", excl_res$cluster[i], ":",
          excl_res$n_condA[i], condition_a, "/",
          excl_res$n_condB[i], condition_b,
          "-> enriched in", excl_res$enrichment[i], "\n")
    }
  }
  
  return(results)
}

###################################################
# ENRICHMENT PLOT (JOOST ET AL. 2020 STYLE)
###################################################

#' Creates a mirrored -log10(p-value) plot where:
#'   - Positive y  = enriched in condition_a (top half)
#'   - Negative y  = enriched in condition_b (bottom half)
#'   - Background shading intensifies past significance threshold
#'   - Dashed horizontal lines at significance cutoff
#'   - Points colored by enrichment call
#'   - Clusters excluded from baseline are marked with a distinct shape
#'
#' @param results      Output data.frame from test_cluster_enrichment()
#' @param group_col    Optional column name in results for x-axis grouping
#'                     with labeled color bars. Add to results before calling.
#' @param group_colors Named vector of colors per group level. If NULL a
#'                     default palette is used.
#' @param alpha        Significance threshold for dashed lines
#' @param max_score    Cap for y-axis in -log10 scale (default 20)
#' @param condA_color  Color for condition_a enrichment (default "green4")
#' @param condB_color  Color for condition_b enrichment (default "red3")
#' @param ns_color     Color for non-significant points (default "black")
#' @param title        Plot title
#' @param base_size    Base font size for theme (default 11)
#'

plot_enrichment <- function(results,
                            group_col    = NULL,
                            group_colors = NULL,
                            alpha        = NULL,
                            max_score    = 20,
                            condA_color  = "green4",
                            condB_color  = "red3",
                            ns_color     = "black",
                            title        = "Condition enrichment of clusters",
                            base_size    = 11) {
  
  suppressPackageStartupMessages({
    require(ggplot2)
    require(dplyr)
  })
  
  condition_a <- attr(results, "condition_a")
  condition_b <- attr(results, "condition_b")
  if (is.null(alpha)) alpha <- attr(results, "alpha") %||% 0.001
  
  thresh <- -log10(alpha)
  
  # Signed enrichment score:
  # positive = more condition_a than expected, negative = fewer
  results <- results %>%
    mutate(
      score = case_when(
        frac_condA >= p_expected ~ -log10(pmax(pval_enriched_condA, 1e-300)),
        TRUE                    ~  log10(pmax(pval_enriched_condB, 1e-300))
      ),
      score_capped = pmax(pmin(score, max_score), -max_score),
      point_color  = case_when(
        enrichment == condition_a ~ "condA",
        enrichment == condition_b ~ "condB",
        TRUE                      ~ "ns"
      ),
      # Mark excluded clusters with different shape
      point_shape = ifelse(excluded_from_baseline, "excluded", "tested")
    )
  
  # Order clusters by group then name
  if (!is.null(group_col) && group_col %in% colnames(results)) {
    results <- results %>% arrange(.data[[group_col]], cluster)
  }
  results$cluster <- factor(results$cluster, levels = results$cluster)
  
  # Y-axis breaks and expression labels
  y_breaks <- sort(unique(c(-max_score, -15, -10, -thresh,
                            0, thresh, 10, 15, max_score)))
  y_breaks <- y_breaks[y_breaks >= -max_score & y_breaks <= max_score]
  
  y_labels <- sapply(y_breaks, function(v) {
    if (v == 0) return(expression(1))
    as.expression(bquote(10^{-.(abs(v))}))
  })
  
  # Base plot
  p <- ggplot(results, aes(x = cluster, y = score_capped)) +
    
    # Background shading
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = 0, ymax = thresh,
             fill = condA_color, alpha = 0.04) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = thresh, ymax = max_score,
             fill = condA_color, alpha = 0.10) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -thresh, ymax = 0,
             fill = condB_color, alpha = 0.04) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -max_score, ymax = -thresh,
             fill = condB_color, alpha = 0.10) +
    
    # Reference lines
    geom_hline(yintercept = thresh,  linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -thresh, linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    
    # Lollipop stems
    geom_segment(aes(x = cluster, xend = cluster,
                     y = 0, yend = score_capped),
                 color = "grey60", linewidth = 0.25) +
    
    # Points: shape distinguishes excluded vs tested clusters
    geom_point(aes(color = point_color, shape = point_shape), size = 2.5) +
    scale_color_manual(
      values = c("condA" = condA_color,
                 "condB" = condB_color,
                 "ns"    = ns_color),
      labels = c("condA" = paste0("Enriched in ", condition_a),
                 "condB" = paste0("Enriched in ", condition_b),
                 "ns"    = "Not significant"),
      name   = NULL
    ) +
    scale_shape_manual(
      values = c("tested" = 16, "excluded" = 17),
      labels = c("tested" = "Tested",
                 "excluded" = "Excluded from baseline"),
      name   = NULL
    ) +
    
    # Y-axis
    scale_y_continuous(
      limits = c(-max_score, max_score),
      breaks = y_breaks,
      labels = y_labels
    ) +
    
    # Direction annotations
    annotate("text", x = 0.5, y = max_score - 0.5,
             label = paste("Enriched in", condition_a),
             hjust = 0, vjust = 1, size = 2.5,
             fontface = "italic", color = "grey30") +
    annotate("text", x = 0.5, y = -max_score + 0.5,
             label = paste("Enriched in", condition_b),
             hjust = 0, vjust = 0, size = 2.5,
             fontface = "italic", color = "grey30") +
    
    labs(x = NULL, y = "P value [binomial]", title = title) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x         = element_text(angle = 45, hjust = 1, size = 9,
                                         face = "bold"),
      axis.text.y         = element_text(size = 8),
      axis.title          = element_text(size = 9),
      legend.position     = "bottom",
      legend.text         = element_text(size = 8),
      panel.grid.major.y  = element_line(color = "grey92", linewidth = 0.2),
      panel.grid.minor    = element_blank(),
      plot.title          = element_text(hjust = 0.5, face = "bold",
                                         size = 12),
      plot.margin         = margin(t = 10, r = 10, b = 30, l = 10)
    )
  
  # Group color bars below x-axis
  if (!is.null(group_col) && group_col %in% colnames(results)) {
    
    group_info <- results %>%
      group_by(.data[[group_col]]) %>%
      summarise(
        xmin = as.numeric(first(cluster)) - 0.4,
        xmax = as.numeric(last(cluster))  + 0.4,
        xmid = mean(c(as.numeric(first(cluster)),
                      as.numeric(last(cluster)))),
        .groups = "drop"
      )
    
    if (is.null(group_colors)) {
      palette <- c("#E8A0BF", "#A0C4E8", "#A8E6CF", "#FFD3B6",
                   "#D5AAFF", "#FFE0AC", "#B5EAD7", "#C7CEEA",
                   "#FFDAC1", "#E2F0CB")
      group_colors <- setNames(
        palette[seq_len(nrow(group_info))],
        group_info[[group_col]]
      )
    }
    
    bar_y <- -max_score - 2.0
    bar_h <- 1.2
    
    p <- p +
      annotate("rect",
               xmin = group_info$xmin, xmax = group_info$xmax,
               ymin = bar_y - bar_h / 2, ymax = bar_y + bar_h / 2,
               fill = group_colors[group_info[[group_col]]],
               color = "grey30", linewidth = 0.3) +
      annotate("text",
               x = group_info$xmid,
               y = bar_y,
               label = group_info[[group_col]],
               fontface = "bold", size = 3) +
      coord_cartesian(
        ylim = c(-max_score - 4, max_score),
        clip = "off"
      )
  }
  
  return(p)
}

###############################################################
# SUBCLUSTER FOR MILO PIPELINE (LIGHTWEIGHT)
# Use when input is already annotated with SCT layer
###############################################################

subcluster_for_milo <- function(obj, output_dir, dims = 1:30, resolution = 0.5, nfeatures = 2000) {
  
  library(Seurat)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(glmGamPoi)
  
  # Cell cycle genes
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  # QC metrics
  obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
  obj <- PercentageFeatureSet(obj, pattern = "^RP[LS]", col.name = "percent.ribo")
  
  # Normalize and score cell cycle
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  obj$CC.Difference <- obj$S.Score - obj$G2M.Score
  
  # SCTransform with regression
  obj <- SCTransform(
    obj,
    n_genes = nfeatures,
    vars.to.regress = c("percent.mt", "percent.ribo", "CC.Difference"),
    variable.features.n = nfeatures,
    conserve.memory = FALSE,
    method = "glmGamPoi",
    seed.use = 123,
    verbose = FALSE
  )
  
  # Remove blacklist genes from variable features
  if (exists("get_blacklist_genes", mode = "function")) {
    blacklist <- get_blacklist_genes(obj)
    var_feats <- VariableFeatures(obj)
    VariableFeatures(obj) <- setdiff(var_feats, blacklist)
    message(paste("Variable features:", length(var_feats), "->",
                  length(VariableFeatures(obj)),
                  "(removed", length(intersect(var_feats, blacklist)), "blacklisted)"))
  }
  
  # PCA, UMAP, clustering
  obj <- RunPCA(obj, features = VariableFeatures(obj), seed.use = 123, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims, seed.use = 123, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, algorithm = 1, verbose = FALSE)
  return(obj)
}
###############################################################
# ADAPTIVE MILO PARAMETERS BASED ON DATASET SIZE
###############################################################

get_milo_params <- function(n_cells, n_samples) {
  
  # Select k, prop for Milo based on dataset size.
  
  # Returns a named list with k and prop.
  # Note: d (PCA dims) is set separately via subcluster_params in main script.
  
  
  if (n_cells < 5000) {
    list(k = 20, prop = 0.2)
  } else if (n_cells < 10000) {
    list(k = 30, prop = 0.1)
  } else if (n_cells < 20000) {
    list(k = 30, prop = 0.1)
  } else {
    list(k = 40, prop = 0.1)
  }
}

###############################################################
# RUN MILO PIPELINE
###############################################################

run_milo_pipeline <- function(obj_sub,
                              k             = 50,
                              d             = 50,
                              prop          = 0.1,
                              sample_col    = "orig.ident",
                              treatment_col = "treatment",
                              pca_name      = "pca",
                              umap_name     = "umap") {
  
  DefaultAssay(obj_sub) <- "RNA"
  obj.sce <- as.SingleCellExperiment(obj_sub)
  
  # MiloR needs a logcounts assay. The Seurat→SCE conversion should
  # produce one from the RNA "data" layer, but if it's missing, add it.
  if (!"logcounts" %in% assayNames(obj.sce)) {
    cat("  logcounts missing — running logNormCounts...\n")
    obj.sce <- scuttle::logNormCounts(obj.sce)
  }
  
  # Inject Seurat reductions
  seurat_pca  <- Embeddings(obj_sub, reduction = pca_name)
  seurat_umap <- Embeddings(obj_sub, reduction = umap_name)
  
  # Cap d to available PCA components
  if (ncol(seurat_pca) < d) {
    d <- ncol(seurat_pca)
    cat("  PCA has", ncol(seurat_pca), "components — using d =", d, "\n")
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
  
  # Count + distances
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
                           norm.method = "RLE",
                           design    = ~ Treatment,
                           design.df = design_matrix)
  
  cat("  Sig neighborhoods (P < 0.05):",
      sum(da_results$PValue < 0.05, na.rm = TRUE),
      "| Up:", sum(da_results$PValue < 0.05 & da_results$logFC > 0, na.rm = TRUE),
      "| Down:", sum(da_results$PValue < 0.05 & da_results$logFC < 0, na.rm = TRUE), "\n")
  
  # Nhood graph + annotate
  obj.milo <- buildNhoodGraph(obj.milo)
  
  if ("fine_clust" %in% colnames(colData(obj.milo))) {
    da_results <- annotateNhoods(obj.milo, da_results,
                                 coldata_col = "fine_clust")
  }
  
  return(list(milo = obj.milo, da_results = da_results))
}

################################################################################
# MILO DA NEIGHBOURHOOD GRAPH
################################################################################

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


################################################################################
# PLOT SPLIT UMAP FOR MILO PIPELINE
################################################################################
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

############################################################################
# PLOT ORDERED DA BEESWARM (CLEANED UP)
############################################################################

plot_da_beeswarm_ordered <- function(da_results, 
                                     group_by = "fine_clust",
                                     title = NULL,
                                     alpha = 0.05,
                                     use_pvalue_as_fdr = FALSE,
                                     pt_size = 1.5) {
  
  library(ggplot2)
  library(ggbeeswarm)
  
  df <- as.data.frame(da_results)
  
  # Determine significance column
  if (use_pvalue_as_fdr) {
    df$sig <- df$PValue < alpha
  } else {
    df$sig <- df$SpatialFDR < alpha
  }
  
  # Reorder: non-significant first, significant last (plotted on top)
  df <- df[order(df$sig, decreasing = FALSE), ]
  
  # Handle missing group column
  if (!group_by %in% colnames(df)) {
    warning(paste0("Column '", group_by, "' not found. Using NhoodGroup if available."))
    if ("NhoodGroup" %in% colnames(df)) {
      group_by <- "NhoodGroup"
    } else {
      df$group <- "All"
      group_by <- "group"
    }
  }
  
  p <- ggplot(df, aes(x = .data[[group_by]], y = logFC, color = sig)) +
    geom_quasirandom(size = pt_size, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "grey70"),
      labels = c("TRUE" = paste0("p < ", alpha), "FALSE" = "NS"),
      name = NULL
    ) +
    labs(
      x = NULL,
      y = "Log Fold Change",
      title = title
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x    = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y    = element_text(size = 10, face = "bold"),
      axis.title.y   = element_text(size = 11, face = "bold"),
      plot.title     = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.position = "top",
      legend.text    = element_text(size = 10, face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  return(p)
}

################################################################################
# PLOT SPLIT UMAP FOR MILO PIPELINE (PRETTIER VERSION)
################################################################################

plot_split_dimred <- function(obj_sub,
                              split_by  = "treatment",
                              color_by  = "fine_clust",
                              reduction = "umap",
                              pt_size   = 0.5,
                              title     = NULL,
                              cmap      = NULL,
                              useRaster = TRUE) {
  
  library(ggplot2)
  library(ggrastr)
  
  # Default stallion palette
  if (is.null(cmap)) {
    cmap <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B",
              "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC","10"="#90D5E4",
              "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8",
              "16"="#6E4B9E","17"="#0C727C","18"="#7E1416","9"="#D8A767","20"="#3D3D3D")
  }
  
  emb <- Embeddings(obj_sub, reduction = reduction)
  df  <- data.frame(
    UMAP_1    = emb[, 1],
    UMAP_2    = emb[, 2],
    split_var = obj_sub@meta.data[[split_by]],
    color_var = obj_sub@meta.data[[color_by]],
    stringsAsFactors = FALSE
  )
  
  # Shuffle for random plotting order
  set.seed(42)
  df <- df[sample(nrow(df)), ]
  
  # Expand color palette if needed
  n_colors <- length(unique(df$color_var))
  if (n_colors > length(cmap)) {
    cmap <- colorRampPalette(cmap)(n_colors)
  } else {
    names(cmap) <- NULL
    cmap <- cmap[1:n_colors]
  }
  
  # Build plot
  if (useRaster) {
    p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = color_var)) +
      geom_point_rast(size = pt_size)
  } else {
    p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = color_var)) +
      geom_point(size = pt_size)
  }
  
  p <- p +
    facet_wrap(~ split_var, ncol = 2) +
    scale_color_manual(values = cmap, name = NULL, na.value = "grey35") +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme_classic(base_size = 11) +
    theme(
      axis.ticks       = element_blank(),
      axis.text        = element_blank(),
      aspect.ratio     = 1,
      strip.text       = element_text(face = "bold", size = 12),
      strip.background = element_blank(),
      legend.position  = "right",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 10, face = "bold"),
      legend.key.size  = unit(0.4, "cm"),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 13)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  # Title: default to cell count, or custom, or NULL
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}

################################################################################
# PLOT SINGLE UMAP (PRETTIER VERSION)
################################################################################

plot_single_dimred <- function(obj_sub,
                               color_by  = "fine_clust",
                               reduction = "umap",
                               pt_size   = 0.5,
                               title     = NULL,
                               cmap      = NULL,
                               useRaster = TRUE) {
  
  library(ggplot2)
  library(ggrastr)
  
  # Default stallion palette
  if (is.null(cmap)) {
    cmap <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B",
              "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","9"="#E6C2DC","10"="#90D5E4",
              "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8",
              "16"="#6E4B9E","17"="#0C727C","18"="#7E1416","19"="#D8A767","20"="#3D3D3D")
  }
  
  emb <- Embeddings(obj_sub, reduction = reduction)
  df  <- data.frame(
    UMAP_1    = emb[, 1],
    UMAP_2    = emb[, 2],
    color_var = obj_sub@meta.data[[color_by]],
    stringsAsFactors = FALSE
  )
  
  # Shuffle for random plotting order
  set.seed(42)
  df <- df[sample(nrow(df)), ]
  
  # Expand color palette if needed
  n_colors <- length(unique(df$color_var))
  if (n_colors > length(cmap)) {
    cmap <- colorRampPalette(cmap)(n_colors)
  } else {
    names(cmap) <- NULL
    cmap <- cmap[1:n_colors]
  }
  
  # Build plot
  if (useRaster) {
    p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = color_var)) +
      geom_point_rast(size = pt_size)
  } else {
    p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = color_var)) +
      geom_point(size = pt_size)
  }
  
  p <- p +
    scale_color_manual(values = cmap, name = NULL, na.value = "grey35") +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme_classic(base_size = 11) +
    theme(
      axis.ticks       = element_blank(),
      axis.text        = element_blank(),
      aspect.ratio     = 1,
      legend.position  = "right",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 10, face = "bold"),
      legend.key.size  = unit(0.4, "cm"),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 13)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}

################################################################################
# MILO WRAPPER PLOT FUNCTIONS
################################################################################

plot_nhood_size_hist <- function(milo_obj, title = "Neighbourhood Sizes") {
  plotNhoodSizeHist(milo_obj) +
    ggtitle(title) +
    theme_classic()
}

plot_pvalue_hist <- function(da_results, title = "P-value Distribution") {
  ggplot(da_results, aes(PValue)) +
    geom_histogram(bins = 50, fill = "#272E6A", color = "white") +
    theme_classic() +
    labs(title = title, x = "P-value", y = "Count")
}

plot_volcano <- function(da_results, pval_thresh = 0.05, title = "Volcano Plot") {
  ggplot(da_results, aes(logFC, -log10(PValue))) +
    geom_point(aes(color = PValue < pval_thresh), size = 2) +
    scale_color_manual(values = c("grey60", "#D51F26"),
                       labels = c("NS", paste0("P < ", pval_thresh)),
                       name = "Significance") +
    geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed", color = "red") +
    theme_classic() +
    labs(title = title, x = "Log Fold Change", y = "-log10(P-value)")
}

summarise_da_by_cluster <- function(da_results, cluster_col = "fine_clust") {
  if (!cluster_col %in% colnames(da_results)) return(NULL)
  
  da_results %>%
    group_by(.data[[cluster_col]]) %>%
    summarise(
      n_nhoods     = n(),
      n_sig_p0.05  = sum(PValue < 0.05, na.rm = TRUE),
      n_up         = sum(PValue < 0.05 & logFC > 0, na.rm = TRUE),
      n_down       = sum(PValue < 0.05 & logFC < 0, na.rm = TRUE),
      mean_logFC   = mean(logFC, na.rm = TRUE),
      median_logFC = median(logFC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_sig_p0.05))
}

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
################################################################################
# PLOT NHOOD GRAPH WITH P-VALUE FILTERING (ADJUSTABLE NODE SIZE)
################################################################################

plot_nhood_umap <- function(milo_obj, 
                                  da_results,
                                  alpha = 0.05,
                                  use_pvalue = TRUE,
                                  layout = "UMAP",
                                  node_size_range = c(0.5, 4),
                                  edge_width = 0.2,
                                  title = NULL) {
  
  library(ggplot2)
  library(igraph)
  library(ggraph)
  
  # Get significance column
  if (use_pvalue) {
    da_results$sig <- da_results$PValue < alpha
  } else {
    da_results$sig <- da_results$SpatialFDR < alpha
  }
  
  # Set non-significant logFC to 0 for coloring
  da_results$logFC_plot <- ifelse(da_results$sig, da_results$logFC, 0)
  
  # Get nhood graph
  nh_graph <- nhoodGraph(milo_obj)
  
  # Get layout coordinates from UMAP
  umap_coords <- reducedDim(milo_obj, layout)
  
  # Calculate nhood centroids
  nh_pos <- do.call(rbind, lapply(seq_len(ncol(nhoods(milo_obj))), function(i) {
    cells_in_nhood <- which(nhoods(milo_obj)[, i] == 1)
    colMeans(umap_coords[cells_in_nhood, , drop = FALSE])
  }))
  
  # Build graph layout
  layout_df <- data.frame(
    x = nh_pos[, 1],
    y = nh_pos[, 2]
  )
  
  # Node data
  node_data <- data.frame(
    nhood = seq_len(ncol(nhoods(milo_obj))),
    logFC = da_results$logFC_plot,
    sig = da_results$sig,
    size = colSums(nhoods(milo_obj))
  )
  
  # Create ggraph
  g <- ggraph(nh_graph, layout = as.matrix(layout_df)) +
    geom_edge_link0(edge_colour = "grey80", edge_width = edge_width) +
    geom_node_point(aes(fill = node_data$logFC, size = node_data$size), 
                    shape = 21, stroke = 0.1, color = "grey30") +
    scale_fill_gradient2(low = "#F8766C", mid = "white", high = "#00BEC4",
                         midpoint = 0, name = "logFC",
                         limits = c(-max(abs(da_results$logFC), na.rm = TRUE),
                                    max(abs(da_results$logFC), na.rm = TRUE))) +
    scale_size_continuous(range = node_size_range, name = "Nhood size") +
    theme_void(base_size = 11) +
    theme(
      legend.position  = "right",
      legend.title     = element_text(size = 10, face = "bold"),
      legend.text      = element_text(size = 9, face = "bold"),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 13)
    )
  
  if (!is.null(title)) {
    g <- g + ggtitle(title)
  }
  
  return(g)
}

################################################################################
# PLOT BROAD CLUSTER AND LABELED UMAP FOR PUBLISHING
################################################################################

plot_umap_hierarchical <- function(
    seurat_obj,
    fine_cluster_col = "fine_clust",
    broad_cluster_col = "mapping_cell_type",
    reduction = "umap",
    broad_labels = NULL,
    colors = NULL,
    point_size = 0.3,
    point_alpha = 0.8,
    legend_title_size = 10,
    legend_text_size = 7,
    legend_width = 0.22
) {
  
  # Default colors
  if (is.null(colors)) {
    colors <- c(
      "Keratinocytes" = "#208A42",
      "Fibroblasts" = "#D51F26",
      "Immune" = "#272E6A",
      "Endothelial" = "#89288F",
      "Remaining" = "#808080"
    )
  }
  
  # Default broad cluster full labels (for legend)
  if (is.null(broad_labels)) {
    broad_labels <- c(
      "Keratinocytes" = "KRT - Keratinocytes",
      "Fibroblasts" = "FIB - Fibroblasts",
      "Immune" = "IMM - Immune cells",
      "Endothelial" = "ENDO - Endothelial cells",
      "Remaining" = "NC - Neural-crest cells"
    )
  }
  
  # Extract embedding coordinates and metadata
  embed_df <- as.data.frame(Embeddings(seurat_obj, reduction))
  colnames(embed_df) <- c("Dim1", "Dim2")
  embed_df$fine_clust <- seurat_obj@meta.data[[fine_cluster_col]]
  embed_df$broad_clust <- seurat_obj@meta.data[[broad_cluster_col]]
  
  # Order by number of cells
  broad_order <- embed_df %>%
    count(broad_clust) %>%
    arrange(desc(n)) %>%
    pull(broad_clust)
  
  embed_df$broad_clust <- factor(embed_df$broad_clust, levels = broad_order)
  
  # Main UMAP plot (no labels)
  p <- ggplot(embed_df, aes(x = Dim1, y = Dim2, color = broad_clust)) +
    geom_point(size = point_size, alpha = point_alpha, stroke = 0) +
    scale_color_manual(values = colors) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 5, 0)
    ) +
    coord_fixed(ratio = 1)
  
  # Build legend data
  legend_df <- embed_df %>%
    group_by(broad_clust, fine_clust) %>%
    tally(name = "n_cells") %>%
    arrange(broad_clust, desc(n_cells)) %>%
    ungroup()
  
  # Create legend using grid graphics
  create_legend <- function() {
    
    n_broad <- length(broad_order)
    n_fine_per_broad <- legend_df %>% 
      group_by(broad_clust) %>% 
      summarise(n = n(), .groups = "drop")
    
    total_lines <- n_broad + sum(n_fine_per_broad$n) + (n_broad - 1) * 0.5
    
    line_height <- 0.9 / total_lines
    bar_width <- 0.02
    left_margin <- 0.03
    text_start <- left_margin + bar_width + 0.02
    
    grobs <- gList()
    y_current <- 0.97
    
    for (bc in broad_order) {
      bc_char <- as.character(bc)
      bc_color <- colors[bc_char]
      bc_label <- broad_labels[bc_char]
      bc_fine <- legend_df %>% filter(broad_clust == bc)
      n_fine <- nrow(bc_fine)
      
      section_height <- (n_fine + 1) * line_height
      
      # Colored vertical bar
      grobs <- gList(grobs, rectGrob(
        x = unit(left_margin, "npc"),
        y = unit(y_current - section_height/2, "npc"),
        width = unit(bar_width, "npc"),
        height = unit(section_height, "npc"),
        hjust = 0,
        vjust = 0.5,
        gp = gpar(fill = bc_color, col = NA)
      ))
      
      # Broad cluster title (colored, bold)
      grobs <- gList(grobs, textGrob(
        label = bc_label,
        x = unit(text_start, "npc"),
        y = unit(y_current, "npc"),
        hjust = 0,
        vjust = 0.5,
        gp = gpar(
          fontsize = legend_title_size,
          fontface = "bold",
          col = bc_color
        )
      ))
      
      y_current <- y_current - line_height
      
      # Fine clusters (black text, bold)
      for (i in seq_len(n_fine)) {
        fine_label <- paste0(
          bc_fine$fine_clust[i],
          " (", bc_fine$n_cells[i], ")"
        )
        
        grobs <- gList(grobs, textGrob(
          label = fine_label,
          x = unit(text_start + 0.02, "npc"),
          y = unit(y_current, "npc"),
          hjust = 0,
          vjust = 0.5,
          gp = gpar(
            fontsize = legend_text_size,
            fontface = "bold",
            col = "black"
          )
        ))
        
        y_current <- y_current - line_height
      }
      
      y_current <- y_current - line_height * 0.3
    }
    
    gTree(children = grobs)
  }
  
  legend_grob <- create_legend()
  
  legend_panel <- ggplot() +
    annotation_custom(legend_grob, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))
  
  # Combine
  final_plot <- plot_grid(
    legend_panel,
    p,
    rel_widths = c(legend_width, 1 - legend_width),
    nrow = 1,
    align = "h"
  )
  
  return(final_plot)
}

#########################################################################
# VIOLIN GENE EXPRESSION PLOT
#########################################################################

plot_custom_expression <- function(seurat_obj,
                                     features,
                                     idents = NULL,
                                     group.by = NULL,
                                     split.by = "datasets",
                                     pt.size = 0) {

  if (!is.null(idents)) {
    seurat_obj <- subset(seurat_obj, idents = idents)
    if (!is.null(group.by)) {
      seurat_obj@meta.data[[group.by]] <- factor(seurat_obj@meta.data[[group.by]], levels = idents)
    } else {
      Seurat::Idents(seurat_obj) <- factor(Seurat::Idents(seurat_obj), levels = idents)
    }
  }

  features <- features[features %in% rownames(seurat_obj)]

  expr_mat <- Seurat::GetAssayData(seurat_obj, layer = "data")[features, ]
  scaled_mat <- t(apply(as.matrix(expr_mat), 1, function(x) {
    if (max(x) == 0) return(x)
    (x - min(x)) / (max(x) - min(x))
  }))

  seurat_obj[["scaled"]] <- Seurat::CreateAssayObject(data = scaled_mat)
  Seurat::DefaultAssay(seurat_obj) <- "scaled"

  p <- Seurat::VlnPlot(seurat_obj,
                        features = features,
                        group.by = group.by,
                        split.by = split.by,
                        pt.size = pt.size,
                        ncol = 1,
                        stack = TRUE,
                        flip = TRUE,
                        cols = c("#F8766C", "#00BEC4")) +
    ggplot2::ylab("Scaled Expression (0-1)")

  Seurat::DefaultAssay(seurat_obj) <- "RNA"

  return(p)
}

#########################################################################
# CLEAN CELLCHAT BUBBLEPLOT
#########################################################################

clean_bubble_xaxis <- function(gg, cond1 = "PBS", cond2 = "sCD83",
                                col1 = "#F8766C", col2 = "#00BEC4") {
  gb <- ggplot_build(gg)
  x_labels <- gb$layout$panel_params[[1]]$x$get_labels()
  n_x <- length(x_labels)
  
  base_names <- gsub(paste0("\\s*\\(", cond1, "\\)|\\s*\\(", cond2, "\\)"), "", x_labels)
  is_cond1 <- grepl(cond1, x_labels, fixed = TRUE)
  
  # Strip x-axis from main plot, increase legend and title fonts
  gg_main <- gg + theme(axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.title.x = element_blank(),
                         axis.text.y = element_text(size = 10),
                         plot.title = element_text(size = 15, face = "bold"),
                         legend.title = element_text(size = 11, face = "bold"),
                         legend.text = element_text(size = 10))
  
  # Build connector data
  seg_list <- list()
  txt_list <- list()
  
  i <- 1
  while (i <= n_x) {
    if (i < n_x && base_names[i] == base_names[i + 1]) {
      mid <- i + 0.5
      seg_list <- c(seg_list, list(
        data.frame(x = i, xend = mid, y = 1, yend = 0.8, clr = col1),
        data.frame(x = i + 1, xend = mid, y = 1, yend = 0.8, clr = col2)
      ))
      txt_list <- c(txt_list, list(
        data.frame(x = mid, y = 0.7, label = base_names[i])
      ))
      i <- i + 2
    } else {
      clr <- ifelse(is_cond1[i], col1, col2)
      seg_list <- c(seg_list, list(
        data.frame(x = i, xend = i, y = 1, yend = 0.8, clr = clr)
      ))
      txt_list <- c(txt_list, list(
        data.frame(x = i, y = 0.7, label = base_names[i])
      ))
      i <- i + 1
    }
  }
  
  seg_df <- do.call(rbind, seg_list)
  txt_df <- do.call(rbind, txt_list)
  
  # Axis panel: continuous x matched to discrete positions
  axis_panel <- ggplot() +
    geom_segment(data = seg_df,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = seg_df$clr, linewidth = 0.8) +
    geom_text(data = txt_df,
              aes(x = x, y = y, label = label),
              angle = 90, hjust = 1, size = 3.5, vjust = 0.5) +
    scale_x_continuous(limits = c(0.4, n_x + 0.6), expand = c(0, 0)) +
    coord_cartesian(clip = "off", ylim = c(-0.2, 1)) +
    theme_void() +
    theme(plot.margin = unit(c(0, 0.5, 5, 0.5), "cm"))
  
  # Stack: main plot on top, axis panel below
  gg_main / axis_panel + plot_layout(heights = c(5, 1.5))
}