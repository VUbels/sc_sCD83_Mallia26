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

###############################################################
# SIMPLE CLUSTERIZATION FOR SUBCLUSTERS
###############################################################

cluster_subcluster <- function(obj, output_dir, n_genes = 2000, features = 2000, dims = 1:40, conserve.memory = FALSE) {
  
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(clustree)
  
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

###############################################################
# BLACKLIST GENE FUNCTION FOR CLUSTERING
###############################################################

get_blacklist_genes <- function(seurat_obj) {
  
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
  
  # Cell cycle genes
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
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


##################################################
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
  plotNhoodGraph(x, colour_by = res_column, layout = layout, size_range = c(1, 5), is.da = TRUE, ...) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
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